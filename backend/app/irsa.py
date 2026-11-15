"""Access to SPHEREx Spectral Images at NASA/IPAC IRSA.

* Discovery uses IRSA's IVOA SIA2 service (ObsCore columns), one query per
  configured collection, so new releases only need a config change.
* Pixels come from IRSA's on-the-fly cutout service: appending
  ``?center=ra,dec&size=deg`` to an image ``access_url`` returns a small
  multi-extension FITS with the same HDUs as the full image.
* Cutouts are immutable, so they are cached on disk by URL.
"""

from __future__ import annotations

import asyncio
import hashlib
import io
import logging
import os
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import httpx
import numpy as np
from astropy.io.votable import parse as parse_votable

from .config import settings

log = logging.getLogger(__name__)

_DETECTOR_RE = re.compile(r"D(\d)\s*$", re.IGNORECASE)


class ArchiveError(RuntimeError):
    """The archive could not be reached or returned something unusable."""


@dataclass(frozen=True)
class Exposure:
    """One SPHEREx spectral image (one detector, one exposure)."""

    obs_id: str
    detector: int
    band: str
    mjd_start: float
    mjd_end: float
    em_min_um: float
    em_max_um: float
    access_url: str
    collection: str

    @property
    def mjd(self) -> float:
        return 0.5 * (self.mjd_start + self.mjd_end)

    def cutout_url(self, ra: float, dec: float, size_deg: float) -> str:
        sep = "&" if "?" in self.access_url else "?"
        return f"{self.access_url}{sep}center={ra:.8f},{dec:.8f}&size={size_deg:.5f}"

    def to_json(self) -> dict:
        data = asdict(self)
        data["mjd"] = self.mjd
        return data


def detector_from_band(band: str) -> int:
    match = _DETECTOR_RE.search(band or "")
    if not match:
        raise ValueError(f"Unrecognised SPHEREx bandpass name {band!r}")
    return int(match.group(1))


def _column(table, name: str, default=None):
    if name in table.colnames:
        return table[name]
    return [default] * len(table)


def _as_float(value, default: float = float("nan")) -> float:
    try:
        if np.ma.is_masked(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_str(value) -> str:
    if np.ma.is_masked(value) or value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def parse_sia_votable(payload: bytes, collection: str) -> list[Exposure]:
    """Turn an SIA2 VOTable response into :class:`Exposure` rows."""
    try:
        votable = parse_votable(io.BytesIO(payload), verify="ignore")
    except Exception as exc:
        raise ArchiveError(f"IRSA returned an unreadable SIA response for {collection}") from exc

    for info in getattr(votable.resources[0], "infos", []) if votable.resources else []:
        if info.name == "QUERY_STATUS" and info.value == "ERROR":
            raise ArchiveError(f"IRSA SIA error for {collection}: {info.content}")

    table = votable.get_first_table().to_table(use_names_over_ids=True)
    exposures: list[Exposure] = []
    meters_to_um = 1e6
    for url, band, t_min, t_max, em_min, em_max, obs_id in zip(
        _column(table, "access_url"),
        _column(table, "energy_bandpassname", ""),
        _column(table, "t_min"),
        _column(table, "t_max"),
        _column(table, "em_min"),
        _column(table, "em_max"),
        _column(table, "obs_id", ""),
    ):
        url_s, band_s = _as_str(url), _as_str(band)
        start = _as_float(t_min)
        if not url_s or not np.isfinite(start):
            continue
        try:
            detector = detector_from_band(band_s)
        except ValueError:
            continue
        end = _as_float(t_max, start)
        exposures.append(
            Exposure(
                obs_id=_as_str(obs_id) or url_s.rsplit("/", 1)[-1],
                detector=detector,
                band=band_s,
                mjd_start=start,
                mjd_end=end if np.isfinite(end) else start,
                em_min_um=_as_float(em_min) * meters_to_um,
                em_max_um=_as_float(em_max) * meters_to_um,
                access_url=url_s,
                collection=collection,
            )
        )
    return exposures


def deduplicate(exposures_by_collection: list[list[Exposure]]) -> list[Exposure]:
    """Merge results from several collections, keeping the first occurrence of
    each (obs_id, detector) in collection-priority order, sorted by time."""
    seen: set[tuple[str, int]] = set()
    merged: list[Exposure] = []
    for exposures in exposures_by_collection:
        for exp in exposures:
            key = (exp.obs_id, exp.detector)
            if key in seen:
                continue
            seen.add(key)
            merged.append(exp)
    merged.sort(key=lambda e: (e.mjd, e.detector))
    return merged


class IrsaClient:
    def __init__(self, client: httpx.AsyncClient):
        self._client = client
        self._sem = asyncio.Semaphore(settings.max_concurrent_downloads)
        self._cache_dir = settings.cache_dir / "cutouts"
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        self._last_prune = 0.0

    # ------------------------------------------------------------------ query
    async def _query_collection(self, collection: str, ra: float, dec: float) -> list[Exposure]:
        params = {
            "COLLECTION": collection,
            # A tiny circle returns only images that actually contain the target.
            "POS": f"CIRCLE {ra:.8f} {dec:.8f} 0.0003",
            "RESPONSEFORMAT": "votable",
            "MAXREC": str(settings.max_images_per_request * 4),
        }
        # IRSA's SIA front end returns transient 502/503s under load; retry
        # with backoff before giving up on a collection.
        last_exc: Exception | None = None
        for attempt in range(4):
            try:
                response = await self._client.get(settings.sia_url, params=params)
                response.raise_for_status()
                return parse_sia_votable(response.content, collection)
            except httpx.HTTPStatusError as exc:
                last_exc = exc
                if exc.response.status_code < 500:
                    break
            except httpx.HTTPError as exc:
                last_exc = exc
            await asyncio.sleep(1.5 * 2**attempt)
        raise ArchiveError(f"IRSA SIA query failed for {collection}: {last_exc}")

    async def find_exposures(self, ra: float, dec: float) -> list[Exposure]:
        results = await asyncio.gather(
            *(self._query_collection(c, ra, dec) for c in settings.collections),
            return_exceptions=True,
        )
        good: list[list[Exposure]] = []
        errors: list[Exception] = []
        for collection, result in zip(settings.collections, results):
            if isinstance(result, Exception):
                log.warning("Collection %s failed: %s", collection, result)
                errors.append(result)
            else:
                good.append(result)
        if not good and errors:
            raise ArchiveError(str(errors[0]))
        return deduplicate(good)

    # ---------------------------------------------------------------- cutouts
    def _cache_path(self, url: str) -> Path:
        return self._cache_dir / (hashlib.sha256(url.encode()).hexdigest() + ".fits")

    async def fetch_cutout(self, url: str) -> bytes:
        path = self._cache_path(url)
        if path.exists():
            os.utime(path)  # LRU bookkeeping
            return path.read_bytes()

        async with self._sem:
            last_exc: Exception | None = None
            for attempt in range(3):
                try:
                    response = await self._client.get(url)
                    response.raise_for_status()
                    payload = response.content
                    if not payload.startswith(b"SIMPLE"):
                        raise ArchiveError("Cutout service did not return a FITS file")
                    break
                except (httpx.HTTPError, ArchiveError) as exc:
                    last_exc = exc
                    await asyncio.sleep(0.5 * 2**attempt)
            else:
                raise ArchiveError(f"Cutout download failed: {last_exc}")

        tmp = path.with_suffix(f".{os.getpid()}.tmp")
        tmp.write_bytes(payload)
        tmp.replace(path)
        self._maybe_prune()
        return payload

    def _maybe_prune(self) -> None:
        now = time.monotonic()
        if now - self._last_prune < 60:
            return
        self._last_prune = now
        files = sorted(self._cache_dir.glob("*.fits"), key=lambda p: p.stat().st_mtime)
        total = sum(p.stat().st_size for p in files)
        while files and total > settings.cache_max_bytes:
            victim = files.pop(0)
            total -= victim.stat().st_size
            victim.unlink(missing_ok=True)
