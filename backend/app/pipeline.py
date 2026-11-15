"""End-to-end request pipeline: find → download → align → stack → measure."""

from __future__ import annotations

import asyncio
import logging
import math
from collections import OrderedDict
from typing import Literal

import numpy as np
from astropy.time import Time

from .config import DETECTORS, settings
from .imaging import (
    Cutout,
    aperture_photometry,
    coadd,
    encode_frame,
    flags_at,
    group_into_exposures,
    group_into_visits,
    load_cutout,
    resample,
    target_grid,
    target_pixel,
    wavelength_at,
)
from .irsa import ArchiveError, Exposure, IrsaClient

log = logging.getLogger(__name__)

Mode = Literal["visit", "exposure"]


class NoDataError(LookupError):
    pass


def _iso(mjd: float) -> str:
    return Time(mjd, format="mjd").isot[:19] + "Z"


def _clean(value: float) -> float | None:
    return None if value is None or not math.isfinite(value) else round(float(value), 6)


class _ResultCache:
    """Tiny LRU for finished responses; identical requests are common when a
    link to a field is shared."""

    def __init__(self, size: int = 64):
        self._size = size
        self._data: OrderedDict[tuple, dict] = OrderedDict()

    def get(self, key: tuple) -> dict | None:
        if key in self._data:
            self._data.move_to_end(key)
            return self._data[key]
        return None

    def put(self, key: tuple, value: dict) -> None:
        self._data[key] = value
        self._data.move_to_end(key)
        while len(self._data) > self._size:
            self._data.popitem(last=False)


_cache = _ResultCache()
_inflight: dict[tuple, asyncio.Future] = {}


async def _load_all(
    irsa: IrsaClient, exposures: list[Exposure], ra: float, dec: float, cut_size: float
) -> tuple[list[Cutout], int]:
    async def one(exp: Exposure) -> Cutout | None:
        try:
            payload = await irsa.fetch_cutout(exp.cutout_url(ra, dec, cut_size))
            return await asyncio.to_thread(load_cutout, exp, payload)
        except Exception as exc:  # one bad exposure must not sink the request
            log.warning("Skipping %s: %s", exp.access_url, exc)
            return None

    results = await asyncio.gather(*(one(e) for e in exposures))
    cutouts = [c for c in results if c is not None]
    return cutouts, len(results) - len(cutouts)


def _stack_groups(
    groups: list[list[Cutout]], grid, shape, ra: float, dec: float, with_layers: bool
) -> tuple[list[dict], dict[int, list[dict]]]:
    frames: list[dict] = []
    layers: dict[int, list[dict]] = {d: [] for d in DETECTORS} if with_layers else {}

    for index, group in enumerate(groups):
        resampled = [(c, *resample(c, grid, shape)) for c in group]
        image = coadd([(d, w) for _, d, w in resampled], shape)

        flagged_at_target = 0
        flag_names: set[str] = set()
        for cutout in group:
            x, y = target_pixel(cutout, ra, dec)
            names = flags_at(cutout, x, y)
            if names:
                flag_names.update(names)
            if any(n != "SOURCE" for n in names):
                flagged_at_target += 1

        mjds = [c.exposure.mjd for c in group]
        frames.append(
            {
                "index": index,
                "mjd": float(np.mean(mjds)),
                "mjd_start": float(min(c.exposure.mjd_start for c in group)),
                "mjd_end": float(max(c.exposure.mjd_end for c in group)),
                "date": _iso(float(np.mean(mjds))),
                "n_exposures": len(group),
                "detectors": sorted({c.exposure.detector for c in group}),
                "flags_at_target": sorted(flag_names),
                "flagged_fraction": round(flagged_at_target / len(group), 3),
                "coverage": round(float(np.isfinite(image).mean()), 3),
                **encode_frame(image),
            }
        )

        if with_layers:
            for det in DETECTORS:
                subset = [(d, w) for c, d, w in resampled if c.exposure.detector == det]
                layer = coadd(subset, shape) if subset else np.full(shape, np.nan)
                layers[det].append({"n_exposures": len(subset), **encode_frame(layer)})

    return frames, layers


def _spectrum(cutouts: list[Cutout], visits: list[list[Cutout]], ra: float, dec: float) -> list[dict]:
    visit_of = {id(c): i for i, group in enumerate(visits) for c in group}
    points = []
    for cutout in cutouts:
        x, y = target_pixel(cutout, ra, dec)
        flux, err = aperture_photometry(cutout, x, y)
        if not math.isfinite(flux):
            continue
        points.append(
            {
                "wavelength_um": _clean(wavelength_at(cutout, x, y)),
                "flux_ujy": _clean(flux),
                "error_ujy": _clean(err),
                "mjd": _clean(cutout.exposure.mjd),
                "detector": cutout.exposure.detector,
                "visit": visit_of.get(id(cutout), 0),
                "flags": flags_at(cutout, x, y),
            }
        )
    points.sort(key=lambda p: (p["wavelength_um"] or 0.0, p["mjd"] or 0.0))
    return points


async def build_sequence(
    irsa: IrsaClient,
    ra: float,
    dec: float,
    fov_deg: float,
    mode: Mode = "visit",
    detectors: tuple[int, ...] = tuple(DETECTORS),
) -> dict:
    fov_deg = min(max(fov_deg, settings.min_fov_deg), settings.max_fov_deg)
    key = (round(ra, 6), round(dec, 6), round(fov_deg, 4), mode, tuple(sorted(detectors)))
    cached = _cache.get(key)
    if cached is not None:
        return cached
    # Collapse concurrent identical requests into one archive round-trip.
    if key in _inflight:
        return await asyncio.shield(_inflight[key])

    future: asyncio.Future = asyncio.get_running_loop().create_future()
    _inflight[key] = future
    try:
        result = await _build(irsa, ra, dec, fov_deg, mode, detectors)
        _cache.put(key, result)
        future.set_result(result)
        return result
    except BaseException as exc:
        future.set_exception(exc)
        future.exception()  # mark retrieved so asyncio does not warn
        raise
    finally:
        _inflight.pop(key, None)


async def _build(irsa, ra, dec, fov_deg, mode, detectors) -> dict:
    warnings: list[str] = []
    exposures = [e for e in await irsa.find_exposures(ra, dec) if e.detector in detectors]
    if not exposures:
        raise NoDataError("SPHEREx has no released images covering these coordinates yet.")

    limit = settings.max_images_per_request
    if len(exposures) > limit:
        # Keep the time sampling even rather than truncating the newest data.
        idx = np.linspace(0, len(exposures) - 1, limit).round().astype(int)
        exposures = [exposures[i] for i in sorted(set(idx))]
        warnings.append(f"Field has many exposures; using an evenly spaced subset of {len(exposures)}.")

    # Cutouts are aligned with the detector, which is rotated on the sky, so
    # request enough margin to fill the north-up output square.
    cut_size = min(fov_deg * math.sqrt(2) + 0.01, 3.4)
    cutouts, failed = await _load_all(irsa, exposures, ra, dec, cut_size)
    if failed:
        warnings.append(f"{failed} exposure(s) could not be downloaded and were skipped.")
    if not cutouts:
        raise ArchiveError("Every cutout download failed; IRSA may be busy, please retry.")

    grid, shape = target_grid(ra, dec, fov_deg)
    by_exp = {id(c.exposure): c for c in cutouts}
    loaded = [e for e in exposures if id(e) in by_exp]

    visit_groups = [[by_exp[id(e)] for e in g] for g in group_into_visits(loaded, settings.visit_gap_days)]
    if mode == "visit":
        groups = visit_groups
    else:
        groups = [[by_exp[id(e)] for e in g] for g in group_into_exposures(loaded)]

    frames, layers = await asyncio.to_thread(
        _stack_groups, groups, grid, shape, ra, dec, mode == "visit"
    )
    spectrum = await asyncio.to_thread(_spectrum, cutouts, visit_groups, ra, dec)

    return {
        "target": {"ra": ra, "dec": dec},
        "fov_deg": fov_deg,
        "mode": mode,
        "detectors": sorted(detectors),
        "grid": {
            "width": shape[1],
            "height": shape[0],
            "pixel_scale_arcsec": settings.pixel_scale_arcsec,
            "orientation": "north up, east left",
        },
        "frames": frames,
        "layers": {str(k): v for k, v in layers.items()},
        "spectrum": spectrum,
        "n_exposures": len(cutouts),
        "n_visits": len(visit_groups),
        "collections": sorted({c.exposure.collection for c in cutouts}),
        "warnings": warnings,
    }
