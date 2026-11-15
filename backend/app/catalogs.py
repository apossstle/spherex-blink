"""External catalogues: Gaia DR3 overlay and CDS Sesame name resolution."""

from __future__ import annotations

import io
import math
import re
from urllib.parse import quote

import httpx
import numpy as np
from astropy.io.votable import parse as parse_votable

from .config import settings
from .coords import now_jyear, propagate

GAIA_EPOCH = 2016.0
SIMBAD_EPOCH = 2000.0


class CatalogError(RuntimeError):
    pass


async def gaia_sources(client: httpx.AsyncClient, ra: float, dec: float, radius_deg: float, max_mag: float = 19.0) -> list[dict]:
    """Gaia DR3 sources in a cone, at their DR3 epoch (J2016.0) with proper
    motions, so the browser can place them at the date of every frame."""
    radius_deg = min(radius_deg, 0.6)
    adql = (
        "SELECT TOP 3000 source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag "
        "FROM gaiadr3.gaia_source "
        f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.8f}, {dec:.8f}, {radius_deg:.6f})) "
        f"AND phot_g_mean_mag < {max_mag:.2f} ORDER BY phot_g_mean_mag"
    )
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable", "QUERY": adql}
    try:
        response = await client.post(settings.gaia_tap_url, data=params)
        response.raise_for_status()
        table = parse_votable(io.BytesIO(response.content), verify="ignore").get_first_table().to_table()
    except Exception as exc:
        raise CatalogError(f"Gaia archive query failed: {exc}") from exc

    def num(v):
        return None if np.ma.is_masked(v) or not math.isfinite(float(v)) else round(float(v), 7)

    return [
        {
            "source_id": str(row["source_id"]),
            "ra": num(row["ra"]),
            "dec": num(row["dec"]),
            "pmra": num(row["pmra"]),
            "pmdec": num(row["pmdec"]),
            "parallax": num(row["parallax"]),
            "g": num(row["phot_g_mean_mag"]),
        }
        for row in table
    ]


_SESAME_J = re.compile(r"^%J\s+([-+\d.]+)\s+([-+\d.]+)", re.M)
_SESAME_P = re.compile(r"^%P\s+([-+\d.]+)\s+([-+\d.]+)", re.M)
_SESAME_NAME = re.compile(r"^%I\.0\s+(.+)$", re.M)


def parse_sesame(text: str) -> dict:
    j = _SESAME_J.search(text)
    if not j:
        raise CatalogError("Object name not recognised by CDS Sesame")
    ra, dec = float(j.group(1)), float(j.group(2))
    result = {"ra_j2000": ra, "dec_j2000": dec, "ra": ra, "dec": dec, "pmra": None, "pmdec": None}
    p = _SESAME_P.search(text)
    if p:
        pmra, pmdec = float(p.group(1)), float(p.group(2))
        result.update(pmra=pmra, pmdec=pmdec)
        # SIMBAD positions are epoch J2000; high proper-motion stars (Barnard's
        # Star moves ~10"/yr) would otherwise fall off a small field.
        result["ra"], result["dec"] = propagate(ra, dec, pmra, pmdec, SIMBAD_EPOCH, now_jyear())
        result["epoch"] = round(now_jyear(), 2)
    name = _SESAME_NAME.search(text)
    if name:
        result["name"] = name.group(1).strip()
    return result


async def resolve_name(client: httpx.AsyncClient, name: str) -> dict:
    try:
        response = await client.get(f"{settings.sesame_url}?{quote(name.strip())}")
        response.raise_for_status()
    except httpx.HTTPError as exc:
        raise CatalogError(f"Name resolver unavailable: {exc}") from exc
    return parse_sesame(response.text)
