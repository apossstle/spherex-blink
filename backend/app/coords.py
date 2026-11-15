"""Coordinate parsing and proper-motion propagation."""

from __future__ import annotations

import math
import re

from astropy import units as u
from astropy.coordinates import Angle
from astropy.time import Time


class CoordinateError(ValueError):
    """Raised when user-supplied coordinates cannot be understood."""


_SEXAGESIMAL_HINT = re.compile(r"[hms:dD°'\"]|\s")


def parse_ra(text: str) -> float:
    """Parse right ascension.

    Plain numbers are degrees, as are values marked with ``d`` or ``°``
    (``134.29d``). Anything else sexagesimal (``08h55m10.8s``, ``08:55:10.8``,
    ``08 55 10.8``) is read as hours.
    """
    value = text.strip()
    if not value:
        raise CoordinateError("Right ascension is empty")
    try:
        lowered = value.lower()
        if ("d" in lowered or "°" in value) and "h" not in lowered:
            deg = Angle(value, unit=u.deg).degree
        elif _SEXAGESIMAL_HINT.search(value):
            deg = Angle(value, unit=u.hourangle).degree
        else:
            deg = float(value)
    except Exception as exc:  # astropy raises a zoo of exception types
        raise CoordinateError(f"Cannot parse right ascension {text!r}") from exc
    if not math.isfinite(deg):
        raise CoordinateError("Right ascension must be finite")
    return deg % 360.0


def parse_dec(text: str) -> float:
    """Parse declination: plain degrees, ``-7d14m42.53s``, ``-07:14:42.5`` …"""
    value = text.strip()
    if not value:
        raise CoordinateError("Declination is empty")
    try:
        deg = Angle(value, unit=u.deg).degree if _SEXAGESIMAL_HINT.search(value) else float(value)
    except Exception as exc:
        raise CoordinateError(f"Cannot parse declination {text!r}") from exc
    if not -90.0 <= deg <= 90.0:
        raise CoordinateError("Declination must be between -90 and +90 degrees")
    return deg


def propagate(
    ra: float,
    dec: float,
    pmra_mas_yr: float,
    pmdec_mas_yr: float,
    from_jyear: float,
    to_jyear: float,
) -> tuple[float, float]:
    """Linear proper-motion propagation; ``pmra`` includes the cos(dec) factor
    (the Gaia / SIMBAD convention). Accurate to well under a SPHEREx pixel over
    decades for anything that is not within a few arcmin of a pole."""
    dt = to_jyear - from_jyear
    cosd = max(math.cos(math.radians(dec)), 1e-6)
    new_dec = dec + pmdec_mas_yr * dt / 3.6e6
    new_ra = ra + pmra_mas_yr * dt / 3.6e6 / cosd
    return new_ra % 360.0, max(-90.0, min(90.0, new_dec))


def mjd_to_jyear(mjd: float) -> float:
    return float(Time(mjd, format="mjd").jyear)


def now_jyear() -> float:
    return float(Time.now().jyear)


def format_sexagesimal(ra: float, dec: float) -> tuple[str, str]:
    ra_s = Angle(ra, unit=u.deg).to_string(unit=u.hourangle, sep="hms", precision=2, pad=True)
    dec_s = Angle(dec, unit=u.deg).to_string(unit=u.deg, sep="dms", precision=1, alwayssign=True, pad=True)
    return ra_s, dec_s
