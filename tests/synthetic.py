"""Synthetic SPHEREx-like data for offline tests and UI development.

Generates multi-extension FITS cutouts with the same HDU layout as IRSA's
SPHEREx Spectral Image cutouts (IMAGE, FLAGS, VARIANCE, ZODI), each rotated on
the sky differently, containing a fixed star, a high proper-motion star and a
variable star. Nothing here is used by the production server.
"""

from __future__ import annotations

import io
import math

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from app.coords import mjd_to_jyear, propagate
from app.irsa import Exposure

SCALE = 6.15 / 3600.0
T0_MJD = 60800.0  # ~May 2025

FIELD_RA, FIELD_DEC = 150.0, 20.0
MOVER_PM = (-4000.0, 12000.0)  # mas/yr, about 2 SPHEREx pixels per year
STATIC_OFFSET = (0.02, -0.015)  # deg from field centre
VARIABLE_OFFSET = (-0.02, 0.01)


def sources_at(mjd: float, ra0: float = FIELD_RA, dec0: float = FIELD_DEC):
    jyear = mjd_to_jyear(mjd)
    mover = propagate(ra0, dec0, *MOVER_PM, mjd_to_jyear(T0_MJD), jyear)
    static = (ra0 + STATIC_OFFSET[0] / math.cos(math.radians(dec0)), dec0 + STATIC_OFFSET[1])
    variable = (ra0 + VARIABLE_OFFSET[0] / math.cos(math.radians(dec0)), dec0 + VARIABLE_OFFSET[1])
    var_amp = 1.0 + 0.6 * math.sin(2 * math.pi * (mjd - T0_MJD) / 300.0)
    return [(mover, 50.0), (static, 80.0), (variable, 40.0 * var_amp)]


def make_exposures(n_visits: int = 3, per_visit: int = 4, ra0: float = FIELD_RA, dec0: float = FIELD_DEC) -> list[Exposure]:
    exps = []
    for v in range(n_visits):
        for k in range(per_visit):
            mjd = T0_MJD + 182.6 * v + 1.5 * k
            for det in range(1, 7):
                exps.append(
                    Exposure(
                        obs_id=f"obs{v}_{k}",
                        detector=det,
                        band=f"SPHEREx-D{det}",
                        mjd_start=mjd,
                        mjd_end=mjd + 0.002,
                        em_min_um=0.75 + 0.7 * (det - 1),
                        em_max_um=0.75 + 0.7 * (det - 1) + 0.02,
                        access_url=f"https://example.invalid/{v}/{k}/D{det}.fits",
                        collection="synthetic",
                    )
                )
    return exps


def make_cutout(exp: Exposure, ra: float, dec: float, size_deg: float, seed: int | None = None) -> bytes:
    rng = np.random.default_rng(seed if seed is not None else abs(hash(exp.access_url)) % 2**32)
    n = max(11, int(size_deg / SCALE))
    angle = math.radians((hash(exp.obs_id) % 360))
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    c, s = math.cos(angle), math.sin(angle)
    wcs.wcs.cd = np.array([[-SCALE * c, SCALE * s], [SCALE * s, SCALE * c]])

    yy, xx = np.mgrid[0:n, 0:n]
    sigma = 0.9 + 0.15 * exp.detector
    image = rng.normal(0.0, 0.05, (n, n))
    for (sra, sdec), amp in sources_at(exp.mjd, FIELD_RA, FIELD_DEC):
        px, py = wcs.all_world2pix([[sra, sdec]], 0)[0]
        image += amp * np.exp(-((xx - px) ** 2 + (yy - py) ** 2) / (2 * sigma**2)) / (2 * math.pi * sigma**2)
    zodi = np.full((n, n), 0.3 + 0.05 * exp.detector)
    image += zodi

    flags = np.zeros((n, n), dtype=np.int32)
    flags[0, 0] = 1 << 10  # a HOT pixel
    flags[1, 1] = 1 << 21  # SOURCE — must not be masked

    flags_hdr = fits.Header()
    flags_hdr["MP_HOT"] = 10
    flags_hdr["MP_SOURCE"] = 21
    flags_hdr["MP_TRANSIENT"] = 0

    hdul = fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(image.astype(np.float32), header=wcs.to_header(), name="IMAGE"),
            fits.ImageHDU(flags, header=flags_hdr, name="FLAGS"),
            fits.ImageHDU(np.full((n, n), 0.0025, dtype=np.float32), name="VARIANCE"),
            fits.ImageHDU(zodi.astype(np.float32), name="ZODI"),
        ]
    )
    buf = io.BytesIO()
    hdul.writeto(buf)
    return buf.getvalue()


class FakeIrsa:
    """Drop-in replacement for :class:`app.irsa.IrsaClient`."""

    def __init__(self, exposures=None):
        self.exposures = exposures if exposures is not None else make_exposures()
        self.fetches = 0

    async def find_exposures(self, ra, dec):
        if abs(dec - FIELD_DEC) > 1 or abs(ra - FIELD_RA) > 1:
            return []
        return list(self.exposures)

    async def fetch_cutout(self, url):
        self.fetches += 1
        from urllib.parse import parse_qs, urlparse

        q = parse_qs(urlparse(url).query)
        ra, dec = map(float, q["center"][0].split(","))
        size = float(q["size"][0])
        base = url.split("?")[0]
        exp = next(e for e in self.exposures if e.access_url == base)
        return make_cutout(exp, ra, dec, size)
