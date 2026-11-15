"""Turning SPHEREx Spectral Image cutouts into aligned, comparable frames.

Each SPHEREx exposure lands on the sky with its own orientation and SIP
distortion, so every cutout is resampled onto one shared north-up tangent-plane
grid centred on the target. Only then can frames be blinked, stacked or
subtracted pixel-for-pixel.
"""

from __future__ import annotations

import base64
import io
import math
import os
import warnings
from dataclasses import dataclass, field

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import FITSFixedWarning, WCS
from astropy.wcs.utils import proj_plane_pixel_area
from reproject import reproject_interp

from .config import settings
from .irsa import Exposure

warnings.simplefilter("ignore", FITSFixedWarning)

# Pixels carrying these flags are excluded from coadds. Names follow the FLAGS
# definitions in the SPHEREx Explanatory Supplement; "SOURCE" is deliberately
# absent because it marks real astrophysical sources, not defects.
BAD_FLAG_NAMES = frozenset(
    os.getenv(
        "BAD_FLAG_NAMES",
        "TRANSIENT,OVERFLOW,SUR_ERROR,PHANTOM,REFERENCE,NONFUNC,MISSING_DATA,HOT,COLD,NONLINEAR,PERSIST,OUTLIER,PHANMISS",
    ).split(",")
)
# Used only if a FLAGS header does not define its own bit names.
FALLBACK_BAD_BITS = tuple(
    int(b) for b in os.getenv("BAD_FLAG_BITS", "0,1,2,6,9,10,11,14,15,17,19").split(",") if b
)


@dataclass
class Cutout:
    exposure: Exposure
    image: np.ndarray  # background-subtracted surface brightness [MJy/sr]
    bad: np.ndarray  # bool mask of pixels to ignore
    variance: np.ndarray
    wcs: WCS
    flag_names: dict[int, str]
    raw_flags: np.ndarray
    hdul_bytes: bytes = field(repr=False)


def target_grid(ra: float, dec: float, fov_deg: float) -> tuple[WCS, tuple[int, int]]:
    """North-up TAN grid at native SPHEREx sampling, odd-sized so the target
    sits exactly on the central pixel."""
    scale_deg = settings.pixel_scale_arcsec / 3600.0
    n = max(9, int(math.ceil(fov_deg / scale_deg)))
    if n % 2 == 0:
        n += 1
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.crpix = [(n + 1) / 2.0, (n + 1) / 2.0]
    wcs.wcs.cdelt = [-scale_deg, scale_deg]
    wcs.wcs.cunit = ["deg", "deg"]
    return wcs, (n, n)


def _flag_definitions(header: fits.Header) -> dict[int, str]:
    """SPHEREx FLAGS headers document their bits as ``MP_<NAME> = <bit>``."""
    names: dict[int, str] = {}
    for key, value in header.items():
        if key.startswith("MP_") and isinstance(value, int) and 0 <= value < 32:
            names[value] = key[3:]
    return names


def _bad_mask(flags: np.ndarray, names: dict[int, str]) -> np.ndarray:
    if names:
        bits = [bit for bit, name in names.items() if name in BAD_FLAG_NAMES]
    else:
        bits = list(FALLBACK_BAD_BITS)
    mask_value = 0
    for bit in bits:
        mask_value |= 1 << bit
    return (flags.astype(np.int64) & mask_value) != 0


def _get_hdu(hdul: fits.HDUList, name: str, fallback_index: int):
    try:
        return hdul[name]
    except KeyError:
        return hdul[fallback_index] if fallback_index < len(hdul) else None


def load_cutout(exposure: Exposure, payload: bytes) -> Cutout:
    with fits.open(io.BytesIO(payload), memmap=False) as hdul:
        image_hdu = _get_hdu(hdul, "IMAGE", 1)
        if image_hdu is None or image_hdu.data is None:
            raise ValueError("Cutout has no IMAGE extension")
        image = np.asarray(image_hdu.data, dtype=np.float64)
        header = image_hdu.header.copy()

        flags_hdu = _get_hdu(hdul, "FLAGS", 2)
        flags = (
            np.asarray(flags_hdu.data)
            if flags_hdu is not None and flags_hdu.data is not None
            else np.zeros(image.shape, dtype=np.int32)
        )
        flag_names = _flag_definitions(flags_hdu.header) if flags_hdu is not None else {}

        var_hdu = _get_hdu(hdul, "VARIANCE", 3)
        variance = (
            np.asarray(var_hdu.data, dtype=np.float64)
            if var_hdu is not None and var_hdu.data is not None
            else np.full(image.shape, np.nan)
        )

        zodi_hdu = _get_hdu(hdul, "ZODI", 4)
        if zodi_hdu is not None and zodi_hdu.data is not None and zodi_hdu.data.shape == image.shape:
            image = image - np.asarray(zodi_hdu.data, dtype=np.float64)

    wcs = WCS(header, naxis=2)
    bad = _bad_mask(flags, flag_names) | ~np.isfinite(image)

    # Remove whatever smooth offset remains after the zodiacal model so that
    # frames from different bands and epochs share a common zero point.
    good = image[~bad]
    if good.size > 10:
        _, median, _ = sigma_clipped_stats(good, sigma=3.0, maxiters=5)
        image = image - median

    return Cutout(exposure, image, bad, variance, wcs, flag_names, flags, payload)


def resample(cutout: Cutout, grid: WCS, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Return (image, weight) on the target grid; weight is 0 where no data."""
    clean = np.where(cutout.bad, np.nan, cutout.image)
    data, footprint = reproject_interp((clean, cutout.wcs), grid, shape_out=shape, order="bilinear")
    valid = np.isfinite(data) & (footprint > 0)

    noise = cutout.variance[~cutout.bad & np.isfinite(cutout.variance)]
    sky_var = float(np.median(noise)) if noise.size else float("nan")
    if not np.isfinite(sky_var) or sky_var <= 0:
        finite = cutout.image[~cutout.bad]
        sky_var = float(np.var(finite)) if finite.size > 1 else 1.0
    # Weight by the image's typical sky variance, not per-pixel variance:
    # per-pixel weighting would down-weight bright sources themselves.
    weight = np.where(valid, 1.0 / max(sky_var, 1e-30), 0.0)
    return np.where(valid, data, 0.0), weight


def coadd(stack: list[tuple[np.ndarray, np.ndarray]], shape: tuple[int, int]) -> np.ndarray:
    num = np.zeros(shape)
    den = np.zeros(shape)
    for data, weight in stack:
        num += data * weight
        den += weight
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out[den == 0] = np.nan
    return out


def target_pixel(cutout: Cutout, ra: float, dec: float) -> tuple[float, float]:
    x, y = cutout.wcs.all_world2pix([[ra, dec]], 0)[0]
    return float(x), float(y)


def flags_at(cutout: Cutout, x: float, y: float) -> list[str]:
    ix, iy = int(round(x)), int(round(y))
    h, w = cutout.raw_flags.shape
    if not (0 <= ix < w and 0 <= iy < h):
        return []
    value = int(cutout.raw_flags[iy, ix])
    names = []
    for bit in range(32):
        if value & (1 << bit):
            names.append(cutout.flag_names.get(bit, f"BIT{bit}"))
    return names


def wavelength_at(cutout: Cutout, x: float, y: float) -> float:
    """Central wavelength [µm] at a pixel from the spectral ("W") WCS, which
    uses a -TAB lookup stored in the WCS-WAVE extension. Falls back to the
    middle of the image's bandpass if that cannot be evaluated."""
    try:
        with fits.open(io.BytesIO(cutout.hdul_bytes), memmap=False) as hdul:
            header = _get_hdu(hdul, "IMAGE", 1).header
            spec_wcs = WCS(header, fobj=hdul, key="W")
            spec_wcs.sip = None
            world = spec_wcs.pixel_to_world(x, y)
            value = world[0] if isinstance(world, (list, tuple)) else world
            if hasattr(value, "to"):
                um = float(value.to(u.um, equivalencies=u.spectral()).value)
            else:
                um = float(np.ravel(spec_wcs.all_pix2world([[x, y]], 0))[0]) * 1e6
            if np.isfinite(um) and 0.5 < um < 5.5:
                return um
    except Exception:
        pass
    exp = cutout.exposure
    return 0.5 * (exp.em_min_um + exp.em_max_um)


def aperture_photometry(cutout: Cutout, x: float, y: float, r_ap: float = 2.0) -> tuple[float, float]:
    """Simple aperture flux [µJy] and 1-sigma error with a local annulus sky.

    This is a quick-look measurement for spotting variability and spectral
    shape, not PSF-fitted science photometry.
    """
    h, w = cutout.image.shape
    yy, xx = np.mgrid[0:h, 0:w]
    rr = np.hypot(xx - x, yy - y)
    good = ~cutout.bad
    ap = (rr <= r_ap) & good
    ann = (rr >= r_ap + 2) & (rr <= r_ap + 6) & good
    if ap.sum() < 3:
        return float("nan"), float("nan")
    sky = float(np.median(cutout.image[ann])) if ann.sum() >= 5 else 0.0
    pix_sr = proj_plane_pixel_area(cutout.wcs) * (math.pi / 180.0) ** 2
    mjy_sr_to_ujy = pix_sr * 1e12
    flux = float(np.sum(cutout.image[ap] - sky)) * mjy_sr_to_ujy
    var = cutout.variance[ap]
    var = var[np.isfinite(var)]
    err = float(np.sqrt(np.sum(var))) * mjy_sr_to_ujy if var.size else float("nan")
    return flux, err


def group_into_visits(exposures: list[Exposure], gap_days: float) -> list[list[Exposure]]:
    """Split a time-sorted list wherever consecutive exposures are separated
    by more than ``gap_days``. Each group is one survey pass over the field."""
    groups: list[list[Exposure]] = []
    for exp in sorted(exposures, key=lambda e: e.mjd):
        if groups and exp.mjd - groups[-1][-1].mjd <= gap_days:
            groups[-1].append(exp)
        else:
            groups.append([exp])
    return groups


def group_into_exposures(exposures: list[Exposure], tol_days: float = 30 / 86400) -> list[list[Exposure]]:
    """All six detectors read out simultaneously; group by observation time."""
    return group_into_visits(exposures, tol_days)


def encode_frame(image: np.ndarray) -> dict:
    """Float32 little-endian, base64. The browser does all display stretching,
    so contrast changes never need another round-trip."""
    arr = np.asarray(image, dtype="<f4")
    return {"data": base64.b64encode(arr.tobytes()).decode("ascii")}
