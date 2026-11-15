"""Runtime configuration.

Everything that changes when SPHEREx publishes a new data release lives here and
can be overridden with environment variables, so a new release is picked up
without touching code (e.g. ``SPHEREX_COLLECTIONS=spherex_dr2,spherex_dr1``).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


def _env_list(name: str, default: str) -> tuple[str, ...]:
    raw = os.getenv(name, default)
    return tuple(item.strip() for item in raw.split(",") if item.strip())


@dataclass(frozen=True)
class Settings:
    # IRSA SIA2 collections, highest priority first. When the same exposure is
    # present in several releases, the one from the earliest collection wins.
    collections: tuple[str, ...] = field(
        default_factory=lambda: _env_list("SPHEREX_COLLECTIONS", "spherex_dr1,spherex_qr3")
    )
    sia_url: str = os.getenv("IRSA_SIA_URL", "https://irsa.ipac.caltech.edu/SIA")
    irsa_base: str = os.getenv("IRSA_BASE_URL", "https://irsa.ipac.caltech.edu")
    gaia_tap_url: str = os.getenv("GAIA_TAP_URL", "https://gea.esac.esa.int/tap-server/tap/sync")
    sesame_url: str = os.getenv("SESAME_URL", "https://cds.unistra.fr/cgi-bin/nph-sesame/-oIp/SNV")

    cache_dir: Path = Path(os.getenv("CACHE_DIR", "/tmp/spherex-blink-cache"))
    cache_max_bytes: int = int(os.getenv("CACHE_MAX_MB", "2048")) * 1024 * 1024

    http_timeout_s: float = float(os.getenv("HTTP_TIMEOUT_S", "60"))
    max_concurrent_downloads: int = int(os.getenv("MAX_CONCURRENT_DOWNLOADS", "6"))
    max_images_per_request: int = int(os.getenv("MAX_IMAGES_PER_REQUEST", "600"))

    # SPHEREx native pixel scale is 6.15 arcsec.
    pixel_scale_arcsec: float = 6.15
    min_fov_deg: float = 0.02
    max_fov_deg: float = float(os.getenv("MAX_FOV_DEG", "0.5"))

    # Epochs: exposures separated by more than this gap start a new "visit".
    # SPHEREx covers a given point over a span of days each survey pass and
    # comes back about six months later.
    visit_gap_days: float = float(os.getenv("VISIT_GAP_DAYS", "20"))


settings = Settings()


# SPHEREx has six detector arrays; each carries a linear variable filter that
# spans 17 of the 102 spectral channels.
DETECTORS: dict[int, tuple[float, float]] = {
    1: (0.75, 1.09),
    2: (1.10, 1.62),
    3: (1.63, 2.41),
    4: (2.42, 3.82),
    5: (3.83, 4.41),
    6: (4.42, 5.00),
}
