"""Run the app against synthetic data — for UI work without network access.

    python tests/dev_server.py   # then open http://127.0.0.1:8765/?ra=150&dec=20&gaia=1
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "backend"), str(ROOT / "tests")]

import uvicorn  # noqa: E402

import app.main as main  # noqa: E402
from app.coords import propagate  # noqa: E402
from synthetic import FIELD_DEC, FIELD_RA, MOVER_PM, T0_MJD, FakeIrsa, sources_at  # noqa: E402


async def fake_gaia(_client, ra, dec, radius, max_mag=19.0):
    from app.coords import mjd_to_jyear

    (mover, _), (static, _), (variable, _) = sources_at(T0_MJD)
    m_ra, m_dec = propagate(*mover, *MOVER_PM, mjd_to_jyear(T0_MJD), 2016.0)
    return [
        {"source_id": "1", "ra": m_ra, "dec": m_dec, "pmra": MOVER_PM[0], "pmdec": MOVER_PM[1], "parallax": 400.0, "g": 12.0},
        {"source_id": "2", "ra": static[0], "dec": static[1], "pmra": 1.0, "pmdec": -2.0, "parallax": 1.0, "g": 11.0},
        {"source_id": "3", "ra": variable[0], "dec": variable[1], "pmra": 0.0, "pmdec": 0.0, "parallax": 1.0, "g": 13.0},
    ]


async def fake_resolve(_client, name):
    return {"ra": FIELD_RA, "dec": FIELD_DEC, "name": name, "pmra": None}


main.IrsaClient = lambda client: FakeIrsa()
main.gaia_sources = fake_gaia
main.resolve_name = fake_resolve

if __name__ == "__main__":
    uvicorn.run(main.app, host="127.0.0.1", port=8765, log_level="warning")
