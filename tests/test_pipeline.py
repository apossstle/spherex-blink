import base64
import math

import numpy as np
import pytest
from fastapi.testclient import TestClient

from app import pipeline
from app.catalogs import parse_sesame
from app.coords import CoordinateError, parse_dec, parse_ra, propagate
from app.imaging import group_into_visits, load_cutout, target_grid
from app.irsa import deduplicate, detector_from_band, parse_sia_votable

from synthetic import FIELD_DEC, FIELD_RA, MOVER_PM, FakeIrsa, make_cutout, make_exposures, sources_at


# ----------------------------------------------------------------- coordinates
@pytest.mark.parametrize(
    "text, expected",
    [("134.0", 134.0), ("08h55m10.831s", 133.795129), ("08:55:10.831", 133.795129), ("08 55 10.831", 133.795129), ("134.5d", 134.5), ("-10", 350.0)],
)
def test_parse_ra(text, expected):
    assert parse_ra(text) == pytest.approx(expected, abs=1e-5)


@pytest.mark.parametrize(
    "text, expected",
    [("-7d14m42.53s", -7.245147), ("-07:14:42.53", -7.245147), ("+30 24 05.1", 30.401417), ("4.5", 4.5)],
)
def test_parse_dec(text, expected):
    assert parse_dec(text) == pytest.approx(expected, abs=1e-5)


@pytest.mark.parametrize("bad", ["", "abc", "95", "-91d"])
def test_parse_dec_rejects(bad):
    with pytest.raises(CoordinateError):
        parse_dec(bad)


def test_propagate_barnards_star():
    # SIMBAD J2000 position and proper motion of Barnard's Star.
    ra, dec = propagate(269.452076, 4.693364, -802.8, 10362.5, 2000.0, 2025.9)
    assert dec == pytest.approx(4.693364 + 10.3625 * 25.9 / 3600, abs=1e-6)
    assert ra < 269.452076


def test_parse_sesame_propagates_proper_motion():
    text = "# Barnard's star\n%I.0 V* V2500 Oph\n%J 269.45207 +04.69336 = 17:57:48.49 +04:41:36.1\n%P -802.803 +10362.542 [0.06 0.04 0] A 2020yCat.1350....0G\n"
    res = parse_sesame(text)
    assert res["ra_j2000"] == pytest.approx(269.45207)
    assert res["dec"] > 4.76  # moved ~4.4 arcmin north since 2000
    assert res["name"] == "V* V2500 Oph"


# --------------------------------------------------------------------- archive
SIA_XML = b"""<?xml version="1.0"?>
<VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
<RESOURCE type="results"><INFO name="QUERY_STATUS" value="OK"/>
<TABLE>
<FIELD name="access_url" datatype="char" arraysize="*"/>
<FIELD name="energy_bandpassname" datatype="char" arraysize="*"/>
<FIELD name="t_min" datatype="double"/><FIELD name="t_max" datatype="double"/>
<FIELD name="em_min" datatype="double"/><FIELD name="em_max" datatype="double"/>
<FIELD name="obs_id" datatype="char" arraysize="*"/>
<DATA><TABLEDATA>
<TR><TD>https://irsa/a.fits</TD><TD>SPHEREx-D2</TD><TD>60800.1</TD><TD>60800.2</TD><TD>1.1e-6</TD><TD>1.15e-6</TD><TD>X1</TD></TR>
<TR><TD>https://irsa/b.fits</TD><TD>SPHEREx-D5</TD><TD>60801.1</TD><TD>60801.2</TD><TD>4.0e-6</TD><TD>4.1e-6</TD><TD>X2</TD></TR>
<TR><TD>https://irsa/c.fits</TD><TD>weird</TD><TD>60801.1</TD><TD>60801.2</TD><TD>4.0e-6</TD><TD>4.1e-6</TD><TD>X3</TD></TR>
</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>"""


def test_parse_sia_votable():
    rows = parse_sia_votable(SIA_XML, "spherex_dr1")
    assert [r.detector for r in rows] == [2, 5]
    assert rows[0].em_min_um == pytest.approx(1.1)
    assert rows[0].cutout_url(1, 2, 0.1).startswith("https://irsa/a.fits?center=1.00000000,2.00000000&size=0.1")


def test_deduplicate_prefers_first_collection():
    a = parse_sia_votable(SIA_XML, "dr1")
    b = parse_sia_votable(SIA_XML, "qr3")
    merged = deduplicate([a, b])
    assert len(merged) == 2 and all(e.collection == "dr1" for e in merged)


def test_detector_from_band():
    assert detector_from_band("SPHEREx-D6") == 6
    with pytest.raises(ValueError):
        detector_from_band("D9x")


# ---------------------------------------------------------------------- images
def test_cutout_loading_masks_defects_not_sources():
    exp = make_exposures(1, 1)[0]
    c = load_cutout(exp, make_cutout(exp, FIELD_RA, FIELD_DEC, 0.1))
    assert c.bad[0, 0] and not c.bad[1, 1]
    assert c.flag_names[21] == "SOURCE"
    assert abs(np.nanmedian(c.image)) < 0.05  # zodi + offset removed


def test_visit_grouping():
    exps = make_exposures(3, 4)
    groups = group_into_visits(exps, 20)
    assert len(groups) == 3 and all(len(g) == 24 for g in groups)


def test_target_grid_is_odd_and_centred():
    wcs, shape = target_grid(10, 10, 0.1)
    assert shape[0] % 2 == 1
    x, y = wcs.all_world2pix([[10, 10]], 0)[0]
    assert x == pytest.approx((shape[1] - 1) / 2) and y == pytest.approx((shape[0] - 1) / 2)


def _decode(frame, shape):
    return np.frombuffer(base64.b64decode(frame["data"]), dtype="<f4").reshape(shape)


def _centroid(img, x0, y0, r=3):
    h, w = img.shape
    yy, xx = np.mgrid[0:h, 0:w]
    m = (np.hypot(xx - x0, yy - y0) <= r) & np.isfinite(img)
    wts = np.clip(np.where(m, img, 0), 0, None)
    return (wts * xx).sum() / wts.sum(), (wts * yy).sum() / wts.sum()


@pytest.mark.asyncio
async def test_sequence_tracks_moving_star_and_holds_static_star():
    pipeline._cache = pipeline._ResultCache()
    irsa = FakeIrsa()
    result = await pipeline.build_sequence(irsa, FIELD_RA, FIELD_DEC, 0.1, "visit")
    assert result["n_visits"] == 3 and len(result["frames"]) == 3
    shape = (result["grid"]["height"], result["grid"]["width"])
    wcs, _ = target_grid(FIELD_RA, FIELD_DEC, result["fov_deg"])

    positions = []
    for frame in result["frames"]:
        img = _decode(frame, shape)
        (mover, _), (static, _), _ = sources_at(frame["mjd"])
        mx, my = wcs.all_world2pix([mover], 0)[0]
        sx, sy = wcs.all_world2pix([static], 0)[0]
        cmx, cmy = _centroid(img, mx, my)
        csx, csy = _centroid(img, sx, sy)
        assert math.hypot(cmx - mx, cmy - my) < 0.35  # star found where it should be
        assert math.hypot(csx - sx, csy - sy) < 0.35
        positions.append((cmx, cmy))

    # North is up: a +Dec proper motion moves the star to larger y.
    assert positions[-1][1] - positions[0][1] > 1.5
    # Six detector layers per visit are available for colour composites.
    assert set(result["layers"]) == {"1", "2", "3", "4", "5", "6"}
    assert len(result["spectrum"]) == len(irsa.exposures)


@pytest.mark.asyncio
async def test_exposure_mode_and_detector_filter():
    pipeline._cache = pipeline._ResultCache()
    result = await pipeline.build_sequence(FakeIrsa(), FIELD_RA, FIELD_DEC, 0.08, "exposure", (2, 3))
    assert len(result["frames"]) == 12
    assert all(f["detectors"] == [2, 3] for f in result["frames"])
    assert result["layers"] == {}


@pytest.mark.asyncio
async def test_repeat_request_is_cached():
    pipeline._cache = pipeline._ResultCache()
    irsa = FakeIrsa()
    await pipeline.build_sequence(irsa, FIELD_RA, FIELD_DEC, 0.05)
    n = irsa.fetches
    await pipeline.build_sequence(irsa, FIELD_RA, FIELD_DEC, 0.05)
    assert irsa.fetches == n


# ------------------------------------------------------------------------- API
@pytest.fixture
def client():
    from app.main import app

    with TestClient(app) as c:
        app.state.irsa = FakeIrsa()
        pipeline._cache = pipeline._ResultCache()
        yield c


def test_api_sequence(client):
    r = client.get("/api/sequence", params={"ra": "10h00m00s", "dec": "+20d00m00s", "fov": 0.06})
    assert r.status_code == 200, r.text
    assert len(r.json()["frames"]) == 3


def test_api_errors(client):
    assert client.get("/api/sequence", params={"ra": "x", "dec": "1"}).status_code == 422
    assert client.get("/api/sequence", params={"ra": "1", "dec": "1", "detectors": "7"}).status_code == 422
    assert client.get("/api/sequence", params={"ra": "1", "dec": "-50"}).status_code == 404


def test_frontend_served(client):
    r = client.get("/")
    assert r.status_code == 200 and "SPHEREx" in r.text
