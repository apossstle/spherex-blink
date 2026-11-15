"""SPHEREx Blink — HTTP API and static frontend."""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager
from pathlib import Path

import httpx
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from . import __version__
from .catalogs import CatalogError, gaia_sources, resolve_name
from .config import DETECTORS, settings
from .coords import CoordinateError, parse_dec, parse_ra
from .irsa import ArchiveError, IrsaClient
from .pipeline import NoDataError, build_sequence

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

FRONTEND_DIR = Path(__file__).resolve().parents[2] / "frontend"


@asynccontextmanager
async def lifespan(app: FastAPI):
    limits = httpx.Limits(max_connections=settings.max_concurrent_downloads * 2)
    async with httpx.AsyncClient(
        timeout=settings.http_timeout_s,
        limits=limits,
        follow_redirects=True,
        headers={"User-Agent": f"spherex-blink/{__version__} (NASA Space Apps 2026)"},
    ) as client:
        app.state.http = client
        app.state.irsa = IrsaClient(client)
        yield


app = FastAPI(title="SPHEREx Blink", version=__version__, lifespan=lifespan)
app.add_middleware(GZipMiddleware, minimum_size=1024)
# Public, read-only API: allow other sites and notebooks to call it directly.
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["GET"], allow_headers=["*"])


def _parse_target(ra: str, dec: str) -> tuple[float, float]:
    try:
        return parse_ra(ra), parse_dec(dec)
    except CoordinateError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


def _parse_detectors(raw: str) -> tuple[int, ...]:
    try:
        dets = tuple(sorted({int(x) for x in raw.split(",") if x.strip()}))
    except ValueError as exc:
        raise HTTPException(422, "detectors must be a comma-separated list of 1-6") from exc
    if not dets or any(d not in DETECTORS for d in dets):
        raise HTTPException(422, "detectors must be a comma-separated list of 1-6")
    return dets


@app.exception_handler(ArchiveError)
async def _archive_error(_: Request, exc: ArchiveError):
    return JSONResponse(status_code=502, content={"detail": str(exc)})


@app.get("/api/health")
async def health():
    return {"status": "ok", "version": __version__}


@app.get("/api/config")
async def config():
    return {
        "version": __version__,
        "collections": list(settings.collections),
        "detectors": {str(k): {"min_um": v[0], "max_um": v[1]} for k, v in DETECTORS.items()},
        "pixel_scale_arcsec": settings.pixel_scale_arcsec,
        "fov_deg": {"min": settings.min_fov_deg, "max": settings.max_fov_deg, "default": 0.1},
    }


@app.get("/api/sequence")
async def sequence(
    request: Request,
    ra: str = Query(..., description="Right ascension: degrees or sexagesimal hours"),
    dec: str = Query(..., description="Declination: degrees or sexagesimal"),
    fov: float = Query(0.1, gt=0, description="Field of view [deg]"),
    mode: str = Query("visit", pattern="^(visit|exposure)$"),
    detectors: str = Query("1,2,3,4,5,6"),
):
    ra_deg, dec_deg = _parse_target(ra, dec)
    try:
        return await build_sequence(
            request.app.state.irsa, ra_deg, dec_deg, fov, mode, _parse_detectors(detectors)  # type: ignore[arg-type]
        )
    except NoDataError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/api/gaia")
async def gaia(request: Request, ra: str, dec: str, radius: float = Query(0.1, gt=0, le=0.6)):
    ra_deg, dec_deg = _parse_target(ra, dec)
    try:
        return {"epoch": 2016.0, "sources": await gaia_sources(request.app.state.http, ra_deg, dec_deg, radius)}
    except CatalogError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


@app.get("/api/resolve")
async def resolve(request: Request, name: str = Query(..., min_length=1, max_length=120)):
    try:
        return await resolve_name(request.app.state.http, name)
    except CatalogError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


if FRONTEND_DIR.exists():

    @app.get("/", include_in_schema=False)
    async def index():
        return FileResponse(FRONTEND_DIR / "index.html")

    app.mount("/", StaticFiles(directory=FRONTEND_DIR), name="frontend")
