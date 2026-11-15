---
title: SPHEREx Blink
emoji: 🔭
colorFrom: indigo
colorTo: yellow
sdk: docker
app_port: 7860
pinned: false
---

# SPHEREx Blink

**A public web tool for watching the infrared sky change, built on NASA's SPHEREx all-sky survey.**

Enter a position on the sky (or just a name). SPHEREx Blink finds every SPHEREx
exposure covering it in the NASA/IPAC Infrared Science Archive (IRSA), aligns
them north-up on one pixel grid and plays them as a short animation. Anything
that jumps between frames is moving: nearby stars, brown dwarfs, asteroids,
comets. Anything that brightens and fades is variable.

Built for the 2026 NASA Space Apps Challenge — *Planet X and SPHEREx*.

## What it does

| Challenge idea | How SPHEREx Blink handles it |
|---|---|
| Enter RA/Dec, see the patch as a video | Degrees or sexagesimal (`08h55m10.831s`, `-7d14m42.53s`), or an object name resolved by CDS Sesame |
| Keep up with new data releases | No local data. Every request goes to IRSA live. Releases are a config value (`SPHEREX_COLLECTIONS`), so a new release is one environment variable |
| Use the 102 bands well | Exposures from each of the six detector arrays are stacked by inverse variance. You choose which arrays to include, or build a false-colour RGB image from any three |
| Spectra on the fly | Aperture photometry at the target in every exposure, plotted against the exact wavelength from the spectral WCS. Points are coloured by survey visit, so changes show up |
| Add / subtract images from different times | **Δ previous**, **Δ mean** (everything static cancels out) and **Deep stack** views |
| Let users choose the field of view | 1.2′ to 30′ slider. Cutouts are fetched with margin, so rotated exposures still fill the square |
| Overlay Gaia / WISE | Gaia DR3 sources are placed at each frame's date using their proper motions; fast movers get a 5-year arrow. An Aladin Lite panel shows WISE, 2MASS, DSS or Pan-STARRS context |
| Coordinates in the URL, embeddable | `/?ra=…&dec=…&fov=…&view=…`, `/?target=R+Dra`, `&embed=1` for a bare player. One-click link and `<iframe>` copy |
| Show metadata like contamination flags | Per-frame table of dates, exposures, arrays, coverage and the SPHEREx FLAGS bits set at the target. Frames where the target is flagged are marked. Flagged points are faded in the plots |
| Frames for slow and fast movers | **Survey visits** (one stacked frame per 6-month pass) for stars and brown dwarfs, **Every exposure** (hours apart) for asteroids and comets |

You can also save the animation as a video, read RA/Dec and the pixel value
under the cursor, double-click to recentre, and use the keyboard: Space, ←/→,
**D** (cycle display), **G** (Gaia).

## Test targets

The challenge's test objects are one-click presets: Barnard's Star, ε Indi B,
WISE 0855−0714, RW Aurigae and R Draconis. Barnard's Star moves ~10″ a year, so
its J2000 catalogue position is now ~4′ off. Name resolution therefore
propagates SIMBAD proper motions to today before querying.

## How it works

```
browser ──► /api/sequence ─► IRSA SIA2 (spherex_dr1, spherex_qr3 …)   find exposures containing the target
                           ─► IRSA cutout service (parallel, disk cache)  small multi-extension FITS per exposure
                           ─► IMAGE − ZODI, mask FLAGS defects, subtract residual sky
                           ─► reproject onto a north-up TAN grid at 6.15″/px
                           ─► group by visit / exposure → inverse-variance coadd
                           ─► aperture photometry + spectral WCS → spectrum points
        ◄── float32 frames (base64) + metadata
browser: stretch, colour map, differencing, RGB, Gaia propagation, plots — all client side
```

* **Backend:** Python, FastAPI, astropy, reproject, httpx (async, up to 16 parallel downloads).
* **Frontend:** plain HTML/CSS/JS with no build step and no framework. Aladin Lite v3 is used for context.
* Flag handling reads the bit names from each FLAGS header (`MP_*` keywords).
  Defect bits (HOT, COLD, NONFUNC, PERSIST …) are masked, while `SOURCE` is kept
  because it marks real astronomical sources. You can override the masked bits
  with `BAD_FLAG_NAMES`.
* Identical concurrent requests are merged into one, and finished results are cached in memory.

The spectrum and light curve use quick-look aperture photometry: a 2-pixel radius
with a local sky annulus. They are good enough to spot variability and spectral
shape, but they are not PSF-fitted science photometry.

## Run it

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --app-dir backend --port 8000
# open http://localhost:8000/?target=Barnard%27s%20Star
```

The first look at a new field takes 30–90 s while IRSA builds the cutouts.
Repeat views come from the cache and are instant.

### Docker / free hosting

```bash
docker build -t spherex-blink . && docker run -p 7860:7860 spherex-blink
```

The repo is ready for a **Hugging Face Space** (SDK: Docker): the front matter at
the top of this README sets that up. Push the repo to a new Space and it builds
on its own. The same Dockerfile works on Render, Fly.io or Railway, which supply `PORT`.

### Configuration

| Variable | Default | Purpose |
|---|---|---|
| `SPHEREX_COLLECTIONS` | `spherex_dr1,spherex_qr3` | IRSA SIA2 collections, highest priority first |
| `MAX_IMAGES_PER_REQUEST` | `600` | Evenly thins very deep fields |
| `MAX_FOV_DEG` | `0.5` | Largest field allowed |
| `VISIT_GAP_DAYS` | `20` | Gap that separates two survey visits |
| `CACHE_DIR`, `CACHE_MAX_MB` | `/tmp/…`, `2048` | Cutout disk cache |
| `BAD_FLAG_NAMES` | defects only | FLAGS bits excluded from stacks |

### Tests

```bash
pip install -r requirements-dev.txt
pytest
python tests/dev_server.py   # UI against synthetic data, no network needed
```

The tests build synthetic SPHEREx-format cutouts (rotated WCS, zodi, flags,
variance) with a fixed star, a high proper-motion star and a variable star. They
check that after alignment the moving star shows up where its proper motion
predicts and the fixed star stays put (both within 0.35 px), that SOURCE pixels
are never masked, and that the API rejects bad input.

## API

| Endpoint | Description |
|---|---|
| `GET /api/sequence?ra&dec&fov&mode=visit\|exposure&detectors=1,…,6` | Aligned frames, per-array layers, spectrum points, metadata |
| `GET /api/gaia?ra&dec&radius` | Gaia DR3 sources (epoch 2016.0) with proper motions |
| `GET /api/resolve?name=` | Name → coordinates, moved forward for proper motion |
| `GET /api/config` | Collections, detector wavelength ranges |

The API is CORS-open and read-only, so notebooks and other sites can use it directly.

## Data and credits

* **NASA SPHEREx** Spectral Images via **NASA/IPAC Infrared Science Archive (IRSA)**:
  SIA2 discovery and the cutout service. When using SPHEREx data, please cite
  the release DOIs listed by IRSA.
* **ESA Gaia DR3** via the ESA Gaia Archive TAP service.
* **CDS Strasbourg:** Sesame name resolver and Aladin Lite (WISE, 2MASS, DSS, Pan-STARRS HiPS).
* Software: astropy, reproject, numpy, FastAPI, httpx.
* Inspired by WiseView (byw.tools) and the Backyard Worlds project.

This is not an official NASA product and uses no NASA logos or mission insignia.

## Use of AI

This project was developed with help from an AI coding assistant (Anthropic Claude),
which generated and refactored much of the source code and documentation under
the team's direction. The team chose the project concept, design decisions and
data choices, and reviewed and tested the code. No AI-generated images, video or
audio are used in the tool. Source files are marked by this README, and the HTML
`generator` meta tag says the same.

## License

MIT, see `LICENSE`.
