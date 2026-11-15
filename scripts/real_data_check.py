"""End-to-end check against the live NASA/IPAC IRSA archive.

Runs the production pipeline on the challenge's test targets and writes, per
target, an animated GIF of the aligned frames plus a JSON summary. Used by the
"Real NASA data" GitHub Actions workflow; can also be run locally:

    python scripts/real_data_check.py out/
"""

from __future__ import annotations

import asyncio
import base64
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import httpx  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

from app.catalogs import resolve_name  # noqa: E402
from app.irsa import IrsaClient  # noqa: E402
from app.pipeline import build_sequence  # noqa: E402

# (label, name for the resolver, fallback RA/Dec, field of view)
TARGETS = [
    ("barnards_star", "Barnard's Star", (269.4463, 4.7676), 0.1),
    ("wise_0855", "WISE J085510.83-071442.5", (133.7951, -7.2451), 0.15),
    ("eps_indi_b", "eps Ind Ba", (331.0959, -56.8009), 0.12),
    ("rw_aur", "RW Aur", (76.95642, 30.40141), 0.08),
    ("r_dra", "R Dra", (248.16742, 66.75497), 0.08),
]


def to_gif(result: dict, path: Path) -> None:
    w, h = result["grid"]["width"], result["grid"]["height"]
    frames = [np.frombuffer(base64.b64decode(f["data"]), dtype="<f4").reshape(h, w) for f in result["frames"]]
    stack = np.concatenate([f[np.isfinite(f)] for f in frames]) if frames else np.array([0.0])
    lo, hi = np.percentile(stack, [5, 99.5]) if stack.size else (0, 1)
    images = []
    for f in frames:
        t = np.clip((np.nan_to_num(f, nan=lo) - lo) / max(hi - lo, 1e-9), 0, 1)
        t = np.arcsinh(10 * t) / np.arcsinh(10)
        img = Image.fromarray((np.flipud(t) * 255).astype(np.uint8)).resize((w * 6, h * 6), Image.NEAREST)
        images.append(img)
    if images:
        images[0].save(path, save_all=True, append_images=images[1:], duration=700, loop=0)


async def main(out: Path) -> int:
    out.mkdir(parents=True, exist_ok=True)
    summary = []
    async with httpx.AsyncClient(timeout=120, follow_redirects=True) as client:
        irsa = IrsaClient(client)
        for label, name, fallback, fov in TARGETS:
            row: dict = {"target": label}
            try:
                try:
                    r = await resolve_name(client, name)
                    ra, dec = r["ra"], r["dec"]
                    row["resolved"] = True
                except Exception as exc:
                    ra, dec = fallback
                    row["resolved"] = f"fallback ({exc})"
                t0 = time.time()
                result = await build_sequence(irsa, ra, dec, fov, "visit")
                row.update(
                    ra=round(ra, 5), dec=round(dec, 5), seconds=round(time.time() - t0, 1),
                    exposures=result["n_exposures"], visits=result["n_visits"],
                    frames=[(f["date"][:10], f["n_exposures"], f["coverage"], f["flagged_fraction"]) for f in result["frames"]],
                    spectrum_points=len(result["spectrum"]),
                    wavelength_range=[
                        min((p["wavelength_um"] for p in result["spectrum"]), default=None),
                        max((p["wavelength_um"] for p in result["spectrum"]), default=None),
                    ],
                    collections=result["collections"], warnings=result["warnings"], ok=True,
                )
                to_gif(result, out / f"{label}.gif")
            except Exception as exc:
                row.update(ok=False, error=f"{type(exc).__name__}: {exc}")
            print(json.dumps(row), flush=True)
            summary.append(row)

    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    lines = ["| Target | OK | Exposures | Visits | Frames (date, exp.) | Spectrum pts | λ range (µm) | Time (s) |", "|---|---|---|---|---|---|---|---|"]
    for r in summary:
        if r.get("ok"):
            frames = "; ".join(f"{d} ({n})" for d, n, *_ in r["frames"])
            lam = r["wavelength_range"]
            lam_s = f"{lam[0]:.2f}–{lam[1]:.2f}" if lam[0] else "—"
            lines.append(f"| {r['target']} | yes | {r['exposures']} | {r['visits']} | {frames} | {r['spectrum_points']} | {lam_s} | {r['seconds']} |")
        else:
            lines.append(f"| {r['target']} | **no** | {r['error'][:120]} | | | | | |")
    lines += ["", "Each animation shows the aligned frames (one per survey visit), north up, east left, target at the centre:", ""]
    lines += [f"**{r['target']}**  \n![{r['target']}]({r['target']}.gif)" for r in summary if r.get("ok")]
    md = "\n".join(lines)
    (out / "README.md").write_text("# SPHEREx Blink on real NASA SPHEREx data\n\nGenerated automatically by GitHub Actions from live NASA/IPAC IRSA data.\n\n" + md + "\n")
    (out / "summary.md").write_text(md)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as fh:
            fh.write("## SPHEREx Blink on real NASA data\n\n" + md + "\n")
    return 0 if any(r.get("ok") for r in summary) else 1


if __name__ == "__main__":
    sys.exit(asyncio.run(main(Path(sys.argv[1] if len(sys.argv) > 1 else "real-data-out"))))
