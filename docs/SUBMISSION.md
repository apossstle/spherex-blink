# Space Apps project page — ready-to-paste text

Replace `[TEAM]`, `[MEMBERS]` and the links before submitting.

## Project Name
SPHEREx Blink

## High-Level Project Summary
SPHEREx Blink is a free public web tool that turns NASA SPHEREx's all-sky infrared survey into a time-lapse of any patch of sky. You type coordinates or a name. The tool fetches every SPHEREx exposure of that spot from NASA/IPAC IRSA, aligns them and plays them as a short animation, so a moving brown dwarf, nearby star, asteroid or comet jumps out at a glance. Difference images, deep stacks, false-colour images from SPHEREx's six detector arrays, a Gaia DR3 overlay, flag metadata and an on-the-fly spectrum and light curve help anyone, expert or not, confirm what they see. SPHEREx returns far more images than professionals can inspect alone. Tools like this let volunteers join in, as they did with WISE and Backyard Worlds, where the public found most known brown dwarfs.

## Link to Project Demo
[YouTube / Google Slides link — ≤30 s video or ≤7 slides, public]

## Link to Final Project
[GitHub repository link] · [live app link, e.g. Hugging Face Space]

## Detailed Project Description
**What it does.** The user enters Right Ascension and Declination in degrees or sexagesimal, or an object name. SPHEREx Blink shows every SPHEREx observation of that field as a blinking animation. There are two time scales. *Survey visits* gives one stacked frame per six-month pass and suits stars and brown dwarfs. *Every exposure* shows frames hours apart and suits asteroids and comets. Display modes are Blink, Δ previous, Δ mean (static objects cancel, movers appear as bright/dark pairs) and Deep stack. Users choose the field of view (1.2′–30′), which of the six detector arrays (0.75–5.0 µm) to stack, and an asinh/log/sqrt/linear stretch. They can also build a false-colour RGB image from any three arrays.

**How it works.** A Python/FastAPI backend queries IRSA's IVOA SIA2 service for all SPHEREx Spectral Images that contain the target, across the configured releases (DR1 and QR3 by default). It downloads small cutouts in parallel from IRSA's cutout service and caches them. For each exposure it subtracts the zodiacal-light model (ZODI extension) and masks defective pixels using the bit names in the FLAGS header. Real-source flags are kept. It then resamples every exposure onto one north-up gnomonic grid at SPHEREx's 6.15″ pixel scale and coadds by inverse variance. For each exposure it also measures aperture photometry at the target and reads the exact wavelength from the spectral WCS. This produces a spectrum and light curve on the fly, coloured by visit so variability is visible. The browser receives raw float frames, so all contrast, colour and differencing changes are instant.

**Checking discoveries.** Gaia DR3 stars are drawn at their predicted position for each frame's date using their proper motions. Fast movers carry a five-year arrow, so an unexplained moving source stands out. An Aladin Lite panel shows WISE, 2MASS, DSS or Pan-STARRS imagery around the field. A frames table lists the date, exposure count, arrays, coverage and the SPHEREx flags set at the target, and flagged frames are marked.

**Sharing.** Every view is a URL (`?ra=…&dec=…&fov=…&view=…`, or `?target=R+Dra`). `&embed=1` gives a bare player for other websites, and the animation can be saved as a video. The API is open (CORS), so notebooks can use it too.

**Ready for new data.** No data is stored locally. New SPHEREx observations appear as soon as IRSA publishes them, and a future release is a one-line configuration change.

**Tested on** the challenge targets: Barnard's Star, ε Indi B, WISE 0855−0714, RW Aurigae and R Draconis. Name lookups apply proper motion from J2000 to today. Otherwise Barnard's Star, which has moved ~4′ since 2000, would fall outside a small field. Automated tests with synthetic SPHEREx-format data confirm that a moving star is recovered at its predicted position and a fixed star stays within 0.35 pixel after alignment.

**Tools:** Python, FastAPI, astropy, reproject, numpy, httpx; HTML/CSS/JavaScript (no framework); Aladin Lite v3; Docker.

**Benefits and goals.** Anyone with a browser can search SPHEREx for new nearby brown dwarfs, high proper-motion stars, moving Solar System objects and variable sources. We hope SPHEREx Blink can serve citizen-science projects the way WiseView served Backyard Worlds.

## Use of Artificial Intelligence (AI)
We used an AI coding assistant (Anthropic Claude) to generate and refactor much of the source code and documentation under our direction. The concept, design choices and data choices are ours, and we reviewed and tested the code. The repository README states this, and the web page's `generator` meta tag says the same. No AI-generated images, video or audio appear in the tool.

## NASA Data
- NASA SPHEREx Spectral Images (Quick Release QR3 and Data Release DR1), from the NASA/IPAC Infrared Science Archive (IRSA): SIA2 service and image cutout service — https://irsa.ipac.caltech.edu/Missions/spherex.html
- IRSA SPHEREx tutorial notebooks (cutouts, spectral WCS) — https://caltech-ipac.github.io/irsa-tutorials/spherex/
- WiseView (design inspiration) — http://byw.tools/wiseview

## Space Agency Partner & Other Data
- ESA Gaia Data Release 3 (Gaia Archive TAP) — positions, proper motions, magnitudes
- CDS Strasbourg: Sesame name resolver; Aladin Lite v3 with WISE, 2MASS, DSS2 and Pan-STARRS HiPS imagery
- Open-source software: astropy, reproject, numpy, FastAPI, uvicorn, httpx (all under permissive licences)

---

## 7-slide demo outline (if you choose slides)
1. **Title:** SPHEREx Blink — [TEAM], [MEMBERS], challenge “Planet X and SPHEREx”.
2. **Problem:** SPHEREx maps the whole sky every 6 months in 102 colours, giving more than a billion objects. Moving objects hide in it. Pluto and most known brown dwarfs were found by comparing images.
3. **Solution:** Screenshot of the player. Type coordinates, see the sky move.
4. **How it works:** IRSA SIA2 → cutouts → align → stack → browser (the diagram from the README).
5. **Demo, Barnard's Star / WISE 0855:** blink frames and the Δ mean image with the Gaia arrow.
6. **Spectra and variability:** R Draconis spectrum by visit, light curve.
7. **Future:** Citizen-science classification, known-asteroid overlay (JPL), automatic mover candidates, more releases.
