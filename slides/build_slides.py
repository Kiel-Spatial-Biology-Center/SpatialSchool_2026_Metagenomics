#!/usr/bin/env python3
"""Inline the study-design SVG and the RLQ figure into the slide deck.

Run from the repository root:
    python3 slides/build_slides.py
"""
import base64, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
TPL = ROOT / "slides" / "slides.template.html"
OUT = ROOT / "slides" / "mssm-practical-slides.html"
SVG = ROOT / "resources" / "images" / "study_design.svg"
PNG = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "slides" / "rlq_axes.png"

svg = SVG.read_text()
# strip the white backdrop so the SVG sits on the slide's own plate
svg = svg.replace('<rect width="920" height="660" fill="#ffffff"/>', "")
# let it scale with its container
svg = re.sub(r'width="920" height="660"', 'width="100%" height="auto"', svg, count=1)

b64 = base64.b64encode(PNG.read_bytes()).decode()
img = (f'<img src="data:image/png;base64,{b64}" '
       f'alt="Caecum cryosection with microsamples coloured by RLQ axis 2 (a coherent '
       f'band across the tissue) and by axis 1 (patchier, locally clustered)">')

html = TPL.read_text().replace("{{DESIGN_SVG}}", svg).replace("{{AXES_IMG}}", img)
OUT.write_text(html)
print(f"wrote {OUT}  ({OUT.stat().st_size/1e6:.2f} MB)")
