#!/usr/bin/env python3
"""Generate the MSSM study-design diagram as a standalone SVG."""
import random, math

W, H = 920, 660
INK, MUTED, RULE = "#1f2933", "#5b6b7a", "#cbd5dd"
CAEC, COLON, BLUE = "#2A9D8F", "#E76F51", "#4a6ab7"
PANEL, HIBG, HIBR = "#f7f9fb", "#fdf4e6", "#e0a458"
FONT = "'Helvetica Neue',Helvetica,Arial,sans-serif"

o = []
A = o.append


def txt(x, y, s, size=12, fill=INK, weight="normal", anchor="start", style=""):
    A(f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" '
      f'fill="{fill}" font-weight="{weight}" text-anchor="{anchor}"{style}>{s}</text>')


def panel(x, y, w, h, fill=PANEL, stroke=RULE, rx=8, sw=1):
    A(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" '
      f'stroke="{stroke}" stroke-width="{sw}"/>')


def step(x, y, n, title):
    A(f'<circle cx="{x + 10}" cy="{y}" r="10" fill="{INK}"/>')
    txt(x + 10, y + 4, str(n), 12, "#ffffff", "bold", "middle")
    txt(x + 28, y + 5, title, 13.5, INK, "bold")


def arrow(x1, y1, x2, y2, colour=MUTED):
    A(f'<path d="M{x1},{y1} L{x2},{y2}" stroke="{colour}" stroke-width="1.6" '
      f'fill="none" marker-end="url(#ah)"/>')


A(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
  f'height="{H}" role="img" aria-label="MSSM study design">')
A(f'<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
  f'markerHeight="6" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{MUTED}"/>'
  f'</marker></defs>')
A(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')

LX, LW = 24, 430
RX, RW = 486, 410

# ---------------------------------------------------------------- 1. tissue
step(LX, 26, 1, "Tissue")
panel(LX, 40, LW, 96)
txt(LX + 16, 62, "Chicken G121", 12, INK, "bold")
for i, (lab, col) in enumerate([("Caecum", CAEC), ("Colon", COLON)]):
    bx = LX + 16 + i * 108
    A(f'<rect x="{bx}" y="{74}" width="96" height="26" rx="13" fill="{col}" '
      f'opacity="0.16" stroke="{col}" stroke-width="1.2"/>')
    txt(bx + 48, 91, lab, 11.5, col, "bold", "middle")
txt(LX + 16, 122, "snap-frozen, cryosectioned and photographed", 11, MUTED)
A(f'<rect x="{LX + 250}" y="{58}" width="70" height="60" rx="4" fill="#e8edf2" '
  f'stroke="{RULE}"/>')
txt(LX + 285, 92, "2 mm", 10, MUTED, "normal", "middle")
txt(LX + 285, 130, "cryosection", 10, MUTED, "normal", "middle")

# --------------------------------------------------- 2. laser microdissection
step(LX, 166, 2, "Laser microdissection")
panel(LX, 180, LW, 158)

random.seed(7)
cx, cy, rr = LX + 92, 262, 66
A(f'<circle cx="{cx}" cy="{cy}" r="{rr}" fill="#eef2f6" stroke="{RULE}"/>')

DOT = 3.0
pts = []
while len(pts) < 62:
    a = random.uniform(0, 2 * math.pi)
    d = (rr - DOT - 3) * math.sqrt(random.random())
    px, py = cx + d * math.cos(a), cy + d * math.sin(a)
    if all((px - qx) ** 2 + (py - qy) ** 2 > 100 for qx, qy in pts):
        pts.append((px, py))
for px, py in pts:
    A(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{DOT}" fill="{CAEC}" '
      f'opacity="0.75"/>')

# distance from one microsample to the gut wall
sx, sy = min(pts, key=lambda q: (q[0] - (cx - 16)) ** 2 + (q[1] - (cy - 26)) ** 2)
ang = math.atan2(sy - cy, sx - cx)
wx, wy = cx + rr * math.cos(ang), cy + rr * math.sin(ang)
A(f'<path d="M{sx:.1f},{sy:.1f} L{wx:.1f},{wy:.1f}" stroke="{INK}" '
  f'stroke-width="1.2" stroke-dasharray="2 2"/>')
A(f'<circle cx="{wx:.1f}" cy="{wy:.1f}" r="2.4" fill="{INK}"/>')
lbx, lby = sx + 8, (sy + wy) / 2 - 13
A(f'<rect x="{lbx:.1f}" y="{lby:.1f}" width="70" height="25" rx="4" '
  f'fill="#ffffff" opacity="0.9"/>')
txt(lbx + 5, lby + 11, "distance to", 9, INK)
txt(lbx + 5, lby + 21, "gut wall", 9, INK)

tx = LX + 178
for i, (lab, val) in enumerate([
        ("microsamples per section", "168"),
        ("area of each", "5,000 µm²"),
        ("recorded per microsample", "X, Y position"),
        ("sequencing depth", "≈ biomass")]):
    y = 216 + i * 30
    txt(tx, y, val, 12, INK, "bold")
    txt(tx, y + 13, lab, 10, MUTED)

# ------------------------------------------- 3. libraries (the unusual step)
step(LX, 368, 3, "Libraries and sequencing")
panel(LX, 382, LW, 240, HIBG, HIBR, 8, 1.5)
txt(LX + 16, 404, "The step that makes MSSM unusual", 12, "#a8641c", "bold")

fills = [0.86, 0.30, 0.62, 0.14, 0.95]
for i, f in enumerate(fills):
    bx = LX + 22 + i * 34
    A(f'<rect x="{bx}" y="{418}" width="20" height="58" rx="3" fill="#ffffff" '
      f'stroke="{MUTED}" stroke-width="1"/>')
    fh = 54 * f
    A(f'<rect x="{bx + 2}" y="{474 - fh:.1f}" width="16" height="{fh:.1f}" rx="2" '
      f'fill="{BLUE}" opacity="0.8"/>')
txt(LX + 22, 490, "same 15 PCR cycles · different yields", 10, MUTED)

arrow(LX + 200, 447, LX + 232, 447)
txt(LX + 216, 436, "pool", 9.5, MUTED, "normal", "middle")
txt(LX + 216, 462, "as-is", 9.5, MUTED, "normal", "middle")

lane_x, lane_w = LX + 244, 168
A(f'<rect x="{lane_x}" y="{418}" width="{lane_w}" height="58" rx="4" '
  f'fill="#ffffff" stroke="{MUTED}"/>')
tot = sum(fills)
xoff = lane_x
for f in fills:
    seg = lane_w * f / tot
    A(f'<rect x="{xoff:.1f}" y="{419}" width="{seg:.1f}" height="56" '
      f'fill="{BLUE}" opacity="0.8" stroke="#ffffff" stroke-width="1.2"/>')
    xoff += seg
txt(lane_x + lane_w / 2, 490, "one run per cryosection", 10, MUTED, "normal", "middle")

A(f'<line x1="{LX + 16}" y1="506" x2="{LX + LW - 16}" y2="506" stroke="{HIBR}" '
  f'stroke-dasharray="3 3"/>')
for i, (bad, good) in enumerate([
        ("Standard metagenomics", "MSSM"),
        ("quantify, pool equimolarly", "no normalisation"),
        ("every sample ≈ equal reads", "reads follow DNA in"),
        ("depth = nuisance", "depth ≈ biomass")]):
    y = 524 + i * 22
    txt(LX + 22, y, bad, 10, MUTED if i else "#8a99a8",
        "normal" if i else "bold")
    txt(LX + 246, y, good, 10, INK if i else "#a8641c", "bold")
txt(LX + 232, 524, "→", 10, MUTED, "normal", "middle")

# -------------------------------------------------- 4. mapping and filtering
step(RX, 26, 4, "Mapping and filtering")
panel(RX, 40, RW, 298)
txt(RX + 16, 64, "223 MAGs", 12, INK, "bold")
txt(RX + 16, 78, "assembled beforehand from macro-scale sequencing", 10, MUTED)
A(f'<rect x="{RX + 16}" y="{90}" width="{RW - 32}" height="1" fill="{RULE}"/>')

boxes = [("reads per microsample", "raw sequencing output", "#eef2f6"),
         ("read counts  ·  covered bases", "one value per genome × microsample", "#eef2f6"),
         ("breadth-of-coverage filter  ≥ 30 %", "removes cross-mapping noise", "#fdf4e6"),
         ("genome counts", "length-normalised community table", "#e9f3f1")]
for i, (t, sub, bg) in enumerate(boxes):
    y = 108 + i * 62
    A(f'<rect x="{RX + 16}" y="{y}" width="{RW - 32}" height="40" rx="5" '
      f'fill="{bg}" stroke="{RULE}"/>')
    txt(RX + 28, y + 18, t, 11.5, INK, "bold")
    txt(RX + 28, y + 32, sub, 9.5, MUTED)
    if i < 3:
        arrow(RX + RW / 2, y + 42, RX + RW / 2, y + 59)

# --------------------------------------------------------- 5. the five tables
step(RX, 368, 5, "Five tables → extended RLQ")
panel(RX, 382, RW, 240)

tables = [("E", "environment", "depth ≈ biomass, richness,\ndistance to gut wall", BLUE),
          ("S", "space", "Gabriel neighbour graph\nof microsample positions", BLUE),
          ("L", "community", "genome counts,\nHellinger-transformed", INK),
          ("T", "traits", "metabolic capacity\n(GIFT), ordinated", CAEC),
          ("P", "phylogeny", "distances on the\nMAG tree", CAEC)]
for i, (k, name, desc, col) in enumerate(tables):
    y = 400 + i * 34
    A(f'<rect x="{RX + 16}" y="{y}" width="24" height="24" rx="4" fill="{col}" '
      f'opacity="0.16" stroke="{col}"/>')
    txt(RX + 28, y + 17, k, 13, col, "bold", "middle")
    txt(RX + 50, y + 11, name, 11, INK, "bold")
    d1, d2 = desc.split("\n")
    txt(RX + 118, y + 10, d1, 9.5, MUTED)
    txt(RX + 118, y + 21, d2, 9.5, MUTED)

A(f'<rect x="{RX + 16}" y="{578}" width="{RW - 32}" height="30" rx="5" '
  f'fill="{INK}"/>')
txt(RX + RW / 2, 597, "axes where space and genome biology covary", 11,
    "#ffffff", "bold", "middle")

A('</svg>')

out = "resources/images/study_design.svg"
open(out, "w").write("\n".join(o))
print("wrote", out)
