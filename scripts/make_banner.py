"""Rebuild docs/banner.png from the project's own logo and wordmark.

    python scripts/make_banner.py

Kept because the banner is generated rather than drawn: replace assets/logo.png or
assets/titulo.png and run this again. Both are dark artwork on transparency, so they are
repainted here -- pasted unchanged onto a dark banner they show nothing at all.
"""
import os
from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSETS = os.path.join(ROOT, "src", "poliscreen", "assets")


def _fonts() -> str:
    """Where DejaVu lives. Asked for rather than written down, so this runs on anyone's machine."""
    try:
        import matplotlib
        found = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
        if os.path.exists(os.path.join(found, "DejaVuSans.ttf")):
            return found
    except Exception:
        pass
    for path in ("/usr/share/fonts/truetype/dejavu",
                 "/usr/share/fonts/dejavu",
                 "C:/Windows/Fonts"):
        if os.path.exists(os.path.join(path, "DejaVuSans.ttf")):
            return path
    raise SystemExit("DejaVu fonts not found: pip install matplotlib, or install fonts-dejavu.")


FONTS = _fonts()

S = 2                                  # drawn at 2x, downsampled at the end
W, H = 1200 * S, 430 * S

BG_TOP, BG_BOT = (11, 42, 46), (8, 24, 31)
CARD = (14, 52, 56)
TEAL = (94, 234, 212)
TEXT = (226, 240, 240)
DIM = (140, 174, 178)

def font(name, size):
    return ImageFont.truetype(os.path.join(FONTS, name), size * S)

BOLD, REG = "DejaVuSans-Bold.ttf", "DejaVuSans.ttf"

img = Image.new("RGB", (W, H), BG_BOT)
d = ImageDraw.Draw(img)
for y in range(H):                     # vertical gradient
    t = y / H
    d.line([(0, y), (W, y)], fill=tuple(int(a + (b - a) * t) for a, b in zip(BG_TOP, BG_BOT)))

def tinted(path, colour, height, alpha=255):
    """The assets are dark artwork on transparency; on a dark banner they have to be repainted."""
    src = Image.open(path).convert("RGBA")
    w = int(height * src.width / src.height)
    src = src.resize((w, height), Image.LANCZOS)
    flat = Image.new("RGBA", src.size, colour + (0,))
    a = src.split()[3].point(lambda v: int(v * alpha / 255))
    flat.putalpha(a)
    return flat

# the logo, large and faint on the right, as a watermark
logo = tinted(os.path.join(ASSETS, "logo.png"), TEAL, 232 * S, alpha=52)
img.paste(logo, (W - logo.width - 56 * S, 8 * S), logo)

# the wordmark
title = tinted(os.path.join(ASSETS, "titulo.png"), (240, 253, 250), 74 * S)
img.paste(title, (60 * S, 46 * S), title)

d.text((62 * S, 138 * S),
       "Reproducible virtual screening  ·  objective interaction-quality scoring",
       font=font(REG, 17), fill=DIM)
d.rounded_rectangle([62 * S, 176 * S, 190 * S, 180 * S], radius=2 * S, fill=TEAL)

STAGES = [
    ("1", "Receptor & cavity", ["PDB / fpocket", "co-crystal control"], (56, 189, 248)),
    ("2", "Design & filter",   ["reactions, R-group", "synthesizability"], (129, 140, 248)),
    ("3", "Docking",           ["AutoDock Vina", "ADCP for peptides"],  (168, 130, 248)),
    ("4", "Interactions",      ["PLIP fingerprint", "per-cavity score"], (52, 211, 153)),
    ("5", "Rank & confidence", ["effectiveness %", "orthogonal metric"], (250, 204, 21)),
]

x0, y0 = 60 * S, 212 * S
bw, bh, gap = 196 * S, 116 * S, 22 * S
for i, (num, name, lines, accent) in enumerate(STAGES):
    x = x0 + i * (bw + gap)
    d.rounded_rectangle([x, y0, x + bw, y0 + bh], radius=10 * S, fill=CARD,
                        outline=(30, 82, 88), width=1 * S)
    d.ellipse([x + 14 * S, y0 + 14 * S, x + 34 * S, y0 + 34 * S], fill=accent)
    d.text((x + 21 * S, y0 + 17 * S), num, font=font(BOLD, 11), fill=(10, 30, 34))
    d.text((x + 14 * S, y0 + 44 * S), name, font=font(BOLD, 13), fill=TEXT)
    for j, line in enumerate(lines):
        d.text((x + 14 * S, y0 + (68 + j * 18) * S), line, font=font(REG, 11), fill=DIM)
    if i < len(STAGES) - 1:
        ax = x + bw + gap // 2
        ay = y0 + bh // 2
        d.line([(ax - 6 * S, ay), (ax + 5 * S, ay)], fill=(70, 130, 136), width=2 * S)
        d.polygon([(ax + 5 * S, ay - 4 * S), (ax + 11 * S, ay), (ax + 5 * S, ay + 4 * S)],
                  fill=(70, 130, 136))

d.text((62 * S, 366 * S),
       "AutoDock Vina · PLIP · RDKit · Open Babel · fpocket · ADMET-AI · OpenMM",
       font=font(REG, 12), fill=(96, 130, 136))
right = "Python 3.11 · Docker · GPL-3.0"
rw = d.textlength(right, font=font(REG, 12))
d.text((W - 60 * S - rw, 366 * S), right, font=font(REG, 12), fill=(96, 130, 136))

out = os.path.join(ROOT, "docs", "banner.png")
img.resize((W // S, H // S), Image.LANCZOS).save(out, optimize=True)
print("wrote", out, os.path.getsize(out) // 1024, "KB")
