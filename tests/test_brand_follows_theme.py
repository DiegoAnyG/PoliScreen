"""The wordmark must change colour with the theme, without waiting for a rerun.

Streamlit does not re-run the script when the theme is switched, so a colour chosen in Python and
written into the tag survives the switch: the wordmark stayed in the previous theme's colour until
something unrelated forced a rerun -- changing the language was what usually did it, which is how
this was reported twice. Painting it as a mask filled with currentColor moves the decision to the
browser, where every repaint resolves it.
"""
import re
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src" / "poliscreen"
APP = (SRC / "ui" / "streamlit_app.py").read_text(encoding="utf-8")


def test_the_wordmark_takes_its_colour_from_the_theme_not_from_python():
    start = APP.index("def _img_tinted")
    block = APP[start:APP.index("</span>\")", start)]
    assert "currentColor" in block
    assert "_tema" not in block, "a Python-side theme decision is stale after a theme switch"


def test_the_wordmark_is_drawn_with_that_function():
    assert re.search(r"_marca\.append\(_img_tinted\(_wm, \d+\)\)", APP)


def test_the_size_comes_from_the_file_so_a_replacement_still_fits():
    """The wordmark is meant to be drop-in: the mask needs a width, and guessing one would crop it."""
    from importlib import util
    spec = util.spec_from_file_location("_wm", SRC / "core" / "viewer.py")
    mod = util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    wm = mod.wordmark_path()
    if wm is None or wm.suffix.lower() != ".png":
        return
    head = wm.read_bytes()[16:24]
    w, h = int.from_bytes(head[:4], "big"), int.from_bytes(head[4:], "big")
    assert w > 0 and h > 0, "PNG header did not parse; the mask would be laid out with no width"
