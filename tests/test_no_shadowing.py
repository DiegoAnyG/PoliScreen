"""No local may shadow the interface helpers.

A local named `t` turns every t("...") in that function into a call on the local:
the 3D complex viewer died with "'DataFrame' object is not callable". It is invisible
until that exact branch runs, so it is checked statically instead.
"""
import ast
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "ui" / "streamlit_app.py"


def _imported_helpers(tree):
    """Names the module imports and then calls, e.g. t, lay, sc, nm."""
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                names.add(alias.asname or alias.name.split(".")[0])
    return names


def test_no_local_shadows_an_imported_helper():
    tree = ast.parse(APP.read_text(encoding="utf-8"))
    helpers = _imported_helpers(tree)
    offenders = []
    for scope in ast.walk(tree):
        if not isinstance(scope, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        # Comprehension variables live in their own scope and cannot shadow the outer name,
        # so those subtrees are skipped whole (ast.walk does not prune by itself).
        inner = set()
        for n in ast.walk(scope):
            if isinstance(n, (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp)):
                inner.update(id(c) for c in ast.walk(n))
        assigned = set()
        for n in ast.walk(scope):
            if id(n) in inner:
                continue
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store):
                assigned.add(n.id)
            elif isinstance(n, ast.arg):
                assigned.add(n.arg)
        clash = assigned & helpers
        if clash:
            offenders.append(f"{scope.name}: {sorted(clash)}")
    assert not offenders, f"locals shadowing imported helpers: {offenders}"
