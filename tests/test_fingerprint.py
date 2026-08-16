"""The diagnostic has to name the stage, not merely report that something differs.

Two machines on the same commit ranked the compounds differently, and every explanation offered
was a hypothesis nobody could check, because "the results differ" says nothing about where they
start to differ. Hashing each stage in pipeline order turns that into a diff whose first differing
line is the answer.
"""
from pathlib import Path

from poliscreen.core import fingerprint as fp
from poliscreen.core import layout as lay


def _project(root: Path) -> Path:
    proj = root / "run"
    (proj / lay.RECEPTORS).mkdir(parents=True)
    (proj / lay.RECEPTORS / "8HTB_ready.pdb").write_text("ATOM      1  N   HIS A  10\nEND\n")
    (proj / "prep").mkdir()
    (proj / "prep" / "8HTB_ready.pdbqt").write_text("ATOM      1  N   UNK A   1\n")
    (proj / "poses").mkdir()
    (proj / "poses" / "8HTB_ready__1-butanol.pdbqt").write_text("REMARK VINA RESULT: -8.1\n")
    return proj


def test_every_stage_of_the_pipeline_is_covered(tmp_path):
    stages = {s for s, _n, _h in fp.entries(_project(tmp_path))}
    assert stages == {"source-receptor", "docking-input", "pose"}


def test_the_same_project_hashes_the_same_twice(tmp_path):
    proj = _project(tmp_path)
    assert fp.entries(proj) == fp.entries(proj)


def test_a_changed_input_changes_only_its_own_stage(tmp_path):
    """This is the whole point: the first differing line says which stage to look at."""
    proj = _project(tmp_path)
    before = dict(((s, n), h) for s, n, h in fp.entries(proj))
    (proj / "prep" / "8HTB_ready.pdbqt").write_text("ATOM      1  N   UNK A   1\nATOM 2\n")
    after = dict(((s, n), h) for s, n, h in fp.entries(proj))
    changed = {k for k in before if before[k] != after[k]}
    assert changed == {("docking-input", "8HTB_ready.pdbqt")}


def test_a_legacy_spanish_project_is_still_read(tmp_path):
    """Older projects use the Spanish folder names, and are exactly the ones being compared."""
    proj = tmp_path / "old"
    (proj / lay.LEGACY[lay.RECEPTORS]).mkdir(parents=True)
    (proj / lay.LEGACY[lay.RECEPTORS] / "r.pdb").write_text("END\n")
    assert [s for s, _n, _h in fp.entries(proj)] == ["source-receptor"]


def test_the_report_is_ordered_so_a_plain_diff_works(tmp_path):
    """Two machines do not share a path, so the order must come from the content, not the disk."""
    rows = fp.entries(_project(tmp_path))
    stage_order = [s for s, _f, _p in fp.STAGES]
    seen = [stage_order.index(s) for s, _n, _h in rows]
    assert seen == sorted(seen), "stages are out of pipeline order"


def test_an_unrun_project_says_so_instead_of_looking_identical(tmp_path):
    """Two empty projects would otherwise match perfectly and prove nothing."""
    empty = tmp_path / "nothing"
    empty.mkdir()
    assert "nothing found" in fp.render(empty)


def test_a_windows_line_ending_is_not_a_divergence(tmp_path):
    """The same file written on Windows and on Linux must hash the same.

    It did not: CRLF made an identical prepared receptor differ at the first stage, which is
    exactly the reading that sends someone to audit a preparation step that in fact agreed.
    """
    lf, crlf = tmp_path / "lf.pdb", tmp_path / "crlf.pdb"
    lf.write_bytes(b"ATOM      1  N   HIS A  10\nEND\n")
    crlf.write_bytes(b"ATOM      1  N   HIS A  10\r\nEND\r\n")
    assert fp.sha256(lf) == fp.sha256(crlf)


def test_real_content_still_changes_the_hash(tmp_path):
    """Normalising line endings must not blunt the diagnostic this exists to serve."""
    a, b = tmp_path / "a.pdb", tmp_path / "b.pdb"
    a.write_bytes(b"ATOM      1  N   HIS A  10\n")
    b.write_bytes(b"ATOM      1  N   HIS A  11\n")
    assert fp.sha256(a) != fp.sha256(b)
