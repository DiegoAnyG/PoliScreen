"""Reading a tunnel calculation someone else ran.

PoliScreen does not run CAVER or CaverDock. It reads what they left behind, which is the half that
costs nothing: no engine to install, no licence to accept, and it works the same in the container,
in a development checkout and in the one-click installer.

The fixtures here are the shape of a real cd-analysis output folder, trimmed to the few lines that
carry the numbers. The energies are the first discs of the reference job -- 8HTB, benzofuroxanic
acid, tunnel 3 -- so a format change shows up as a wrong energy and not just a missing file.
"""
import pytest

from poliscreen.core import tunnels

pytestmark = pytest.mark.skipif(not tunnels.available(),
                                reason="caver-translate is not installed in this environment")

# One MODEL per disc. The first CAVERDOCK TUNNEL remark of a MODEL is that model's disc; a run can
# append a second one repeating disc 0, which is why counting remarks reads the reference as data.
DISCS = [(0, -3.2, 1.5), (1, -2.7, 1.6), (2, -2.3, 1.7), (3, -0.5, 1.7), (4, -6.7, 2.7)]

# x y z, the disc normal, then the radius: a straight run one angstrom apart, widening inwards.
DSD = "\n".join(f"{i}.0 0.0 0.0 1.0 0.0 0.0 {1.5 + i * 0.1}" for i in range(5)) + "\n"


def _model(disc, energy, radius):
    return "\n".join([
        f"MODEL {disc + 1}",
        f"REMARK CAVERDOCK RESULT:      {energy}      0.000      0.000",
        f"REMARK CAVERDOCK TUNNEL: {disc}     {energy}      {radius}      1.3",
        "ENDMDL",
    ])


def _run(folder, direction="in"):
    """One finished CaverDock job, named the way caverdock-run names them."""
    folder = folder / f"r8HTB-lbenzo-ttun_cl_003-d{direction}-lowerbound"
    folder.mkdir(parents=True)
    (folder / "analysis-lb.pdbqt").write_text(
        "\n".join(_model(d, e, r) for d, e, r in DISCS) + "\n")
    (folder / "tunnel.dsd").write_text(DSD)
    return folder


def test_a_finished_run_becomes_a_table(tmp_path):
    _run(tmp_path)
    table, _cov = tunnels.read(tmp_path)
    assert len(table) == 1
    row = table.iloc[0]
    assert row["receptor"] == "8HTB"
    assert row["ligand"] == "benzo"
    assert row["tunnel"] == 3
    assert row["direction"] == "in"
    assert row["Ea"] == pytest.approx(2.7)
    assert row["dE_BS"] == pytest.approx(-3.5)


def test_the_missing_direction_is_counted_not_hidden(tmp_path):
    """A calculation that was never run is a result: the gap is the only evidence of the gap."""
    _run(tmp_path, "in")
    _table, cov = tunnels.read(tmp_path)
    assert cov["present"] == 1
    assert cov["expected"] == 2
    assert cov["missing"] == [("8HTB", "benzo", 3, "out")]


def test_an_empty_folder_reads_as_empty_not_as_an_error(tmp_path):
    """The Results tab asks before it knows there is anything to show."""
    table, cov = tunnels.read(tmp_path)
    assert table.empty
    assert cov["present"] == 0


def test_a_caverweb_download_reads_the_same_way(tmp_path):
    """The same folder argument, whether the job ran here or came back from the server."""
    import json
    import zipfile

    receptor = tmp_path / "8HTB"
    receptor.mkdir()
    profile = [{"distance": float(i), "disc": i, "radius": r, "energyLb": e,
                "energyUbMin": e, "energyUbMax": e} for i, e, r in DISCS]
    with zipfile.ZipFile(receptor / "benzo3in0a1b_results.zip", "w") as zf:
        zf.writestr("results.json", json.dumps(
            [{"name": "ligand", "hasUb": True, "profile": profile}]))

    table, _cov = tunnels.read(tmp_path)
    assert len(table) == 1
    assert table.iloc[0]["receptor"] == "8HTB"
    assert table.iloc[0]["tunnel"] == 3


def test_what_would_mislead_the_numbers_travels_with_the_row(tmp_path):
    """Ea without its caveats is how a two-angstrom dent in the surface wins a ranking."""
    _run(tmp_path)
    table, _cov = tunnels.read(tmp_path)
    assert "lower_bound_only" in table.iloc[0]["flags"]
    assert tunnels.FLAG_TEXT["lower_bound_only"]


def test_the_export_writes_the_page_and_the_tables(tmp_path):
    _run(tmp_path)
    out = tunnels.export(tmp_path, tmp_path / "report")
    assert (out / "transport.csv").exists()
    assert (out / "report.html").exists()
