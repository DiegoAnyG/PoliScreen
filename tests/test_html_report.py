"""Tests for the Zero-Server Standalone Interactive HTML Report generator."""
from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from poliscreen.core import html_report as hr
from poliscreen.core import layout as lay
from poliscreen.core import session as ss


@pytest.fixture
def mock_project(tmp_path: Path) -> Path:
    """Create a populated PoliScreen mock project for report testing."""
    proj = tmp_path / "mock_project_01"
    proj.mkdir(parents=True, exist_ok=True)

    # 1. run.json
    meta = {
        "project": "MockScreen_Test",
        "receptors": ["target_ready.pdb"],
        "controls": ["ctrl_lig.sdf"],
        "control_assign": {"ctrllig": "target_ready"},
        "catalytic": {"target_ready": ["HIS57", "ASP102"]},
        "secondary": {"target_ready": ["SER195"]},
        "crystal_feats": {"target_ready": ["HIS57_hbond", "TRP215_hydrophobic"]},
        "exhaustiveness": 8,
        "n_poses": 9,
        "energy_range": 3.0,
        "ph": 7.4,
        "seed": 42,
    }
    (proj / "run.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    # 2. ranking.csv
    rk_df = pd.DataFrame([
        {
            "compound": "ctrllig",
            "is_control": 1,
            "best_dock": -8.5,
            "best_inter": 1.0,
            "effectiveness_pct": 100.0,
            "confidence": 0.85,
            "pKi": 6.2,
            "LE": 0.32,
            "pareto_rank": 1,
            "is_pareto": True,
            "receptor": "target_ready",
            "pose": 1,
        },
        {
            "compound": "cand_alpha",
            "is_control": 0,
            "best_dock": -9.2,
            "best_inter": 0.95,
            "effectiveness_pct": 115.0,
            "confidence": 0.78,
            "pKi": 6.8,
            "LE": 0.36,
            "pareto_rank": 1,
            "is_pareto": True,
            "receptor": "target_ready",
            "pose": 1,
        },
        {
            "compound": "cand_beta",
            "is_control": 0,
            "best_dock": -7.1,
            "best_inter": 0.40,
            "effectiveness_pct": 75.0,
            "confidence": 0.45,
            "pKi": 5.2,
            "LE": 0.25,
            "pareto_rank": 2,
            "is_pareto": False,
            "receptor": "target_ready",
            "pose": 1,
        },
    ])
    rk_df.to_csv(proj / "ranking.csv", index=False)

    # 3. interactions.csv
    inter_df = pd.DataFrame([
        {
            "name": "Complex_docking_target_ready_compounds_a_ctrllig-model1",
            "compound": "ctrllig",
            "HIS57_hbond": 1,
            "TRP215_hydrophobic": 1,
        },
        {
            "name": "Complex_docking_target_ready_compounds_a_cand_alpha-model1",
            "compound": "cand_alpha",
            "HIS57_hbond": 1,
            "ASP102_saltbridge": 1,
            "TRP215_hydrophobic": 1,
        },
        {
            "name": "Complex_docking_target_ready_compounds_a_cand_beta-model1",
            "compound": "cand_beta",
            "GLY193_hbond": 1,
        },
    ])
    inter_df.to_csv(proj / "interactions.csv", index=False)

    # 4. ligands_meta.csv
    lm_df = pd.DataFrame([
        {"name": "ctrllig", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
        {"name": "cand_alpha", "smiles": "c1ccccc1NC(=O)C"},
        {"name": "cand_beta", "smiles": "c1ccncc1"},
    ])
    lm_df.to_csv(proj / "ligands_meta.csv", index=False)

    # 5. Receptors PDB
    rec_dir = proj / "receptors"
    rec_dir.mkdir()
    sample_pdb = (
        "ATOM      1  N   HIS A  57      20.154  13.442  12.331  1.00 20.00           N\n"
        "ATOM      2  CA  HIS A  57      21.234  14.221  13.001  1.00 20.00           C\n"
        "ATOM      3  N   ASP A 102      18.112  11.230  10.111  1.00 20.00           N\n"
        "END\n"
    )
    (rec_dir / "target_ready.pdb").write_text(sample_pdb, encoding="utf-8")

    # 6. Poses PDB
    poses_dir = proj / "poses"
    poses_dir.mkdir()
    lig_pdb = (
        "HETATM    1  C1  LIG A   1      20.500  13.800  12.500  1.00 20.00           C\n"
        "HETATM    2  O1  LIG A   1      21.100  14.100  12.800  1.00 20.00           O\n"
        "END\n"
    )
    (poses_dir / "docking_target_ready_compounds_a_ctrllig-model1.pdb").write_text(lig_pdb, encoding="utf-8")
    (poses_dir / "docking_target_ready_compounds_a_cand_alpha-model1.pdb").write_text(lig_pdb, encoding="utf-8")

    return proj


def test_build_interactive_report_full(mock_project: Path, tmp_path: Path):
    """Test generating a full standalone interactive HTML dossier from a project."""
    out_html = tmp_path / "test_dossier.html"
    html_str = hr.build_interactive_report(mock_project, out_path=out_html)

    assert len(html_str) > 1000
    assert out_html.is_file()
    assert out_html.stat().st_size > 1000

    # Key HTML components present
    assert "<!DOCTYPE html>" in html_str
    assert "MockScreen_Test" in html_str
    assert "PoliScreen Dossier" in html_str
    assert "cand_alpha" in html_str
    assert "ctrllig" in html_str
    assert "Pareto Frontier" in html_str
    assert "PLIP" in html_str
    assert "3Dmol" in html_str
    assert "init3DViewer" in html_str
    assert "initPlotlyChart" in html_str

    # Embedded PDB data
    assert "ATOM      1  N   HIS" in html_str
    assert "HETATM    1  C1  LIG" in html_str

    # Methods section
    assert "Methodology" in html_str or "Reproducibility" in html_str
    assert "exhaustiveness" in html_str

    # Spanish language test
    html_es = hr.build_interactive_report(mock_project, lang="es")
    assert "Frontera de Pareto" in html_es


def test_build_interactive_report_minimal(tmp_path: Path):
    """Test generating report gracefully when minimal data exists."""
    min_proj = tmp_path / "min_proj"
    min_proj.mkdir()
    rk_df = pd.DataFrame([
        {"compound": "aspirin", "best_dock": -7.2, "best_inter": 0.8, "effectiveness_pct": 95.0}
    ])
    rk_df.to_csv(min_proj / "ranking.csv", index=False)

    html_str = hr.build_interactive_report(min_proj)
    assert "<!DOCTYPE html>" in html_str
    assert "aspirin" in html_str
    assert "min_proj" in html_str


def test_session_catalog_and_package_bytes(mock_project: Path):
    """Test that html_report is in catalog and properly bundled in session packages."""
    cat = ss.catalog(mock_project)
    assert "html_report" in cat
    assert cat["html_report"]["has"] is True
    assert "Single-file" in cat["html_report"]["reason"] or "Self-contained" in cat["html_report"]["reason"]

    # Test package_bytes
    data, items = ss.package_bytes(mock_project, ["html_report", "results_csv"])
    assert "reporte_interactivo.html" in items
    assert len(data) > 2000

    # Verify zip content
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        names = zf.namelist()
        assert "reporte_interactivo.html" in names
        assert "CONTENTS.txt" in names
        report_content = zf.read("reporte_interactivo.html").decode("utf-8")
        assert "<!DOCTYPE html>" in report_content
        assert "MockScreen_Test" in report_content


def test_multi_pocket_isolation(tmp_path: Path):
    """Ensure that multi-pocket targets (e.g. rec~principal vs rec~Pk16) isolate only the selected pocket."""
    proj = tmp_path / "multi_pocket_proj"
    proj.mkdir()
    meta = {
        "project": "MultiPocketTest",
        "receptors": ["8HTB_ready.pdb"],
        "controls": ["control_zi9.sdf"],
        "control_assign": {"control_zi9": "8HTB_ready"},
    }
    (proj / "run.json").write_text(json.dumps(meta), encoding="utf-8")

    # 1 control + 2 candidates in principal pocket, and 1 control + 2 candidates in Pk16 pocket
    rk_df = pd.DataFrame([
        {"compound": "control_zi9", "receptor": "8HTB_ready~principal", "is_control": 1, "best_dock": -10.0, "best_inter": 1.0, "LE": 0.22, "effectiveness_pct": 100.0, "pose": 1},
        {"compound": "benzofuroxan", "receptor": "8HTB_ready~principal", "is_control": 0, "best_dock": -7.0, "best_inter": 0.66, "LE": 0.44, "effectiveness_pct": 66.0, "pose": 1},
        {"compound": "cand_x", "receptor": "8HTB_ready~principal", "is_control": 0, "best_dock": -6.5, "best_inter": 0.50, "LE": 0.30, "effectiveness_pct": 55.0, "pose": 1},
        {"compound": "control_zi9", "receptor": "8HTB_ready~Pk16", "is_control": 1, "best_dock": -9.0, "best_inter": 0.90, "LE": 0.20, "effectiveness_pct": 100.0, "pose": 1},
        {"compound": "benzofuroxan", "receptor": "8HTB_ready~Pk16", "is_control": 0, "best_dock": -6.0, "best_inter": 0.40, "LE": 0.38, "effectiveness_pct": 50.0, "pose": 1},
        {"compound": "cand_y", "receptor": "8HTB_ready~Pk16", "is_control": 0, "best_dock": -5.5, "best_inter": 0.35, "LE": 0.25, "effectiveness_pct": 45.0, "pose": 1},
    ])
    rk_df.to_csv(proj / "ranking.csv", index=False)

    # 1. Single target isolation mode
    html_isolated = hr.build_interactive_report(proj, target="8HTB_ready~principal", all_targets=False)
    assert "cand_x" in html_isolated
    assert "cand_y" not in html_isolated
    assert '<select id="target-selector"' not in html_isolated

    # 2. Multi-target dashboard with interactive pocket switcher
    html_dashboard = hr.build_interactive_report(proj, target="8HTB_ready~principal", all_targets=True)
    assert '<select id="target-selector"' in html_dashboard
    assert "TARGETS_DATA" in html_dashboard
    assert "cand_x" in html_dashboard
    assert "cand_y" in html_dashboard

    # Verify LE effectiveness is rendered
    assert "Eff (LE %)" in html_dashboard or "Effectiveness (LE)" in html_dashboard or "th_eff_le" in html_dashboard


def test_fpocket_cavity_extraction(tmp_path: Path):
    """Test extracting fpocket alpha spheres and cavity-lining residues for HTML report."""
    proj = tmp_path / "fpocket_proj"
    proj.mkdir()

    # Create a mock .poliscreen session archive with fpocket data in ui_state.json
    zf_path = proj / "mock_session.poliscreen"
    with zipfile.ZipFile(zf_path, "w") as zf:
        ui_state = {
            "pockets": {
                "{project}/receptors/test_rec.pdb": [
                    {
                        "n": 1,
                        "druggability": 0.88,
                        "volume": 950.0,
                        "spheres": 3,
                        "alpha_xyz": [[10.0, 20.0, 30.0, 3.5], [11.0, 21.0, 31.0, 3.8], [12.0, 22.0, 32.0, 4.0]],
                        "residues": ["His57", "Asp102", "Ser195"],
                    },
                    {
                        "n": 16,
                        "druggability": 0.12,
                        "volume": 500.0,
                        "spheres": 2,
                        "alpha_xyz": [[50.0, 60.0, 70.0, 3.2], [51.0, 61.0, 71.0, 3.4]],
                        "residues": ["Ala10", "Val12"],
                    }
                ]
            }
        }
        zf.writestr("ui_state.json", json.dumps(ui_state))

    # Test extraction for principal (matched via catalytic His57, Asp102)
    sph1, nums1, info1 = hr._extract_fpocket_data_for_target(
        proj, "test_rec~principal", "", [57, 102]
    )
    assert len(sph1) == 3
    assert sph1[0]["x"] == 10.0
    assert 57 in nums1
    assert 102 in nums1
    assert info1["n"] == 1
    assert info1["druggability"] == 0.88

    # Test extraction for Pk16 (matched via pocket index 16)
    sph16, nums16, info16 = hr._extract_fpocket_data_for_target(
        proj, "test_rec~Pk16", "", []
    )
    assert len(sph16) == 2
    assert sph16[0]["x"] == 50.0
    assert 10 in nums16
    assert 12 in nums16
    assert info16["n"] == 16


def test_pymol_tunnel_zip_compilation(tmp_path: Path):
    """Test generating a self-contained PyMOL visualizer bundle ZIP for Caver tunnels."""
    from poliscreen.ui.views.results import _build_pymol_tunnel_zip

    proj_dir = tmp_path / "tunnels_proj"
    tunnels_dir = proj_dir / "tunnels"
    tunnels_dir.mkdir(parents=True)

    # Mock caver directory and a caverdock run directory
    caver_out = tunnels_dir / "caver_test_rec" / "out"
    caver_out.mkdir(parents=True)
    (caver_out / "summary.txt").write_text("ID Avg_BR ...\n1 2.5 ...\n")
    (caver_out / "tun_cl_001.pdb").write_text("ATOM      1  N   TUN A   1       0.0   0.0   0.0  1.00  1.50\nEND\n")

    run_dir = tunnels_dir / "rTestRec-lBenzofuroxan-ttun_cl_001-din-upperbound"
    run_dir.mkdir()
    (run_dir / "analysis-lb.pdbqt").write_text("MODEL 1\nHETATM    1  C1  LIG A   1       1.0   1.0   1.0\nENDMDL\n")
    (run_dir / "figure.pml").write_text("# PyMOL figure script\nload analysis-lb.pdbqt\n")

    zip_name, zip_bytes = _build_pymol_tunnel_zip(proj_dir, proj_dir)

    assert zip_bytes is not None
    assert len(zip_bytes) > 500

    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        names = zf.namelist()
        assert "README.txt" in names
        readme = zf.read("README.txt").decode("utf-8")
        assert "Visualization Package" in readme
        # Check that the visualizer .pml script is present
        assert any(n.endswith(".pml") for n in names)


