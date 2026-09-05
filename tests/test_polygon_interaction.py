"""Tests for polygon_interaction module."""
import pandas as pd
import pytest
import matplotlib.pyplot as plt

from poliscreen.core.polygon_interaction import (
    draw_interaction_polygon,
    draw_interaction_polygon_bytes,
    extract_pocket_residues,
    generate_top_interactions_report,
    parse_contacts_from_features,
    parse_contacts_from_row,
    residue_sort_key,
)


def test_residue_sort_key():
    assert residue_sort_key("Val7") < residue_sort_key("Ala8")
    assert residue_sort_key("Ala8") < residue_sort_key("Lys78")
    assert residue_sort_key("Lys78") < residue_sort_key("Asp150")
    # Non-numbered fallback
    assert residue_sort_key("LIG")[0] == 9999


def test_extract_pocket_residues():
    df = pd.DataFrame([
        {"receptor": "3ACX", "Val7_hbond": 1, "Lys78_saltbridge": 1, "Asp150_hydrophobic": 0},
        {"receptor": "3ACX", "Val7_hbond": 0, "Lys78_saltbridge": 0, "Asp150_hydrophobic": 2},
        {"receptor": "OtherRec", "His200_hbond": 1},
    ])
    res = extract_pocket_residues(df, "3ACX", extra_residues=["Met152"])
    assert res == ["Val7", "Lys78", "Asp150", "Met152"]


def test_parse_contacts_from_row():
    row = pd.Series({
        "name": "pose_1",
        "Lys78_saltbridge": 1,
        "Asp150_hbond": 2,
        "Asp150_hydrophobic": 0,
        "unrelated_col": 99,
    })
    contacts = parse_contacts_from_row(row)
    assert "Lys78" in contacts
    assert contacts["Lys78"] == [("saltbridge", 1)]
    assert "Asp150" in contacts
    assert contacts["Asp150"] == [("hbond", 2)]
    assert "unrelated_col" not in contacts


def test_parse_contacts_from_features():
    feats = ["Lys78_saltbridge", "Asp150_hbond", "Trp92_pistack"]
    contacts = parse_contacts_from_features(feats)
    assert contacts["Lys78"] == [("saltbridge", 1)]
    assert contacts["Asp150"] == [("hbond", 1)]
    assert contacts["Trp92"] == [("pistack", 1)]


def test_draw_interaction_polygon():
    pocket = ["Val7", "Ala8", "Lys78", "Trp92", "Tyr110", "Leu127", "Asp150"]
    cmp_contacts = {
        "Lys78": [("saltbridge", 1)],
        "Asp150": [("hbond", 2)],
        "Trp92": [("pistack", 1)],
        "Leu127": [("hydrophobic", 1)],
    }
    ctrl_contacts = {
        "Lys78": [("saltbridge", 1)],
        "Asp150": [("hbond", 1)],
        "Tyr110": [("hbond", 1)],
    }
    fig = draw_interaction_polygon(
        all_pocket_residues=pocket,
        compound_contacts=cmp_contacts,
        control_contacts=ctrl_contacts,
        catalytic_residues=["Asp150"],
        secondary_residues=["Trp92"],
        title="Candidate 1",
    )
    assert fig is not None
    plt.close(fig)


def test_draw_interaction_polygon_bytes():
    pocket = ["Val7", "Ala8", "Lys78", "Asp150"]
    cmp_contacts = {"Lys78": [("saltbridge", 1)], "Asp150": [("hbond", 1)]}
    data = draw_interaction_polygon_bytes(
        all_pocket_residues=pocket,
        compound_contacts=cmp_contacts,
        title="Test Bytes",
    )
    assert isinstance(data, bytes)
    assert len(data) > 1000


def test_generate_top_interactions_report():
    pocket = ["Val7", "Lys78", "Asp150"]
    compounds = [
        {
            "rank": 1,
            "name": "Bf-1",
            "eff": 85.0,
            "dock": -8.1,
            "quality": 0.8,
            "contacts": {"Lys78": [("saltbridge", 1)], "Asp150": [("hbond", 1)]},
        },
        {
            "rank": 2,
            "name": "Bf-2",
            "eff": 78.0,
            "dock": -7.5,
            "quality": 0.7,
            "contacts": {"Val7": [("hydrophobic", 1)]},
        },
    ]
    pdf_bytes = generate_top_interactions_report(
        compounds_data=compounds,
        all_pocket_residues=pocket,
        control_contacts={"Lys78": [("saltbridge", 1)]},
        catalytic_residues=["Lys78"],
        secondary_residues=[],
        target_name="3ACX",
        format="pdf",
    )
    assert isinstance(pdf_bytes, bytes)
    assert pdf_bytes.startswith(b"%PDF")


def test_get_compound_pocket_residues():
    from poliscreen.core.polygon_interaction import get_compound_pocket_residues

    cmp_contacts = {"Leu91": [("hydrophobic", 1)], "Ile377": [("hbond", 1)]}
    ctrl_contacts = {"Tyr126": [("pistack", 1)], "Phe230": [("hydrophobic", 1)]}
    cat = ["Tyr126"]
    sec = ["Leu125"]

    res = get_compound_pocket_residues(
        compound_contacts=cmp_contacts,
        control_contacts=ctrl_contacts,
        catalytic_residues=cat,
        secondary_residues=sec,
    )
    # Total residues is strictly the union of these key anchors and contacts
    assert res == ["Leu91", "Leu125", "Tyr126", "Phe230", "Ile377"]
    assert len(res) == 5


