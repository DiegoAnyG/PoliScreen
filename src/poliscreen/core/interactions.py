"""Receptor-pose complexes and interactions with PLIP.

A single interaction engine. PLIP protonates internally and keeps the PDB author numbering, which is
the one that appears in the literature and in the viewers. The docked ligand is written as HETATM/LIG
so PLIP identifies it and does not confuse cofactors with the ligand.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from . import layout as lay
from typing import Optional, Sequence

import pandas as pd

from .screening import compound_from_pose_name, receptor_from_name

LIGAND_RESNAME = "LIG"
LIGAND_CHAIN = "Z"
LIGAND_RESSEQ = 900
WATERS = {"HOH", "WAT", "H2O", "DOD"}

PLIP_TAGS = {
    "hbond": ["hydrogen_bond"],
    "hydrophobic": ["hydrophobic_interaction"],
    "saltbridge": ["salt_bridge"],
    "water": ["water_bridge"],
    "halogen": ["halogen_bond"],
    "pistack": ["pi_stack"],
    "pication": ["pi_cation_interaction"],
    "metal": ["metal_complex"],
}


class InteractionsError(RuntimeError):
    pass


def receptor_text(pdb, cache_dir=None) -> str:
    """Receptor body without the closing lines, converting to PDB if needed."""
    pdb = Path(pdb)
    if pdb.suffix.lower() != ".pdb":
        cache = Path(cache_dir or pdb.parent)
        cache.mkdir(parents=True, exist_ok=True)
        conv = cache / f"{pdb.stem}_rec.pdb"
        if not conv.exists():
            subprocess.run(["obabel", str(pdb), "-O", str(conv)], capture_output=True)
        pdb = conv
    return "".join(l for l in pdb.read_text(errors="ignore").splitlines(keepends=True)
                   if not l.strip().startswith(("END", "ENDMDL")))


def max_serial(pdb_text: str) -> int:
    """Highest atom serial in a PDB body, so whatever is appended to it can start above."""
    top = 0
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                top = max(top, int(line[6:11]))
            except ValueError:
                continue
    return top


def ligand_text(pose_pdb, serial_offset: int = 0) -> str:
    """Rewrites the pose as HETATM/LIG in its own chain, with its bonds still meaning something.

    A pose numbers its atoms from 1 and so does the receptor, so fusing the two produced a file in
    which CONECT 1 named two different atoms. Those records are the ligand's real bonds -- a PDB
    has no other way to state them -- and a complex handed to PLIP without them leaves openbabel to
    infer the connectivity from distances, which costs hydrogen bonds that were known all along.
    Offsetting the serials makes the records unambiguous, so they can be kept.
    """
    text = Path(pose_pdb).read_text(errors="ignore")
    out, renumber = [], {}
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            old = int(line[6:11])
            name = line[12:16]
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        el = (line[76:78].strip() or name.strip()[:1])
        renumber[old] = old + serial_offset
        out.append(f"{'HETATM':<6}{renumber[old]:>5} {name:<4} {LIGAND_RESNAME:>3} "
                   f"{LIGAND_CHAIN}{LIGAND_RESSEQ:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
                   f"  1.00  0.00          {el:>2}\n")
    for line in text.splitlines():
        if not line.startswith("CONECT"):
            continue
        fields = [line[i:i + 5] for i in range(6, len(line), 5)]
        mapped = [renumber[int(f)] for f in fields if f.strip().isdigit() and int(f) in renumber]
        if len(mapped) > 1:
            out.append("CONECT" + "".join(f"{m:>5}" for m in mapped) + "\n")
    return "".join(out)


def fuse(receptor, pose_pdb, out_pdb, cache_dir=None) -> Path:
    out_pdb = Path(out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    body = receptor_text(receptor, cache_dir)
    if not body.endswith("\n"):
        body += "\n"
    out_pdb.write_text(body + ligand_text(pose_pdb, serial_offset=max_serial(body))\
                       + "END\n")
    return out_pdb


def fuse_batch(receptors: Sequence, poses_dir, out_dir, cache_dir=None, on_progress=None,
               stem_to_file=None) -> list:
    """Fuses each pose with its receptor. Skips the ones that already exist.

    The receptor is inferred from the pose name. `stem_to_file` (site_id -> file) lets several docking
    sites point to the same PDB (hybrid docking); if not passed, each receptor is its own site.
    """
    by_stem = dict(stem_to_file) if stem_to_file else {Path(r).stem: r for r in receptors}
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    models = sorted(Path(poses_dir).glob("*-model*.pdb"))
    made = []
    for i, pose in enumerate(models, 1):
        out = out_dir / f"{lay.COMPLEX_PREFIX}{pose.stem}.pdb"
        if not out.exists():
            rec = by_stem.get(receptor_from_name(pose.name))
            if rec is None:
                continue
            fuse(rec, pose, out, cache_dir)
            made.append(out)
        if on_progress:
            on_progress(i, len(models), out.name)
    return made


def sanitize_pdb(src, dst) -> bool:
    """Rewrites the complex with clean columns; PLIP is strict about the format."""
    clean = []
    for line in Path(src).read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                name = line[12:16]
                res = line[17:20].strip()
                chain = line[21] if len(line) > 21 else "A"
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                serial, resseq = line[6:11].strip(), line[22:26].strip()
                el = (line[76:78].strip() or name.strip()[0])
                clean.append(f"{line[0:6]:<6}{serial:>5} {name:<4} {res:>3} {chain}{resseq:>4}"
                             f"    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 20.00          {el:>2}\n")
            except Exception:
                continue
        elif line.startswith("TER"):
            clean.append("TER\n")
        elif line.startswith("CONECT"):
            # The ligand's bonds. Dropping them left openbabel re-inferring the connectivity
            # from distances inside PLIP, and hydrogen bonds went with it. The column rewriting
            # above is what this function is for -- a misaligned residue number is read as the
            # neighbouring residue -- and it never required throwing the bonds away.
            clean.append(line.rstrip() + "\n")
    if not clean:
        return False
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    Path(dst).write_text("".join(clean) + "END\n")
    return True


def parse_plip_xml(xml_path, ligand_resname: str = LIGAND_RESNAME) -> dict:
    """Returns {Residue+number_type: count} only of the docked ligand; ignores cofactors."""
    feats = {}
    try:
        root = ET.parse(xml_path).getroot()
    except Exception:
        return feats
    for site in root.iter("bindingsite"):
        hetid = (site.findtext(".//hetid") or "").strip().upper()
        if hetid and hetid != ligand_resname.upper():
            continue
        for itype, tags in PLIP_TAGS.items():
            for tag in tags:
                for it in site.iter(tag):
                    rn = (it.findtext("resnr") or "").strip()
                    rt = (it.findtext("restype") or "").strip()
                    if rn and rt:
                        k = f"{rt.capitalize()}{rn}_{itype}"
                        feats[k] = feats.get(k, 0) + 1
    return feats


def plip_available() -> bool:
    return shutil.which("plip") is not None


def run_plip(complex_pdb, xml_dir, san_dir) -> Optional[Path]:
    """Runs PLIP on a complex and leaves its XML. Reuses the one that already exists."""
    complex_pdb = Path(complex_pdb)
    xml_dir, san_dir = Path(xml_dir), Path(san_dir)
    xml_dir.mkdir(parents=True, exist_ok=True)
    san_dir.mkdir(parents=True, exist_ok=True)
    xml = xml_dir / f"{complex_pdb.stem}.xml"
    if xml.exists():
        return xml
    san = san_dir / f"fx_{complex_pdb.stem}.pdb"
    if not sanitize_pdb(complex_pdb, san):
        return None
    tmp = xml_dir / f"_tmp_{complex_pdb.stem}"
    tmp.mkdir(exist_ok=True)
    subprocess.run(["plip", "-f", str(san), "-x", "--nopdb", "-o", str(tmp)], capture_output=True, text=True)
    found = list(tmp.glob("*.xml"))
    if found:
        found[0].rename(xml)
    shutil.rmtree(tmp, ignore_errors=True)
    return xml if xml.exists() else None


def crystal_fingerprint(receptor, control_file, work_dir, cache_dir=None) -> dict:
    """PLIP fingerprint of the crystallographic ligand in its real pose (the one from the control
    file, not a docked pose). Fuses receptor + control and profiles with PLIP. It is the objective
    baseline: the compounds' interactions are compared against this, not against the control's docking."""
    work = Path(work_dir); work.mkdir(parents=True, exist_ok=True)
    cf = Path(control_file)
    if cf.suffix.lower() != ".pdb":
        pdbc = work / f"{cf.stem}_xtal.pdb"
        if not pdbc.exists():
            # Largest fragment: a fragmented control would put spurious remnants in the reference fingerprint.
            subprocess.run(["obabel", str(cf), "-O", str(pdbc), "-r"], capture_output=True)
        cf = pdbc
    if not cf.exists() or cf.stat().st_size == 0:
        return {}
    comp = work / f"{lay.COMPLEX_PREFIX}xtal_{Path(receptor).stem}.pdb"
    fuse(receptor, cf, comp, cache_dir)
    xml = run_plip(comp, work / "xml", work / "san")
    return parse_plip_xml(xml) if xml else {}


def hetero_fingerprint(receptor_pdb, box, work_dir, min_atoms: int = 6) -> tuple:
    """PLIP fingerprint of a co-crystallized cofactor/ligand that falls INSIDE a site's box. Serves as
    an objective reference for a secondary pocket of the hybrid docking. Returns (feats, label) or
    ({}, None) if there is none of sufficient size (1-atom ions do not count).

    PLIP is run on the receptor as is (which already contains the cofactor) and its binding site is
    read, instead of re-attaching it, so the cofactor is not counted twice."""
    receptor_pdb = Path(receptor_pdb)
    hx, hy, hz = box.sx / 2.0, box.sy / 2.0, box.sz / 2.0
    groups_ = {}
    for l in receptor_pdb.read_text(errors="ignore").splitlines():
        if not l.startswith("HETATM"):
            continue
        rn = l[17:20].strip()
        if rn in WATERS:
            continue
        try:
            x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        except ValueError:
            continue
        if abs(x - box.cx) <= hx and abs(y - box.cy) <= hy and abs(z - box.cz) <= hz:
            groups_.setdefault((rn, (l[21].strip() or "_"), l[22:26].strip()), []).append(l)
    if not groups_:
        return {}, None
    (rn, ch, seq), atoms = max(groups_.items(), key=lambda kv: len(kv[1]))
    if len(atoms) < min_atoms:
        return {}, None
    work = Path(work_dir); work.mkdir(parents=True, exist_ok=True)
    xml = run_plip(receptor_pdb, work / "xml", work / "san")
    if not xml:
        return {}, None
    return parse_plip_xml(xml, ligand_resname=rn), f"{rn} {ch}:{seq}"


def crystal_fingerprints(receptors: Sequence, controls: Sequence, control_assign: dict,
                         work_dir, cache_dir=None) -> dict:
    """{receptor_stem: feats} taking, for each receptor, the control assigned to it.
    With a single control it is used for the only receptor. Serves as the objective ranking reference."""
    from .screening import normalize_key
    by_ck = {normalize_key(Path(c).stem): c for c in controls}
    out = {}
    for r in receptors:
        stem = Path(r).stem
        ck = next((k for k, rc in control_assign.items() if rc == stem), None)
        cfile = by_ck.get(ck) if ck else None
        if cfile is None and len(controls) == 1:
            cfile = controls[0]
        if cfile:
            fp = crystal_fingerprint(r, cfile, work_dir, cache_dir)
            if fp:
                out[stem] = fp
    return out


def plip_batch(complexes: Sequence, work_dir, workers: int = 0, force: bool = False,
               on_progress=None) -> pd.DataFrame:
    """Profiles all complexes and returns the pose x interaction matrix.

    Reuses interacciones.csv unless force=True. Writes the CSV in work_dir.
    """
    work = Path(work_dir)
    csv = lay.artifact(work, lay.INTERACTIONS_CSV)
    xml_dir, san_dir = work / "xml_plip", work / "san"
    expected_items = {Path(c).stem for c in complexes}
    # The cache is only reused if it covers every current complex; otherwise the ranking is falsified.
    if csv.exists() and not force:
        previous = pd.read_csv(csv)
        if expected_items and expected_items.issubset(set(previous.get("name", []))):
            return previous
    if force:
        for d in (xml_dir, san_dir):
            shutil.rmtree(d, ignore_errors=True)
    if not plip_available():
        raise InteractionsError(
            "PLIP is not installed. Install it without compiling OpenBabel: "
            "conda install -c conda-forge -c bioconda plip"
        )
    complexes = [Path(c) for c in complexes]
    if not complexes:
        raise InteractionsError("No complexes to analyze. Run the fusion first.")
    n = workers if workers > 0 else max(1, (os.cpu_count() or 2) - 2)
    with ThreadPoolExecutor(max_workers=n) as ex:
        for i, _ in enumerate(ex.map(lambda c: run_plip(c, xml_dir, san_dir), complexes), 1):
            if on_progress:
                on_progress(i, len(complexes))
    rows = []
    for c in complexes:
        x = xml_dir / f"{c.stem}.xml"
        if not x.exists():
            continue
        row = parse_plip_xml(x)
        row["name"] = x.stem
        rows.append(row)
    if not rows:
        raise InteractionsError("PLIP produced no results. Check that the complexes have the ligand as LIG.")
    df = pd.DataFrame(rows).fillna(0)
    cols = ["name"] + sorted(c for c in df.columns if c != "name")
    df = df[cols]
    for c in cols[1:]:
        df[c] = df[c].astype(int)
    work.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv, index=False)
    return df
