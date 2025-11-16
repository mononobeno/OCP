#!/usr/bin/env python

import os
import sys
import sqlite3
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

def log(msg: str):
    print(f"[INFO] {msg}", flush=True)

def warn(msg: str):
    print(f"[WARN] {msg}", flush=True)

def die(msg: str, code: int = 1):
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)
    sys.exit(code)

def main():
    if len(sys.argv) < 2:
        print("Usage: rdkit_gen3d_to_pdbqt.py INPUT.smi [LIB_CODE] [MAX_COUNT]", file=sys.stderr)
        sys.exit(1)

    script_dir = Path(__file__).resolve().parent
    root_dir = script_dir.parent
    db_path = root_dir / "catalog" / "db" / "ocp_results.sqlite"

    input_smi = Path(sys.argv[1]).resolve()
    lib_code  = sys.argv[2] if len(sys.argv) > 2 else "zinc_2d_smi_v1"
    max_count = int(sys.argv[3]) if len(sys.argv) > 3 else 100

    if not input_smi.is_file():
        die(f"INPUT.smi not found: {input_smi}")

    if not db_path.is_file():
        die(f"DB not found: {db_path}")

    lib_dir   = root_dir / "catalog" / "libraries" / lib_code
    rdkit_dir = lib_dir / "processed" / "rdkit_3d"
    pdbqt_dir = lib_dir / "processed" / "pdbqt"

    rdkit_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    log(f"ROOT_DIR  = {root_dir}")
    log(f"DB        = {db_path}")
    log(f"INPUT_SMI = {input_smi}")
    log(f"LIB_CODE  = {lib_code}")
    log(f"MAX_COUNT = {max_count}")
    log(f"RDKIT_DIR = {rdkit_dir}")
    log(f"PDBQT_DIR = {pdbqt_dir}")

    # DB 接続
    conn = sqlite3.connect(str(db_path))
    cur  = conn.cursor()

    # library_id
    cur.execute("SELECT id FROM libraries WHERE code = ?", (lib_code,))
    row = cur.fetchone()
    if not row:
        die(f"libraries.code='{lib_code}' not found in DB.")
    lib_id = row[0]
    log(f"LIB_ID = {lib_id}")

    # obabel チェック
    if not shutil_which("obabel"):
        die("obabel not found in PATH. Please install Open Babel.")

    # 実処理
    total = 0
    success = 0
    failed = 0

    with input_smi.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            smiles = parts[0]
            zinc_id = parts[1]

            total += 1
            if max_count != 0 and total > max_count:
                break

            zinc_id = zinc_id.strip()
            smiles  = smiles.strip()

            if not smiles or not zinc_id:
                continue

            pdb_file   = rdkit_dir / f"{zinc_id}.pdb"
            pdbqt_file = pdbqt_dir / f"{zinc_id}.pdbqt"

            # 既に pdbqt があればスキップ（DB メタ情報だけ補完）
            if pdbqt_file.is_file():
                log(f"[SKIP] exists: {pdbqt_file.name}")
                update_ligand_meta(cur, lib_id, zinc_id, smiles, input_smi.name, has_3d=1)
                continue

            log(f"[{total}] {zinc_id} -> 3D PDB + PDBQT")

            ok = gen3d_rdkit_to_pdb(smiles, pdb_file)
            if not ok:
                warn(f"RDKit 3D generation failed for {zinc_id}")
                failed += 1
                if pdb_file.is_file():
                    pdb_file.unlink()
                continue

            ok = pdb_to_pdbqt_with_obabel(pdb_file, pdbqt_file)
            if not ok:
                warn(f"obabel PDB->PDBQT failed for {zinc_id}")
                failed += 1
                if pdbqt_file.is_file():
                    pdbqt_file.unlink()
                continue

            # DB 更新
            update_ligand_meta(cur, lib_id, zinc_id, smiles, input_smi.name, has_3d=1)

            success += 1

    conn.commit()
    conn.close()

    log(f"processed ligands = {total}")
    log(f"success           = {success}")
    log(f"failed            = {failed}")

def shutil_which(cmd: str):
    from shutil import which
    return which(cmd)

def gen3d_rdkit_to_pdb(smiles: str, out_pdb: Path) -> bool:
    """RDKit で SMILES→3D を作り、PDB に書き出す。"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xf00d
        if AllChem.EmbedMolecule(mol, params) != 0:
            return False
        # 軽く最適化（回数は控えめ）
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        Chem.MolToPDBFile(mol, str(out_pdb))
        return True
    except Exception as e:
        warn(f"RDKit exception: {e}")
        return False

def pdb_to_pdbqt_with_obabel(in_pdb: Path, out_pdbqt: Path) -> bool:
    """PDB を PDBQT へ変換（ここでは --gen3d 不要）。"""
    try:
        cmd = ["obabel", "-ipdb", str(in_pdb), "-opdbqt", "-O", str(out_pdbqt)]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if res.returncode != 0:
            warn(f"obabel error: {res.stderr.decode(errors='ignore')}")
            return False
        return True
    except Exception as e:
        warn(f"obabel exception: {e}")
        return False

def update_ligand_meta(cur, lib_id: int, zinc_id: str, smiles: str,
                       source_file: str, has_3d: int = 0):
    cur.execute(
        """
        INSERT OR IGNORE INTO ligands (zinc_id, library_id)
        VALUES (?, ?)
        """,
        (zinc_id, lib_id),
    )
    cur.execute(
        """
        UPDATE ligands
           SET smiles           = COALESCE(smiles, ?),
               source_file      = COALESCE(source_file, ?),
               has_3d           = COALESCE(has_3d, ?),
               conformer_method = COALESCE(conformer_method, 'rdkit_etkdg')
         WHERE zinc_id = ? AND library_id = ?
        """,
        (smiles, source_file, has_3d, zinc_id, lib_id),
    )

if __name__ == "__main__":
    main()
