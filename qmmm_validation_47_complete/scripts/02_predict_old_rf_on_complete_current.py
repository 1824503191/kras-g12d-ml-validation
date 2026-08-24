from pathlib import Path
import csv
import joblib
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances

WORK = Path("/data/student/mingsong/qmmm_validation_47_complete")
COMPLETE = Path("/home/edina/shared/qmmm_validation/complete_structures")

MODEL_PATH = WORK / "model_features/rf_model_dist_below_10_new.pkl"
FEATURE_LIST_PATH = WORK / "model_features/dist_features_below_10_ang_list_new.npy"

FEATURE_DIR = WORK / "features"
PRED_DIR = WORK / "predictions"
FEATURE_DIR.mkdir(exist_ok=True)
PRED_DIR.mkdir(exist_ok=True)

def compute_pp_dists_from_pdb(pdb_path, ligands=("GTP", "LIG"), radius=25.0, molid="Molecule"):
    u = mda.Universe(str(pdb_path))

    lig_expr = " or ".join(f"resname {rn}" for rn in ligands)
    shell = u.select_atoms(f"protein and around {radius} ({lig_expr})").residues

    res_exprs = []
    seen = set()

    for res in shell:
        seg = (res.segid or "").strip()
        expr = f"segid {seg} and resid {res.resid}" if seg else f"protein and resid {res.resid}"
        if expr not in seen:
            res_exprs.append(expr)
            seen.add(expr)

    lig_exprs = []
    for rn in ligands:
        sel = u.select_atoms(f"resname {rn}")
        if sel.n_atoms > 0:
            lig_exprs.append(f"resname {rn}")

    sel_exprs = res_exprs + lig_exprs
    if len(sel_exprs) < 2:
        raise ValueError(f"Not enough selections. Ligand selections found: {lig_exprs}")

    pairs = [(a, b) for i, a in enumerate(sel_exprs) for b in sel_exprs[i+1:]]

    def _label(expr):
        if expr.startswith("resname "):
            return expr.replace("resname ", "")
        if expr.startswith("segid "):
            parts = expr.split()
            seg = parts[1]
            resid = parts[-1]
            return f"PRO{seg}_{resid}"
        return "PRO" + expr.split()[-1]

    col_names = [f"{_label(a)}-{_label(b)}" for a, b in pairs]
    sel_objs = {expr: u.select_atoms(expr) for expr in sel_exprs}

    dims = getattr(u, "dimensions", None)
    box = None if dims is None or np.allclose(dims[:3], 0) else dims

    row = np.empty(len(pairs), dtype=np.float32)

    for j, (a, b) in enumerate(pairs):
        A, B = sel_objs[a], sel_objs[b]
        if A.n_atoms == 0 or B.n_atoms == 0:
            row[j] = np.nan
        else:
            row[j] = distances.distance_array(A.positions, B.positions, box=box).min()

    df = pd.DataFrame([np.round(row, 3)], columns=col_names)
    df.insert(0, "Molecule", molid)
    return df

def clean_protein_labels(df):
    def _clean_label(label):
        if not isinstance(label, str):
            return label
        parts = label.split("-")
        cleaned_parts = []
        for p in parts:
            if p.startswith("PRO"):
                cleaned_parts.append("PRO" + "".join(ch for ch in p if ch.isdigit()))
            else:
                cleaned_parts.append(p)
        return "-".join(cleaned_parts)

    df = df.copy()
    df.columns = [_clean_label(c) for c in df.columns]
    return df

def find_columns(rows):
    cols = list(rows[0].keys())
    id_col = "ligand_id" if "ligand_id" in cols else cols[0]

    if "barrier_kcal_per_mol" in cols:
        barrier_col = "barrier_kcal_per_mol"
    else:
        candidates = [c for c in cols if "barrier" in c.lower()]
        if not candidates:
            raise RuntimeError(f"Cannot find barrier column. Columns: {cols}")
        barrier_col = candidates[0]

    return id_col, barrier_col

def main():
    rows = list(csv.DictReader(open(WORK / "barriers.csv")))
    id_col, barrier_col = find_columns(rows)

    required_features = np.load(FEATURE_LIST_PATH, allow_pickle=True).astype(str).tolist()
    model = joblib.load(MODEL_PATH)

    print("Rows in barriers.csv:", len(rows))
    print("ID column:", id_col)
    print("Barrier column:", barrier_col)
    print("Required model features:", len(required_features))

    out_rows = []
    failures = []

    for i, r in enumerate(rows, start=1):
        molid = r[id_col]
        computed = float(r[barrier_col])
        pdb = COMPLETE / f"{molid}.pdb"

        print(f"\n[{i}/{len(rows)}] Processing {molid}")

        try:
            df = compute_pp_dists_from_pdb(
                pdb,
                ligands=("GTP", "LIG"),
                radius=25.0,
                molid=molid
            )
            df = clean_protein_labels(df)

            duplicate_cols = int(df.columns.duplicated().sum())
            missing = [c for c in required_features if c not in df.columns]

            print("Generated columns:", len(df.columns))
            print("Duplicate columns:", duplicate_cols)
            print("Missing required features:", len(missing))

            if duplicate_cols != 0:
                raise RuntimeError(f"Duplicate columns: {duplicate_cols}")

            if missing:
                raise RuntimeError(f"Missing required features: {len(missing)}; first missing: {missing[:20]}")

            X = df[required_features]
            n_na = int(X.isna().sum().sum())
            print("NA values in required features:", n_na)

            if n_na != 0:
                raise RuntimeError(f"NA values in required features: {n_na}")

            feature_out = FEATURE_DIR / f"{molid}_required_features.csv.gz"
            feature_df = pd.concat([df[["Molecule"]], X], axis=1)
            feature_df.to_csv(feature_out, index=False, compression="gzip")

            pred = float(model.predict(X)[0])

            out_rows.append({
                "ligand_id": molid,
                "computed_barrier_kcal_per_mol": computed,
                "old_rf_predicted_barrier": pred,
                "offset_pred_minus_computed": pred - computed,
                "source_batch": r.get("source_batch", ""),
                "n_sp_windows": r.get("n_sp_windows", ""),
                "converged_scfs": r.get("converged_scfs", "")
            })

            print(f"PASS: computed={computed:.4f}, predicted={pred:.4f}, offset={pred-computed:+.4f}")

        except Exception as e:
            failures.append((molid, str(e)))
            print("FAILED:", e)

    pred_df = pd.DataFrame(out_rows)
    pred_df.to_csv(PRED_DIR / "old_rf_vs_qmmm_complete_current.csv", index=False)

    if len(pred_df) > 0:
        pred_df[["ligand_id", "old_rf_predicted_barrier"]].to_csv(
            PRED_DIR / "old_rf_predictions_complete_current.tsv",
            sep="\t",
            index=False
        )

    with open(PRED_DIR / "failures_complete_current.tsv", "w") as f:
        f.write("ligand_id\terror\n")
        for molid, err in failures:
            f.write(f"{molid}\t{err}\n")

    print("\n=== SUMMARY ===")
    print("Successful predictions:", len(out_rows))
    print("Failed predictions:", len(failures))
    print("Saved:", PRED_DIR / "old_rf_vs_qmmm_complete_current.csv")
    print("Saved failures:", PRED_DIR / "failures_complete_current.tsv")

if __name__ == "__main__":
    main()
