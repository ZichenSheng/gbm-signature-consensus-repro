"""End-to-end demonstration on SYNTHETIC data (examples/).

Runs every released procedure once: RNA-seq processing, membership and weighted scores, reference
coordinates, the frozen-loading consensus, pairwise-structure reproducibility, ecological share,
leave-one-cohort-out recovery, driver-associated variance with a cluster bootstrap, fixed-effect
pooling, ICC(1) contrasts, and the expression-matched random collection (N1) regenerated from the
supplied pool features.

The synthetic numbers demonstrate execution only. They are NOT evidence for, and must not be compared
with, any result in the manuscript.

Usage:  python scripts/run_example.py [--out example_output] [--n-boot 200]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import signature_analysis as sa  # noqa: E402
import statistical_validation as sv  # noqa: E402

DATA, EX = ROOT / "data", ROOT / "examples"


def load_package_data():
    ids, sigs = sa.load_signatures(DATA / "signatures_201.tsv")
    ref = json.load(open(DATA / "reference_gene_sets.json"))
    return {
        "ids": ids, "sigs": sigs, "ref": ref,
        "formulas": sa.load_formulas(DATA / "weighted_formulas_104.tsv"),
        "pcg": [g.strip() for g in open(DATA / "protein_coding_genes.txt") if g.strip()],
        "loadings": pd.read_csv(DATA / "consensus_loadings.tsv", sep="\t").set_index("signature_id").loading,
        "sd_ref": pd.read_csv(DATA / "formula_gene_reference_sd.tsv", sep="\t").set_index("gene_symbol").sd_ref,
        "pool_features": pd.read_csv(DATA / "matching_pool_features.tsv", sep="\t").set_index("gene_symbol"),
    }


def score_cohort(A: pd.DataFrame, P: dict):
    """Linear-scale RNA-seq (genes x samples) -> Z, coordinates (+ DEV_rel, AXIS), membership and weighted scores."""
    X = sa.normalize_rnaseq(A, P["pcg"])
    Z = sa.zscore_genes(X)
    C, _, _ = sa.score_reference_coordinates(Z, P["ref"]["coordinates"], X.mean(axis=1))
    C["DEV_rel"] = sa.relative_programme_scores(Z, {"DEV": P["ref"]["programme_sets"]["DEV"]}, P["ref"]["pan_malignant_reference"])["DEV_rel"]
    C["AXIS"], used = sa.consensus_score(Z, P["sigs"], P["loadings"])
    S, cov = sa.membership_scores(Z, P["sigs"], P["ids"])
    W = sa.weighted_formula_scores(Z, P["formulas"], P["sd_ref"])
    return Z, C, S, W, used


def main(out: Path, n_boot: int = 200):
    P = load_package_data()
    expr = pd.read_csv(EX / "synthetic_example.csv", index_col=0)
    samples = pd.read_csv(EX / "synthetic_samples.csv").set_index("sample_id")
    expr = sa.harmonize_symbols(expr, how="sum")
    out.mkdir(parents=True, exist_ok=True)
    summary = {"WARNING": "SYNTHETIC DATA - demonstrates execution only; not evidence for any manuscript result"}

    R = {}
    for coh, ss in samples.groupby("cohort"):
        Z, C, S, W, used = score_cohort(expr[ss.index], P)
        R[coh] = dict(Z=Z, C=C, S=S, W=W)
        C.to_csv(out / f"{coh}_coordinates_and_consensus.tsv", sep="\t")
        S.to_csv(out / f"{coh}_membership_scores.tsv", sep="\t")
        summary[coh] = {"samples": int(Z.shape[1]), "genes_after_processing": int(Z.shape[0]),
                        "membership_scores_computable": int(S.notna().all().sum()), "consensus_signatures_used": len(used),
                        "weighted_formulas_computable": int(W.shape[1]),
                        "ecological_share_of_consensus": round(sv.ecological_share(C["AXIS"], C), 3)}

    # pairwise correlation structure between two cohorts (signatures scorable in both)
    keep = [r for r in P["ids"] if R["SYN_A"]["S"][r].notna().all() and R["SYN_B"]["S"][r].notna().all()]
    summary["structure_reproducibility_A_vs_B"] = round(sv.structure_reproducibility(R["SYN_A"]["S"][keep], R["SYN_B"]["S"][keep]), 3)

    # leave-one-cohort-out recovery of the AC-developmental balance
    common = [r for r in P["ids"] if all(R[c]["S"][r].notna().all() for c in R)]
    X = {c: R[c]["S"][common].to_numpy() for c in R}
    y = {c: R[c]["C"]["MAL_AC"].to_numpy() for c in R}
    summary["heldout_R2_balance_all_scores"] = sv.loco_ridge_recovery(X, y).round(3).to_dict(orient="records")
    summary["heldout_R2_balance_consensus_alone"] = sv.loco_one_variable({c: R[c]["C"]["AXIS"].to_numpy() for c in R}, y).round(3).to_dict(orient="records")

    # driver-associated variance (one sample per patient cohorts pooled with first region of SYN_C)
    first = samples.reset_index().groupby("patient_id").head(1).set_index("sample_id")
    D = pd.concat([R[c]["C"] for c in R]).loc[first.index]
    G = first[sv.usable_drivers(first)]
    res, _ = sv.driver_variance_bootstrap(D, G.to_numpy(float), first.patient_id.to_numpy(), n_boot=n_boot)
    summary["driver_associated_variance"] = {"drivers": list(G.columns), **json.loads(json.dumps(res, default=float))}
    summary["driver_variance_beyond_purity"] = {c: round(sv.purity_adjusted_driver_variance(D[c], G.to_numpy(float), first.dna_purity.to_numpy()), 3)
                                                for c in ["MAL_AC", "AXIS"]}
    d = res["differences"]["MAL_AC-AXIS"]
    summary["fixed_effect_meta_example"] = sv.fixed_effect_meta([d["point"], d["point"] / 2], [d["se"], d["se"] * 1.5])

    # between-patient share in the multi-region cohort
    mr = samples[samples.cohort == "SYN_C"]
    Cm = R["SYN_C"]["C"].loc[mr.index].reset_index(drop=True)
    g = mr.patient_id.to_numpy()
    summary["ICC1"] = {c: round(sv.icc1(Cm[c], g), 3) for c in ["MAL_AC"] + sv.ECO7}
    con, _ = sv.bootstrap_icc_contrast(Cm, g, ["MAL_AC"], sv.ECO7, n_boot=n_boot)
    summary["ICC_balance_minus_mean_ECO7"] = con

    # expression-matched collection (N1) from the supplied pool features (no sample-level data needed)
    pool = list(P["pool_features"].index)
    lit = sv.literature_sets_on_pool(P["sigs"], P["ids"], pool)
    F = sv.standardized_matching_features(P["pool_features"])
    n1 = sv.expression_matched_collection(lit, F, np.random.default_rng(sv.SEED_COLLECTIONS + 1000 + 0))
    pd.DataFrame([(k, g_) for k, v in n1.items() for g_ in v], columns=["set_id", "gene_symbol"]).to_csv(out / "N1_collection_0.tsv", sep="\t", index=False)
    summary["N1_collection_0"] = {"sets": len(n1), "gene_occurrences": int(sum(len(v) for v in n1.values()))}

    json.dump(summary, open(out / "example_summary.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in summary.items() if not k.startswith("heldout")}, indent=1, default=float)[:3000])
    print(f"\nOutputs written to {out}/  (SYNTHETIC; execution demonstration only)")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "example_output"))
    ap.add_argument("--n-boot", type=int, default=200)
    a = ap.parse_args()
    main(Path(a.out), a.n_boot)
