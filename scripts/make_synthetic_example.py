"""Generate the SYNTHETIC example data in examples/ (deterministic, seed 1).

The data are simulated from a simple latent-factor count model. They exist only to demonstrate that the
code runs end to end. They contain no patient information and carry no biological evidence; numbers
computed from them must not be compared with the manuscript.

Usage:  python scripts/make_synthetic_example.py
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA, EX = ROOT / "data", ROOT / "examples"


def main(seed: int = 1, n_per_cohort: int = 24, n_filler: int = 300):
    rng = np.random.default_rng(seed)
    sig = pd.read_csv(DATA / "signatures_201.tsv", sep="\t")
    ref = json.load(open(DATA / "reference_gene_sets.json"))
    pcg = [g.strip() for g in open(DATA / "protein_coding_genes.txt") if g.strip()]
    panel = ref["coordinates"]
    panel_genes = sorted(set(g for c in panel.values() for k in ("A", "B") for g in c[k]) | set(ref["programme_sets"]["OPC"] + ref["programme_sets"]["NPC"]))
    sig_genes = sorted(set(sig.gene_symbol))
    used = set(panel_genes) | set(sig_genes)
    filler = [g for g in rng.permutation(pcg) if g not in used][:n_filler]
    genes = sorted(used | set(filler))

    cohorts = ["SYN_A", "SYN_B", "SYN_C"]
    n = n_per_cohort * len(cohorts)
    samples = pd.DataFrame({"sample_id": [f"S{i + 1:03d}" for i in range(n)],
                            "cohort": np.repeat(cohorts, n_per_cohort)})
    # cohort SYN_C is multi-region: 8 patients x 3 regions; others one sample per patient
    pid = [f"P{i + 1:03d}" for i in range(2 * n_per_cohort)] + [f"P{100 + i // 3:03d}" for i in range(n_per_cohort)]
    samples["patient_id"] = pid
    npat = samples.patient_id.nunique()
    pat_index = {p: i for i, p in enumerate(samples.patient_id.unique())}
    # patient-level latent factors and drivers; sample-level ecology varies within patients
    drivers = pd.DataFrame(rng.random((npat, 4)) < [0.45, 0.15, 0.15, 0.2], columns=["EGFR_amp", "PDGFRA_amp", "CDK4_amp", "NF1_loss"]).astype(int)
    balance_p = 0.9 * drivers.EGFR_amp - 0.8 * drivers.PDGFRA_amp + 0.4 * drivers.NF1_loss + rng.normal(0, 0.7, npat)
    eco_s = rng.normal(0, 1, n)
    bal_s = balance_p.to_numpy()[[pat_index[p] for p in samples.patient_id]] + rng.normal(0, 0.3, n)
    purity = 1 / (1 + np.exp(0.9 * eco_s + rng.normal(0, 0.3, n)))

    gidx = {g: i for i, g in enumerate(genes)}
    load_eco = rng.normal(0, 0.15, len(genes))
    load_bal = np.zeros(len(genes))
    for c in ["MYE", "LYM", "VAS"]:
        for g in panel[c]["A"]:
            load_eco[gidx[g]] = 0.8
    for g in panel["MYE"]["B"]:
        load_eco[gidx[g]] = -0.5
    for g in panel["HYPOX"]["A"] + panel["MAL_MES"]["A"]:
        load_eco[gidx[g]] = 0.5
    for g in panel["MAL_AC"]["A"]:
        load_bal[gidx[g]] = 0.7
    for g in panel["MAL_AC"]["B"]:
        load_bal[gidx[g]] = -0.7
    sig_g = [gidx[g] for g in sig_genes if g in gidx]
    load_eco[sig_g] += rng.normal(0.35, 0.3, len(sig_g))           # published genes lean ecological (by construction)
    mu = rng.normal(3.0, 1.5, len(genes))
    batch = {c: rng.normal(0, 0.2, len(genes)) for c in cohorts}
    lib = rng.lognormal(0, 0.25, n)
    eta = (mu[:, None] + np.outer(load_eco, eco_s) + np.outer(load_bal, bal_s)
           + np.column_stack([batch[c] for c in samples.cohort]) + rng.normal(0, 0.25, (len(genes), n)))
    counts = rng.poisson(np.exp(eta) * lib)
    X = pd.DataFrame(counts, index=genes, columns=samples.sample_id)
    X.index.name = "gene_symbol"
    X.to_csv(EX / "synthetic_example.csv")
    s = samples.copy()
    for d in drivers.columns:
        s[d] = drivers[d].to_numpy()[[pat_index[p] for p in samples.patient_id]]
    s["dna_purity"] = np.round(purity, 4)
    s.to_csv(EX / "synthetic_samples.csv", index=False)
    print(f"wrote {X.shape[0]} genes x {X.shape[1]} samples (SYNTHETIC)")


if __name__ == "__main__":
    main()
