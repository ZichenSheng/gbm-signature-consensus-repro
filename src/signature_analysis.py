"""Signature scoring, reference coordinates and the literature consensus.

Every function in this module is extracted from the validated analysis code that produced the
manuscript's results. Formulas, filters, degrees of freedom, tie handling and random-number usage are
kept as they were, so that the functions reproduce the frozen per-cohort scores numerically
(see README, "Verification against frozen outputs").

Conventions
-----------
* Expression matrices are pandas DataFrames with genes as rows (upper-case HGNC symbols) and samples
  as columns.
* ``Z`` always denotes a per-gene standardized matrix (z-scores across the samples of one
  cohort-platform, sample standard deviation, ddof = 1).
* Score tables are samples x scores.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SEED_CONTROL_GENES = 20261107      # seed of the expression-matched hypoxia control genes

COORDINATES = ["MAL_MES", "MAL_AC", "PROL", "MYE_ORIG", "MYE_HYP", "MYE", "LYM", "VAS", "GLN", "HYPOX"]
ECOLOGICAL = ["MYE", "LYM", "VAS", "GLN", "HYPOX"]
COORDINATE_NAMES = {
    "MYE": "Myeloid balance", "LYM": "Lymphoid balance", "VAS": "Vascular balance",
    "GLN": "Non-malignant glial balance", "HYPOX": "Hypoxia",
    "MAL_AC": "AC-developmental balance", "MAL_MES": "MES-like vs developmental",
    "PROL": "Proliferation", "MYE_ORIG": "Myeloid origin", "MYE_HYP": "Hypoxic myeloid state",
    "DEV_rel": "Relative developmental programme", "AXIS": "Literature consensus",
}


# ----------------------------------------------------------------------------- inputs
def harmonize_symbols(expr: pd.DataFrame, how: str = "sum") -> pd.DataFrame:
    """Strip and upper-case gene symbols and collapse duplicated symbols.

    Production used ``sum`` for count-like units (CPTAC and CPTAC-3 counts, GBmap and GBM-Space
    pseudobulks) and ``mean`` for the distributed CGGA, GLASS, TCGA and array matrices.
    """
    x = expr.copy()
    x.index = x.index.astype(str).str.strip().str.upper()
    g = x.groupby(level=0)
    return g.sum() if how == "sum" else g.mean()


def load_signatures(path: str) -> tuple[list[str], dict[str, list[str]]]:
    """Read the long-format membership table (signature_id, gene_symbol).

    Returns signature identifiers in file order and, per signature, its de-duplicated gene list in
    file order (both orders matter for exact reproduction of seeded procedures).
    """
    c = pd.read_csv(path, sep="\t", dtype=str)
    c["gene_symbol"] = c.gene_symbol.str.strip().str.upper()
    ids = list(dict.fromkeys(c.signature_id))
    sets = {m: list(dict.fromkeys(c.loc[c.signature_id == m, "gene_symbol"])) for m in ids}
    return ids, sets


def load_formulas(path: str) -> dict[str, dict[str, float]]:
    """Read published weighted formulas (signature_id, gene_symbol, coefficient)."""
    co = pd.read_csv(path, sep="\t", dtype={"signature_id": str, "gene_symbol": str})
    co["gene_symbol"] = co.gene_symbol.str.upper()
    return {m: dict(zip(g.gene_symbol, g.coefficient.astype(float))) for m, g in co.groupby("signature_id")}


# ----------------------------------------------------------------------------- normalization
def detected_protein_coding(A: pd.DataFrame, protein_coding, min_detected: float = 0.5) -> pd.Index:
    """Nuclear protein-coding genes detected (> 0) in at least ``min_detected`` of samples.

    ``A`` is on a linear scale. ``protein_coding`` is the supplied symbol list (MT- genes excluded).
    """
    B = A.loc[A.index.isin(set(protein_coding))]
    return B.index[(B > 0).mean(axis=1) >= min_detected]


def median_of_ratios_size_factors(B: pd.DataFrame) -> pd.Series:
    """Median-of-ratios size factors computed over the upper half (by mean log expression) of the
    genes expressed in every sample; rescaled to geometric mean one."""
    lg = np.log(B.where(B > 0))
    full = lg.dropna()
    full = full.loc[full.mean(axis=1) >= full.mean(axis=1).median()]
    ref = full.mean(axis=1)
    sf = np.exp((full.sub(ref, axis=0)).median(axis=0))
    return sf / np.exp(np.log(sf).mean())


def normalize_rnaseq(A: pd.DataFrame, protein_coding, min_detected: float = 0.5) -> pd.DataFrame:
    """RNA-sequencing processing used for every RNA-seq cohort-platform.

    ``A``: genes x samples on a linear scale (counts, FPKM or TPM as distributed; log-scale
    distributions must first be returned to linear scale with 2**x - 1). Keeps detected protein-coding
    genes, divides by median-of-ratios size factors and returns log2(x / size factor + 1).
    """
    genes = detected_protein_coding(A, protein_coding, min_detected)
    B = A.loc[genes].astype(float)
    sf = median_of_ratios_size_factors(B)
    return np.log2(B / sf + 1)


def zscore_genes(X: pd.DataFrame) -> pd.DataFrame:
    """Per-gene standardization across samples (genes with zero variance removed; ddof = 1)."""
    X = X.loc[X.std(axis=1) > 0]
    return X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1), axis=0)


def zs(x) -> np.ndarray:
    """Standardize a vector (population SD, NaN-aware)."""
    x = np.asarray(x, float)
    s = np.nanstd(x)
    return (x - np.nanmean(x)) / s if s > 0 else x * 0


def _std_cols(M) -> np.ndarray:
    M = np.asarray(M, float)
    return (M - M.mean(axis=0)) / M.std(axis=0)


# ----------------------------------------------------------------------------- published signatures
def membership_scores(Z: pd.DataFrame, signatures: dict, ids=None, min_genes: int = 3,
                      min_coverage: float = 0.5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Unsigned gene-membership score: mean z-score of a signature's measured genes.

    A score is computed only when at least ``min_genes`` genes and at least ``min_coverage`` of the
    signature's genes are measured; otherwise the column is NaN. Returns (scores, coverage table).
    """
    ids = list(signatures) if ids is None else list(ids)
    gi = {g: i for i, g in enumerate(Z.index)}
    A = Z.to_numpy(float)
    out = np.full((A.shape[1], len(ids)), np.nan)
    cov = []
    for j, m in enumerate(ids):
        g = signatures[m]
        idx = [gi[x] for x in g if x in gi]
        c = len(idx) / max(len(g), 1)
        cov.append((len(idx), c))
        if len(idx) >= min_genes and c >= min_coverage:
            out[:, j] = A[idx, :].mean(axis=0)
    return (pd.DataFrame(out, index=Z.columns, columns=ids),
            pd.DataFrame(cov, index=ids, columns=["n_mapped", "coverage"]))


def scorable_signatures(Zs, signatures: dict, ids, min_coverage: float = 0.5) -> list[str]:
    """Signatures (in ``ids`` order) whose membership score is computable in every matrix of ``Zs``."""
    out = []
    for r in ids:
        if all(membership_scores(Z, signatures, [r], min_coverage=min_coverage)[0][r].notna().all() for Z in Zs):
            out.append(r)
    return out


def weighted_formula_scores(Z: pd.DataFrame, formulas: dict, sd_ref: pd.Series,
                            min_coverage: float = 0.8) -> pd.DataFrame:
    """Weighted score of a published formula on a platform-comparable scale.

    score = sum_g beta_g * sd_ref_g * z_g, where sd_ref_g is the median standard deviation of gene g's
    log2 expression across the five development cohorts. Computed when at least max(2, 80%) of the
    formula's genes are measured.
    """
    out = {}
    for m, beta in formulas.items():
        gg = [g for g in beta if g in Z.index and g in sd_ref.index and np.isfinite(sd_ref[g])]
        if len(gg) < max(2, min_coverage * len(beta)):
            continue
        b = np.array([beta[g] for g in gg])
        out[m] = (b * sd_ref[gg].to_numpy()) @ Z.loc[gg].to_numpy(float)
    return pd.DataFrame(out, index=Z.columns)


# ----------------------------------------------------------------------------- literature consensus
def fit_consensus_loadings(S: pd.DataFrame, orient_to=None) -> pd.Series:
    """First principal component of standardized membership scores (loading cohort).

    ``S``: samples x signatures, columns containing NaN are dropped. Scores are standardized with
    ddof = 1 and the first right singular vector is returned. PCA is sign-invariant; if ``orient_to``
    is given the sign is chosen so that the component correlates positively with it.

    Production fitted the frozen loadings once in CGGA693 (see README, "Literature consensus", for the
    loading-cohort preprocessing); the frozen vector is supplied in ``data/consensus_loadings.tsv`` and
    is applied unchanged everywhere else.
    """
    S = S.dropna(axis=1)
    A = S.to_numpy(float)
    A = (A - A.mean(axis=0)) / A.std(axis=0, ddof=1)
    u = np.linalg.svd(A, full_matrices=False)[2][0]
    if orient_to is not None and np.corrcoef(A @ u, np.asarray(orient_to, float))[0, 1] < 0:
        u = -u
    return pd.Series(u, index=S.columns, name="loading")


def consensus_score(Z: pd.DataFrame, signatures: dict, loadings: pd.Series,
                    min_coverage: float = 0.5) -> tuple[pd.Series, list[str]]:
    """Literature consensus: frozen loadings applied unchanged to standardized membership scores.

    Uses the loading signatures that are scorable in this cohort-platform (>= 3 genes and >= 50% of
    genes measured), standardizes each score within the cohort, projects on the frozen loadings and
    standardizes the result. Returns (consensus, signatures used).
    """
    reps = [r for r in loadings.index]
    S, _ = membership_scores(Z, signatures, reps, min_coverage=min_coverage)
    used = [r for r in reps if S[r].notna().all()]
    ax = _std_cols(S[used].to_numpy()) @ loadings.loc[used].to_numpy()
    return pd.Series(zs(ax), index=Z.columns, name="AXIS"), used


# ----------------------------------------------------------------------------- reference coordinates
def set_mean(Z: pd.DataFrame, genes, min_frac: float = 0.5, min_n: int = 5):
    """Mean z-score of a gene set if >= ``min_n`` and >= ``min_frac`` of its genes are measured."""
    gg = [g for g in dict.fromkeys(genes) if g in Z.index]
    if len(gg) < min_n or len(gg) < min_frac * len(genes):
        return None, gg
    return Z.loc[gg].mean(axis=0), gg


def expression_matched_controls(mean_expr: pd.Series, genes, rng, n_bins: int = 25, per_gene: int = 10,
                                pool=None) -> list[str]:
    """Expression-bin matched control genes (control-matched construction of the hypoxia coordinate).

    Genes are split into ``n_bins`` bins of mean log expression (ties broken by order); for each set
    gene, ``per_gene`` genes from its bin (excluding the set) are drawn without replacement.
    """
    m = mean_expr.dropna()
    if pool is not None:
        m = m.loc[m.index.intersection(pool)]
    bins = pd.qcut(m.rank(method="first"), n_bins, labels=False)
    gs = set(genes)
    out = []
    for g in genes:
        if g not in bins.index:
            continue
        cand = bins.index[(bins == bins[g]).to_numpy()]
        cand = [c for c in cand if c not in gs]
        out += list(rng.choice(cand, size=min(per_gene, len(cand)), replace=False))
    return out


def score_reference_coordinates(Z: pd.DataFrame, gene_sets: dict, mean_expr: pd.Series,
                                seed: int = SEED_CONTROL_GENES, min_frac: float = 0.5, min_n: int = 5,
                                standardize: bool = True):
    """Score the frozen reference coordinates.

    Each coordinate is mean z(set A) - mean z(set B); hypoxia is its canonical 16-gene set minus
    expression-matched control genes drawn from ``mean_expr`` (per-gene mean log expression of the same
    cohort-platform) with a fresh generator seeded by ``seed``. Contrasts are standardized within the
    cohort-platform (ddof = 1). Returns (contrasts, set-A levels, coverage).
    """
    rng = np.random.default_rng(seed)
    con, lev, cov = {}, {}, {}
    for c in COORDINATES:
        A = gene_sets[c]["A"]
        B = gene_sets[c]["B"]
        a, ga = set_mean(Z, A, min_frac, min_n)
        if c == "HYPOX":
            if a is None or mean_expr is None:
                continue
            ctl = expression_matched_controls(mean_expr.loc[Z.index.intersection(mean_expr.index)], ga, rng)
            b = Z.loc[[g for g in ctl if g in Z.index]].mean(axis=0)
            gb = ctl
        else:
            b, gb = set_mean(Z, B, min_frac, min_n)
        cov[c] = (len(ga), len(A), len(gb), len(B))
        if a is None or b is None:
            continue
        con[c] = a - b
        lev[c] = a
    C = pd.DataFrame(con)
    Lv = pd.DataFrame(lev)
    if standardize:
        C = (C - C.mean()) / C.std()
        Lv = (Lv - Lv.mean()) / Lv.std()
    return C, Lv, cov


def relative_programme_scores(Z: pd.DataFrame, programme_sets: dict, reference_set) -> pd.DataFrame:
    """Relative programme score: zs(mean z(programme) - mean z(pan-malignant reference set)).

    With ``programme_sets = {"DEV": developmental genes}`` and the pan-malignant set this gives the
    relative developmental programme (DEV_rel) used in the manuscript.
    """
    b, _ = set_mean(Z, reference_set)
    out = {}
    for nm, gset in programme_sets.items():
        a, _ = set_mean(Z, gset)
        if a is None or b is None:
            continue
        out[f"{nm}_rel"] = zs((a - b).to_numpy())
    return pd.DataFrame(out, index=Z.columns)


def load_reference_gene_sets(path: str) -> dict:
    """Read data/reference_gene_sets.json into {coordinate: {"A": [...], "B": [...]}} plus extras."""
    import json
    d = json.load(open(path))
    return d


# ----------------------------------------------------------------------------- reference-set derivation
# GBmap level-4 annotation groups used by the pre-specified derivation rules.
_MES = ["MES-like hypoxia independent", "MES-like hypoxia/MHC"]
_AC = ["AC-like"]
_DEV = ["OPC-like", "NPC-like OPC", "NPC-like neural"]
_CYC = ["AC-like Prolif", "NPC-like Prolif", "OPC-like Prolif"]
_BDM = ["TAM-BDM INF", "TAM-BDM MHC", "TAM-BDM anti-infl", "Mono anti-infl", "Mono naive"]
_HYPM = ["Mono hypoxia", "TAM-BDM hypoxia/MES"]
_MG = ["TAM-MG aging sig", "TAM-MG pro-infl I", "TAM-MG pro-infl II"]
COMPARTMENTS = {"Neoplastic": "MAL", "Myeloid": "MYE", "Lymphoid": "LYM", "Vascular": "VAS", "Glial-Neuronal": "GLN"}
GROUP_MAP = {**{s: "MES" for s in _MES}, **{s: "AC" for s in _AC}, **{s: "DEV" for s in _DEV},
             **{s: "CYC" for s in _CYC}, **{s: "BDM" for s in _BDM}, **{s: "HYPM" for s in _HYPM},
             **{s: "MG" for s in _MG}, **{s: "DC" for s in ["DC1", "DC2", "DC3", "cDC1", "cDC2", "pDC"]},
             "Mast": "MAST"}
SET_ORDER = ["COMP_MAL", "COMP_MYE", "COMP_LYM", "COMP_VAS", "COMP_GLN", "MAL_CYC", "MAL_MES", "MAL_AC",
             "MAL_DEV", "MYE_HYPM", "MYE_BDMALL", "MYE_MG"]
PANEL_DEFINITION = {"MAL_MES": ("MAL_MES", "MAL_DEV"), "MAL_AC": ("MAL_AC", "MAL_DEV"), "PROL": ("MAL_CYC", "COMP_MAL"),
                    "MYE_ORIG": ("MYE_BDMALL", "MYE_MG"), "MYE_HYP": ("MYE_HYPM", "COMP_MYE"),
                    "MYE": ("COMP_MYE", "COMP_MAL"), "LYM": ("COMP_LYM", "COMP_MAL"), "VAS": ("COMP_VAS", "COMP_MAL"),
                    "GLN": ("COMP_GLN", "COMP_MAL")}
HYPOXIA_GENES = ("VEGFA ADM CA9 NDRG1 ANGPTL4 PGK1 SLC2A1 BNIP3 P4HA1 LDHA ERO1A ANKRD37 EGLN3 PDK1 "
                 "BHLHE40 DDIT4").split()


def _profiles(S, gm, mask, key, min_cells=20):
    sub = gm[mask & gm[key].notna()]
    lab = sub.donor.astype(str) + "||" + sub[key].astype(str)
    sums = S.loc[sub.index].groupby(lab.to_numpy()).sum()
    ncell = sub.groupby(lab.to_numpy()).n.sum()
    sums = sums.loc[ncell[ncell >= min_cells].index]
    lc = np.log2(sums.div(sums.sum(axis=1), axis=0) * 1e6 + 1)
    grp = pd.Index(lc.index).str.split("||", regex=False).str[1]
    return lc.groupby(grp.to_numpy()).median()


def derive_reference_gene_sets(S: pd.DataFrame, meta: pd.DataFrame, protein_coding, top: int = 20,
                               d_comp: float = 2.0, d_ref: float = 1.0, d_rc: float = 1.0,
                               min_expr: float = np.log2(11)) -> dict:
    """Pre-specified marker rules applied to the GBmap core atlas.

    ``S``: rows = (donor x level-4 annotation) groups, columns = genes, values = summed counts.
    ``meta``: indexed like ``S`` with columns donor, l1, l2, l4 and n (cells). Groups with >= 20 cells
    are converted to log2 CPM + 1 per donor and summarized by the median over donors. Compartment
    markers need >= 2 log2 units above every other compartment, malignant-state and myeloid-state
    markers >= 1 unit above every other state of their compartment and >= 1 unit above every other
    compartment, and every marker >= log2(11); the 20 most specific genes are kept, gene-disjoint
    across sets (canonical hypoxia genes reserved).
    """
    pcg = set(protein_coding)
    idx = pd.Index(S.columns)
    bad = idx.str.startswith(("MT-", "RPL", "RPS", "ENSG", "LINC", "MIR")) | ~idx.isin(pcg)
    S = S[idx[~bad]]
    gm = meta.copy()
    gm["comp"] = np.where(gm.l1 == "Neoplastic", "MAL", gm.l2.map(COMPARTMENTS).fillna("OTHER"))
    gm["grp"] = gm.l4.map(GROUP_MAP)
    gm["bdm_all"] = np.where(gm.grp.isin(["BDM", "HYPM"]), "BDMALL", gm.grp)
    mask = pd.Series(True, index=gm.index)
    pc = _profiles(S, gm, mask, "comp").loc[["MAL", "MYE", "LYM", "VAS", "GLN"]]
    pm = _profiles(S, gm, mask & (gm.comp == "MAL"), "grp")
    py = _profiles(S, gm, mask & (gm.comp == "MYE"), "grp")
    pb = _profiles(S, gm, mask & (gm.comp == "MYE"), "bdm_all")

    def sel(s, refs, own, dref, dcomp):
        spec_comp = s - pc.drop(index=own).max(axis=0)
        if refs is None:
            ok = (spec_comp >= dcomp) & (s >= min_expr)
            score = spec_comp[ok]
        else:
            spec_ref = s - refs.max(axis=0)
            ok = (spec_ref >= dref) & (spec_comp >= dcomp) & (s >= min_expr)
            score = spec_ref[ok]
        return list(score.sort_values(ascending=False).index)

    sets = {f"COMP_{c}": sel(pc.loc[c], None, c, None, d_comp) for c in ["MAL", "MYE", "LYM", "VAS", "GLN"]}
    mal = ["MES", "AC", "DEV", "CYC"]
    for s_ in mal:
        sets[f"MAL_{s_}"] = sel(pm.loc[s_], pm.loc[[r for r in mal if r != s_]], "MAL", d_ref, d_rc)
    mye = ["BDM", "HYPM", "MG", "DC", "MAST"]
    for s_ in ["HYPM", "MG"]:
        sets[f"MYE_{s_}"] = sel(py.loc[s_], py.loc[[r for r in mye if r != s_ and r in py.index]], "MYE", d_ref, d_rc)
    sets["MYE_BDMALL"] = sel(pb.loc["BDMALL"], pb.loc[[r for r in ["MG", "DC", "MAST"] if r in pb.index]], "MYE", d_ref, d_rc)
    used = set(HYPOXIA_GENES)
    out = {}
    for k in SET_ORDER:
        pick = [g for g in sets[k] if g not in used][:top]
        used |= set(pick)
        out[k] = pick
    return out


def panel_from_sets(sets: dict) -> dict:
    """Arrange derived gene sets into the coordinate panel ({coordinate: {"A": ..., "B": ...}})."""
    out = {c: {"A": sets[a], "B": sets[b]} for c, (a, b) in PANEL_DEFINITION.items()}
    out["HYPOX"] = {"A": list(HYPOXIA_GENES), "B": []}
    return out
