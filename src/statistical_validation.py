"""Statistical procedures behind the manuscript's central quantitative claims.

Extracted from the validated analysis code. Resampling functions consume a numpy Generator in exactly
the order used in production, so that, with the production seeds listed in the README, they reproduce
the frozen confidence intervals and null distributions draw for draw.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SEED_DRIVER_BOOTSTRAP = 20261108       # driver-associated variance bootstrap and random contrasts
SEED_COLLECTIONS = 20261108            # base seed of the matched random collections
SEED_ICC_COORDINATE = 20261108         # per-coordinate ICC intervals
SEED_ICC_CONTRAST = 20261110           # single-coordinate ICC contrast (balance vs ecological mean)
SEED_STRUCTURE_NULL = 20261107         # matched-panel null of pairwise correlation structure

CORE_DRIVERS = ["EGFR_amp", "PDGFRA_amp", "CDK4_amp", "NF1_loss"]
ECO7 = ["AXIS", "MYE", "LYM", "VAS", "GLN", "HYPOX", "MAL_MES"]


# ----------------------------------------------------------------------------- small helpers
def zs(x) -> np.ndarray:
    x = np.asarray(x, float)
    s = np.nanstd(x)
    return (x - np.nanmean(x)) / s if s > 0 else x * 0


def std_cols(M) -> np.ndarray:
    """Column standardization (population SD; constant columns left centred)."""
    M = np.asarray(M, float)
    s = M.std(axis=0)
    s[s == 0] = 1
    return (M - M.mean(axis=0)) / s


def rank(x) -> np.ndarray:
    """Average ranks (ties averaged)."""
    return pd.Series(np.asarray(x, float)).rank().to_numpy()


def percentile_ci(values, q=(2.5, 97.5)):
    v = np.asarray([x for x in values if x is not None and np.isfinite(x)], float)
    return [float(np.percentile(v, q[0])), float(np.percentile(v, q[1]))] if len(v) else [np.nan, np.nan]


def r2(y, X) -> float:
    """In-sample R^2 of OLS with intercept."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    Xd = np.column_stack([np.ones(len(y)), X])
    b = np.linalg.lstsq(Xd, y, rcond=None)[0]
    return float(1 - ((y - Xd @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def adjusted_r2(y, X) -> float:
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    n, k = X.shape
    if k == 0:
        return np.nan
    return 1 - (1 - r2(y, X)) * (n - 1) / (n - k - 1)


def residualize(y, X) -> np.ndarray:
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    Xd = np.column_stack([np.ones(len(y)), X])
    return y - Xd @ np.linalg.lstsq(Xd, y, rcond=None)[0]


# ----------------------------------------------------------------------------- pairwise structure
def spearman_corr_matrix(S) -> np.ndarray:
    """Spearman correlation matrix of the columns of S (samples x scores)."""
    return np.corrcoef(pd.DataFrame(np.asarray(S, float)).rank(axis=0).to_numpy(), rowvar=False)


def structure_reproducibility(S1, S2) -> float:
    """Cross-cohort reproducibility of pairwise score correlations: Spearman correlation between the
    upper triangles of the two cohorts' Spearman correlation matrices (same signatures, same order)."""
    C1, C2 = spearman_corr_matrix(S1), spearman_corr_matrix(S2)
    iu = np.triu_indices(np.asarray(S1).shape[1], 1)
    return float(np.corrcoef(rank(C1[iu]), rank(C2[iu]))[0, 1])


def gene_features(x: pd.DataFrame, transform_is_log: bool) -> pd.DataFrame:
    """Per-gene matching features: mean and SD of log2 expression and detected fraction."""
    lx = x if transform_is_log else np.log2(x + 1.0)
    return pd.DataFrame({"mean": lx.mean(axis=1), "sd": lx.std(axis=1, ddof=1), "det": (x > 0).mean(axis=1)})


class FeatureMatcher:
    """Matched random gene panels for the pairwise-structure null.

    Each member gene is replaced by a uniform draw from its k nearest genes of the same biotype in
    standardized (mean, sd, detection) space, excluding the panel's own genes; set size is preserved.
    """

    def __init__(self, feats: pd.DataFrame, biotype: dict, k: int = 30, exclude=()):
        f = feats.copy()
        f = f.loc[~f.index.isin(set(exclude))]
        f["bt"] = [biotype.get(g, "unknown") for g in f.index]
        self.F = f
        self.k = k
        cols = ["mean", "sd", "det"]
        self.std = (f[cols] - f[cols].mean()) / f[cols].std()
        self.by_bt = {b: grp.index.to_numpy() for b, grp in f.groupby("bt")}
        self.cache = {}

    def neighbours(self, g, featsrc):
        if g in self.cache:
            return self.cache[g]
        if g in self.std.index:
            v = self.std.loc[g].to_numpy()
            bt = self.F.loc[g, "bt"]
        else:
            v = featsrc(g)
            bt = "protein-coding"
            if v is None:
                return None
        pool = self.by_bt.get(bt, self.F.index.to_numpy())
        if len(pool) < self.k + 5:
            pool = self.F.index.to_numpy()
        P = self.std.loc[pool].to_numpy()
        d = np.sqrt(((P - v) ** 2).sum(axis=1))
        nn = pool[np.argsort(d)[: self.k + 10]]
        self.cache[g] = nn
        return nn

    def draw(self, genes, rng, featsrc=lambda g: None):
        own = set(genes)
        out = []
        for g in genes:
            nn = self.neighbours(g, featsrc)
            if nn is None:
                continue
            cand = [c for c in nn if c not in own and c not in out][: self.k]
            if not cand:
                continue
            out.append(cand[rng.integers(len(cand))])
        return out


# ----------------------------------------------------------------------------- matched random collections
def literature_sets_on_pool(signatures: dict, ids, pool, min_genes: int = 3, min_frac: float = 0.5) -> dict:
    """Published gene sets restricted to the common gene pool (kept if >= 3 genes and >= 50%)."""
    pset = set(pool)
    out = {}
    for r in ids:
        g0 = list(dict.fromkeys(signatures[r]))
        gg = [g for g in g0 if g in pset]
        if len(gg) >= min_genes and len(gg) >= min_frac * len(g0):
            out[r] = gg
    return out


def standardized_matching_features(features: pd.DataFrame) -> pd.DataFrame:
    """Standardize the pool's (mean, sd, det) features (ddof = 1) for nearest-neighbour matching."""
    M = features[["mean", "sd", "det"]]
    return (M - M.mean()) / M.std()


def pooled_standardized_expression(Z_by_cohort, pool) -> np.ndarray:
    """Genes x samples: per-cohort z matrices on the pool, concatenated over cohorts, missing set to 0,
    then re-standardized per gene (population SD). Used for within-set co-expression."""
    Zp = np.hstack([Z.reindex(pool).to_numpy(float) for Z in Z_by_cohort])
    Zp = np.nan_to_num(Zp)
    return (Zp - Zp.mean(1, keepdims=True)) / Zp.std(1, keepdims=True)


def within_set_coherence(genes, Zp: np.ndarray, pool_index: dict) -> float:
    """Mean pairwise correlation among a set's genes in the pooled standardized expression."""
    ii = [pool_index[g] for g in genes]
    X = Zp[ii]
    C = X @ X.T / Zp.shape[1]
    k = len(ii)
    return float((C.sum() - k) / (k * (k - 1))) if k > 1 else np.nan


def expression_matched_collection(lit_sets: dict, features: pd.DataFrame, rng, k_nn: int = 10) -> dict:
    """Expression-matched collection (N1).

    Every published gene is mapped, through one map shared by all signatures that contain it, to one of
    its ``k_nn`` nearest non-published pool genes in standardized (mean, sd, det) space, without reuse;
    set sizes and overlaps are preserved. ``features`` = standardized_matching_features(pool features).
    """
    litgenes = sorted(set(g for v in lit_sets.values() for g in v))
    cand = np.array([g for g in features.index if g not in set(litgenes)])
    A = features.loc[cand].to_numpy()
    used = np.zeros(len(cand), bool)
    mp = {}
    for g in rng.permutation(litgenes):
        d = ((A - features.loc[g].to_numpy()) ** 2).sum(axis=1)
        d[used] = np.inf
        nn = np.argpartition(d, k_nn)[:k_nn]
        j = rng.choice(nn)
        used[j] = True
        mp[g] = cand[j]
    return {r: [mp[g] for g in v] for r, v in lit_sets.items()}


def coexpression_matched_collection(lit_sets: dict, features: pd.DataFrame, Zp: np.ndarray, pool, rng,
                                    n_cand: int = 300, n_best: int = 3) -> dict:
    """Co-expression-matched collection (N2).

    Built set by set from the ``n_cand`` nearest feature-matched non-published candidates of each
    published gene: the first gene is a random candidate; each further gene is drawn at random from the
    ``n_best`` candidates that keep the mean within-set correlation closest to the published set's.
    """
    pool_arr = np.array(pool)
    pidx = {g: i for i, g in enumerate(pool)}
    litgenes = sorted(set(g for v in lit_sets.values() for g in v))
    litmask = np.isin(pool_arr, litgenes)
    A1 = features.loc[list(pool)].to_numpy()
    nS = Zp.shape[1]
    out = {}
    for r, v in lit_sets.items():
        target = within_set_coherence(v, Zp, pidx)
        chosen = []
        for g in v:
            d = ((A1 - features.loc[g].to_numpy()) ** 2).sum(axis=1)
            d[litmask] = np.inf
            for gg in chosen:
                d[pidx[gg]] = np.inf
            nn = np.argpartition(d, n_cand)[:n_cand]
            if not chosen:
                out_g = pool_arr[rng.choice(nn)]
            else:
                X = Zp[[pidx[x] for x in chosen]]
                c = (Zp[nn] @ X.T / nS).sum(axis=1)
                k = len(chosen)
                cur = (Zp[[pidx[x] for x in chosen]] @ X.T / nS)
                cur_sum = (cur.sum() - k) / 2
                new_mean = (cur_sum + c) / (k * (k + 1) / 2)
                best = np.argsort(np.abs(new_mean - target))[:n_best]
                out_g = pool_arr[nn[rng.choice(best)]]
            chosen.append(out_g)
        out[r] = chosen
    return out


def collection_scores(sets: dict, Z_pool: np.ndarray, pool) -> np.ndarray:
    """Samples x sets membership scores of a collection on a pool-indexed z matrix (genes x samples)."""
    pidx = {g: i for i, g in enumerate(pool)}
    return np.vstack([Z_pool[[pidx[g] for g in v]].mean(axis=0) for v in sets.values()]).T


def leading_components(scores, n: int = 20):
    """Principal-component scores of a collection (columns standardized first) and explained variance."""
    S = std_cols(scores)
    U, s, _ = np.linalg.svd(S - S.mean(axis=0), full_matrices=False)
    m = min(n, len(s))
    return U[:, :m] * s[:m], s ** 2 / (s ** 2).sum()


# ----------------------------------------------------------------------------- ecological alignment and recovery
def ecological_share(score, coordinates: pd.DataFrame, ecological=("MYE", "LYM", "VAS", "GLN", "HYPOX")) -> float:
    """In-sample R^2 of a score on the five ecological coordinates (OLS with intercept)."""
    return r2(np.asarray(score, float), std_cols(coordinates[list(ecological)].to_numpy()))


def loco_ridge_recovery(X_by_cohort: dict, y_by_cohort: dict, alphas=np.logspace(-1, 5, 25)) -> pd.DataFrame:
    """Leave-one-cohort-out ridge prediction of a reference coordinate.

    Predictors are standardized within each cohort (missing -> 0) and the target is standardized within
    each cohort. The penalty is chosen by RidgeCV (efficient leave-one-out) on the pooled training
    cohorts only; held-out R^2 is computed around the held-out cohort's mean.
    """
    from sklearn.linear_model import RidgeCV
    rows = []
    coh = list(X_by_cohort)
    for h in coh:
        tr = [c for c in coh if c != h]
        xtr = np.nan_to_num(np.vstack([_std_nan(X_by_cohort[c]) for c in tr]))
        xte = np.nan_to_num(_std_nan(X_by_cohort[h]))
        ytr = np.concatenate([zs(y_by_cohort[c]) for c in tr])
        yte = zs(y_by_cohort[h])
        m = RidgeCV(alphas=alphas).fit(xtr, ytr)
        yhat = m.predict(xte)
        rows.append({"heldout": h, "R2": float(1 - ((yte - yhat) ** 2).sum() / ((yte - yte.mean()) ** 2).sum()),
                     "alpha": float(m.alpha_)})
    return pd.DataFrame(rows)


def loco_one_variable(x_by_cohort: dict, y_by_cohort: dict) -> pd.DataFrame:
    """Leave-one-cohort-out one-variable (linear) prediction, e.g. from the consensus alone."""
    rows = []
    coh = list(x_by_cohort)
    for h in coh:
        tr = [c for c in coh if c != h]
        xtr = np.concatenate([zs(x_by_cohort[c]) for c in tr])
        ytr = np.concatenate([zs(y_by_cohort[c]) for c in tr])
        yte = zs(y_by_cohort[h])
        yhat = np.polyval(np.polyfit(xtr, ytr, 1), zs(x_by_cohort[h]))
        rows.append({"heldout": h, "R2": float(1 - ((yte - yhat) ** 2).sum() / ((yte - yte.mean()) ** 2).sum())})
    return pd.DataFrame(rows)


def _std_nan(M):
    M = np.asarray(M, float)
    return (M - M.mean(axis=0)) / M.std(axis=0)


# ----------------------------------------------------------------------------- driver associations
def copy_number_driver_calls(cn: pd.DataFrame, autosomal: np.ndarray,
                             drivers=(("EGFR", "amp"), ("PDGFRA", "amp"), ("CDK4", "amp"), ("NF1", "loss"),
                                      ("CDKN2A", "homdel"), ("PTEN", "loss"), ("MDM2", "amp"), ("MET", "amp"))):
    """Whole-genome driver calls (CPTAC-3 rule): ploidy = median autosomal gene-level copy number;
    amplification if CN >= 2.5 x ploidy, loss if CN <= ploidy - 1, homozygous deletion if CN == 0.
    ``cn``: genes x samples (absolute gene-level copy number). Returns (calls, ploidy)."""
    ploidy = cn[autosomal].median(axis=0)
    out = {}
    for gene, kind in drivers:
        v = cn.loc[gene]
        v = v.iloc[0] if isinstance(v, pd.DataFrame) else v
        if kind == "amp":
            out[f"{gene}_{kind}"] = (v >= 2.5 * ploidy).astype(float)
        elif kind == "loss":
            out[f"{gene}_{kind}"] = (v <= ploidy - 1).astype(float)
        else:
            out[f"{gene}_{kind}"] = (v == 0).astype(float)
    return pd.DataFrame(out), ploidy


def gistic_driver_calls(gistic: pd.DataFrame) -> pd.DataFrame:
    """TCGA GISTIC2 thresholded calls (patients x genes): amplification >= 2, loss <= -1, homdel <= -2."""
    g = gistic
    out = {"EGFR_amp": g.EGFR >= 2, "PDGFRA_amp": g.PDGFRA >= 2, "CDK4_amp": g.CDK4 >= 2, "NF1_loss": g.NF1 <= -1,
           "CDKN2A_homdel": g.CDKN2A <= -2, "PTEN_loss": g.PTEN <= -1, "MDM2_amp": g.MDM2 >= 2, "MET_amp": g.MET >= 2}
    return pd.DataFrame(out).astype(float).where(g.notna().all(1))


def vaf_purity(mutations: pd.DataFrame, min_vaf: float = 0.05, min_snvs: int = 10) -> pd.Series:
    """DNA purity = 2 x median VAF of somatic SNVs (VAF >= 0.05), capped at 1, if >= 10 SNVs.
    ``mutations`` columns: sample, dna_vaf, ref, alt."""
    m = mutations[(mutations.dna_vaf >= min_vaf) & (mutations.ref.str.len() == 1) & (mutations.alt.str.len() == 1)]
    v = m.groupby("sample").dna_vaf.agg(["median", "size"])
    return (2 * v["median"]).clip(upper=1).where(v["size"] >= min_snvs)


def usable_drivers(calls: pd.DataFrame, drivers=CORE_DRIVERS, min_events: int = 5) -> list[str]:
    """Drivers modelled in a population: >= 5 altered and >= 5 unaltered tumours."""
    return [c for c in drivers if c in calls and calls[c].sum() >= min_events and (calls[c] == 0).sum() >= min_events]


def driver_associated_variance(score, indicators) -> float:
    """Adjusted R^2 of the rank-transformed score on the driver alteration indicators."""
    return adjusted_r2(zs(rank(score)), np.asarray(indicators, float))


def driver_variance_bootstrap(scores: pd.DataFrame, indicators, patients, contrasts=(("MAL_AC", "AXIS"),),
                              nuisance=None, n_boot: int = 2000, seed: int = SEED_DRIVER_BOOTSTRAP,
                              min_events_boot: int = 2, rng=None):
    """Driver-associated variance of each score and differences between scores, with a
    participant-cluster bootstrap.

    Scores are ranked once in the full sample (optionally residualized on ``nuisance``) and standardized;
    resamples draw patients with replacement and are skipped when any modelled driver has fewer than two
    events. Returns dict with point estimates, bootstrap SE (ddof = 1) and 95% percentile intervals, and
    the generator (positioned after the bootstrap, as in production).
    """
    G = np.asarray(indicators, float)
    pat = np.asarray(patients)
    cols = sorted(set(c for pair in contrasts for c in pair), key=list(scores.columns).index)

    def yv(c):
        y = rank(scores[c])
        if nuisance is not None:
            y = residualize(y, np.asarray(nuisance, float))
        return zs(y)
    Y = {c: yv(c) for c in cols}
    point = {c: adjusted_r2(Y[c], G) for c in cols}
    up = np.unique(pat)
    gi = {u: np.where(pat == u)[0] for u in up}
    rng = np.random.default_rng(seed) if rng is None else rng
    bt = []
    for _ in range(n_boot):
        ii = np.concatenate([gi[u] for u in rng.choice(up, len(up))])
        if any(G[ii, j].sum() < min_events_boot for j in range(G.shape[1])):
            continue
        r = {c: adjusted_r2(zs(Y[c][ii]), G[ii]) for c in cols}
        rec = dict(r)
        for a, b in contrasts:
            rec[f"{a}-{b}"] = r[a] - r[b]
        bt.append(rec)
    BT = pd.DataFrame(bt)
    res = {"n": int(len(pat)), "n_patients": int(len(up)), "n_boot_used": int(len(BT)), "R2adj": point,
           "R2adj_CI": {c: np.nanpercentile(BT[c], [2.5, 97.5]).tolist() for c in cols}, "differences": {}}
    for a, b in contrasts:
        k = f"{a}-{b}"
        res["differences"][k] = {"point": point[a] - point[b], "se": float(BT[k].std()),
                                 "CI": np.nanpercentile(BT[k], [2.5, 97.5]).tolist()}
    return res, rng


def purity_adjusted_driver_variance(score, indicators, purity) -> float:
    """Driver-associated variance beyond DNA purity: adjR2(drivers + purity) - adjR2(purity)."""
    y = zs(rank(score))
    pur = np.asarray(purity, float)
    ok = np.isfinite(pur)
    G = np.asarray(indicators, float)
    return adjusted_r2(y[ok], np.column_stack([G[ok], pur[ok]])) - adjusted_r2(y[ok], pur[ok, None])


def matched_random_contrasts(Z: pd.DataFrame, mean_expr: pd.Series, set_a, set_b, excluded, indicators,
                             sample_index, rng, n: int = 1000, tol: float = 0.02) -> np.ndarray:
    """Null of driver-associated variance for 1,000 expression-matched random contrasts.

    Each balance gene is replaced by a random gene whose mean-expression percentile lies within
    +-``tol`` of its own, drawn from genes outside the reference panel and dosage regions (``excluded``).
    ``sample_index`` orders the samples as the rows of ``indicators``.
    """
    G = np.asarray(indicators, float)
    out = []
    pool = [g for g in Z.index if g not in set(excluded)]
    mr = mean_expr.reindex(pool).rank(pct=True)
    pr_all = mean_expr.rank(pct=True)
    for _ in range(n):
        pick = []
        for gset in (set_a, set_b):
            sel = []
            for g in gset:
                if g not in mean_expr.index:
                    continue
                t = pr_all.get(g, np.nan)
                cand = mr.index[(mr - t).abs() <= tol].difference(sel)
                if len(cand):
                    sel.append(rng.choice(cand))
            pick.append(sel)
        x = (Z.loc[pick[0]].mean(axis=0) - Z.loc[pick[1]].mean(axis=0)).reindex(sample_index)
        out.append(adjusted_r2(zs(rank(x)), G))
    return np.asarray(out)


def dosage_region_genes(genes, positions: pd.DataFrame, loci=("EGFR", "PDGFRA", "CDK4", "NF1", "CDKN2A", "PTEN", "MDM2", "MET"),
                        window: int = 10_000_000) -> set:
    """Genes on chromosomes 7 or 10 or within +-10 Mb of any of eight driver loci.
    ``positions``: indexed by gene symbol with columns chrom ('chr7', ...) and chromStart (GENCODE v36)."""
    lp = {g: (positions.loc[g, "chrom"], positions.loc[g, "chromStart"]) for g in loci}
    out = set()
    for g in genes:
        if g not in positions.index:
            continue
        c, s = positions.loc[g, "chrom"], positions.loc[g, "chromStart"]
        if c in ("chr7", "chr10") or any(c == lc and abs(s - ls) <= window for lc, ls in lp.values()):
            out.add(g)
    return out


def fixed_effect_meta(estimates, ses) -> dict:
    """Inverse-variance fixed-effect pooling with Cochran's Q."""
    from scipy.stats import chi2
    e = np.asarray(estimates, float)
    s = np.asarray(ses, float)
    w = 1 / s ** 2
    m = (w * e).sum() / w.sum()
    se = 1 / np.sqrt(w.sum())
    q = float((w * (e - m) ** 2).sum())
    return {"estimate": float(m), "se": float(se), "CI": [float(m - 1.96 * se), float(m + 1.96 * se)], "Q": q,
            "Q_P": float(1 - chi2.cdf(q, len(e) - 1)) if len(e) > 1 else np.nan}


# ----------------------------------------------------------------------------- between-patient share
def icc1(y, groups) -> float:
    """One-way ANOVA ICC(1) for unequal group sizes (k0 correction); negative variance set to zero."""
    y = np.asarray(y, float)
    g = np.asarray(groups)
    ok = np.isfinite(y)
    y, g = y[ok], g[ok]
    ug = np.unique(g)
    n = len(y)
    a = len(ug)
    ni = np.array([(g == u).sum() for u in ug])
    mi = np.array([y[g == u].mean() for u in ug])
    msb = (ni * (mi - y.mean()) ** 2).sum() / (a - 1)
    msw = sum(((y[g == u] - y[g == u].mean()) ** 2).sum() for u in ug) / (n - a)
    k0 = (n - (ni ** 2).sum() / n) / (a - 1)
    s2b = max((msb - msw) / k0, 0)
    return s2b / (s2b + msw)


def icc_contrast(S: pd.DataFrame, groups, a, b) -> float:
    """Mean ICC(1) of coordinates ``a`` minus mean ICC(1) of coordinates ``b``."""
    return float(np.mean([icc1(S[c], groups) for c in a]) - np.mean([icc1(S[c], groups) for c in b]))


def _cluster_resample(g, ug, gi, rng):
    pick = rng.choice(ug, len(ug))
    ii = np.concatenate([gi[u] for u in pick])
    gg = np.concatenate([[str(j)] * len(gi[u]) for j, u in enumerate(pick)])
    return ii, gg


def bootstrap_icc_contrast(S: pd.DataFrame, groups, a, b, n_boot: int = 2000, rng=None, seed: int = SEED_ICC_CONTRAST):
    """Patient-cluster bootstrap of an ICC contrast (resampled patients relabelled as distinct)."""
    rng = np.random.default_rng(seed) if rng is None else rng
    g = np.asarray(groups)
    ug = np.unique(g)
    gi = {u: np.where(g == u)[0] for u in ug}
    bs = []
    for _ in range(n_boot):
        ii, gg = _cluster_resample(g, ug, gi, rng)
        bs.append(icc_contrast(S.iloc[ii].reset_index(drop=True), gg, a, b))
    bs = np.array(bs)
    return {"estimate": icc_contrast(S, g, a, b), "CI": percentile_ci(bs), "boot_frac_le0": float((bs <= 0).mean())}, rng


def bootstrap_icc(y, groups, n_boot: int = 1000, rng=None, seed: int = SEED_ICC_COORDINATE):
    """Patient-cluster bootstrap interval of one coordinate's ICC(1)."""
    rng = np.random.default_rng(seed) if rng is None else rng
    y = np.asarray(y, float)
    g = np.asarray(groups)
    ug = np.unique(g)
    gi = {u: np.where(g == u)[0] for u in ug}
    bs = []
    for _ in range(n_boot):
        ii, gg = _cluster_resample(g, ug, gi, rng)
        bs.append(icc1(y[ii], gg))
    return {"estimate": icc1(y, g), "CI": percentile_ci(bs)}, rng
