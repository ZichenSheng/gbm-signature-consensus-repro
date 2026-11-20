"""Numerical and smoke tests.  Run:  python -m unittest discover -s tests -v"""
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import signature_analysis as sa  # noqa: E402
import statistical_validation as sv  # noqa: E402

DATA = ROOT / "data"
# SHA-256 of the frozen production collections (JSON, insertion order, compact separators).
FROZEN_LIT_DIGEST = "754007dce025f7ce27303a62a1540ed416d57c5911378b3b7d4857fd2407ba1e"
FROZEN_N1_DIGESTS = {0: "f1b28e0e39f9ee72b51edb428d90504c30af38629cce2353b579a2918b3fcc18",
                     1: "d2b91c067e53713098f5c2553c289975c3f60dd3a6e2c0b926823324ebbee26e"}


def digest(obj):
    return hashlib.sha256(json.dumps(obj, separators=(",", ":")).encode()).hexdigest()


class DataFiles(unittest.TestCase):
    def test_checksums(self):
        for line in open(DATA / "SHA256SUMS"):
            h, name = line.split()
            self.assertEqual(hashlib.sha256(open(DATA / name, "rb").read()).hexdigest(), h, name)

    def test_registry_counts(self):
        ids, sigs = sa.load_signatures(DATA / "signatures_201.tsv")
        self.assertEqual(len(ids), 201)
        self.assertEqual(len({g for v in sigs.values() for g in v}), 1277)
        src = pd.read_csv(DATA / "signature_sources.tsv", sep="\t")
        self.assertEqual(src.doi.str.lower().nunique(), 197)
        self.assertEqual(len(sa.load_formulas(DATA / "weighted_formulas_104.tsv")), 104)
        L = pd.read_csv(DATA / "consensus_loadings.tsv", sep="\t")
        self.assertEqual(len(L), 201)
        self.assertAlmostEqual(float((L.loading ** 2).sum()), 1.0, places=12)

    def test_reference_sets(self):
        ref = json.load(open(DATA / "reference_gene_sets.json"))
        C = ref["coordinates"]
        self.assertEqual(sorted(C), sorted(sa.COORDINATES))
        self.assertEqual({k: (len(v["A"]), len(v["B"])) for k, v in C.items()}["MAL_AC"], (20, 20))
        self.assertEqual(C["HYPOX"]["A"], sa.HYPOXIA_GENES)
        self.assertEqual(C["MAL_AC"]["B"], C["MAL_MES"]["B"])          # shared developmental set
        self.assertEqual(ref["pan_malignant_reference"], C["MYE"]["B"])


class Normalization(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.genes = [f"G{i}" for i in range(400)]
        self.A = pd.DataFrame(rng.poisson(np.exp(rng.normal(3, 1, (400, 1))) * np.ones((1, 12))).astype(float),
                              index=self.genes, columns=[f"s{j}" for j in range(12)])

    def test_size_factors(self):
        B = self.A.copy()
        B["s3"] *= 3.0                                     # scaling one library scales its factor ratio exactly
        sf1 = sa.median_of_ratios_size_factors(self.A)
        sf2 = sa.median_of_ratios_size_factors(B)
        self.assertAlmostEqual(float(sf2["s3"] / sf2["s0"]), 3.0 * float(sf1["s3"] / sf1["s0"]), places=10)
        self.assertAlmostEqual(float(np.exp(np.log(sf1).mean())), 1.0, places=12)
        X = sa.normalize_rnaseq(self.A, self.genes)
        np.testing.assert_allclose(X.to_numpy(), np.log2(self.A.loc[X.index] / sf1 + 1).to_numpy(), atol=1e-12)

    def test_zscore(self):
        Z = sa.zscore_genes(np.log2(self.A + 1))
        np.testing.assert_allclose(Z.mean(axis=1), 0, atol=1e-12)
        np.testing.assert_allclose(Z.std(axis=1), 1, atol=1e-12)


class Scores(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        self.Z = sa.zscore_genes(pd.DataFrame(rng.normal(size=(60, 30)), index=[f"G{i}" for i in range(60)]))

    def test_membership_coverage_rule(self):
        sigs = {"a": ["G1", "G2", "G3"], "b": ["G1", "G2", "X1", "X2", "X3"], "c": ["G1", "X1"]}
        S, cov = sa.membership_scores(self.Z, sigs)
        np.testing.assert_allclose(S["a"], self.Z.loc[["G1", "G2", "G3"]].mean(axis=0))
        self.assertTrue(S["b"].isna().all())                             # 2 of 5 genes: < 50%
        self.assertTrue(S["c"].isna().all())                             # fewer than 3 genes

    def test_consensus_projection(self):
        sigs = {f"s{i}": [f"G{3 * i}", f"G{3 * i + 1}", f"G{3 * i + 2}"] for i in range(20)}
        u = pd.Series(np.linspace(-1, 1, 20), index=list(sigs))
        ax, used = sa.consensus_score(self.Z, sigs, u)
        S, _ = sa.membership_scores(self.Z, sigs)
        A = S.to_numpy()
        manual = ((A - A.mean(0)) / A.std(0)) @ u.to_numpy()
        np.testing.assert_allclose(ax, (manual - manual.mean()) / manual.std(), atol=1e-12)
        self.assertEqual(len(used), 20)

    def test_contrast_invariant_to_sample_shift(self):
        ref = json.load(open(DATA / "reference_gene_sets.json"))["coordinates"]
        genes = sorted({g for c in ref.values() for k in ("A", "B") for g in c[k]})
        rng = np.random.default_rng(2)
        Z = sa.zscore_genes(pd.DataFrame(rng.normal(size=(len(genes) + 300, 25)), index=genes + [f"F{i}" for i in range(300)]))
        me = pd.Series(rng.normal(size=len(Z)), index=Z.index)
        C1, _, _ = sa.score_reference_coordinates(Z, ref, me)
        C2, _, _ = sa.score_reference_coordinates(Z + rng.normal(0, 2, Z.shape[1]), ref, me)
        np.testing.assert_allclose(C1.to_numpy(), C2.to_numpy(), atol=1e-10)   # zero-sum contrasts


class Statistics(unittest.TestCase):
    def test_icc_balanced_closed_form(self):
        rng = np.random.default_rng(3)
        g = np.repeat(np.arange(8), 4)
        y = rng.normal(size=8)[g] + rng.normal(0, 0.7, 32)
        a, k = 8, 4
        mi = np.array([y[g == u].mean() for u in range(a)])
        msb = k * ((mi - y.mean()) ** 2).sum() / (a - 1)
        msw = sum(((y[g == u] - mi[u]) ** 2).sum() for u in range(a)) / (a * (k - 1))
        expect = max((msb - msw) / k, 0) / (max((msb - msw) / k, 0) + msw)
        self.assertAlmostEqual(sv.icc1(y, g), expect, places=12)
        self.assertEqual(sv.icc1(np.tile([1.0, 2.0, 3.0], 5), np.repeat(np.arange(5), 3)), 0.0)

    def test_adjusted_r2(self):
        rng = np.random.default_rng(4)
        X = rng.normal(size=(50, 3))
        y = X @ [1, 0.5, 0] + rng.normal(size=50)
        r2 = sv.r2(y, X)
        self.assertAlmostEqual(sv.adjusted_r2(y, X), 1 - (1 - r2) * 49 / 46, places=12)

    def test_fixed_effect_meta(self):
        m = sv.fixed_effect_meta([0.2, 0.1], [0.1, 0.05])
        w = np.array([100.0, 400.0])
        est = (w * [0.2, 0.1]).sum() / w.sum()
        self.assertAlmostEqual(m["estimate"], est, places=12)
        self.assertAlmostEqual(m["se"], 1 / np.sqrt(w.sum()), places=12)

    def test_bootstrap_reproducible(self):
        rng = np.random.default_rng(5)
        D = pd.DataFrame({"MAL_AC": rng.normal(size=60), "AXIS": rng.normal(size=60)})
        G = (rng.random((60, 2)) < 0.3).astype(float)
        D["MAL_AC"] += 2 * G[:, 0]
        r1, _ = sv.driver_variance_bootstrap(D, G, np.arange(60), n_boot=100, seed=7)
        r2, _ = sv.driver_variance_bootstrap(D, G, np.arange(60), n_boot=100, seed=7)
        self.assertEqual(r1, r2)
        self.assertGreater(r1["R2adj"]["MAL_AC"], r1["R2adj"]["AXIS"])

    def test_structure_reproducibility_identity(self):
        rng = np.random.default_rng(6)
        S = rng.normal(size=(40, 10))
        self.assertAlmostEqual(sv.structure_reproducibility(S, S), 1.0, places=12)


class FrozenCollections(unittest.TestCase):
    """The expression-matched collections need only the supplied pool features, so a reviewer can
    regenerate them exactly; digests are those of the frozen production collections."""

    def test_n1_regenerated_exactly(self):
        ids, sigs = sa.load_signatures(DATA / "signatures_201.tsv")
        pf = pd.read_csv(DATA / "matching_pool_features.tsv", sep="\t").set_index("gene_symbol")
        lit = sv.literature_sets_on_pool(sigs, ids, list(pf.index))
        self.assertEqual(digest(lit), FROZEN_LIT_DIGEST)
        F = sv.standardized_matching_features(pf)
        for i, h in FROZEN_N1_DIGESTS.items():
            n1 = sv.expression_matched_collection(lit, F, np.random.default_rng(sv.SEED_COLLECTIONS + 1000 + i))
            self.assertEqual(digest(n1), h)


class Example(unittest.TestCase):
    def test_example_runs(self):
        import run_example
        with tempfile.TemporaryDirectory() as d:
            s = run_example.main(Path(d), n_boot=50)
            self.assertTrue((Path(d) / "example_summary.json").exists())
        self.assertEqual(s["SYN_A"]["membership_scores_computable"], 201)
        self.assertTrue(np.isfinite(s["driver_associated_variance"]["R2adj"]["MAL_AC"]))
        self.assertIn("SYNTHETIC", s["WARNING"])


if __name__ == "__main__":
    unittest.main()
