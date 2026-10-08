import csv
import itertools
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from extract_unitigs_5 import iter_fasta
from fdr_5 import benjamini_hochberg, storey_qvalues
from unitig_fisher_5 import fisher_two_sided


class FisherTests(unittest.TestCase):
    def test_reference_tables(self):
        # Two-sided reference examples from scipy.stats.fisher_exact.
        for table, expected in [
            ((6, 2, 1, 4), 0.10256410256410257),
            ((8, 2, 1, 5), 0.034965034965034975),
            ((5, 0, 0, 5), 1 / 126),
        ]:
            with self.subTest(table=table):
                self.assertTrue(math.isclose(
                    fisher_two_sided(*table), expected, rel_tol=1e-12
                ))

    def test_extreme_tables_preserve_tiny_probabilities(self):
        # Only the two perfectly separated tables are as improbable.
        for n in (10, 30, 50, 100):
            with self.subTest(n=n):
                expected = 2 / math.comb(2 * n, n)
                actual = fisher_two_sided(n, 0, 0, n)
                self.assertTrue(math.isclose(actual, expected, rel_tol=1e-12),
                                f"{actual} != {expected}")

    def test_degenerate_tables(self):
        for table in ((0, 0, 0, 0), (0, 0, 3, 7), (5, 0, 8, 0)):
            with self.subTest(table=table):
                self.assertEqual(fisher_two_sided(*table), 1.0)

    def test_small_tables_match_exact_label_permutations(self):
        # Enumerate assignments of carriers to labelled individuals, without
        # using the production function's binomial-coefficient calculation.
        for n1, n2 in ((3, 4), (4, 4)):
            for carriers in range(n1 + n2 + 1):
                frequencies = {}
                for assignment in itertools.combinations(range(n1 + n2), carriers):
                    a = sum(i < n1 for i in assignment)
                    frequencies[a] = frequencies.get(a, 0) + 1
                for a, frequency in frequencies.items():
                    c = carriers - a
                    table = (a, n1 - a, c, n2 - c)
                    expected = sum(v for v in frequencies.values() if v <= frequency)
                    expected /= sum(frequencies.values())
                    with self.subTest(table=table):
                        self.assertTrue(math.isclose(
                            fisher_two_sided(*table), expected, rel_tol=1e-12
                        ))


class FDRTests(unittest.TestCase):
    def test_bh_preserves_order_and_handles_ties(self):
        actual = benjamini_hochberg([0.04, 0.01, 0.01, 0.2])
        for value, expected in zip(actual, [0.04 * 4 / 3, 0.02, 0.02, 0.2]):
            self.assertAlmostEqual(value, expected)

    def test_empty_inputs(self):
        self.assertEqual(benjamini_hochberg([]), [])
        self.assertEqual(storey_qvalues([]), ([], 1.0))

    def test_storey_with_positive_pi0(self):
        values, pi0 = storey_qvalues([0.01, 0.02, 0.2, 0.3, 0.9], [0.5])
        self.assertAlmostEqual(pi0, 0.4)
        for value, expected in zip(values, [0.02, 0.02, 2 / 15, 0.15, 0.36]):
            self.assertAlmostEqual(value, expected)

    def test_storey_uses_bh_when_pi0_estimate_is_zero(self):
        for p_values, expected in [
            ([0.1, 0.2, 0.3], [0.3, 0.3, 0.3]),
            ([0.1, 0.2, 0.6], [0.3, 0.3, 0.6]),
            ([0.0, 0.1, 0.3], [0.0, 0.15, 0.3]),
        ]:
            with self.subTest(p_values=p_values):
                values, pi0 = storey_qvalues(p_values)
                self.assertEqual(pi0, 1.0)
                for value, reference in zip(values, expected):
                    self.assertAlmostEqual(value, reference)


class SyntheticPassVTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="bituga-stats-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.results = self.root / "results"
        self.results.mkdir()
        self.metadata = self.root / "metadata.tsv"
        # Two read files per individual must still count as one individual.
        self.metadata.write_text("path\tsample_id\ttrait\n" + "".join(
            f"/unused/{trait}{i}_{mate}.fq\t{trait}{i}\t{trait}\n"
            for trait in ("A", "B") for i in range(5) for mate in (1, 2)
        ))

    def run_pass(self, counts, mode):
        sequences = {uid: "ACGT" * (i + 8) for i, (uid, _, _) in enumerate(counts)}
        (self.results / "unitigs.fa").write_text("".join(
            f">{uid}\n{seq}\n" for uid, seq in sequences.items()
        ))
        (self.results / "unitig_presence.tsv").write_text(
            "unitig\tA\tB\n" + "".join(f"{uid}\t{a}\t{b}\n" for uid, a, b in counts)
        )
        proc = subprocess.run(
            [sys.executable, "-B", str(ROOT / "bin/orchestrator_5.py"),
             "--outdir", str(self.root), "--metadata", str(self.metadata),
             "--fdr", mode, "--split-by-bias"],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        stats = self.results / "statistics"
        summary = dict(line.split("\t", 1) for line in
                       (stats / "unitig_association_stats.txt").read_text().splitlines()
                       if "\t" in line)
        with (stats / "unitig_association_details.tsv").open() as fh:
            rows = {r["unitig"]: r for r in csv.DictReader(
                itertools.takewhile(lambda line: not line.startswith("#"), fh),
                delimiter="\t",
            )}
        exported = {h[1:].split()[0]: seq for h, seq in
                    iter_fasta(self.results / "significant_unitigs.fa")}
        for uid, seq in exported.items():
            self.assertEqual(seq, sequences[uid])
        return summary, rows, exported

    def test_bh_counts_probabilities_and_fasta_export(self):
        summary, rows, exported = self.run_pass([
            ("a_only", 5, 0), ("b_only", 0, 5), ("weak", 4, 1),
            ("neutral", 2, 2), ("absent", 0, 0),
        ], "bh")
        self.assertEqual(summary["tested_unitigs"], "4")
        self.assertEqual(summary["significant_unitigs"], "2")
        self.assertEqual(summary["multiple_testing_correction"], "BH")
        self.assertEqual(set(exported), {"a_only", "b_only"})
        for uid, trait in (("a_only", "A"), ("b_only", "B")):
            self.assertAlmostEqual(float(rows[uid]["p_value"]), 1 / 126, places=7)
            self.assertAlmostEqual(float(rows[uid]["p_adj"]), 1 / 63, places=7)
            self.assertEqual(rows[uid]["bias"], trait)
            self.assertEqual(float(rows[uid][f"prev_{trait}"]), 1.0)
        self.assertEqual(float(rows["neutral"]["p_value"]), 1.0)
        self.assertEqual(rows["absent"]["p_value"], "NA")

    def test_trait_specific_fasta_contains_hits(self):
        _, _, exported = self.run_pass([("a_only", 5, 0), ("b_only", 0, 5)], "bh")
        for uid, trait in (("a_only", "A"), ("b_only", "B")):
            split = list(iter_fasta(self.results / f"significant_unitigs_{trait}.fa"))
            self.assertEqual([h[1:].split()[0] for h, _ in split], [uid])
            self.assertEqual(split[0][1], exported[uid])

    def test_trait_specific_fasta_leaves_other_group_empty(self):
        _, _, exported = self.run_pass([("a_only", 5, 0), ("neutral", 2, 2)], "bh")
        split = list(iter_fasta(self.results / "significant_unitigs_A.fa"))
        self.assertEqual([h[1:].split()[0] for h, _ in split], ["a_only"])
        self.assertEqual(split[0][1], exported["a_only"])
        self.assertEqual((self.results / "significant_unitigs_B.fa").read_text(), "")

    def test_storey_does_not_export_weak_signals(self):
        summary, rows, exported = self.run_pass(
            [("weak_a", 4, 1), ("weak_b", 1, 4)], "storey"
        )
        self.assertEqual(summary["multiple_testing_correction"], "storey")
        self.assertEqual(float(summary["storey_pi0"]), 1.0)
        self.assertEqual(summary["significant_unitigs"], "0")
        for uid in ("weak_a", "weak_b"):
            self.assertAlmostEqual(float(rows[uid]["q_value"]), 13 / 63, places=6)
        self.assertEqual(exported, {})

    def sparse_signal(self):
        return [("signal", 5, 0)] + [(f"null{i}", 2, 2) for i in range(99)]

    def test_bh_rejects_sparse_signal_after_correction(self):
        summary, rows, exported = self.run_pass(self.sparse_signal(), "bh")
        self.assertEqual(summary["significant_unitigs"], "0")
        self.assertAlmostEqual(float(rows["signal"]["p_adj"]), 100 / 126, places=6)
        self.assertEqual(exported, {})

    def test_explicit_uncorrected_mode_exports_sparse_signal(self):
        summary, _, exported = self.run_pass(self.sparse_signal(), "none")
        self.assertEqual(summary["multiple_testing_correction"], "none")
        self.assertEqual(set(exported), {"signal"})

    def test_auto_preserves_uncorrected_fallback(self):
        summary, _, exported = self.run_pass(self.sparse_signal(), "auto")
        self.assertEqual(summary["multiple_testing_correction"], "none")
        self.assertEqual(set(exported), {"signal"})

    def test_auto_preserves_bh_when_feasible(self):
        summary, rows, exported = self.run_pass([("signal", 5, 0)], "auto")
        self.assertEqual(summary["multiple_testing_correction"], "BH")
        self.assertAlmostEqual(float(rows["signal"]["p_adj"]), 1 / 126, places=7)
        self.assertEqual(set(exported), {"signal"})


if __name__ == "__main__":
    unittest.main()
