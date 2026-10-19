import contextlib
import gzip
import io
import os
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from make_metadata import build_trait_mapping
from orchestrator_4 import compute_length_params_heuristic


class TraitTests(unittest.TestCase):
    def mapping(self, labels):
        with contextlib.redirect_stderr(io.StringIO()):
            return build_trait_mapping([(str(i), label) for i, label in enumerate(labels)])

    def test_numeric_labels_remain_distinct(self):
        self.assertEqual(self.mapping(["0", "1"]), {"0": "0", "1": "1"})

    def test_letter_labels_remain_distinct(self):
        self.assertEqual(self.mapping(["A", "B"]), {"A": "A", "B": "B"})

    def test_anagram_labels_remain_distinct(self):
        self.assertEqual(self.mapping(["east", "seat"]), {"east": "east", "seat": "seat"})

    def test_sex_aliases_still_normalize(self):
        self.assertEqual(self.mapping(["M", "male", "F", "female"]),
                         {"M": "male", "male": "male", "F": "female", "female": "female"})


class PoolThreadTests(unittest.TestCase):
    def test_kmc_tools_union_receives_requested_threads(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            bindir = work / "bin"
            bindir.mkdir()
            args_log = work / "kmc_tools.args"
            fake_kmc_tools = bindir / "kmc_tools"
            fake_kmc_tools.write_text(
                "#!/usr/bin/env bash\n"
                "printf '%s\\n' \"$@\" > \"$KMC_TOOLS_ARGS_LOG\"\n"
                "for output in \"$@\"; do :; done\n"
                ": > \"${output}.kmc_pre\"\n"
                ": > \"${output}.kmc_suf\"\n"
            )
            fake_kmc_tools.chmod(0o755)
            sample = work / "sample"
            results = work / "results"
            pool = results / "pooled_kmc_files" / "pooled_case"
            pool.parent.mkdir(parents=True)
            for base in (sample, pool):
                Path(f"{base}.kmc_pre").touch()
                Path(f"{base}.kmc_suf").touch()
            env = os.environ.copy()
            env["KMC_TOOLS_ARGS_LOG"] = str(args_log)
            result = subprocess.run(
                [str(ROOT / "bin/update_pool_1.sh"), "--sample-db", str(sample),
                 "--trait", "case", "--outdir", str(results),
                 "--kmc-bin", str(bindir), "--threads", "7"],
                capture_output=True, text=True, env=env, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(args_log.read_text().splitlines()[:2], ["-t7", "simple"])
            self.assertIn('--threads "${THREADS}"',
                          (ROOT / "bin/build_groups_1.sh").read_text())


class NativeInputTests(unittest.TestCase):
    @unittest.skipUnless((ROOT / "bin/nt_mini_matcher.x").exists(), "Build MiniMatcher first")
    def test_formats_threads_and_exact_batch_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            rng = random.Random(42)
            seq = "".join(rng.choices("ACGT", k=300))
            ref = work / "ref.fa"
            ref.write_text(f">reference\n{seq}\n")
            fasta, fastq, zipped = (work / name for name in ("reads.fa", "reads.fq", "reads.fq.gz"))
            reads = [seq[start:start + 150] for start in (0, 50, 100, 150)]
            fasta.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(reads)))
            fastq.write_text("".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(reads)))
            with gzip.open(zipped, "wt") as handle:
                handle.write(fastq.read_text())
            expected = None
            for query in (fasta, fastq, zipped):
                for threads in (1, 2):
                    with self.subTest(query=query.name, threads=threads):
                        # Four reads exactly fill two batches; EOF must not
                        # yield a third, empty batch to the hashing code.
                        command = [str(ROOT / "bin/nt_mini_matcher.x"),
                                   "--threads", str(threads), "--query_split_size", "2",
                                   "--kmer_length", "19", "--minimum_mem_length", "32",
                                   "--reverse_complement", "for_reference", str(ref), str(query)]
                        result = subprocess.run(command, cwd=work, capture_output=True,
                                                text=True, timeout=30)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        rows = sorted(line for line in result.stdout.splitlines()
                                      if line.strip() and not line.startswith("#"))
                        self.assertTrue(rows)
                        if expected is None:
                            expected = rows
                        else:
                            self.assertEqual(rows, expected)

    @unittest.skipUnless((ROOT / "bin/postprocess_minimatcher").exists(), "Build postprocessing first")
    def test_short_reads_cover_long_unitig_but_partial_coverage_fails(self):
        params = compute_length_params_heuristic({900: 2}, 31, 900, 900)
        option_names = ("short-max", "medium-max", "short-min-bp", "medium-min-bp",
                        "long-min-bp", "short-min-segment", "medium-min-segment", "long-min-segment")
        # Total unitig-coverage requirements must not be lowered by the fix.
        self.assertEqual(params[2:5], (810, 811, 812))
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            unitigs = work / "unitigs.fa"
            unitigs.write_text(">covered\n" + "ACGT" * 225 + "\n")
            for last_start, expected_count in ((750, 1), (600, 0)):
                with self.subTest(last_start=last_start):
                    # Independent 150-bp reads tile a 900-bp unitig; the second
                    # case omits the final 150 bp and must still be rejected.
                    matches = "".join(f"covered\t{start}\t150\tread{i}\n"
                                      for i, start in enumerate(range(0, last_start + 1, 150)))
                    output = work / "present.bin"
                    command = [str(ROOT / "bin/postprocess_minimatcher"), str(unitigs), "-", str(output)]
                    for name, value in zip(option_names, params):
                        command += [f"--{name}", str(value)]
                    result = subprocess.run(command, input=matches, capture_output=True,
                                            text=True, timeout=10)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    count = struct.unpack("=I", output.read_bytes()[:4])[0]
                    self.assertEqual(count, expected_count)
