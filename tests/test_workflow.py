import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "run_workflow.sh"


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.genome = self.root / "genome"
        self.jum = self.root / "jum"
        self.genome.mkdir()
        self.jum.mkdir()
        for name in ("JUM_A.sh", "JUM_B.sh", "JUM_C.sh", "R_script_JUM.R"):
            (self.jum / name).touch()
        self.ref_flat = self.root / "refFlat.txt"
        self.design = self.root / "design.txt"
        self.ref_flat.touch()
        self.design.touch()

    def tearDown(self):
        self.tempdir.cleanup()

    def write_sheet(self, rows):
        sheet = self.root / "samples.tsv"
        text = "sample\tcondition\tread1\tread2\n"
        for sample, condition in rows:
            read1 = self.root / f"{sample}_R1.fastq.gz"
            read2 = self.root / f"{sample}_R2.fastq.gz"
            read1.touch()
            read2.touch()
            text += f"{sample}\t{condition}\t{read1}\t{read2}\n"
        sheet.write_text(text, encoding="utf-8")
        return sheet

    def run_script(self, sheet):
        output = self.root / "output"
        result = subprocess.run(
            [
                "bash", str(SCRIPT),
                "--samples", str(sheet),
                "--genome-dir", str(self.genome),
                "--ref-flat", str(self.ref_flat),
                "--jum-dir", str(self.jum),
                "--design", str(self.design),
                "--output-dir", str(output),
            ],
            text=True,
            capture_output=True,
        )
        return result, output

    def test_dry_run_builds_two_pass_plan_without_writing_output(self):
        sheet = self.write_sheet([
            ("control_a", "control"),
            ("control_b", "control"),
            ("treated_a", "treated"),
            ("treated_b", "treated"),
        ])
        result, output = self.run_script(sheet)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("STAR", result.stdout)
        self.assertIn("samtools", result.stdout)
        self.assertIn("JUM_A.sh", result.stdout)
        self.assertIn("JUM_B.sh", result.stdout)
        self.assertIn("JUM_C.sh", result.stdout)
        self.assertIn(r"control_a_\,control_b_", result.stdout)
        self.assertIn(r"treated_a_\,treated_b_", result.stdout)
        self.assertIn("Dry run only", result.stdout)
        self.assertFalse(output.exists())

    def test_rejects_more_than_two_conditions(self):
        sheet = self.write_sheet([
            ("sample_a", "first"),
            ("sample_b", "second"),
            ("sample_c", "third"),
        ])
        result, _output = self.run_script(sheet)
        self.assertEqual(result.returncode, 2)
        self.assertIn("exactly two conditions", result.stderr)

    def test_rejects_duplicate_samples(self):
        sheet = self.write_sheet([
            ("sample_a", "control"),
            ("sample_a", "treated"),
        ])
        result, _output = self.run_script(sheet)
        self.assertEqual(result.returncode, 2)
        self.assertIn("duplicate sample", result.stderr)


if __name__ == "__main__":
    unittest.main()
