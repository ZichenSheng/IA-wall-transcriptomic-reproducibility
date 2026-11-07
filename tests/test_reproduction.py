"""Clean run and fail-closed checks for data substitutions."""
from pathlib import Path
import shutil
import tempfile
import unittest

from scripts.common import ROOT, validate_inputs
from scripts.verify_source_data import run


class ReproductionTests(unittest.TestCase):
    def test_complete_standalone_run(self):
        result = run()
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['fig2_level_B']['n_claims'], 239)
        self.assertEqual(result['fig2_level_B']['sign_reversals'], 43)
        self.assertEqual(result['frozen_checks_level_C']['strict_claims'], 58)
        self.assertGreater(result['fig3']['drop']['archived_draws_level_C']['p'], 0.05)

    def test_reordered_claims_are_rejected(self):
        self.reject_corruption('data/q1/claims.tsv', self.reorder)

    def test_cohort_substitution_is_rejected(self):
        self.reject_corruption('data/rupture_counts.tsv', lambda s: s.replace('GSE54083', 'GSE75436'))

    def test_mapping_substitution_is_rejected(self):
        self.reject_corruption('data/q1/claims.tsv', lambda s: s.replace('HGNC:29', 'HGNC:92', 1))

    def reject_corruption(self, rel, mutate):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            shutil.copytree(ROOT / 'data', root / 'data')
            shutil.copytree(ROOT / 'expected', root / 'expected')
            file = root / rel
            source = file.read_text()
            modified = mutate(source)
            self.assertNotEqual(source, modified)
            file.write_text(modified)
            with self.assertRaisesRegex(ValueError, 'Input identity failed'):
                validate_inputs(root)

    @staticmethod
    def reorder(source):
        header, first, second, rest = source.split('\n', 3)
        return '\n'.join((header, second, first, rest))


if __name__ == '__main__':
    unittest.main()
