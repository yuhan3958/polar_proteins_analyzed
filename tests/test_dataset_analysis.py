import unittest

from dataset_analysis import compare_sequences, sequence_features, group_summaries


class DatasetAnalysisTests(unittest.TestCase):
    def test_gap_position_mapping_and_denominators(self):
        metrics, variants, windows, columns, changes = compare_sequences('AAC', 'AC', 50)
        self.assertEqual(metrics['identical_sites'], 2)
        self.assertEqual(metrics['gap_columns'], 1)
        self.assertEqual(metrics['aligned_non_gap_sites'], 2)
        self.assertEqual(metrics['identity_on_comparable_sites_percent'], 100)
        self.assertAlmostEqual(metrics['identity_on_all_alignment_columns_percent'], 200 / 3)
        self.assertEqual(variants[0]['b_position'], '')
        self.assertEqual(''.join(a for a, _, _ in columns).replace('-', ''), 'AAC')
        self.assertEqual(''.join(b for _, b, _ in columns).replace('-', ''), 'AC')

    def test_gpr_related_sites_count_once(self):
        metrics, variants, _, _, changes = compare_sequences('AGPRA', 'APGKA')
        self.assertEqual(metrics['GPR_substitution_sites'], 3)
        self.assertEqual(sum(changes.values()), 3)
        self.assertEqual(sum(v['focus_GPR'] for v in variants), 3)
        self.assertEqual(metrics['positive_score_substitutions'] + metrics['nonpositive_score_substitutions'], metrics['substitutions'])

    def test_per_protein_mean_is_distinct_from_pooled_residue_fraction(self):
        a = sequence_features({'sequence': 'AG', 'accession': 'a', 'group': 'cold'})
        b = sequence_features({'sequence': 'AAAA', 'accession': 'b', 'group': 'cold'})
        c = sequence_features({'sequence': 'GG', 'accession': 'c', 'group': 'meso'})
        self.assertEqual(sum(a[aa + '_count'] for aa in 'ACDEFGHIKLMNPQRSTVWY'), 2)
        row = next(row for row in group_summaries([a, b, c], 'test') if row['amino_acid'] == 'G')
        self.assertEqual(row['cold_mean_percent'], 25)
        self.assertAlmostEqual(row['cold_pooled_residue_percent'], 100 / 6)
        self.assertEqual(row['difference_pp_cold_minus_meso'], -75)


if __name__ == '__main__':
    unittest.main()
