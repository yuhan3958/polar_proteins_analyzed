import unittest
from scripts.extended_research import aligner, map_reference_sites, comparable_identity, audit_thermal_metadata

class ExtendedResearchTests(unittest.TestCase):
    def test_identity(self):
        e = aligner()
        self.assertAlmostEqual(comparable_identity("MABCDEFG", "MABCDEFG", e), 1.)
    def test_active_site_self_mapping(self):
        seq = "A" * 600
        rows = map_reference_sites(seq, seq, aligner())
        self.assertEqual(len(rows), 5)
        self.assertTrue(all(row["conserved"] and row["target_position"] == row["reference_position"]
                            for row in rows))
    def test_thermal_distinction(self):
        rows = audit_thermal_metadata([{
            "accession":"TEST", "group":"cold", "organism":"Alpha species",
            "thermal_class":"psychrotolerant", "growth_optimum_c":"",
            "culture_temperature_c":"20", "enzyme_activity_optimum_c":"45",
            "temperature_source_url":"https://example.org/source",
        }])
        self.assertFalse(rows[0]["growth_optimum_recorded"])
        self.assertTrue(rows[0]["source_link_recorded"])
        self.assertTrue(rows[0]["cold_tolerant_label"])

if __name__ == "__main__":
    unittest.main()
