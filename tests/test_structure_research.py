import unittest
from scripts.structure_distances import map_positions, pair_distances
from scripts.tree_support_audit import audit

class StructureResearchTests(unittest.TestCase):
    def test_alignment_self(self):
        mapping, identity=map_positions("ACDEFGHIK","ACDEFGHIK")
        self.assertEqual(mapping[5],5)
        self.assertEqual(identity,1.0)
    def test_geometry(self):
        left={"atoms":{"CA":(0.,0.,0.),"CB":(1.,0.,0.)}}
        right={"atoms":{"CA":(3.,0.,0.),"CB":(2.,0.,0.)}}
        ca,heavy=pair_distances(left,right)
        self.assertEqual((ca,heavy),(3.,1.))

if __name__=="__main__":
    unittest.main()
