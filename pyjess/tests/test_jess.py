import unittest

from .._jess import Template, Molecule, Jess
from .utils import files
from . import data


class TestAtom(unittest.TestCase):

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_query(self):

        with files(data).joinpath("template_01.qry").open() as f:
            template = Template.load(f) 
            jess = Jess([template])

        with files(data).joinpath("pdb1lnb.pdb").open() as f:
            molecule = Molecule.load(f)

        hits = list(jess.query(molecule, 1, 2, 2))
        self.assertEqual(len(hits), 1)
        self.assertAlmostEqual(hits[0].rmsd, 0.555, places=3)
        self.assertAlmostEqual(hits[0].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[0].log_evalue, -2.04, places=1)

        hits = list(jess.query(molecule, 2, 5, 3))
        self.assertEqual(len(hits), 2)
        self.assertAlmostEqual(hits[0].rmsd, 0.555, places=3)
        self.assertAlmostEqual(hits[0].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[0].log_evalue, -2.04, places=1)
        self.assertAlmostEqual(hits[1].rmsd, 1.440, places=3)
        self.assertAlmostEqual(hits[1].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[1].log_evalue, 0.17, places=1)

        hits = list(jess.query(molecule, 2, 5, 5))
        self.assertEqual(len(hits), 3)
        self.assertAlmostEqual(hits[0].rmsd, 0.555, places=3)
        self.assertAlmostEqual(hits[0].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[0].log_evalue, -2.04, places=1)
        self.assertAlmostEqual(hits[1].rmsd, 1.440, places=3)
        self.assertAlmostEqual(hits[1].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[1].log_evalue, 0.17, places=1)
        self.assertAlmostEqual(hits[2].rmsd, 1.644, places=3)
        self.assertAlmostEqual(hits[2].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[2].log_evalue, 0.68, places=1)

        hits = list(jess.query(molecule.conserved(10.0), 1, 2, 2))
        self.assertEqual(len(hits), 1)
        self.assertAlmostEqual(hits[0].rmsd, 0.555, places=3)
        self.assertAlmostEqual(hits[0].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[0].log_evalue, -2.10, places=1)