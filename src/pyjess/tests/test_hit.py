import math
import unittest
import sys
import pickle
import textwrap

from .._jess import Template, Molecule, Jess
from .utils import files
from . import data


class TestHit(unittest.TestCase):
    
    maxDiff = None

    @unittest.skipUnless(files, "importlib.resources not available")
    @classmethod
    def setUpClass(cls):
        with files(data).joinpath("template_01.qry").open() as f:
            template = Template.load(f, id="T1")
            jess = Jess([template])
        with files(data).joinpath("pdb1lnb.pdb").open() as f:
            molecule = Molecule.load(f)

        cls.hit = next(jess.query(molecule, 1, 2, 2))
        
    def test_pickle(self):
        hit = pickle.loads(pickle.dumps(self.hit))
        self.assertIsNot(self.hit, hit)
        self.assertEqual(self.hit.rmsd, hit.rmsd)
        self.assertEqual(self.hit.determinant, hit.determinant)
        self.assertEqual(self.hit.evalue, hit.evalue)
        self.assertEqual(self.hit.template(transform=False), hit.template(transform=False))
        self.assertListEqual(self.hit.atoms(transform=True), hit.atoms(transform=True))
        self.assertListEqual(self.hit.atoms(transform=False), hit.atoms(transform=False))
        self.assertEqual(self.hit.molecule(transform=False), hit.molecule(transform=False))

    def test_dumps(self):
        expected = textwrap.dedent(
            """
            REMARK 1LNB 0.555 T1 Det= 1.0 log(E)~ -2.04
            ATOM   1112  ND1 HIS E 142       4.157  -2.574   4.091  1.00 10.22           N
            ATOM   1115  NE2 HIS E 142       3.805  -2.316   1.890  1.00 16.00           N
            ATOM   1123  OE1 GLU E 143       0.480  -4.544  -0.363  1.00 14.72           O
            ATOM   1124  OE2 GLU E 143       1.870  -6.097   0.337  1.00 22.48           O
            ATOM   1145  CG  HIS E 146       0.031   0.017   0.019  1.00 12.45           C
            ATOM   1146  ND1 HIS E 146       0.823   1.160   0.000  1.00 15.79           N
            ATOM   1147  CD2 HIS E 146       0.878  -1.064   0.006  1.00 20.85           C
            ATOM   1148  CE1 HIS E 146       2.115   0.777  -0.013  1.00 13.50           C
            ATOM   1149  NE2 HIS E 146       2.181  -0.553  -0.032  1.00 12.44           N
            ATOM   1298  OE1 GLU E 166       6.432  -0.605   1.056  1.00 16.03           O
            ATOM   1299  OE2 GLU E 166       4.839   0.070  -0.343  1.00 19.31           O
            ENDMDL
            """
        ).strip()
        actual = self.hit.dumps().strip()
        self.assertMultiLineEqual(actual, expected)