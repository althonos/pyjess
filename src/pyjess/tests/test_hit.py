import math
import unittest
import sys
import pickle

from .._jess import Template, Molecule, Jess
from .utils import files
from . import data


class TestHit(unittest.TestCase):

    @unittest.skipUnless(files, "importlib.resources not available")
    @classmethod
    def setUpClass(cls):
        with files(data).joinpath("template_01.qry").open() as f:
            template = Template.load(f)
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
        self.assertEqual(self.hit.template, hit.template)
        self.assertListEqual(self.hit.atoms(transform=True), hit.atoms(transform=True))
        self.assertListEqual(self.hit.atoms(transform=False), hit.atoms(transform=False))
        self.assertEqual(self.hit.molecule(transform=False), hit.molecule(transform=False))