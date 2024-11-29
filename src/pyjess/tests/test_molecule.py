import os
import pickle
import unittest
import tempfile
import textwrap
import sys

from .._jess import Atom, Molecule

MOLECULE = textwrap.dedent(
    """
    HEADER    DNA RECOMBINATION                       05-DEC-97   1A0P              
    COMPND    MOL_ID: 1;                                                            
    ATOM      1  N   GLN A   3       8.171 -51.403  42.886  1.00 55.63           N
    ATOM      2  CA  GLN A   3       9.475 -50.697  42.743  1.00 56.29           C
    ATOM      3  C   GLN A   3      10.215 -51.213  41.516  1.00 55.54           C
    ATOM      4  O   GLN A   3      10.401 -50.398  40.585  1.00 56.57           O
    ATOM      5  CB  GLN A   3      10.267 -50.747  44.040  1.00 72.29           C
    """
).strip()


class TestMolecule(unittest.TestCase):

    def _create_atom(self, **kwargs):
        default = dict(serial=1, name='N', altloc=' ', residue_name='GLN', chain_id='A', residue_number=3, x=8.171, y=-51.403, z=42.886, segment='', insertion_code=' ', occupancy=1.0, temperature_factor=55.63, charge=0, element='N')
        default.update(kwargs)
        return Atom(**default)

    def test_loads(self):
        molecule = Molecule.loads(MOLECULE)
        self.assertEqual(len(molecule), 5)
        self.assertEqual(molecule.id, '1A0P')

    @unittest.skipUnless(sys.implementation.name == "cpython", "only available on CPython")
    def test_sizeof(self):
        molecule = Molecule.loads(MOLECULE)
        self.assertGreater(sys.getsizeof(molecule), 0)

    @unittest.skipIf(os.name == "nt", "permission errors on Windows")
    def test_load_filename(self):
        with tempfile.NamedTemporaryFile("w", suffix=".pdb") as f:
            f.write(MOLECULE)
            f.flush()
            molecule = Molecule.load(f.name)
        self.assertEqual(len(molecule), 5)
        self.assertEqual(molecule.id, '1A0P')

    @unittest.skipIf(os.name == "nt", "permission errors on Windows")
    def test_load_file(self):
        with tempfile.NamedTemporaryFile("r+", suffix=".pdb") as f:
            f.write(MOLECULE)
            f.flush()
            f.seek(0)
            molecule = Molecule.load(f)
        self.assertEqual(len(molecule), 5)
        self.assertEqual(molecule.id, '1A0P')

    @unittest.skipIf(os.name == "nt", "permission errors on Windows")
    def test_load_error(self):
        self.assertRaises(FileNotFoundError, Molecule.load, "/some/nonsensical/file")
        self.assertRaises(IsADirectoryError, Molecule.load, os.path.dirname(__file__))

    def test_init(self):
        atoms = [
            self._create_atom(serial=1, name='N'),
            self._create_atom(serial=2, name='CA'),
            self._create_atom(serial=3, name='C'),
            self._create_atom(serial=4, name='O'),
        ]
        molecule = Molecule(atoms)
        self.assertEqual(molecule[0].name, 'N')
        self.assertEqual(molecule[1].name, 'CA')
        self.assertEqual(molecule[2].name, 'C')
        self.assertEqual(molecule[3].name, 'O')

    def test_init_long_id(self):
        mol = Molecule.loads(MOLECULE, id="long identifier")
        self.assertEqual(mol.id, "long identifier")

    def test_getitem_slicing(self):
        mol = Molecule.loads(MOLECULE)
        mol2 = mol[1:3]
        self.assertEqual(mol2.id, mol.id)
        self.assertEqual(len(mol2), 2)
        self.assertEqual(mol2[0].name, "CA")
        self.assertEqual(mol2[1].name, "C")

    def test_hash(self):
        atoms = [
            self._create_atom(serial=1, name='N'),
            self._create_atom(serial=2, name='CA'),
            self._create_atom(serial=3, name='C'),
            self._create_atom(serial=4, name='O'),
        ]
        mol1 = Molecule(atoms)
        mol2 = Molecule(atoms)
        self.assertEqual(hash(mol1), hash(mol2))
        self.assertIsNot(mol1, mol2)
        mol3 = Molecule(atoms[:-1])
        self.assertNotEqual(hash(mol1), hash(mol3))

    def test_eq(self):
        atoms = [
            self._create_atom(serial=1, name='N'),
            self._create_atom(serial=2, name='CA'),
            self._create_atom(serial=3, name='C'),
            self._create_atom(serial=4, name='O'),
        ]
        mol1 = Molecule(atoms)
        mol2 = Molecule(atoms)
        self.assertEqual(mol1, mol2)
        self.assertIsNot(mol1, mol2)
        mol3 = Molecule(atoms[:-1])
        self.assertNotEqual(mol1, mol3)

    def test_copy(self):
        atoms = [
            self._create_atom(serial=1, name='N'),
            self._create_atom(serial=2, name='CA'),
            self._create_atom(serial=3, name='C'),
            self._create_atom(serial=4, name='O'),
        ]
        mol1 = Molecule(atoms)
        mol2 = mol1.copy()
        self.assertEqual(list(mol1), list(mol2))
        self.assertEqual(mol1.id, mol2.id)
        self.assertEqual(mol1, mol2)

    def test_pickle_roundtrip(self):
        atoms = [
            self._create_atom(serial=1, name='N'),
            self._create_atom(serial=2, name='CA'),
            self._create_atom(serial=3, name='C'),
            self._create_atom(serial=4, name='O'),
        ]
        mol1 = Molecule(atoms)
        mol2 = pickle.loads(pickle.dumps(mol1))
        self.assertEqual(list(mol1), list(mol2))
        self.assertEqual(mol1.id, mol2.id)
        self.assertEqual(mol1, mol2)
