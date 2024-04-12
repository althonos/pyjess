import unittest
import textwrap

from .._jess import Atom, Molecule


class TestMolecule(unittest.TestCase):

    def _create_atom(self, **kwargs):
        default = dict(serial=1, name='N', altloc=' ', residue_name='GLN', chain_id='A', residue_number=3, x=8.171, y=-51.403, z=42.886, segment='', insertion_code=' ', occupancy=1.0, temperature_factor=55.63, charge=0, element='N')
        default.update(kwargs)
        return Atom(**default)

    def test_load(self):
        molecule = Molecule.loads(
            textwrap.dedent(
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
        )
        self.assertEqual(len(molecule), 5)
        self.assertEqual(molecule.id, '1A0P')

    def test_init_invalid_id(self):
        self.assertRaises(ValueError, Molecule, atoms=[], id="too long")
        self.assertRaises(ValueError, Molecule, atoms=[], id="x")

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

