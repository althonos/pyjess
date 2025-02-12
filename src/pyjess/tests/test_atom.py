import ctypes
import unittest
import pickle

from .._jess import Atom


class TestAtom(unittest.TestCase):
    
    def _create_atom(self, **kwargs):
        default = dict(serial=1, name='N', altloc=' ', residue_name='GLN', chain_id='A', residue_number=3, x=8.171, y=-51.403, z=42.886, segment='', insertion_code=' ', occupancy=1.0, temperature_factor=55.63, charge=0, element='N')
        default.update(kwargs)
        return Atom(**default)

    def test_load(self):
        atom = Atom.loads("ATOM     39  CA  PRO A 469     -14.948   2.091  10.228  1.00 27.71           C")
        self.assertEqual(atom.serial, 39)
        self.assertEqual(atom.name, 'CA')
        self.assertEqual(atom.residue_name, 'PRO')
        self.assertEqual(atom.chain_id, 'A')
        self.assertEqual(atom.residue_number, 469)
        self.assertAlmostEqual(atom.x, -14.948, places=3)
        self.assertAlmostEqual(atom.y, 2.091, places=3)
        self.assertAlmostEqual(atom.z, 10.228, places=3)
        self.assertAlmostEqual(atom.occupancy, 1.00, places=2)
        self.assertAlmostEqual(atom.temperature_factor, 27.71, places=2)
        self.assertEqual(atom.segment, '')
        self.assertEqual(atom.element, 'C')
        self.assertEqual(atom.charge, 0)

    def test_init_consistency(self):
        loaded = Atom.loads("ATOM     39  CA  PRO A 469     -14.948   2.091  10.228  1.00 27.71           C")
        created = Atom(
            serial=loaded.serial,
            name=loaded.name,
            residue_name=loaded.residue_name,
            chain_id=loaded.chain_id,
            altloc=loaded.altloc,
            insertion_code=loaded.insertion_code,
            residue_number=loaded.residue_number,
            x=loaded.x,
            y=loaded.y,
            z=loaded.z,
            occupancy=loaded.occupancy,
            temperature_factor=loaded.temperature_factor,
            element=loaded.element,
        )
        self.assertEqual(loaded, created)
        
    def test_init_invalid_chain_id(self):
        self.assertRaises(ValueError, self._create_atom, chain_id="too long")

    def test_init_invalid_name(self):
        self.assertRaises(ValueError, self._create_atom, name="too long")

    def test_init_invalid_residue_name(self):
        self.assertRaises(ValueError, self._create_atom, residue_name="too long")

    def test_hash(self):
        a1 = self._create_atom()
        a2 = self._create_atom()
        self.assertEqual(hash(a1), hash(a2))
        self.assertIsNot(a1, a2)
        a3 = self._create_atom(x=1.0)
        self.assertNotEqual(hash(a1), hash(a3))

    def test_eq(self):
        a1 = self._create_atom()
        a2 = self._create_atom()
        self.assertEqual(a1, a2)
        self.assertIsNot(a1, a2)
        a3 = self._create_atom(x=1.0)
        self.assertNotEqual(a1, a3)

    def test_repr_roundtrip(self):
        atom = self._create_atom()
        copy = eval(repr(atom))
        self.assertEqual(atom.serial, copy.serial)
        self.assertEqual(atom.altloc, copy.altloc)
        self.assertEqual(atom.name, copy.name)
        self.assertEqual(atom.residue_name, copy.residue_name)
        self.assertEqual(atom.residue_number, copy.residue_number)
        self.assertEqual(atom.element, copy.element)
        self.assertEqual(atom.insertion_code, copy.insertion_code)
        self.assertEqual(atom.chain_id, copy.chain_id)
        self.assertEqual(atom.occupancy, copy.occupancy)
        self.assertEqual(atom.temperature_factor, copy.temperature_factor)
        self.assertEqual(atom.charge, copy.charge)
        self.assertEqual(atom.x, copy.x)
        self.assertEqual(atom.y, copy.y)
        self.assertEqual(atom.z, copy.z)
        self.assertEqual(atom, copy)

    def test_pickle_roundtrip(self):
        atom = self._create_atom()
        copy = pickle.loads(pickle.dumps(atom))
        self.assertEqual(atom.serial, copy.serial)
        self.assertEqual(atom.altloc, copy.altloc)
        self.assertEqual(atom.name, copy.name)
        self.assertEqual(atom.residue_name, copy.residue_name)
        self.assertEqual(atom.residue_number, copy.residue_number)
        self.assertEqual(atom.element, copy.element)
        self.assertEqual(atom.insertion_code, copy.insertion_code)
        self.assertEqual(atom.chain_id, copy.chain_id)
        self.assertEqual(atom.occupancy, copy.occupancy)
        self.assertEqual(atom.temperature_factor, copy.temperature_factor)
        self.assertEqual(atom.charge, copy.charge)
        self.assertEqual(atom.x, copy.x)
        self.assertEqual(atom.y, copy.y)
        self.assertEqual(atom.z, copy.z)
        self.assertEqual(atom, copy)