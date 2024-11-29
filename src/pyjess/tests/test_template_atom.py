import unittest
import sys
import pickle

from .._jess import TemplateAtom


class TestTemplateAtom(unittest.TestCase):

    def test_load(self):
        atom = TemplateAtom.loads("ATOM      1  NE  ARG A1136       3.953   0.597  -1.721 K")
        self.assertEqual(atom.match_mode, 1)
        self.assertEqual(atom.atom_names, ("NE",))
        self.assertEqual(atom.residue_names, ("ARG", "LYS",))
        self.assertEqual(atom.chain_id, "A")
        self.assertEqual(atom.residue_number, 1136)
        self.assertEqual(atom.x, 3.953)
        self.assertEqual(atom.y, 0.597)
        self.assertEqual(atom.z, -1.721)

    def _create_atom(self, **kwargs):
        default = {
            "atom_names": ["NE"],
            "residue_names": ["ARG"],
            "chain_id": "A",
            "residue_number": 1,
            "x": 0.0,
            "y": 0.0,
            "z": 0.0,
            "match_mode": 0,
        }
        default.update(kwargs)
        return TemplateAtom(**default)

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

    @unittest.skipUnless(sys.implementation.name == "cpython", "only available on CPython")
    def test_sizeof(self):
        atom = self._create_atom()
        self.assertGreater(sys.getsizeof(atom), 0)

    def test_init_invalid_match_code(self):
        self.assertRaises(ValueError, self._create_atom, match_mode=2000)
        self.assertRaises(ValueError, self._create_atom, match_mode=-10)

    def test_init_invalid_residue_name(self):
        self.assertRaises(ValueError, self._create_atom, residue_names=["something"])

    def test_init_consistency(self):
        loaded = TemplateAtom.loads("ATOM      1  NE  ARG A1136       3.953   0.597  -1.721 K")
        created = TemplateAtom(
            chain_id=loaded.chain_id, 
            residue_number=loaded.residue_number,
            x=loaded.x,
            y=loaded.y,
            z=loaded.z,
            residue_names=loaded.residue_names,
            atom_names=loaded.atom_names,
            distance_weight=loaded.distance_weight,
            match_mode=loaded.match_mode
        )
        for attribute in ("atom_names", "residue_names", "chain_id", "x", "y", "z", "match_mode"):
            self.assertEqual(getattr(loaded, attribute), getattr(created, attribute))
        self.assertEqual(loaded, created)

    def test_repr_roundtrip(self):
        atom = self._create_atom()
        copy = eval(repr(atom))
        for attribute in ("atom_names", "residue_names", "chain_id", "x", "y", "z", "match_mode"):
            self.assertEqual(getattr(copy, attribute), getattr(atom, attribute))
        self.assertEqual(atom, copy)

    def test_pickle_roundtrip(self):
        atom = self._create_atom()
        copy = pickle.loads(pickle.dumps(atom))
        for attribute in ("atom_names", "residue_names", "chain_id", "x", "y", "z", "match_mode"):
            self.assertEqual(getattr(copy, attribute), getattr(atom, attribute))
        self.assertEqual(atom, copy)