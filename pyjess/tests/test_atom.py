import unittest

from .._jess import Atom


class TestAtom(unittest.TestCase):
    
    def _create_atom(self, **kwargs):
        default = dict(serial=1, name='N', altloc=' ', residue_name='GLN', chain_id='A', residue_number=3, x=8.171, y=-51.403, z=42.886, segment='', insertion_code=' ', occupancy=1.0, temperature_factor=55.63, charge=0, element='N')
        default.update(kwargs)
        return Atom(**default)

    def test_init_invalid_chain_id(self):
        self.assertRaises(ValueError, self._create_atom, chain_id="hello")