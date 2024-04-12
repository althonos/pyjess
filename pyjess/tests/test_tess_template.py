import unittest
import tempfile
import textwrap

from .._jess import Template

TEMPLATE = textwrap.dedent(
    """
    ATOM      0  NZ  LYS A1030       1.431  -9.717  -1.722
    ATOM      0  CG  ASP A1132       0.000   0.000   0.000
    ATOM      0  OD1 ASP A1132       0.595   1.085   0.000
    ATOM      0  OD2 ASP A1132       0.609  -1.111   0.000
    ATOM      1  NE  ARG A1136       3.953   0.597  -1.721 K
    ATOM      0  CG  ASN A1137       0.736  -1.869  -4.261
    ATOM      0  OD1 ASN A1137       1.760  -2.435  -4.344
    ATOM      0  ND2 ASN A1137       0.264  -1.477  -3.134
    ATOM      0  CG  ASP A1150      -0.010  -5.408  -2.070
    ATOM      0  OD1 ASP A1150      -0.587  -5.769  -1.153
    ATOM      0  OD2 ASP A1150       1.085  -5.212  -2.011
    """
).strip()


class TestTessTemplate(unittest.TestCase):

    def test_load(self):

        with tempfile.NamedTemporaryFile(mode="w") as dst:
            dst.write(TEMPLATE)
            dst.flush()
            template = Template("template", dst.name)
        
        self.assertEqual(len(template), 11)
        self.assertEqual(template.dimension, 5)
        self.assertEqual(template[0].residue_names, ["LYS"])
        self.assertEqual(template[0].atom_names, ["NZ"])
        self.assertEqual(template[1].atom_names, ["CG"])
        self.assertEqual(template[2].residue_number, 1132)
        self.assertEqual(template[-1].residue_number, 1150)