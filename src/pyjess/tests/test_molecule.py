import datetime
import itertools
import io
import os
import pickle
import unittest
import tempfile
import textwrap
import sys
import shutil
import warnings

from .utils import files
from . import data
from .._jess import Atom, Molecule

try:
    import gemmi
except ImportError:
    gemmi = None

try:
    import Bio.PDB
except ImportError:
    Bio = None

try:
    from biotite.structure.io.pdb import PDBFile
except ImportError:
    PDBFile = None

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

    def assertHeteroAtomEqual(self, a1, a2):
        # self.assertEqual(a1.serial, a2.serial)
        self.assertEqual(a1.name, a2.name)
        self.assertEqual(a1.altloc, a2.altloc)
        self.assertEqual(a1.residue_name, a2.residue_name)
        # self.assertEqual(a1.chain_id, a2.chain_id)
        # self.assertEqual(a1.residue_number, a2.residue_number)
        self.assertEqual(a1.insertion_code, a2.insertion_code)
        self.assertEqual(a1.element, a2.element)
        self.assertAlmostEqual(a1.x, a2.x, places=5)
        self.assertAlmostEqual(a1.y, a2.y, places=5)
        self.assertAlmostEqual(a1.z, a2.z, places=5)
        self.assertAlmostEqual(a1.occupancy, a2.occupancy, places=5)
        self.assertAlmostEqual(a1.temperature_factor, a2.temperature_factor, places=5)
        self.assertAlmostEqual(a1.charge, a2.charge, places=5)

    def assertAtomEqual(self, a1, a2):
        self.assertEqual(a1.serial, a2.serial)
        self.assertEqual(a1.name, a2.name)
        self.assertEqual(a1.altloc, a2.altloc)
        self.assertEqual(a1.residue_name, a2.residue_name)
        self.assertEqual(a1.chain_id, a2.chain_id)
        self.assertEqual(a1.residue_number, a2.residue_number)
        self.assertEqual(a1.insertion_code, a2.insertion_code)
        self.assertEqual(a1.element, a2.element)
        self.assertAlmostEqual(a1.x, a2.x, places=5)
        self.assertAlmostEqual(a1.y, a2.y, places=5)
        self.assertAlmostEqual(a1.z, a2.z, places=5)
        self.assertAlmostEqual(a1.occupancy, a2.occupancy, places=5)
        self.assertAlmostEqual(a1.temperature_factor, a2.temperature_factor, places=5)
        self.assertAlmostEqual(a1.charge, a2.charge, places=5)

    def _create_atom(self, **kwargs):
        default = dict(serial=1, name='N', altloc=' ', residue_name='GLN', chain_id='A', residue_number=3, x=8.171, y=-51.403, z=42.886, segment='', insertion_code=' ', occupancy=1.0, temperature_factor=55.63, charge=0, element='N')
        default.update(kwargs)
        return Atom(**default)

    def test_loads(self):
        molecule = Molecule.loads(MOLECULE)
        self.assertEqual(len(molecule), 5)
        self.assertEqual(molecule.id, '1A0P')
        self.assertEqual(molecule.date, datetime.date(1997, 12, 5))
        self.assertEqual(molecule.name, 'DNA RECOMBINATION')

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
        self.assertEqual(molecule.date, datetime.date(1997, 12, 5))
        self.assertEqual(molecule.name, 'DNA RECOMBINATION')

    @unittest.skipIf(os.name == "nt", "permission errors on Windows")
    def test_load_file(self):
        with tempfile.NamedTemporaryFile("r+", suffix=".pdb") as f:
            f.write(MOLECULE)
            f.flush()
            f.seek(0)
            molecule = Molecule.load(f)
        self.assertEqual(len(molecule), 5)
        self.assertEqual(molecule.id, '1A0P')
        self.assertEqual(molecule.date, datetime.date(1997, 12, 5))
        self.assertEqual(molecule.name, 'DNA RECOMBINATION')

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
        self.assertIs(molecule.id, None)
        self.assertIs(molecule.name, None)
        self.assertIs(molecule.date, None)
        self.assertEqual(molecule[0].name, 'N')
        self.assertEqual(molecule[1].name, 'CA')
        self.assertEqual(molecule[2].name, 'C')
        self.assertEqual(molecule[3].name, 'O')

    def test_loads_with_kwargs(self):
        mol = Molecule.loads(
            MOLECULE,
            id="long identifier",
            name='struct',
            date=datetime.date.today(),
        )
        self.assertEqual(mol.id, "long identifier")
        self.assertEqual(mol.name, "struct")
        self.assertEqual(mol.date, datetime.date.today())

    def test_getitem_slicing(self):
        mol = Molecule.loads(MOLECULE)
        mol2 = mol[1:3]
        self.assertEqual(mol2.id, mol.id)
        self.assertEqual(mol2.date, mol.date)
        self.assertEqual(mol2.name, mol.name)
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
        self.assertEqual(mol1.date, mol2.date)
        self.assertEqual(mol1.name, mol2.name)
        self.assertEqual(mol1, mol2)

    def test_pickle_metadata_roundtrip(self):
        atoms = [
            self._create_atom(serial=1, name='N'),
            self._create_atom(serial=2, name='CA'),
            self._create_atom(serial=3, name='C'),
            self._create_atom(serial=4, name='O'),
        ]
        mol1 = Molecule(atoms, id="ABC1", name="TEST MOLECULE", date=datetime.date.today())
        mol2 = pickle.loads(pickle.dumps(mol1))
        self.assertEqual(list(mol1), list(mol2))
        self.assertEqual(mol1.id, mol2.id)
        self.assertEqual(mol1.date, mol2.date)
        self.assertEqual(mol1.name, mol2.name)
        self.assertEqual(mol1, mol2)

    def test_dumps_roundtrip(self):
        molecule = Molecule.loads(MOLECULE)
        dump = Molecule.loads(molecule.dumps().strip())
        self.assertEqual(molecule, dump)

    def test_dumps(self):
        molecule = Molecule.loads(MOLECULE)
        expected = textwrap.dedent(
            """
            HEADER    DNA RECOMBINATION                       05-DEC-97   1A0P
            ATOM      1  N   GLN A   3       8.171 -51.403  42.886  1.00 55.63           N
            ATOM      2  CA  GLN A   3       9.475 -50.697  42.743  1.00 56.29           C
            ATOM      3  C   GLN A   3      10.215 -51.213  41.516  1.00 55.54           C
            ATOM      4  O   GLN A   3      10.401 -50.398  40.585  1.00 56.57           O
            ATOM      5  CB  GLN A   3      10.267 -50.747  44.040  1.00 72.29           C
            """
        ).strip()
        self.maxDiff = None
        actual = molecule.dumps().strip()
        self.assertMultiLineEqual(actual, expected)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_cif_dump(self):
        with files(data).joinpath("1AMY.cif").open() as f:
            text= f.read()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cif_molecule = Molecule.loads(text, format="cif")

        first_5 = cif_molecule[0:5]
        expected = textwrap.dedent(
            """
            HEADER    HYDROLASE (O-GLYCOSYL)                  13-MAY-95   1AMY
            ATOM      1  N   GLN A   1       6.240  48.686  17.460  1.00 27.79           N
            ATOM      2  CA  GLN A   1       5.440  49.851  17.773  1.00 16.62           C
            ATOM      3  C   GLN A   1       6.628  50.721  18.086  1.00 14.24           C
            ATOM      4  O   GLN A   1       7.313  50.396  19.052  1.00 11.61           O
            ATOM      5  CB  GLN A   1       4.588  49.715  19.030  1.00 18.54           C
            """
        ).strip()

        actual = first_5.dumps().strip()
        self.maxDiff = None
        self.assertMultiLineEqual(actual, expected)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_cif_vs_pdb(self):
        with files(data).joinpath("1AMY.cif").open() as f:
            text = f.read()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cif_molecule = Molecule.loads(text, format="cif")

        with files(data).joinpath("1AMY.pdb").open() as f:
            text = f.read()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pdb_molecule = Molecule.loads(text, format="pdb")

        first5_pdb = pdb_molecule[0:5]
        first5_cif = cif_molecule[0:5]

        self.maxDiff = None
        self.assertEqual(first5_pdb, first5_cif)
        # I set write_header to False since the dates in the cif and pdb files dont match
        self.assertMultiLineEqual(
            first5_pdb.dumps(write_header=False).strip(), 
            first5_cif.dumps(write_header=False).strip()
        )

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_load_mmcif_comment(self):
        # load mmCIF file into a buffer and add a comment at the top
        buffer = io.StringIO()
        buffer.write("# this is a comment line \n")
        buffer.write("# and a second comment line \n")
        with files(data).joinpath("1AMY.cif").open() as f:
            shutil.copyfileobj(f, buffer)
            
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # should work when explicitly given "CIF" formula
            buffer.seek(0)
            cif_molecule = Molecule.load(buffer, format="cif")
            self.assertEqual(len(cif_molecule), 3339)
            # should work when not given format
            buffer.seek(0)
            cif_molecule = Molecule.load(buffer, format="detect")
            self.assertEqual(len(cif_molecule), 3339)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_loads_mmcif_comment(self):
        # load mmCIF file into a buffer and add a comment at the top
        buffer = io.StringIO()
        buffer.write("# this is a comment line \n")
        buffer.write("# and a second comment line \n")
        with files(data).joinpath("1AMY.cif").open() as f:
            shutil.copyfileobj(f, buffer)
        text = buffer.getvalue()
            
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # should work when explicitly given "CIF" formula
            cif_molecule = Molecule.loads(text, format="cif")
            self.assertEqual(len(cif_molecule), 3339)
            # should work when not given format
            cif_molecule = Molecule.loads(text, format="detect")
            self.assertEqual(len(cif_molecule), 3339)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_load_consistency_no_skip_hetatm(self):
        with files(data).joinpath("1AMY.pdb").open() as f:
            pdb_molecule = Molecule.load(f, format="pdb")
        with files(data).joinpath("1AMY.cif").open() as f:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                cif_molecule = Molecule.load(f, format="cif")
        self.assertEqual(len(cif_molecule), 3339)
        self.assertEqual(len(pdb_molecule), 3339)
        for pdb_atom, cif_atom in itertools.islice(zip(pdb_molecule, cif_molecule), 3184):
            self.assertAtomEqual(pdb_atom, cif_atom)
        for pdb_atom, cif_atom in itertools.islice(zip(pdb_molecule, cif_molecule), 3184, 3339):
            self.assertHeteroAtomEqual(pdb_atom, cif_atom)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_load_consistency_no_skip_hetatm_use_author(self):
        with files(data).joinpath("1AMY.pdb").open() as f:
            pdb_molecule = Molecule.load(f, format="pdb")
        with files(data).joinpath("1AMY.cif").open() as f:
            cif_molecule = Molecule.load(f, format="cif", use_author=True)
        self.assertEqual(len(cif_molecule), 3339)
        self.assertEqual(len(pdb_molecule), 3339)
        for pdb_atom, cif_atom in itertools.islice(zip(pdb_molecule, cif_molecule), 3184):
            self.assertAtomEqual(pdb_atom, cif_atom)
        for pdb_atom, cif_atom in itertools.islice(zip(pdb_molecule, cif_molecule), 3184, 3339):
            self.assertHeteroAtomEqual(pdb_atom, cif_atom)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_load_consistency_skip_hetatm(self):
        with files(data).joinpath("1AMY.pdb").open() as f:
            pdb_molecule = Molecule.load(f, format="pdb", skip_hetatm=True)
        with files(data).joinpath("1AMY.cif").open() as f:
            cif_molecule = Molecule.load(f, format="cif", skip_hetatm=True)
        self.assertEqual(len(cif_molecule), 3184)
        self.assertEqual(len(pdb_molecule), 3184)
        for pdb_atom, cif_atom in itertools.islice(zip(pdb_molecule, cif_molecule), 3184):
            self.assertAtomEqual(pdb_atom, cif_atom)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_load_consistency_skip_hetatm_use_author(self):
        with files(data).joinpath("1AMY.pdb").open() as f:
            pdb_molecule = Molecule.load(f, format="pdb", skip_hetatm=True)
        with files(data).joinpath("1AMY.cif").open() as f:
            cif_molecule = Molecule.load(f, format="cif", skip_hetatm=True, use_author=True)
        self.assertEqual(len(cif_molecule), 3184)
        self.assertEqual(len(pdb_molecule), 3184)
        for pdb_atom, cif_atom in itertools.islice(zip(pdb_molecule, cif_molecule), 3184):
            self.assertAtomEqual(pdb_atom, cif_atom)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(Bio, "biopython not available")
    def test_from_biopython(self):
        parser = Bio.PDB.PDBParser()
        with files(data).joinpath("1AMY.pdb").open() as f:
            structure = parser.get_structure('1amy', f)
            bio_mol = Molecule.from_biopython(structure)
        with files(data).joinpath("1AMY.pdb").open() as f:
            pdb_mol = Molecule.load(f, format="pdb")
        self.assertEqual(len(pdb_mol), 3339)
        self.assertEqual(len(bio_mol), 3339)
        for pdb_atom, bio_atom in zip(pdb_mol, bio_mol):
            self.assertAtomEqual(pdb_atom, bio_atom)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(gemmi, "gemmi not available")
    def test_from_gemmi(self):
        with files(data).joinpath("1AMY.pdb").open() as f:
            structure = gemmi.read_pdb_string(f.read())
            gemmi_mol = Molecule.from_gemmi(structure[0])
        with files(data).joinpath("1AMY.pdb").open() as f:
            pdb_mol = Molecule.load(f, format="pdb")
        self.assertEqual(len(pdb_mol), 3339)
        self.assertEqual(len(gemmi_mol), 3339)
        for pdb_atom, gemmi_atom in zip(pdb_mol, gemmi_mol):
            self.assertAtomEqual(pdb_atom, gemmi_atom)

    @unittest.skipUnless(files, "importlib.resources not available")
    @unittest.skipUnless(PDBFile, "biotite not available")
    def test_from_biotite(self):
        with files(data).joinpath("1AMY.pdb").open() as f:
            pdb_file = PDBFile.read(f)
            structure = pdb_file.get_structure(altloc="all", extra_fields=["atom_id", "b_factor", "occupancy", "charge"])
            biotite_mol = Molecule.from_biotite(structure[0])
        with files(data).joinpath("1AMY.pdb").open() as f:
            pdb_mol = Molecule.load(f, format="pdb")
        self.assertEqual(len(pdb_mol), 3339)
        self.assertEqual(len(biotite_mol), 3339)
        for pdb_atom, biotite_atom in zip(pdb_mol, biotite_mol):
            self.assertAtomEqual(pdb_atom, biotite_atom)