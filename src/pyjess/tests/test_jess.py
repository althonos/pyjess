import math
import unittest
import sys
import pickle

from .._jess import Template, Molecule, Jess
from .utils import files
from . import data


class TestJess(unittest.TestCase):

    def test_init_empty(self):
        jess = Jess()
        self.assertEqual(len(jess), 0)
        self.assertFalse(bool(jess))

    def test_query_empty(self):
        mol = Molecule()
        jess = Jess()
        hits = jess.query(mol, 2, 2, 4)
        self.assertRaises(StopIteration, next, hits)

    def test_hash_empty(self):
        j1 = Jess()
        j2 = Jess()
        self.assertEqual(j1, j2)
        self.assertEqual(hash(j1), hash(j2))
        self.assertIsNot(j1, j2)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_copy(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template1 = Template.load(f)
        with files(data).joinpath("template_02.qry").open() as f:
            template2 = Template.load(f)
        jess = Jess([template1, template2])
        copy = jess.copy()
        self.assertEqual(jess, copy)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_hash(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template1 = Template.load(f)
        with files(data).joinpath("template_02.qry").open() as f:
            template2 = Template.load(f)
        j1 = Jess([template1, template2])
        j2 = Jess([template1, template2])
        self.assertEqual(hash(j1), hash(j2))
        self.assertIsNot(j1, j2)
        j3 = Jess([template1])
        self.assertNotEqual(hash(j1), hash(j3))

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_eq(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template1 = Template.load(f)
        with files(data).joinpath("template_02.qry").open() as f:
            template2 = Template.load(f)
        j1 = Jess([template1, template2])
        j2 = Jess([template1, template2])
        self.assertEqual(j1, j2)
        self.assertIsNot(j1, j2)
        j3 = Jess([template1])
        self.assertNotEqual(j1, j3)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_pickle_roundtrip(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template1 = Template.load(f)
        with files(data).joinpath("template_02.qry").open() as f:
            template2 = Template.load(f)
        jess = Jess([template1, template2])
        copy = pickle.loads(pickle.dumps(jess))
        self.assertEqual(jess, copy)

    @unittest.skipUnless(sys.implementation.name == "cpython", "only available on CPython")
    @unittest.skipUnless(files, "importlib.resources not available")
    def test_sizeof(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template = Template.load(f)
            jess = Jess([template])
        self.assertGreater(sys.getsizeof(jess), 0)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_getitem(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template = Template.load(f)
            jess = Jess([template])

        self.assertEqual(jess[0], template)
        self.assertEqual(len(jess[:1]), 1)
        self.assertEqual(len(jess[1:]), 0)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_multiple_query_split(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template1 = Template.load(f)
        with files(data).joinpath("template_02.qry").open() as f:
            template2 = Template.load(f)
        j0 = Jess([template1, template2])
        j1 = Jess([template1])
        j2 = Jess([template2])

        with files(data).joinpath("1AMY.pdb").open() as f:
            molecule = Molecule.load(f)

        h0 = list(j0.query(molecule, 2, 5, 5))
        h1 = list(j1.query(molecule, 2, 5, 5))
        h2 = list(j2.query(molecule, 2, 5, 5))

        self.assertEqual(len(h0), len(h1) + len(h2))
        for hit0, hit1 in zip(h0, (*h1, *h2)):
            self.assertEqual(hit0.atoms(), hit1.atoms())
            self.assertEqual(hit0.evalue, hit1.evalue)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_query_template_subclass(self):

        class MyTemplate(Template):
            pass
        
        with files(data).joinpath("template_01.qry").open() as f:
            template = MyTemplate.load(f)
            jess = Jess([template])
        with files(data).joinpath("pdb1lnb.pdb").open() as f:
            molecule = Molecule.load(f)

        self.assertEqual(jess[0], template)
        self.assertIsInstance(jess[0], MyTemplate)
        hits = list(jess.query(molecule, 1, 2, 2))
        self.assertIsInstance(hits[0].template, MyTemplate)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_query_max_candidates(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template = Template.load(f)
            jess = Jess([template])
        with files(data).joinpath("pdb1lnb.pdb").open() as f:
            molecule = Molecule.load(f)

        hits = list(jess.query(molecule, 2, 5, 5))
        self.assertEqual(len(hits), 3)

        hits = list(jess.query(molecule, 2, 5, 5, max_candidates=3))
        self.assertEqual(len(hits), 3)

        hits = list(jess.query(molecule, 2, 5, 5, max_candidates=2))
        self.assertEqual(len(hits), 2)

        hits = list(jess.query(molecule, 2, 5, 5, max_candidates=1))
        self.assertEqual(len(hits), 1)

        with self.assertRaises(ValueError):
            hits = list(jess.query(molecule, 2, 5, 5, max_candidates=-1))
        
    @unittest.skipUnless(files, "importlib.resources not available")
    def test_query(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template = Template.load(f)
            jess = Jess([template])
        with files(data).joinpath("pdb1lnb.pdb").open() as f:
            molecule = Molecule.load(f)

        hits = list(jess.query(molecule, 1, 2, 2))
        self.assertEqual(len(hits), 1)
        self.assertIs(hits[0].template, template)
        self.assertAlmostEqual(hits[0].rmsd, 0.555, places=3)
        self.assertAlmostEqual(hits[0].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[0].log_evalue, -2.04, places=1)
        self.assertAlmostEqual(hits[0].evalue, math.exp(-2.04), places=1)

        hits = list(jess.query(molecule, 2, 5, 3))
        self.assertEqual(len(hits), 2)
        self.assertIs(hits[0].template, template)
        self.assertAlmostEqual(hits[0].rmsd, 0.555, places=3)
        self.assertAlmostEqual(hits[0].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[0].log_evalue, -2.04, places=1)
        self.assertAlmostEqual(hits[0].evalue, math.exp(-2.04), places=1)
        self.assertIs(hits[1].template, template)
        self.assertAlmostEqual(hits[1].rmsd, 1.440, places=3)
        self.assertAlmostEqual(hits[1].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[1].log_evalue, 0.17, places=1)
        self.assertAlmostEqual(hits[1].evalue, math.exp(0.17), places=1)

        hits = list(jess.query(molecule, 2, 5, 5))
        self.assertEqual(len(hits), 3)
        self.assertIs(hits[0].template, template)
        self.assertAlmostEqual(hits[0].rmsd, 0.555, places=3)
        self.assertAlmostEqual(hits[0].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[0].log_evalue, -2.04, places=1)
        self.assertAlmostEqual(hits[0].evalue, math.exp(-2.04), places=1)
        self.assertIs(hits[1].template, template)
        self.assertAlmostEqual(hits[1].rmsd, 1.440, places=3)
        self.assertAlmostEqual(hits[1].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[1].log_evalue, 0.17, places=1)
        self.assertAlmostEqual(hits[1].evalue, math.exp(0.17), places=1)
        self.assertIs(hits[2].template, template)
        self.assertAlmostEqual(hits[2].rmsd, 1.644, places=3)
        self.assertAlmostEqual(hits[2].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[2].log_evalue, 0.68, places=1)
        self.assertAlmostEqual(hits[2].evalue, math.exp(0.68), places=1)

        hits = list(jess.query(molecule.conserved(10.0), 1, 2, 2))
        self.assertEqual(len(hits), 1)
        self.assertIs(hits[0].template, template)
        self.assertAlmostEqual(hits[0].rmsd, 0.555, places=3)
        self.assertAlmostEqual(hits[0].determinant, 1.0, places=3)
        self.assertAlmostEqual(hits[0].log_evalue, -2.10, places=1)
        self.assertAlmostEqual(hits[0].evalue, math.exp(-2.10), places=1)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_query_best_match(self):
        with files(data).joinpath("template_01.qry").open() as f:
            template = Template.load(f)
        jess = Jess([template, template, template])
        with files(data).joinpath("pdb1lnb.pdb").open() as f:
            molecule = Molecule.load(f)

        hits = list(jess.query(molecule, 2, 5, 5, best_match=True))
        self.assertEqual(len(hits), 3)
        for hit in hits:
            self.assertIs(hit.template, template)
            self.assertAlmostEqual(hit.rmsd, 0.555, places=3)
            self.assertAlmostEqual(hit.determinant, 1.0, places=3)
            self.assertAlmostEqual(hit.log_evalue, -2.04, places=1)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_mcsa_query_no_reorder(self):
        with files(data).joinpath("1.3.3.tpl").open() as f:
            template = Template.load(f)
            jess = Jess([template])
        with files(data).joinpath("1AMY.pdb").open() as f:
            molecule = Molecule.load(f)
        with files(data).joinpath("1AMY+1.3.3.txt").open() as f:
            results = list(filter(None, f.read().split("REMARK")))

        hits = list(jess.query(molecule, 2, 4, 4, reorder=False))
        self.assertEqual(len(hits), len(results))
        for hit, block in zip(hits, results):
            self.assertIs(hit.template, template)

            lines = block.strip().splitlines()
            query_id, rmsd, template_id, _, determinant, _, logE = lines[0].split()
            self.assertEqual(query_id, "1AMY")
            self.assertAlmostEqual(float(rmsd), hit.rmsd, places=3)
            self.assertAlmostEqual(float(determinant), hit.determinant, places=1)
            self.assertAlmostEqual(float(logE), hit.log_evalue, places=1)

            atom_lines = lines[1:-1]
            atoms = hit.atoms()
            self.assertEqual(len(atoms), len(atom_lines))
            for atom, atom_line in zip(atoms, atom_lines):
                self.assertEqual(atom.serial, int(atom_line[7:12]))
                self.assertEqual(atom.name, atom_line[13:17].strip())
                self.assertEqual(atom.residue_name, atom_line[17:21].strip())
                self.assertEqual(atom.chain_id, atom_line[21:23].strip())
                self.assertEqual(atom.residue_number, int(atom_line[23:27]))
                self.assertAlmostEqual(atom.x, float(atom_line[31:39]), places=3)
                self.assertAlmostEqual(atom.y, float(atom_line[39:47]), places=3)
                self.assertAlmostEqual(atom.z, float(atom_line[47:55]), places=3)
                self.assertAlmostEqual(atom.occupancy, float(atom_line[55:61]), places=3)
                self.assertAlmostEqual(atom.temperature_factor, float(atom_line[61:67]), places=3)

            atoms = hit.atoms(transform=False)
            self.assertEqual(len(atoms), len(atom_lines))
            for atom, atom_line in zip(atoms, atom_lines):
                self.assertEqual(atom.serial, int(atom_line[7:12]))
                self.assertEqual(atom.name, atom_line[13:17].strip())
                self.assertEqual(atom.residue_name, atom_line[17:21].strip())
                self.assertEqual(atom.chain_id, atom_line[21:23].strip())
                self.assertEqual(atom.residue_number, int(atom_line[23:27]))
                self.assertAlmostEqual(atom.occupancy, float(atom_line[55:61]), places=3)
                self.assertAlmostEqual(atom.temperature_factor, float(atom_line[61:67]), places=3)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_mcsa_query_reorder(self):
        with files(data).joinpath("1.3.3.tpl").open() as f:
            template = Template.load(f)
            jess = Jess([template])
        with files(data).joinpath("1AMY.pdb").open() as f:
            molecule = Molecule.load(f)
        with files(data).joinpath("1AMY+1.3.3.txt").open() as f:
            results = list(filter(None, f.read().split("REMARK")))

        hits = list(jess.query(molecule, 2, 4, 4, reorder=True))
        self.assertEqual(len(hits), len(results))

        # `reorder=True` means that we may get results in a different order
        # to Jess, so we need to match the hits in the file by residue number
        # to make sure we compare them consistently.
        
        results_by_serials = {}
        for block in results:
            lines = block.strip().splitlines()
            serials = tuple([ int(line.split()[1]) for line in lines[1:-1] ])
            results_by_serials[serials] = block

        for hit in hits:
            self.assertIs(hit.template, template)
            block = results_by_serials[tuple(atom.serial for atom in hit.atoms(False))]

            lines = block.strip().splitlines()
            query_id, rmsd, template_id, _, determinant, _, logE = lines[0].split()
            self.assertEqual(query_id, "1AMY")
            self.assertAlmostEqual(float(rmsd), hit.rmsd, places=3)
            self.assertAlmostEqual(float(determinant), hit.determinant, places=1)
            self.assertAlmostEqual(float(logE), hit.log_evalue, places=1)

            atom_lines = lines[1:-1]
            atoms = hit.atoms()
            self.assertEqual(len(atoms), len(atom_lines))
            for atom, atom_line in zip(atoms, atom_lines):
                self.assertEqual(atom.serial, int(atom_line[7:12]))
                self.assertEqual(atom.name, atom_line[13:17].strip())
                self.assertEqual(atom.residue_name, atom_line[17:21].strip())
                self.assertEqual(atom.chain_id, atom_line[21:23].strip())
                self.assertEqual(atom.residue_number, int(atom_line[23:27]))
                self.assertAlmostEqual(atom.x, float(atom_line[31:39]), places=3)
                self.assertAlmostEqual(atom.y, float(atom_line[39:47]), places=3)
                self.assertAlmostEqual(atom.z, float(atom_line[47:55]), places=3)
                self.assertAlmostEqual(atom.occupancy, float(atom_line[55:61]), places=3)
                self.assertAlmostEqual(atom.temperature_factor, float(atom_line[61:67]), places=3)

            atoms = hit.atoms(transform=False)
            self.assertEqual(len(atoms), len(atom_lines))
            for atom, atom_line in zip(atoms, atom_lines):
                self.assertEqual(atom.serial, int(atom_line[7:12]))
                self.assertEqual(atom.name, atom_line[13:17].strip())
                self.assertEqual(atom.residue_name, atom_line[17:21].strip())
                self.assertEqual(atom.chain_id, atom_line[21:23].strip())
                self.assertEqual(atom.residue_number, int(atom_line[23:27]))
                self.assertAlmostEqual(atom.occupancy, float(atom_line[55:61]), places=3)
                self.assertAlmostEqual(atom.temperature_factor, float(atom_line[61:67]), places=3)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_ignore_chain_none(self):
        with files(data).joinpath("4.1.2.tpl").open() as f:
            template = Template.load(f)
            jess = Jess([template])
        with files(data).joinpath("5ayx.EF.pdb").open() as f:
            molecule = Molecule.load(f)
        
        hits = list(jess.query(molecule, 3, 3, 3, ignore_chain=None))
        self.assertEqual(len(hits), 0)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_ignore_chain_residues(self):
        with files(data).joinpath("4.1.2.tpl").open() as f:
            template = Template.load(f)
            jess = Jess([template])
        with files(data).joinpath("5ayx.EF.pdb").open() as f:
            molecule = Molecule.load(f)
        
        hits = list(jess.query(molecule, 3, 3, 3, ignore_chain="residues"))
        self.assertEqual(len(hits), 2)
        for hit in hits:
            atoms = hit.atoms(transform=False)
            for i in range(1, len(atoms)):
                if atoms[i-1].residue_number == atoms[i].residue_number:
                    self.assertEqual(atoms[i-1].chain_id, atoms[i].chain_id)

    @unittest.skipUnless(files, "importlib.resources not available")
    def test_ignore_chain_atoms(self):
        with files(data).joinpath("4.1.2.tpl").open() as f:
            template = Template.load(f)
            jess = Jess([template])
        with files(data).joinpath("5ayx.EF.pdb").open() as f:
            molecule = Molecule.load(f)
        
        hits = list(jess.query(molecule, 3, 3, 3, ignore_chain="atoms"))
        self.assertEqual(len(hits), 7)
        
        