import abc
import datetime
import os
import typing
import warnings
from typing import Optional, Type, TextIO, Union

from ._jess import Atom, Molecule

_M = typing.TypeVar("_M", bound=Molecule)


class nullcontext:
    def __init__(self, return_value=None):
        self.retval = return_value
    def __enter__(self):
        return self.retval
    def __exit__(self, exc_type, exc_value, traceback):
        return False


class MoleculeParser(abc.ABC):
    id: Optional[str]
    date: Optional[datetime.date]
    name: Optional[str]

    def __init__(
        self, 
        id: Optional[str] = None, 
        date: Optional[datetime.date] = None, 
        name: Optional[str] = None,
    ):
        self.id = id
        self.date = date
        self.name = name

    def loads(self, text: str, molecule_type: Type[_M]) -> _M:
        return self.load(io.StringIO(text), molecule_type)

    @abc.abstractmethod
    def load(self, file: Union[TextIO, "os.PathLike[str]"], molecule_type: Type[_M]) -> _M:
        raise NotImplementedError


class PDBMoleculeParser(MoleculeParser):
    ignore_endmdl: bool
    skip_hetatm: bool

    def __init__(
        self,
        id: Optional[str] = None,
        date: Optional[datetime.date] = None,
        name: Optional[str] = None,
        ignore_endmdl: bool = False,
        skip_hetatm: bool = False
    ):
        super().__init__(id=id, date=date, name=name)
        self.ignore_endmdl = ignore_endmdl
        self.skip_hetatm = skip_hetatm

    def load(self, file: Union[TextIO, "os.PathLike[str]"], molecule_type: Type[_M]) -> _M:
        # cdef str    line
        id_ = self.id
        date = self.date
        name = self.name
        atoms = []

        try:
            handle = open(file)
        except TypeError:
            handle = nullcontext(file)
        with handle as f:
            for line in f:
                if line.startswith("HEADER"):
                    if id_ is None:
                        id_ = line[62:66].strip() or None
                    if date is None:
                        date_text = line[50:59].strip()
                        if date_text:
                            date = datetime.datetime.strptime(date_text, "%d-%b-%y").date()
                    if name is None:
                        name = line[10:50].strip() or None
                elif line.startswith("ATOM"):
                    atoms.append(Atom.loads(line))
                elif line.startswith("HETATM") and not self.skip_hetatm:
                    atoms.append(Atom.loads(line))
                elif line.startswith("ENDMDL"):
                    if not self.ignore_endmdl:
                        break
                elif line.lower().startswith(("data_", "loop_")):
                    raise ValueError("mmCIF data tags found, file is not in PDB format")
        return molecule_type(atoms, id=id_, date=date, name=name)


class CIFMoleculeParser(MoleculeParser):
    use_author: bool
    skip_hetatm: bool
    ignore_endmdl: bool

    _PRIMARY_COLUMNS = [
        'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id',
        'label_asym_id', 'label_seq_id', '?pdbx_PDB_ins_code', 'Cartn_x',
        'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv',
        '?pdbx_formal_charge', '?group_PDB', 'pdbx_PDB_model_num',
    ]

    _AUTH_COLUMNS = [
        'id', 'type_symbol', 'auth_atom_id', 'label_alt_id', 'auth_comp_id',
        'auth_asym_id', 'auth_seq_id', '?pdbx_PDB_ins_code', 'Cartn_x',
        'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv',
        '?pdbx_formal_charge', '?group_PDB', 'pdbx_PDB_model_num',
    ]

    def __init__(
        self,
        id: Optional[str] = None,
        date: Optional[datetime.date] = None,
        name: Optional[str] = None,
        use_author: bool = False,
        skip_hetatm: bool = False,
        ignore_endmdl: bool = False,
    ):
        super().__init__(id=id, date=date, name=name)
        self.gemmi = __import__('gemmi')
        self.use_author = use_author
        self.skip_hetatm = skip_hetatm
        self.ignore_endmdl = ignore_endmdl

    def _load_block(self, document, molecule_type):
        block = document.sole_block()
        cols = self._AUTH_COLUMNS if self.use_author else self._PRIMARY_COLUMNS
        table = block.find('_atom_site.', cols)
        max_residue_number = 0

        if not table:
            raise ValueError("missing columns in CIF files")

        atoms = []
        for row in table:

            # row[15] contains _atom_site.pdbx_PDB_model_num
            # by default (if ignore_endmdl is False) we break on model number
            if not self.ignore_endmdl and row[15] != '1':
                break

            if row[14] != "ATOM" and (row[14] != "HETATM" or self.skip_hetatm):
                continue

            if row[6] == "." and row[14] == "HETATM":
                warnings.warn(
                    "HETATM line found without residue number. Consider "
                    "parsing with use_author=True to use author-defined "
                    "residue numbers, or skip_hetatm=True to disable "
                    "parsing of HETATM altogether.",
                    UserWarning,
                    stacklevel=3,
                )
                residue_number = max_residue_number
                max_residue_number += 1
            else:
                residue_number = int(row[6])
                max_residue_number = max(residue_number, max_residue_number)

            atom = Atom(
                serial=int(row[0]),
                element=row[1],
                name=row[2].strip('"'),
                altloc=' ' if row[3] == "." else row[3], # FIXME: replace with None?
                residue_name=row[4],
                chain_id=row[5],
                residue_number=residue_number,
                insertion_code=' ' if not row.has(7) or row[7] == "?" else row[7],
                x=float(row[8]),
                y=float(row[9]),
                z=float(row[10]),
                occupancy=0.0 if row[11] == '.' else float(row[11]),
                temperature_factor=0.0 if row[12] == '.' else float(row[12]),
                charge=0 if not row.has(13) or row[13] == "?" else int(row[13]),
            )
            atoms.append(atom)

        id = block.name if self.id is None else self.id

        entry_id = block.find_value("_entry.id")
        pdb_kwds = block.find_value("_struct_keywords.pdbx_keywords")
        title    = block.find_value("_struct.title")

        if pdb_kwds:
            name = pdb_kwds.strip("'")
        elif title:
            name = title.strip("'")
        elif entry_id:
            name = entry_id.strip("'")
        else:
            name = None

        date_tbl = block.find('_pdbx_audit_revision_history.', ["revision_date"])
        if not date_tbl:
            date=None
        else:
            # take earliest date as deposition date
            date = min(datetime.datetime.strptime(row[0], "%Y-%m-%d") for row in date_tbl)

        return molecule_type(atoms, id=id, date=date, name=name)

    def loads(self, text: str, molecule_type: Type[_M]) -> _M:
        document = self.gemmi.cif.read_string(text)
        return self._load_block(document, molecule_type)

    def load(self, file: Union[TextIO, "os.PathLike[str]"], molecule_type: Type[_M]) -> _M:
        if hasattr(file, "read"):
            document = self.gemmi.cif.read_string(file.read())
        else:
            document = self.gemmi.cif.read_file(file)
        return self._load_block(document, molecule_type)

