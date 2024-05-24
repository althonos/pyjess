import os
from typing import Any, Union, Optional, Sequence, Iterator, Iterable, List, TextIO, Sized

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

MATCH_MODE = Literal[
    -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 100, 101, 102, 103, 104, 105, 106, 107
]

class Molecule(Sequence[Atom]):
    @classmethod
    def load(cls, file: TextIO, id: Optional[str] = None, ignore_endmdl: bool = False) -> Molecule: ...
    @classmethod
    def loads(cls, text: str, id: Optional[str] = None, ignore_endmdl: bool = False) -> Molecule: ...
    def __init__(self, atoms: Sequence[Atom] = (), id: Optional[str] = None): ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> Atom: ...
    @property
    def id(self) -> Optional[str]: ...
    def conserved(self, cutoff: float = 0.0) -> Molecule: ...

class Atom:
    @classmethod
    def load(cls, file: Union[str, bytes, os.PathLike[Any], TextIO]) -> Atom: ...
    @classmethod
    def loads(cls, text: str) -> Atom: ...
    def __init__(
        self,
        *,
        serial: int,
        name: str,
        altloc: str,
        residue_name: str,
        chain_id: str,
        residue_number: int,
        insertion_code: str,
        x: float,
        y: float,
        z: float,
        occupancy: float = 0.0,
        temperature_factor: float = 0.0,
        segment: str = "",
        element: str = "",
        charge: int = 0,
    ): ...
    def __repr__(self) -> str: ...
    @property
    def serial(self) -> int: ...
    @property
    def altloc(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def residue_name(self) -> str: ...
    @property
    def residue_number(self) -> int: ...
    @property
    def segment(self) -> str: ...
    @property
    def element(self) -> str: ...
    @property
    def insertion_code(self) -> str: ...
    @property
    def chain_id(self) -> str: ...
    @property
    def occupancy(self) -> float: ...
    @property
    def temperature_factor(self) -> float: ...
    @property
    def charge(self) -> int: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def z(self) -> float: ...

class TemplateAtom:
    @classmethod
    def load(cls, file: Union[str, bytes, os.PathLike[Any], TextIO]) -> Atom: ...
    @classmethod
    def loads(cls, text: str) -> TemplateAtom: ...
    def __init__(
        self,
        *,
        chain_id: str,
        residue_number: int,
        x: float,
        y: float,
        z: float,
        residue_names: Sequence[str],
        atom_names: Sequence[str],
        distance_weight: float = 0.0,
        match_mode: MATCH_MODE = 0,
    ): ...
    def __repr__(self) -> str: ...
    @property
    def match_mode(self) -> MATCH_MODE: ...
    @property
    def residue_number(self) -> int: ...
    @property
    def chain_id(self) -> int: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def z(self) -> float: ...
    @property
    def atom_names(self) -> List[str]: ...
    @property
    def residue_names(self) -> List[str]: ...
    @property
    def distance_weight(self) -> float: ...

class Template(Sequence[TemplateAtom]):
    @classmethod
    def load(cls, file: TextIO, id: Optional[str] = None) -> Template: ...
    @classmethod
    def loads(cls, text: str, id: Optional[str] = None) -> Template: ...
    def __init__(self, atoms: Sequence[TemplateAtom] = (), id: Optional[str] = None): ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> TemplateAtom: ...
    @property
    def id(self) -> str: ...
    @property
    def dimension(self) -> int: ...

class Query(Iterator[Hit]):
    @property
    def jess(self) -> Jess: ...
    @property
    def molecule(self) -> Molecule: ...
    @property
    def ignore_chain(self) -> bool: ...
    @property
    def rmsd_threshold(self) -> float: ...
    @property
    def max_candidates(self) -> int: ...
    @property
    def best_match(self) -> bool: ...
    def __iter__(self) -> Query: ...
    def __next__(self) -> Hit: ...

class Hit:
    @property
    def rmsd(self) -> float: ...
    @property
    def template(self) -> Template: ...
    @property
    def molecule(self) -> Molecule: ...
    @property
    def determinant(self) -> float: ...
    @property
    def log_evalue(self) -> float: ...
    @property
    def evalue(self) -> float: ...
    def atoms(self, transform: bool = True) -> List[Atom]: ...

class Jess(Sized):
    def __init__(self, templates: Iterable[Template] = ()): ...
    def __len__(self) -> int: ...
    def query(
        self,
        molecule: Molecule,
        rmsd_threshold: float,
        distance_cutoff: float,
        max_dynamic_distance: float,
        *,
        max_candidates: int = 1000,
        ignore_chain: bool = False,
        best_match: bool = False,
    ) -> Query: ...
