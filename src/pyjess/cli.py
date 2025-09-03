# coding: utf-8

"""The PyJess CLI.
"""
import argparse
import contextlib
import functools
import io
import itertools
import os
import operator
import sys
import warnings
import typing
from pathlib import Path

try:
    import multiprocessing.pool
except ImportError:  # multiprocessing.pool may be missing, e.g. on AWS
    multiprocessing = None

from . import __name__ as prog, __author__, __version__
from ._jess import Template, Molecule, Jess


_BZ2_MAGIC = b"BZh"
_GZIP_MAGIC = b"\x1f\x8b"
_XZ_MAGIC = b"\xfd7zXZ"
_LZ4_MAGIC = b"\x04\x22\x4d\x18"
_ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"


class SmartFormatter(argparse.ArgumentDefaultsHelpFormatter):

    def _format_text(self, text):
        # if '%(prog)' in text:
        #     text = text % dict(prog=self._prog)
        # text_width = max(self._width - self._current_indent, 11)
        # indent = ' ' * self._current_indent
        # return self._fill_text(text, text_width, indent) + '\n\n'
        return text + '\n\n'

    # def _split_lines(self, text, width):
    #     print(text)
    #     if text.startswith('Copyright'):
    #         return text.splitlines()
    #     return super()._split_lines(text, width)
    #     # this is the RawTextHelpFormatter._split_lines
    #     return argparse.HelpFormatter._split_lines(self, text, width)


@contextlib.contextmanager
def zopen(path, mode='r', encoding=None, errors=None, newline=None) -> typing.Iterator[typing.BinaryIO]:
    with contextlib.ExitStack() as ctx:
        file = ctx.enter_context(open(path, "rb"))
        peek = file.peek()
        if peek.startswith(_GZIP_MAGIC):
            import gzip
            file = ctx.enter_context(gzip.open(file, mode="rb"))
        elif peek.startswith(_BZ2_MAGIC):
            import bz2
            file = ctx.enter_context(bz2.open(file, mode="rb"))
        elif peek.startswith(_XZ_MAGIC):
            import lzma
            file = ctx.enter_context(lzma.open(file, mode="rb"))
        elif peek.startswith(_LZ4_MAGIC):
            try:
                import lz4.frame
            except ImportError as err:
                raise RuntimeError("File compression is LZ4 but lz4 is not installed") from err
            file = ctx.enter_context(lz4.frame.open(file))
        elif peek.startswith(_ZSTD_MAGIC):
            try:
                import zstandard
            except ImportError as err:
                raise RuntimeError("File compression is ZSTD but zstandard is not installed") from err
            decompressor = zstandard.ZstdDecompressor()
            file = decompressor.stream_reader(file)
        if mode == "r":
            file = io.TextIOWrapper(file, encoding=encoding, errors=errors, newline=newline)
        yield file


def argument_parser(
    prog: str = prog,
    version: str = __version__,
    # formatter_class: argparse.HelpFormatter = argparse.ArgumentDefaultsHelpFormatter,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog, 
        add_help=False,
        formatter_class=SmartFormatter,
        description=(
            "PyJess - Optimized Python bindings to Jess, a 3D template matching software.\n\n"
            "MIT License\n\n"
            "Copyright (c) 2025 Martin Larralde <martin.larralde@embl.de>\n"
            "Copyright (c) 2002 Jonathan Barker <jbarker@ebi.ac.uk>\n\n"
        ),
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit."
    )
    parser.add_argument(
        "-V",
        "--version",
        help="Show the version number and exit.",
        action="version",
        version=f"PyJess {__version__}",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        help="The number of jobs to use for multithreading.",
        type=int,
        default=os.cpu_count() or 1,
    )

    group = parser.add_argument_group("Mandatory Parameters")
    group.add_argument(
        "-T",
        "--templates",
        help="The path to the template list file.",
        type=Path,
        required=True,
    )
    group.add_argument(
        "-Q",
        "--queries",
        help="The path to the query list file.",
        type=Path,
        required=True,
    )
    group.add_argument(
        "-R",
        "--rmsd",
        help="The RMSD threshold.",
        type=float,
        required=True,
    )
    group.add_argument(
        "-D",
        "--distance-cutoff",
        help="The distance-cutoff.",
        type=float,
        required=True,
    )
    group.add_argument(
        "-M",
        "--maximum-distance",
        help=(
            "The maximum allowed template/query atom distance after adding the "
            "global distance cutoff and the individual atom distance cutoff "
            "defined in the temperature field of the ATOM record in the "
            "template file."
        ),
        type=float,
        required=True,
    )
    
    group = parser.add_argument_group("Flags")
    group.add_argument(
        "-n",
        "--no-transform",
        help="Do not transform coordinates of hit into the template coordinate frame",
        action="store_false",
        dest="transform"
    )
    group.add_argument(
        "-f",
        "--filenames",
        help="Show PDB filenames in progress on stderr",
        action="store_true"
    )
    group.add_argument(
        "-i",
        "--ignore-chain",
        help=(
            "Include matches composed of residues belonging to multiple chains "
            "(if template is single-chain) or matches with residues from a single "
            "chain (if template has residues from multiple chains)."
        ),
        action="store_true",
    )
    group.add_argument(
        "--ignore-res-chain",
        help=(
            "Include matches composed of residues belonging to multiple chains "
            "but still enforce all atoms of a residue to be part of the same chain."
        ),
        action="store_true",
    )
    group.add_argument(
        "-q",
        "--query-filename",
        help="Write filename of query instead of PDB ID from HEADER",
        action="store_true",
    )
    group.add_argument(
        "-e",
        "--ignore-endmdl",
        help="Parse atoms from all models separated by ENDMDL (use with care).",
        action="store_true",
    )
    group.add_argument(
        "-c",
        "--max-candidates",
        help="Set a maximum number of candidates to return by template.",
        type=int,
        default=None,
    )
    group.add_argument(
        "--no-reorder",
        help=(
            "Disable template atom reordering in the matching process, useful "
            "to enforce results to be returned exactly in the same order as "
            "the original Jess, at the cost of longer runtimes."
        ),
        action="store_false",
        dest="reorder",
    )
    group.add_argument(
        "-b",
        "--best-match",
        help="Return only the best match for each template/query pair.",
        action="store_true",
    )

    return parser


def _process(gene_finder, sequence):
    if not sequence.id:
        warnings.warn("Input file contains a sequence without identifier", stacklevel=2)
    return sequence.id, gene_finder.find_genes(sequence.seq)


def main(
    argv: typing.Optional[typing.List[str]] = None,
    stdout: typing.TextIO = sys.stdout,
    stderr: typing.TextIO = sys.stderr,
) -> int:
    parser = argument_parser()
    args = parser.parse_args(argv)

    ignore_chain = "all" if args.ignore_chain else "residues" if args.ignore_res_chain else None

    with contextlib.ExitStack() as ctx:
        try:
            
            with args.templates.open() as f:
                templates = [Template.load(n, id=n) for n in map(str.strip, f)]
                jess = Jess(templates)

            with args.queries.open() as f:
                for filename in map(str.strip, f):
                    id_ = filename if args.query_filename else None
                    mol = Molecule.load(filename, id=id_, ignore_endmdl=args.ignore_endmdl)
                    if args.filenames:
                        print(filename, file=stderr)
                    query = jess.query(
                        mol, 
                        args.rmsd, 
                        args.distance_cutoff, 
                        args.maximum_distance, 
                        max_candidates=args.max_candidates,
                        ignore_chain=ignore_chain,
                        reorder=args.reorder,
                        best_match=args.best_match,
                    )
                    for hit in query:
                        print(f"REMARK {mol.id} {hit.rmsd:5.3f} {hit.template.id} Det={hit.determinant:4,.1f} log(E)~ {hit.log_evalue:4.2f}", file=stdout)
                        for atom in hit.atoms(transform=args.transform):
                            name = atom.name
                            name = f"  {name:<3}" if len(name) < 4 else name
                            print(
                                f"ATOM  {atom.serial:5}{name}{atom.altloc}{atom.residue_name:<3}{atom.chain_id:>2}{atom.residue_number:4}{atom.insertion_code:4}{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}{atom.occupancy:6.2f}{atom.temperature_factor:6.2f}",
                                file=stdout,
                            )
                        print("ENDMDL\n", file=stdout)

        except Exception as err:
            print("Error: {}".format(err))
            return getattr(err, "errno", 1)
        else:
            return 0