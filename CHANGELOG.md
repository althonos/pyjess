# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyjess/compare/v0.7.0...HEAD


## [v0.8.0] - 2025-11-12
[v0.8.0]: https://github.com/althonos/pyjess/compare/v0.7.0...v0.8.0

### Added
- Testing for Python 3.14.
- `ignore_endmdl` parameter value in mmCIF parser to ignore models past the first one.
- `dump` and `dumps` methods to `Template` and `TemplateAtom` ([#5](https://github.com/althonos/pyjess/issues/5)).

### Fixed
- Parsing issues with CIF files (by [@RayHackett](https://github.com/RayHackett/)).
- Incorrect example use of the `pyjess.API` in `README.md` ([#6](https://github.com/althonos/pyjess/issues/6)).
- Broken format detection for PDB or mmCIF files starting with comment lines (and fix handling of non-seekable files).


## [v0.7.0] - 2025-09-13
[v0.7.0]: https://github.com/althonos/pyjess/compare/v0.6.0...v0.7.0

### Added
- Command Line Interface with the same features as the Jess CLI.
- `Molecule.from_biopython`, `Molecule.from_biotite` and `Molecule.from_gemmi` to create a `Molecule` with objects from other libraries.
- `Hit.dump` and `Hit.dumps` methods to write a hit to PDB format.
- `Molecule.load` option to skip parsing of `HETATM` from PDB files.
- Support for [Selenocysteine](https://en.wikipedia.org/wiki/Selenocysteine) and [Pyrrolysine](https://en.wikipedia.org/wiki/Pyrrolysine) residues in `TemplateAtom`.
- Support for parsing `Molecule` objects from CIF files using [`gemmi`](https://gemmi.readthedocs.io/).

### Fixed
- **breaking**: Incorrect handling of `max_candidates` in `Jess.query`, causing PyJess to erroneously ignore some templates.
- `Template.dimension` reporting incorrect numbers for identical residues across different residues.
- Issue with alignment of atom names in `TemplateAtom`.
- Validate data in PDB parser to ensure it is not parsing a mmCIF file.
- Parsing of mmCIF files with atoms missing occupancy values.
- Incorrect alignment of `TemplateAtom` names for names of less than 4 characters.
- `max_candidates` causing `Query` to stop before reaching the actual number of maximum candidates.

### Changed
- **breaking**: Set the `max_candidates` default value to `None` in `Jess.query`, disabling max candidates filtering by default.
- **breaking**: Use string variants instead of `bool` to control the behaviour of `ignore_chain` argument.
- **breaking**: Make `id` and `ignore_endmdl` arguments of `Molecule.load` and `Molecule.loads` keyword-only.
- Use unrolled string comparison to compare atom names instead of `strcasecmp` in Jess code.
- Make `altloc` and `insertion_code` arguments of `Atom` optional.
- Setup deployment of Limited API wheels for CPython 3.8 onwards.
- Implement format detection between PDB and mmCIF in `Molecule.load`.
- Implement `Template.copy` and `TemplateAtom.copy` using C-level API rather than Python copy.


## [v0.6.0] - 2025-09-01
[v0.6.0]: https://github.com/althonos/pyjess/compare/v0.5.2...v0.6.0

### Added
- Several [algorithmic optimizations](https://pyjess.readthedocs.io/en/v0.6.0/guide/optimizations.html) to Jess, greatly improving runtime:
  - Use QuickSelect ($O(n)$) instead of QuickSort ($O(nlog(n))$) to select the medians used for space partitioning on k-d tree initialization.
  - Use approximate intersection on bounding boxes rather than exact code involving multiple annuli when querying the k-d tree for atoms.
  - **breaking**: Reorder the matching order of template atoms to reduce the amount of backtracking and k-d tree querying performed in the average and worst case.
- `reorder` argument to `Jess.query` to support disabling template atom reordering if needed for Jess 1-to-1 consistency.

### Changed
- Hardcode space dimensions to 3 to encourage compilers to unroll loops over dimensions.
- Recycle memory between templates within a query to reduce total amount of allocation/deallocation in hot paths.


## [v0.5.2] - 2025-08-26
[v0.5.2]: https://github.com/althonos/pyjess/compare/v0.5.1...v0.5.2

### Added
- `pickle` protocol support to `Hit` objects.

### Fixed
- Improve performance of residue index using binary search on name queries.
- Avoid copying atoms in `Hit.atoms` if no transformation is required.


## [v0.5.1] - 2025-08-20
[v0.5.1]: https://github.com/althonos/pyjess/compare/v0.5.0...v0.5.1

### Fixed
- Typo in `Hit.molecule` returning incorrect coordinates when transformation is required.

### Changed
- Implement an index to quickly access atoms by residue names in a `Molecule`, and accelerate initial candidate set creation for compatible match modes.


## [v0.5.0] - 2025-04-09
[v0.5.0]: https://github.com/althonos/pyjess/compare/v0.4.1...v0.5.0

### Fixed
- Remove invalid `platform` key from `pyproject.toml`.

### Changed
- `Hit.molecule` is now a method and can optionally rotate the hit molecule.


## [v0.4.1] - 2024-11-29
[v0.4.1]: https://github.com/althonos/pyjess/compare/v0.4.0...v0.4.1

### Fixed
- `TemplateAtom.__hash__` implementation.

### Changed
- `TemplateAtom.residue_names` and `TemplateAtom.atom_names` are now tuples.


## [v0.4.0] - 2024-11-29
[v0.4.0]: https://github.com/althonos/pyjess/compare/v0.3.3...v0.4.0

### Added
- `__eq__`, `__hash__`, `__copy__`, and `__reduce__` to `Molecule`, `Atom`, `Template`, `TemplateAtom` and `Jess`.

### Changed
- Move project development stage from *Alpha* to *Beta*. 
- List `conda` support in the [Install page](https://pyjess.readthedocs.io/en/latest/guide/install.html) of the documentation.


## [v0.3.3] - 2024-11-05
[v0.3.3]: https://github.com/althonos/pyjess/compare/v0.3.2...v0.3.3

### Fixed
- Tests crashing on missing `importlib.resources.files` function.


## [v0.3.2] - 2024-11-05
[v0.3.2]: https://github.com/althonos/pyjess/compare/v0.3.1...v0.3.2

### Added
- Support for Python 3.13.

### Changed
- Use `scikit-build-core` to build the project.
- Use CMake to detect availability of `qsort_r` / `qsort_s` on the current platform.
- Use thread-local storage to support parallelism on platforms without `qsort_r`.

### Removed
- Support for Python 3.6.


## [v0.3.1] - 2024-07-18
[v0.3.1]: https://github.com/althonos/pyjess/compare/v0.3.0...v0.3.1

### Changed
- Migrate documentation to `pydata-sphinx-theme`.
- Make `Hit.evalue` and `Hit.log_evalue` release the GIL when possible.

### Fixed
- Signatures of `__init__` methods missing from all Cython types after the `v3.0` update.


## [v0.3.0] - 2024-06-11
[v0.3.0]: https://github.com/althonos/pyjess/compare/v0.2.1...v0.3.0

### Added
- Slicing for `Template`, `Molecule` and `Jess` objects.

### Fixed
- Typing for `chain_id` property of `TemplateAtom`.
- `Query.__next__` ignores and raises a warning on `NaN` matrices caused by planar molecule coordinates.

### Changed
- Make `Jess` generic over the internal template class to allow `Template` subclasses as inputs and `Hit` outputs attribute.
- Make `Molecule.conserved` return an instance of the caller class rather than a `Molecule` object.


## [v0.2.1] - 2024-05-30
[v0.2.1]: https://github.com/althonos/pyjess/compare/v0.2.0...v0.2.1

### Fixed
- Type hints for `Molecule.load` and `Template.load` not marked as accepting paths.
- Thread-unsafe use of `qsort` in Jess code, replaced with `qsort_r` or `qsort_s` to allow multithreading.


## [v0.2.0] - 2024-05-24
[v0.2.0]: https://github.com/althonos/pyjess/compare/v0.1.1...v0.2.0

### Added
- `best_match` argument to `Jess.query` to only return the best match to each template for a query molecule.
- `id` argument to `Molecule.load` and `Molecule.loads` to allow overriding the PDB ID stored in the file.
- `ignore_endmdl` argument to `Molecule.load` and `Molecule.loads` to control whether the parser should stop at the first model in a file.

### Changed
- Make `Molecule.load` stop at the first model when parsing a PDB file.

### Fixed
- Invalid pointer assignments in `TessAtom.c` causing compilation errors on stricter compilers.


## [v0.1.1] - 2024-04-18
[v0.1.1]: https://github.com/althonos/pyjess/compare/v0.1.0...v0.1.1

### Added
- Support for passing filenames to `Template.load` and `Molecule.load` directly.
- Default initialization of empty `Jess` objects from an empty iterable.
- Make `Jess` implement the `Sized` abstract base class interface.

### Changed
- Skip displaying default attribute values in `repr` implementation of `Atom` and `TemplateAtom`.
- Make `Jess.query` optional parameters keyword-only.


## [v0.1.0] - 2024-04-18
[v0.1.0]: https://github.com/althonos/pyjess/compare/3f2a7e9...v0.1.0

Initial release.
