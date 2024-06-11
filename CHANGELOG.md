# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyjess/compare/v0.3.0...HEAD


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
