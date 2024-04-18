# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyjess/compare/v0.1.1...HEAD


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
