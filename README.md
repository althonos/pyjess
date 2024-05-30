# 🐍🔍 PyJess [![Stars](https://img.shields.io/github/stars/althonos/pyjess.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyjess/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [Jess](https://github.com/iriziotis/jess), a 3D template matching software.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyjess/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyjess/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyjess?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pyjess/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![PyPI](https://img.shields.io/pypi/v/pyjess.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyjess)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyjess?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyjess)
[![AUR](https://img.shields.io/aur/version/python-pyjess?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyjess)
[![Wheel](https://img.shields.io/pypi/wheel/pyjess.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyjess/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyjess.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyjess/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyjess.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyjess/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyjess/)
[![Mirror](https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=2678400&style=flat-square)](https://git.lumc.nl/mflarralde/pyjess/)
[![Issues](https://img.shields.io/github/issues/althonos/pyjess.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyjess/issues)
[![Docs](https://img.shields.io/readthedocs/pyjess/latest?style=flat-square&maxAge=600)](https://pyjess.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyjess/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyjess?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyjess)


## 🗺️ Overview

Jess is an algorithm for constraint-based structural template matching
proposed by Jonathan Barker *et al.*[\[1\]](#ref1). It can be used to identify
catalytic residues from a known template inside a protein structure. Jess
is an evolution of TESS, a geometric hashing algorithm developed by
Andrew Wallace *et al.*[\[2\]](#ref2), removing some pre-computation and 
structural requirements from the original algorithm. Jess was further 
updated and maintained by [Ioannis Riziotis](https://github.com/iriziotis)
during his PhD in the [Thornton group](https://www.ebi.ac.uk/research/thornton/).

PyJess is a Python module that provides bindings to Jess using
[Cython](https://cython.org/). It allows creating templates, querying them
with protein structures, and retrieving the hits using a Python API without
performing any external I/O.


## 🔧 Installing

PyJess is available for all modern Python versions (3.6+).

It can be installed directly from [PyPI](https://pypi.org/project/pyjess/),
which hosts some pre-built x86-64 wheels for Linux, MacOS, and Windows,
as well as the code required to compile from source with Cython:
```console
$ pip install pyjess
```

<!-- Otherwise, PyJess is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda pyjess
``` -->

Check the [*install* page](https://pyjess.readthedocs.io/en/stable/install.html)
of the documentation for other ways to install PyJess on your machine.

## 💡 Example

Load templates to be used as references from different template files:

```python
import glob
import pyjess

templates = []
for path in sorted(glob.iglob("vendor/jess/examples/template_*.qry")):
    templates.append(Template.load(path, id=os.path.basename(path)))
```

Create a `Jess` instance and use it to query a molecule (a PDB structure)
against the stored templates:

```python
jess = Jess(templates)
mol = Molecule("vendor/jess/examples/test_pdbs/pdb1a0p.ent")
query = jess.query(mol, rmsd_threshold=2.0, distance_cutoff=3.0, max_dynamic_distance=3.0)
```

The hits are computed iteratively, and the different output statistics are
computed on-the-fly when requested:

```python
for hit in query:
    print(hit.molecule.id, hit.template.id, hit.rmsd, hit.log_evalue)
    for atom in hit.atoms():
        print(atom.name, atom.x, atom.y, atom.z)
```


## 🧶 Thread-safety

Once a `Jess` instance has been created, the templates cannot be edited anymore,
making the `Jess.query` method re-entrant. This allows querying several
molecules against the same templates in parallel using a thread pool:

```python
molecules = []
for path in glob.glob("vendor/jess/examples/test_pdbs/*.ent"):
    molecules.append(Molecule.load(path))

with multiprocessing.ThreadPool() as pool:
    hits = pool.map(jess.query, molecules)
```

*⚠️ Prior to PyJess `v0.2.1`, the Jess code was running some thread-unsafe operations which have now been patched. 
If running Jess in parallel, make sure to use `v0.2.1` or later to use the code patched with re-entrant functions*.

<!-- ## ⏱️ Benchmarks -->


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue tracker](https://github.com/althonos/[pyjess]/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.


### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyjess/blob/main/CONTRIBUTING.md)
for more details.


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pyjess/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/). The JESS code is distributed under the [MIT License](https://choosealicense.com/licenses/mit/) as well.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the JESS authors. It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*


## 📚 References

- <a id="ref1">\[1\]</a> Barker, J. A., & Thornton, J. M. (2003). An algorithm for constraint-based structural template matching: application to 3D templates with statistical analysis. Bioinformatics (Oxford, England), 19(13), 1644–1649. [doi:10.1093/bioinformatics/btg226](https://doi.org/10.1093/bioinformatics/btg226).
- <a id="ref2">\[2\]</a> Wallace, A. C., Borkakoti, N., & Thornton, J. M. (1997). TESS: a geometric hashing algorithm for deriving 3D coordinate templates for searching structural databases. Application to enzyme active sites. Protein science : a publication of the Protein Society, 6(11), 2308–2323. [doi:10.1002/pro.5560061104](https://doi.org/10.1002/pro.5560061104).
