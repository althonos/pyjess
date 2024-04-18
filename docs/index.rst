PyJess |Stars|
==============

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyjess.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyjess/stargazers

*A Python interface to Jess, a 3D template matching software.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyjess/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyjess/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyjess?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyjess/

.. |PyPI| image:: https://img.shields.io/pypi/v/pyjess.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyjess

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyjess?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyjess

.. |AUR| image:: https://img.shields.io/aur/version/python-pyjess?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyjess

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyjess?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyjess/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyjess.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyjess/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyjess.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyjess/#files

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyjess/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyjess.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyjess/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pyjess?style=flat-square&maxAge=3600
   :target: http://pyjess.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyjess/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyjess?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pyjess


Overview
--------

Jess is an algorithm for constraint-based structural template matching
proposed by Jonathan Barker & Janet Thornton. It can be used to identify
catalytic residues from a known template inside a protein structure. Jess
is an evolution of TESS, a geometric hashing algorithm developed by
Andrew Wallace *et al.*, removing some pre-computation and 
structural requirements from the original algorithm. Jess was further 
updated and maintained by `Ioannis Riziotis <https://github.com/iriziotis>`_
during his PhD in the `Thornton group <https://www.ebi.ac.uk/research/thornton/>`_.

PyJess is a Python module that provides bindings to Jess using
`Cython <https://cython.org/>`_. It allows creating templates, querying them
with protein structures, and retrieving the hits using a Python API without
performing any external I/O.

PyJess is available for all modern Python versions (3.6+).


Setup
-----

Run ``pip install pyjess`` in a shell to download the latest release from PyPI,
or have a look at the :doc:`Installation page <install>` to find other ways 
to install ``pyjess`.


Library
-------

.. toctree::
   :maxdepth: 2

   Installation <install>
   Contributing <contributing>
   API Reference <api/index>
   Changelog <changes>


License
-------

This library is provided under the `MIT License <https://choosealicense.com/licenses/mit/>`_.
The JESS code is distributed under the 
`MIT License <https://choosealicense.com/licenses/mit/>`_ as well.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the original JESS authors. It was developed
by* `Martin Larralde <https://github.com/althonos/>`_ *during his PhD project
at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_ *in
the* `Zeller team <https://github.com/zellerlab>`_.
