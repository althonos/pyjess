PyJess |Stars|
==============

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyjess.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyjess/stargazers
   :class: dark-light

*A Python interface to Jess, a 3D template matching software.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyjess/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyjess/actions
   :class: dark-light

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyjess?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyjess/
   :class: dark-light

.. |PyPI| image:: https://img.shields.io/pypi/v/pyjess.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyjess
   :class: dark-light

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyjess?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyjess
   :class: dark-light

.. |AUR| image:: https://img.shields.io/aur/version/python-pyjess?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyjess
   :class: dark-light

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyjess?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyjess/#files
   :class: dark-light

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyjess.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyjess/#files
   :class: dark-light

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyjess.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyjess/#files
   :class: dark-light

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/
   :class: dark-light

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pyjess/
   :class: dark-light

.. |Mirror| image:: https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=3600&style=flat-square
   :target:https://git.lumc.nl/mflarralde/pyjess/
   :class: dark-light

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyjess.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyjess/issues
   :class: dark-light

.. |Docs| image:: https://img.shields.io/readthedocs/pyjess?style=flat-square&maxAge=3600
   :target: http://pyjess.readthedocs.io/en/stable/?badge=stable
   :class: dark-light

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pyjess/blob/main/CHANGELOG.md
   :class: dark-light

.. |Downloads| image:: https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fpepy.tech%2Fprojects%2Fpyjess&search=%5B0-9%5D%2B.%5B0-9%5D%2B(k%7CM)&style=flat-square&label=downloads&color=303f9f&cacheSeconds=86400
   :target: https://pepy.tech/project/pyjess
   :class: dark-light


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
`Cython <https://cython.org/>`_:

.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Just add ``pyjess`` as a ``pip`` or ``conda`` dependency, no need
      for the Jess binary or any external dependency.

   .. grid-item-card:: :fas:`screwdriver-wrench` Flexible

      Load a `~pyjess.Molecule` from a :wiki:`.pdb` file or create 
      it programmatically through the :doc:`Python API <api/index>`.

   .. grid-item-card:: :fas:`gears` Practical

      Retrieve results as they become availble as dedicated 
      `~pyjess.Hit` objects, and compute statistics on-the-fly.

   .. grid-item-card:: :fas:`gauge-high` Fast

      Compute matches about 10x faster, thanks to
      several :doc:`algorithmic optimizations <guide/optimizations>`
      made to the original Jess code.

   .. grid-item-card:: :fas:`server` Parallel

      Easily run computations in parallel querying thread-safe 
      `~pyjess.Jess` with several `~pyjess.Molecule` in parallel.

   .. grid-item-card:: :fas:`toolbox` Feature-complete

      Access all the features of the original CLI through the 
      :doc:`Python API <api/index>`, including atom coordinate 
      transformation.


Setup
-----

PyJess is available for all modern Python versions (3.7+).

Run ``pip install pyjess`` in a shell to download the latest release from PyPI,
or have a look at the :doc:`Installation page <guide/install>` to find other ways 
to install PyJess.

Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst


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
