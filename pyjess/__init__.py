# noqa: D104

__version__ = "0.3.0"
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__all__ = ["Atom", "Hit", "Jess", "Query", "Molecule", "Template", "TemplateAtom"]

from . import _jess
from ._jess import Atom, Hit, Jess, Query, Molecule, Template, TemplateAtom

__doc__ = _jess.__doc__

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pyjess.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )
