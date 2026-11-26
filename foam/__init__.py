"""FOAM helper library.

Everything in this package is plain standard-library Python and must stay
compatible with Python 3.6, because the annotation step runs inside the
MITOS 2.1.3 conda environment (which pins Python 3.6).  In practice that means:
no dataclasses, no walrus operator, no ``subprocess.run(capture_output=...)``,
no builtin-generic type hints (``list[int]``) and no ``from __future__ import
annotations``.
"""

__version__ = "0.8.0"
