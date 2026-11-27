#!/usr/bin/env python3
"""FOAM command-line tools (see ``foam.py --help``)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from foam.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
