"""Make the repository root importable so ``import qudit_qec`` works under pytest
regardless of the working directory or invocation style."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
