"""``python -m nurbforge`` prints the compact API reference for language models and agents."""
import sys

from ._reference import reference

sys.stdout.write(reference())
