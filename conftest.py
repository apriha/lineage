"""
Pytest configuration for testing code examples in README.md using Sybil.

This conftest.py enables Sybil to parse and test Python code blocks in the
README.md file as part of the pytest test suite. The PythonCodeBlockParser
evaluates fenced Python code blocks (```python), while SkipParser allows
selective skipping of examples using Markdown comments when needed.
"""

from sybil import Sybil
from sybil.parsers.markdown import PythonCodeBlockParser, SkipParser

pytest_collect_file = Sybil(
    parsers=[
        PythonCodeBlockParser(),
        SkipParser(),
    ],
    patterns=["README.md"],
).pytest()
