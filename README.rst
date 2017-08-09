lineage
=======
``lineage`` provides a framework for analyzing raw DNA data provided by
consumer DNA testing companies (e.g., `23andMe <https://www.23andme.com>`_,
`Family Tree DNA <https://www.familytreedna.com>`_, and
`Ancestry <http://www.ancestry.com>`_), primarily for the purposes of genetic
genealogy.

Capabilities
------------
- Merge raw data files from different DNA testing companies, identifying discrepant SNPs in the process
- Compute centimorgans (cMs) of shared DNA between individuals using HapMap tables
- Plot shared DNA between individuals
- Find discordant SNPs between child and parent(s)
- Remap SNPs between assemblies / builds (e.g., convert SNPs from build 36 to build 37, etc.)

Dependencies
------------
``lineage`` requires Python 3.4+, `pandas <http://pandas.pydata.org>`_, and
`matplotlib <http://matplotlib.org>`_.

Installation
------------
``lineage`` is `available <https://pypi.python.org/pypi/lineage/>`_ on the
`Python Package Index <https://pypi.python.org/pypi>`_. Install ``lineage`` via
``pip``::

    pip install lineage

Documentation
-------------
Documentation is available `here <https://apriha.github.io/lineage/>`_.

Acknowledgements
----------------
Thanks to Whit Athey, Ryan Dale, Binh Bui, Gopal Vashishtha, and
`CS50 <https://cs50.harvard.edu>`_.

License
-------
Copyright (C) 2016 Andrew Riha

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
