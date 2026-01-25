.. image:: https://raw.githubusercontent.com/apriha/lineage/main/docs/images/lineage_banner.png

|ci| |codecov| |docs| |pypi| |python| |downloads| |license| |ruff|

lineage
=======
tools for analyzing and exploring genetic relationships 🧬

``lineage`` *strives to be an easy-to-use and accessible open-source library for genetic genealogy*

Features
--------
Shared DNA Analysis
```````````````````
- Find shared DNA between individuals using genetic maps from the
  `HapMap Project <https://www.genome.gov/10001688/international-hapmap-project/>`_ and
  `1000 Genomes Project <https://www.internationalgenome.org>`_
- Compute centiMorgans (cMs) of shared DNA with configurable thresholds
- Detect IBD1 (half-identical) and IBD2 (fully identical) regions
- Visualize shared DNA segments across all chromosomes
- Find discordant SNPs inconsistent with Mendelian inheritance patterns

Shared Gene Analysis
````````````````````
- Identify genes shared between individuals with the same genetic variations
- Determine which genes produce the same proteins across related individuals

Synthetic Data Generation
`````````````````````````
- Generate realistic synthetic genotype data for related individuals
- Create parent-child pairs with proper single-allele inheritance
- Create sibling pairs with realistic meiotic recombination patterns

Supported Genotype Files
------------------------
``lineage`` supports all genotype files supported by `snps <https://github.com/apriha/snps>`_.

Installation
------------
``lineage`` is `available <https://pypi.org/project/lineage/>`_ on the
`Python Package Index <https://pypi.org>`_. Install ``lineage`` (and its required
Python dependencies) via ``pip``::

    $ pip install lineage

Examples
--------
The examples below demonstrate the core features of ``lineage``. For detailed explanations of
genetic concepts like IBD (Identity By Descent), centiMorgans, and how to interpret results,
see the `Concepts Guide <https://lineage.readthedocs.io/en/stable/concepts.html>`_.

**Optional:** To see file save notifications, configure logging before running examples::

    import logging
    logging.basicConfig(level=logging.INFO, format='%(message)s')

To try these examples, first generate some sample data:

>>> from lineage import Lineage
>>> l = Lineage()
>>> paths = l.create_example_datasets()

Load Individuals
````````````````
Load genotype files and create ``Individual`` objects:

>>> parent = l.create_individual('Parent', paths['parent'])
>>> child = l.create_individual('Child', paths['child'])

Each ``Individual`` inherits from `snps.SNPs <https://snps.readthedocs.io/en/stable/snps.html>`_,
providing access to all SNPs properties and methods:

>>> parent.build
37
>>> parent.assembly
'GRCh37'
>>> parent.count
899992

Find Discordant SNPs
````````````````````
Identify SNPs inconsistent with Mendelian inheritance between parent and child:

>>> discordant_snps = l.find_discordant_snps(parent, child, save_output=True)  # doctest: +SKIP

The example datasets include a small number of simulated genotyping errors (~0.01%) to
demonstrate discordant SNP detection.

Find Shared DNA
```````````````
Parent-Child Example
~~~~~~~~~~~~~~~~~~~~
Compute shared DNA segments and generate a visualization:

>>> results = l.find_shared_dna([parent, child])  # doctest: +SKIP

For parent-child relationships, all shared DNA appears on one chromosome only (IBD1),
representing the chromosome inherited from that parent (~3400-3700 cM total):

.. image:: https://raw.githubusercontent.com/apriha/lineage/main/docs/images/shared_dna_Parent_Child_0p75cM_1100snps_GRCh37_HapMap2.png

Sibling Example
~~~~~~~~~~~~~~~
Analyze siblings using a population-specific genetic map:

>>> sibling1 = l.create_individual('Sibling1', paths['sibling1'])  # doctest: +SKIP
>>> sibling2 = l.create_individual('Sibling2', paths['sibling2'])  # doctest: +SKIP
>>> results = l.find_shared_dna([sibling1, sibling2], shared_genes=True, genetic_map="CEU")  # doctest: +SKIP

Siblings share DNA on one chromosome (IBD1) and both chromosomes (IBD2), reflecting
segments where they inherited the same DNA from one or both parents:

.. image:: https://raw.githubusercontent.com/apriha/lineage/main/docs/images/shared_dna_Sibling1_Sibling2_0p75cM_1100snps_GRCh37_CEU.png

Documentation
-------------
Documentation is available `here <https://lineage.readthedocs.io/>`_.

Acknowledgements
----------------
Thanks to Whit Athey, Ryan Dale, Binh Bui, Jeff Gill, Gopal Vashishtha,
`CS50 <https://cs50.harvard.edu>`_. This project was historically validated using
data from `openSNP <https://opensnp.org>`_.

``lineage`` incorporates code and concepts generated with the assistance of various
generative AI tools (including but not limited to `ChatGPT <https://chatgpt.com>`_,
`Grok <https://grok.com>`_, and `Claude <https://claude.ai>`_). ✨

License
-------
``lineage`` is licensed under the `MIT License <https://github.com/apriha/lineage/blob/main/LICENSE.txt>`_.

.. https://github.com/rtfd/readthedocs.org/blob/master/docs/badges.rst
.. |ci| image:: https://github.com/apriha/lineage/actions/workflows/ci.yml/badge.svg?branch=main
   :target: https://github.com/apriha/lineage/actions/workflows/ci.yml
.. |codecov| image:: https://codecov.io/gh/apriha/lineage/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/apriha/lineage
.. |docs| image:: https://readthedocs.org/projects/lineage/badge/?version=stable
   :target: https://lineage.readthedocs.io/
.. |pypi| image:: https://img.shields.io/pypi/v/lineage.svg
   :target: https://pypi.python.org/pypi/lineage
.. |python| image:: https://img.shields.io/pypi/pyversions/lineage.svg
   :target: https://www.python.org
.. |downloads| image:: https://pepy.tech/badge/lineage
   :target: https://pepy.tech/project/lineage
.. |license| image:: https://img.shields.io/pypi/l/lineage.svg
   :target: https://github.com/apriha/lineage/blob/main/LICENSE.txt
.. |ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
   :target: https://github.com/astral-sh/ruff
   :alt: Ruff
