.. image:: https://raw.githubusercontent.com/apriha/lineage/master/docs/images/lineage_banner.png

|ci| |codecov| |docs| |pypi| |python| |downloads| |license| |black|

lineage
=======
``lineage`` provides a framework for analyzing genotype (raw data) files from direct-to-consumer
(DTC) DNA testing companies, primarily for the purposes of genetic genealogy.

Capabilities
------------
- Find shared DNA and genes between individuals
- Compute centiMorgans (cMs) of shared DNA using a variety of genetic maps (e.g., HapMap Phase II, 1000 Genomes Project)
- Plot shared DNA between individuals
- Find discordant SNPs between child and parent(s)
- Read, write, merge, and remap SNPs for an individual via the `snps <https://github.com/apriha/snps>`_ package

Supported Genotype Files
------------------------
``lineage`` supports all genotype files supported by `snps <https://github.com/apriha/snps>`_.

Installation
------------
``lineage`` is `available <https://pypi.org/project/lineage/>`_ on the
`Python Package Index <https://pypi.org>`_. Install ``lineage`` (and its required
Python dependencies) via ``pip``::

    $ pip install lineage

Also see the `installation documentation <https://lineage.readthedocs.io/en/stable/installation.html>`_.

Dependencies
------------
``lineage`` requires `Python <https://www.python.org>`_ 3.8+ and the following Python packages:

- `numpy <https://numpy.org>`_
- `pandas <https://pandas.pydata.org>`_
- `matplotlib <https://matplotlib.org>`_
- `atomicwrites <https://github.com/untitaker/python-atomicwrites>`_
- `snps <https://github.com/apriha/snps>`_

Examples
--------
Initialize the lineage Framework
````````````````````````````````
Import ``Lineage`` and instantiate a ``Lineage`` object:

>>> from lineage import Lineage
>>> l = Lineage()

Download Example Data
`````````````````````
First, let's setup logging to get some helpful output:

>>> import logging, sys
>>> logger = logging.getLogger()
>>> logger.setLevel(logging.INFO)
>>> logger.addHandler(logging.StreamHandler(sys.stdout))

Now we're ready to download some example data from `openSNP <https://opensnp.org>`_:

>>> paths = l.download_example_datasets()
Downloading resources/662.23andme.340.txt.gz
Downloading resources/662.ftdna-illumina.341.csv.gz
Downloading resources/663.23andme.305.txt.gz
Downloading resources/4583.ftdna-illumina.3482.csv.gz
Downloading resources/4584.ftdna-illumina.3483.csv.gz

We'll call these datasets ``User662``, ``User663``, ``User4583``, and ``User4584``.

Load Raw Data
`````````````
Create an ``Individual`` in the context of the ``lineage`` framework to interact with the
``User662`` dataset:

>>> user662 = l.create_individual('User662', ['resources/662.23andme.340.txt.gz', 'resources/662.ftdna-illumina.341.csv.gz'])
Loading SNPs('662.23andme.340.txt.gz')
Merging SNPs('662.ftdna-illumina.341.csv.gz')
SNPs('662.ftdna-illumina.341.csv.gz') has Build 36; remapping to Build 37
Downloading resources/NCBI36_GRCh37.tar.gz
27 SNP positions were discrepant; keeping original positions
151 SNP genotypes were discrepant; marking those as null

Here we created ``user662`` with the name ``User662``. In the process, we merged two raw data
files for this individual. Specifically:

- ``662.23andme.340.txt.gz`` was loaded.
- Then, ``662.ftdna-illumina.341.csv.gz`` was merged. In the process, it was found to have Build
  36. So, it was automatically remapped to Build 37 (downloading the remapping data in the
  process) to match the build of the SNPs already loaded. After this merge, 27 SNP positions and
  151 SNP genotypes were found to be discrepant.

``user662`` is represented by an ``Individual`` object, which inherits from ``snps.SNPs``.
Therefore, all of the `properties and methods <https://snps.readthedocs.io/en/stable/snps.html>`_
available to a ``SNPs`` object are available here; for example:

>>> len(user662.discrepant_merge_genotypes)
151
>>> user662.build
37
>>> user662.build_detected
True
>>> user662.assembly
'GRCh37'
>>> user662.count
1006960

As such, SNPs can be saved, remapped, merged, etc. See the
`snps <https://github.com/apriha/snps>`_ package for further examples.

Compare Individuals
```````````````````
Let's create another ``Individual`` for the ``User663`` dataset:

>>> user663 = l.create_individual('User663', 'resources/663.23andme.305.txt.gz')
Loading SNPs('663.23andme.305.txt.gz')

Now we can perform some analysis between the ``User662`` and ``User663`` datasets.

`Find Discordant SNPs <https://lineage.readthedocs.io/en/stable/lineage.html#lineage.Lineage.find_discordant_snps>`_
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
First, let's find discordant SNPs (i.e., SNP data that is not consistent with Mendelian
inheritance):

>>> discordant_snps = l.find_discordant_snps(user662, user663, save_output=True)
Saving output/discordant_snps_User662_User663_GRCh37.csv

All `output files <https://lineage.readthedocs.io/en/stable/output_files.html>`_ are saved to
the output directory (a parameter to ``Lineage``).

This method also returns a ``pandas.DataFrame``, and it can be inspected interactively at
the prompt, although the same output is available in the CSV file.

>>> len(discordant_snps.loc[discordant_snps['chrom'] != 'MT'])
37

Not counting mtDNA SNPs, there are 37 discordant SNPs between these two datasets.

`Find Shared DNA <https://lineage.readthedocs.io/en/stable/lineage.html#lineage.Lineage.find_shared_dna>`_
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
``lineage`` uses the probabilistic recombination rates throughout the human genome from the
`International HapMap Project <https://www.genome.gov/10001688/international-hapmap-project/>`_
and the `1000 Genomes Project <https://www.internationalgenome.org>`_ to compute the shared DNA
(in centiMorgans) between two individuals. Additionally, ``lineage`` denotes when the shared DNA
is shared on either one or both chromosomes in a pair. For example, when siblings share a segment
of DNA on both chromosomes, they inherited the same DNA from their mother and father for that
segment.

With that background, let's find the shared DNA between the ``User662`` and ``User663`` datasets,
calculating the centiMorgans of shared DNA and plotting the results:

>>> results = l.find_shared_dna([user662, user663], cM_threshold=0.75, snp_threshold=1100)
Downloading resources/genetic_map_HapMapII_GRCh37.tar.gz
Downloading resources/cytoBand_hg19.txt.gz
Saving output/shared_dna_User662_User663_0p75cM_1100snps_GRCh37_HapMap2.png
Saving output/shared_dna_one_chrom_User662_User663_0p75cM_1100snps_GRCh37_HapMap2.csv

Notice that the centiMorgan and SNP thresholds for each DNA segment can be tuned. Additionally,
notice that two files were downloaded to facilitate the analysis and plotting - future analyses
will use the downloaded files instead of downloading the files again. Finally, notice that a list
of individuals is passed to ``find_shared_dna``... This list can contain an arbitrary number of
individuals, and ``lineage`` will find shared DNA across all individuals in the list (i.e.,
where all individuals share segments of DNA on either one or both chromosomes).

Output is returned as a dictionary with the following keys (``pandas.DataFrame`` and
``pandas.Index`` items):

>>> sorted(results.keys())
['one_chrom_discrepant_snps', 'one_chrom_shared_dna', 'one_chrom_shared_genes', 'two_chrom_discrepant_snps', 'two_chrom_shared_dna', 'two_chrom_shared_genes']

In this example, there are 27 segments of shared DNA:

>>> len(results['one_chrom_shared_dna'])
27

Also, `output files <https://lineage.readthedocs.io/en/stable/output_files.html>`_ are
created; these files are detailed in the documentation and their generation can be disabled with a
``save_output=False`` argument. In this example, the output files consist of a CSV file that
details the shared segments of DNA on one chromosome and a plot that illustrates the shared DNA:

.. image:: https://raw.githubusercontent.com/apriha/lineage/master/docs/images/shared_dna_User662_User663_0p75cM_1100snps_GRCh37_HapMap2.png

`Find Shared Genes <https://lineage.readthedocs.io/en/stable/lineage.html#lineage.Lineage.find_shared_dna>`_
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
The `Central Dogma of Molecular Biology <https://en.wikipedia.org/wiki/Central_dogma_of_molecular_biology>`_
states that genetic information flows from DNA to mRNA to proteins: DNA is transcribed into
mRNA, and mRNA is translated into a protein. It's more complicated than this (it's biology
after all), but generally, one mRNA produces one protein, and the mRNA / protein is considered a
gene.

Therefore, it would be interesting to understand not just what DNA is shared between individuals,
but what *genes* are shared between individuals *with the same variations*. In other words,
what genes are producing the *same* proteins? [*]_ Since ``lineage`` can determine the shared DNA
between individuals, it can use that information to determine what genes are also shared on
either one or both chromosomes.

.. [*] In theory, shared segments of DNA should be producing the same proteins, but there are many
 complexities, such as copy number variation (CNV), gene expression, etc.

For this example, let's create two more ``Individuals`` for the ``User4583`` and ``User4584``
datasets:

>>> user4583 = l.create_individual('User4583', 'resources/4583.ftdna-illumina.3482.csv.gz')
Loading SNPs('4583.ftdna-illumina.3482.csv.gz')

>>> user4584 = l.create_individual('User4584', 'resources/4584.ftdna-illumina.3483.csv.gz')
Loading SNPs('4584.ftdna-illumina.3483.csv.gz')

Now let's find the shared genes, specifying a
`population-specific <https://www.internationalgenome.org/faq/which-populations-are-part-your-study/>`_
1000 Genomes Project genetic map (e.g., as predicted by `ezancestry <https://github.com/arvkevi/ezancestry>`_!):

>>> results = l.find_shared_dna([user4583, user4584], shared_genes=True, genetic_map="CEU")
Downloading resources/CEU_omni_recombination_20130507.tar
Downloading resources/knownGene_hg19.txt.gz
Downloading resources/kgXref_hg19.txt.gz
Saving output/shared_dna_User4583_User4584_0p75cM_1100snps_GRCh37_CEU.png
Saving output/shared_dna_one_chrom_User4583_User4584_0p75cM_1100snps_GRCh37_CEU.csv
Saving output/shared_dna_two_chroms_User4583_User4584_0p75cM_1100snps_GRCh37_CEU.csv
Saving output/shared_genes_one_chrom_User4583_User4584_0p75cM_1100snps_GRCh37_CEU.csv
Saving output/shared_genes_two_chroms_User4583_User4584_0p75cM_1100snps_GRCh37_CEU.csv

The plot that illustrates the shared DNA is shown below. Note that in addition to outputting the
shared DNA segments on either one or both chromosomes, the shared genes on either one or both
chromosomes are also output.

.. note:: Shared DNA is not computed on the X chromosome with the 1000 Genomes Project genetic
          maps since the X chromosome is not included in these genetic maps.

In this example, there are 15,976 shared genes on both chromosomes transcribed from 36 segments
of shared DNA:

>>> len(results['two_chrom_shared_genes'])
15976
>>> len(results['two_chrom_shared_dna'])
36

.. image:: https://raw.githubusercontent.com/apriha/lineage/master/docs/images/shared_dna_User4583_User4584_0p75cM_1100snps_GRCh37_CEU.png

Documentation
-------------
Documentation is available `here <https://lineage.readthedocs.io/>`_.

Acknowledgements
----------------
Thanks to Whit Athey, Ryan Dale, Binh Bui, Jeff Gill, Gopal Vashishtha,
`CS50 <https://cs50.harvard.edu>`_, and `openSNP <https://opensnp.org>`_.

``lineage`` incorporates code and concepts generated with the assistance of
`OpenAI's <https://openai.com>`_ `ChatGPT <https://chatgpt.com>`_ . ✨

.. https://github.com/rtfd/readthedocs.org/blob/master/docs/badges.rst
.. |ci| image:: https://github.com/apriha/lineage/actions/workflows/ci.yml/badge.svg?branch=master
   :target: https://github.com/apriha/lineage/actions/workflows/ci.yml
.. |codecov| image:: https://codecov.io/gh/apriha/lineage/branch/master/graph/badge.svg
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
   :target: https://github.com/apriha/lineage/blob/master/LICENSE.txt
.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
