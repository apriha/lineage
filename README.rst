.. image:: https://raw.githubusercontent.com/apriha/lineage/master/docs/images/lineage_banner.png

|build| |codecov| |docs| |pypi| |python| |downloads| |license|

lineage
=======
``lineage`` provides a framework for analyzing genotype (raw data) files from direct-to-consumer
(DTC) DNA testing companies, primarily for the purposes of genetic genealogy.

Capabilities
------------
- Compute centiMorgans (cMs) of shared DNA between individuals using the HapMap Phase II genetic map
- Plot shared DNA between individuals
- Determine genes shared between individuals (i.e., genes transcribed from shared DNA segments)
- Find discordant SNPs between child and parent(s)
- Read, write, merge, and remaps SNPs for an individual via the  `snps <https://github.com/apriha/snps>`_ package

Supported Genotype Files
------------------------
``lineage`` supports `VCF <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137218/>`_ files and
genotype files from the following DTC DNA testing companies:

- `23andMe <https://www.23andme.com>`_
- `Ancestry <https://www.ancestry.com>`_
- `Family Tree DNA <https://www.familytreedna.com>`_
- `MyHeritage <https://www.myheritage.com>`_

Dependencies
------------
``lineage`` requires `Python <https://www.python.org>`_ 3.5+ and the following Python packages:

- `numpy <http://www.numpy.org>`_
- `pandas <http://pandas.pydata.org>`_
- `matplotlib <http://matplotlib.org>`_
- `atomicwrites <https://github.com/untitaker/python-atomicwrites>`_
- `snps <https://github.com/apriha/snps>`_

On Linux systems, the following system-level installs may also be required::

    $ sudo apt-get install python3-tk
    $ sudo apt-get install gfortran
    $ sudo apt-get install python-dev
    $ sudo apt-get install python-devel
    $ sudo apt-get install python3.X-dev # (where X == Python minor version)

Installation
------------
``lineage`` is `available <https://pypi.org/project/lineage/>`_ on the
`Python Package Index <https://pypi.org>`_. Install ``lineage`` (and its required
Python dependencies) via ``pip``::

    $ pip install lineage

Examples
--------
Initialize the lineage Framework
````````````````````````````````
Import ``Lineage`` and instantiate a ``Lineage`` object:

>>> from lineage import Lineage
>>> l = Lineage()

Download Example Data
`````````````````````
Let's download some example data from `openSNP <https://opensnp.org>`_:

>>> paths = l.download_example_datasets()
Downloading resources/662.23andme.304.txt.gz
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

>>> user662 = l.create_individual('User662', 'resources/662.ftdna-illumina.341.csv.gz')
Loading resources/662.ftdna-illumina.341.csv.gz

Here we created ``user662`` with the name ``User662`` and loaded a raw data file.

Remap SNPs
``````````
Oops! The data we just loaded is Build 36, but we want Build 37 since the other files in the
datasets are Build 37... Let's remap the SNPs:

>>> user662.build
36
>>> chromosomes_remapped, chromosomes_not_remapped = user662.remap_snps(37)
Downloading resources/NCBI36_GRCh37.tar.gz
>>> user662.build
37
>>> user662.assembly
'GRCh37'

SNPs can be re-mapped between Build 36 (``NCBI36``), Build 37 (``GRCh37``), and Build 38
(``GRCh38``).

Merge Raw Data Files
````````````````````
The dataset for ``User662`` consists of three raw data files from two different DNA testing
companies. Let's load the remaining two files.

As the data gets added, it's compared to the existing data, and SNP position and genotype
discrepancies are identified. (The discrepancy thresholds can be tuned via parameters.)

>>> user662.snp_count
708092
>>> user662.load_snps(['resources/662.23andme.304.txt.gz', 'resources/662.23andme.340.txt.gz'],
...                   discrepant_genotypes_threshold=300)
Loading resources/662.23andme.304.txt.gz
3 SNP positions were discrepant; keeping original positions
8 SNP genotypes were discrepant; marking those as null
Loading resources/662.23andme.340.txt.gz
27 SNP positions were discrepant; keeping original positions
156 SNP genotypes were discrepant; marking those as null
>>> len(user662.discrepant_positions)
30
>>> user662.snp_count
1006960

Save SNPs
`````````
Ok, so far we've remapped the SNPs to the same build and merged the SNPs from three files,
identifying discrepancies along the way. Let's save the merged dataset consisting of over 1M+
SNPs to a CSV file:

>>> saved_snps = user662.save_snps()
Saving output/User662_GRCh37.csv

All `output files <https://lineage.readthedocs.io/en/latest/output_files.html>`_ are saved to the output
directory.

Compare Individuals
```````````````````
Let's create another ``Individual`` for the ``User663`` dataset:

>>> user663 = l.create_individual('User663', 'resources/663.23andme.305.txt.gz')
Loading resources/663.23andme.305.txt.gz

Now we can perform some analysis between the ``User662`` and ``User663`` datasets.

Find Discordant SNPs
''''''''''''''''''''
First, let's find discordant SNPs (i.e., SNP data that is not consistent with Mendelian
inheritance):

>>> discordant_snps = l.find_discordant_snps(user662, user663, save_output=True)
Saving output/discordant_snps_User662_User663_GRCh37.csv

This method also returns a ``pandas.DataFrame``, and it can be inspected interactively at
the prompt, although the same output is available in the CSV file.

>>> len(discordant_snps.loc[discordant_snps['chrom'] != 'MT'])
37

Not counting mtDNA SNPs, there are 37 discordant SNPs between these two datasets.

Find Shared DNA
'''''''''''''''
``lineage`` uses the probabilistic recombination rates throughout the human genome from the
`International HapMap Project <https://www.genome.gov/10001688/international-hapmap-project/>`_ to
compute the shared DNA (in centiMorgans) between two individuals. Additionally, ``lineage``
denotes when the shared DNA is shared on either one or both chromosomes in a pair. For example,
when siblings share a segment of DNA on both chromosomes, they inherited the same DNA from their
mother and father for that segment.

With that background, let's find the shared DNA between the ``User662`` and ``User663`` datasets,
calculating the centiMorgans of shared DNA and plotting the results:

>>> results = l.find_shared_dna([user662, user663], cM_threshold=0.75, snp_threshold=1100)
Downloading resources/genetic_map_HapMapII_GRCh37.tar.gz
Downloading resources/cytoBand_hg19.txt.gz
Saving output/shared_dna_User662_User663.png
Saving output/shared_dna_one_chrom_User662_User663_GRCh37.csv

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

Also, `output files <https://lineage.readthedocs.io/en/latest/output_files.html>`_ are
created; these files are detailed in the documentation and their generation can be disabled with a
``save_output=False`` argument. In this example, the output files consist of a CSV file that
details the shared segments of DNA on one chromosome and a plot that illustrates the shared DNA:

.. image:: https://raw.githubusercontent.com/apriha/lineage/master/docs/images/shared_dna_User662_User663.png

Find Shared Genes
'''''''''''''''''
The `Central Dogma of Molecular Biology <https://www.nature.com/nature/focus/crick/pdf/crick227.pdf>`_
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
Loading resources/4583.ftdna-illumina.3482.csv.gz

>>> user4584 = l.create_individual('User4584', 'resources/4584.ftdna-illumina.3483.csv.gz')
Loading resources/4584.ftdna-illumina.3483.csv.gz

Now let's find the shared genes:

>>> results = l.find_shared_dna([user4583, user4584], shared_genes=True)
Downloading resources/knownGene_hg19.txt.gz
Downloading resources/kgXref_hg19.txt.gz
Saving output/shared_dna_User4583_User4584.png
Saving output/shared_dna_one_chrom_User4583_User4584_GRCh37.csv
Saving output/shared_dna_two_chroms_User4583_User4584_GRCh37.csv
Saving output/shared_genes_one_chrom_User4583_User4584_GRCh37.csv
Saving output/shared_genes_two_chroms_User4583_User4584_GRCh37.csv

The plot that illustrates the shared DNA is shown below. Note that in addition to outputting the
shared DNA segments on either one or both chromosomes, the shared genes on either one or both
chromosomes are also output.

In this example, there are 15,976 shared genes on both chromosomes transcribed from 36 segments
of shared DNA:

>>> len(results['two_chrom_shared_genes'])
15976
>>> len(results['two_chrom_shared_dna'])
36

.. image:: https://raw.githubusercontent.com/apriha/lineage/master/docs/images/shared_dna_User4583_User4584.png

Documentation
-------------
Documentation is available `here <https://lineage.readthedocs.io/>`_.

Acknowledgements
----------------
Thanks to Whit Athey, Ryan Dale, Binh Bui, Jeff Gill, Gopal Vashishtha,
`CS50 <https://cs50.harvard.edu>`_, and `openSNP <https://opensnp.org>`_.

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

.. https://github.com/rtfd/readthedocs.org/blob/master/docs/badges.rst
.. |build| image:: https://travis-ci.org/apriha/lineage.svg?branch=master
   :target: https://travis-ci.org/apriha/lineage
.. |codecov| image:: https://codecov.io/gh/apriha/lineage/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/apriha/lineage
.. |docs| image:: https://readthedocs.org/projects/lineage/badge/?version=latest
   :target: https://lineage.readthedocs.io/
.. |pypi| image:: https://img.shields.io/pypi/v/lineage.svg
   :target: https://pypi.python.org/pypi/lineage
.. |python| image:: https://img.shields.io/pypi/pyversions/lineage.svg
   :target: https://www.python.org
.. |downloads| image:: https://pepy.tech/badge/lineage
   :target: https://pepy.tech/project/lineage
.. |license| image:: https://img.shields.io/pypi/l/lineage.svg
   :target: https://github.com/apriha/lineage/blob/master/LICENSE.txt
