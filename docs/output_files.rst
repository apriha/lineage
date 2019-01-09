Output Files
============
The various output files produced by ``lineage`` are detailed below. Output files are saved in
the output directory, which is defined at the instantiation of the :class:`~lineage.Lineage`
object. Generation of output files can usually be enabled or disabled via a ``save_output``
argument to the associated method.

Load SNPs
---------
Multiple raw data files can be loaded when an :class:`~lineage.individual.Individual`
is instantiated. Alternatively, additional files can be loaded after instantiation using
:meth:`~lineage.individual.Individual.load_snps`.

When loading multiple raw data files, if there are any discrepancies between the existing data
and the new data, those are noted. Specifically, discrepancies in SNP positions and genotypes
are output as CSV files. Output of these files is enabled via the ``save_output=True`` argument to
:meth:`~lineage.individual.Individual.load_snps`.

<name>_discrepant_positions_<num>.csv
`````````````````````````````````````
Where ``name`` is the name of the :class:`~lineage.individual.Individual` and ``num``
indicates the file count for discrepant positions files.

==============  ===========
Column          Description
==============  ===========
rsid            SNP ID
chrom           Chromosome of existing SNP
pos             Position of existing SNP
genotype        Genotype of existing SNP
chrom_added     Chromosome of added SNP
pos_added       Position of added SNP (discrepant with pos)
genotype_added  Genotype of added SNP
==============  ===========

A large number of discrepant positions could indicate that the files contain SNPs from different
builds / assemblies.

<name>_discrepant_genotypes_<num>.csv
`````````````````````````````````````
Where ``name`` is the name of the :class:`~lineage.individual.Individual` and ``num``
indicates the file count for discrepant genotypes files.

===============  ===========
Column           Description
===============  ===========
rsid             SNP ID
chrom            Chromosome of existing SNP
pos              Position of existing SNP
genotype         Genotype of existing SNP
chrom_added      Chromosome of added SNP
pos_added        Position of added SNP
genotype_added   Genotype of added SNP (discrepant with genotype)
===============  ===========

A large number of discrepant genotypes could indicate that the files contain SNPs for different
individuals.

Discrepant SNPs
---------------
Summaries can be saved of the discrepant SNPs found while loading files.

<name>_discrepant_positions.csv
```````````````````````````````
Where ``name`` is the name of the :class:`~lineage.individual.Individual`. SNPs with discrepant
positions can be saved with :meth:`~lineage.individual.Individual.save_discrepant_positions`.

==============  ===========
Column          Description
==============  ===========
rsid            SNP ID
chrom           Chromosome of existing SNP
pos             Position of existing SNP
genotype        Genotype of existing SNP
chrom_added     Chromosome of added SNP
pos_added       Position of added SNP (discrepant with pos)
genotype_added  Genotype of added SNP
==============  ===========

<name>_discrepant_genotypes.csv
```````````````````````````````
Where ``name`` is the name of the :class:`~lineage.individual.Individual`. SNPs with discrepant
genotypes can be saved with :meth:`~lineage.individual.Individual.save_discrepant_genotypes`.

===============  ===========
Column           Description
===============  ===========
rsid             SNP ID
chrom            Chromosome of existing SNP
pos              Position of existing SNP
genotype         Genotype of existing SNP
chrom_added      Chromosome of added SNP
pos_added        Position of added SNP
genotype_added   Genotype of added SNP (discrepant with genotype)
===============  ===========

<name>_discrepant_snps.csv
``````````````````````````
Where ``name`` is the name of the :class:`~lineage.individual.Individual`. SNPs with discrepant
positions and / or genotypes can be saved with
:meth:`~lineage.individual.Individual.save_discrepant_snps`.

===============  ===========
Column           Description
===============  ===========
rsid             SNP ID
chrom            Chromosome of existing SNP
pos              Position of existing SNP
genotype         Genotype of existing SNP
chrom_added      Chromosome of added SNP
pos_added        Position of added SNP (possibly discrepant with pos)
genotype_added   Genotype of added SNP (possibly discrepant with genotype)
===============  ===========

Save SNPs
---------
The SNPs for an :class:`~lineage.individual.Individual` can be saved with
:meth:`~lineage.individual.Individual.save_snps`. One CSV file is output when SNPs are saved.

<name>_lineage_<assembly>.csv
`````````````````````````````
Where ``name`` is the name of the :class:`~lineage.individual.Individual` and ``assembly`` is the
assembly of the SNPs being saved.

==========  ===========
Column      Description
==========  ===========
rsid        SNP ID
chromosome  Chromosome of SNP
position    Position of SNP
genotype    Genotype of SNP
==========  ===========

Find Discordant SNPs
--------------------
Discordant SNPs between two or three individuals can be identified with
:meth:`~lineage.Lineage.find_discordant_snps`. One CSV file is optionally output when
``save_output=True``.

discordant_snps_<name1>_<name2>_GRCh37.csv
``````````````````````````````````````````
Where ``name1`` is the name of the first :class:`~lineage.individual.Individual` and
``name2`` is the name of the second :class:`~lineage.individual.Individual`.

================  ===========
Column            Description
================  ===========
rsid              SNP ID
chrom             Chromosome of SNP
pos               Position of SNP
genotype_<name1>  Genotype of first individual
genotype_<name2>  Genotype of second individual
================  ===========

discordant_snps_<name1>_<name2>_<name3>_GRCh37.csv
``````````````````````````````````````````````````
Where ``name1`` is the name of the first :class:`~lineage.individual.Individual`,
``name2`` is the name of the second :class:`~lineage.individual.Individual`, and ``name3`` is
the name of the third :class:`~lineage.individual.Individual`.

================  ===========
Column            Description
================  ===========
rsid              SNP ID
chrom             Chromosome of SNP
pos               Position of SNP
genotype_<name1>  Genotype of first individual
genotype_<name2>  Genotype of second individual
genotype_<name3>  Genotype of third individual
================  ===========

Find Shared DNA
---------------
Shared DNA between two individuals can be identified with
:meth:`~lineage.Lineage.find_shared_dna`. One PNG file and up to two CSV files are output when
``save_output=True``.

In the filenames below, ``name1`` is the name of the first
:class:`~lineage.individual.Individual` and ``name2`` is the name of the second
:class:`~lineage.individual.Individual`.

Note that shared DNA will not be shown on the Y chromosome since the Y chromosome does not
recombine; therefore, genetic maps do not have recombination rates for the Y chromosome.

shared_dna_<name1>_<name2>.png
``````````````````````````````
The plot that illustrates shared DNA (i.e., no shared DNA, shared DNA on one chromosome, and
shared DNA on both chromosomes). The centromere for each chromosome is also detailed. Two examples
of this plot are shown below.

.. image:: https://raw.githubusercontent.com/apriha/lineage/master/docs/images/shared_dna_User662_User663.png

In the above plot, note that the two individuals only share DNA on one chromosome. In this plot,
the larger regions where "No shared DNA" is indicated are due to SNPs not being available in
those regions (i.e., SNPs were not tested in those regions).

.. image:: https://raw.githubusercontent.com/apriha/lineage/master/docs/images/shared_dna_User4583_User4584.png

In the above plot, the areas where "No shared DNA" is indicated are the regions where SNPs were
not tested or where DNA is not shared. The areas where "One chromosome shared" is indicated are
regions where the individuals share DNA on one chromosome. The areas where "Two chromosomes
shared" is indicated are regions where the individuals share DNA on both chromosomes in the pair
(i.e., the individuals inherited the same DNA from their father and mother for those regions).
Note that the regions where DNA is shared on both chromosomes is a subset of the regions where
one chromosome is shared.

shared_dna_one_chrom_<name1>_<name2>_GRCh37.csv
```````````````````````````````````````````````
If DNA is shared on one chromosome, a CSV file details the shared segments of DNA.

=======  ===========
Column   Description
=======  ===========
segment  Shared DNA segment number
chrom    Chromosome with matching DNA segment
start    Start position of matching DNA segment
end      End position of matching DNA segment
cMs      CentiMorgans of matching DNA segment
snps     Number of SNPs in matching DNA segment
=======  ===========

shared_dna_two_chroms_<name1>_<name2>_GRCh37.csv
````````````````````````````````````````````````
If DNA is shared on two chromosomes, a CSV file details the shared segments of DNA.

=======  ===========
Column   Description
=======  ===========
segment  Shared DNA segment number
chrom    Pair of chromosomes with matching DNA segment
start    Start position of matching DNA segment on each chromosome
end      End position of matching DNA segment on each chromosome
cMs      CentiMorgans of matching DNA segment on each chromosome
snps     Number of SNPs in matching DNA segment on each chromosome
=======  ===========

Find Shared Genes
-----------------
Shared genes (with the *same genetic variations*) between two individuals can be identified with
:meth:`~lineage.Lineage.find_shared_dna`, with the parameter ``shared_genes=True``.
In addition to the outputs produced by `Find Shared DNA`_, up to two additional CSV files are
output that detail the shared genes when ``save_output=True``.

In the filenames below, ``name1`` is the name of the first
:class:`~lineage.individual.Individual` and ``name2`` is the name of the second
:class:`~lineage.individual.Individual`.

shared_genes_one_chrom_<name1>_<name2>_GRCh37.csv
`````````````````````````````````````````````````
If DNA is shared on one chromosome, this file details the genes shared between the two
individuals on at least one chromosome; these genes are located in the shared DNA segments
specified in `shared_dna_one_chrom_<name1>_<name2>.csv`_.

===========  ============
Column*      Description*
===========  ============
name         Name of gene
geneSymbol   Gene symbol
chrom        Reference sequence chromosome or scaffold
strand       \+ or - for strand
txStart      Transcription start position (or end position for minus strand item)
txEnd        Transcription end position (or start position for minus strand item)
refseq       RefSeq ID
proteinID    UniProt display ID, UniProt accession, or RefSeq protein ID
description  Description
===========  ============

\* `UCSC Genome Browser <http://genome.ucsc.edu>`_ /
`UCSC Table Browser <http://genome.ucsc.edu/cgi-bin/hgTables>`_

shared_genes_two_chroms_<name1>_<name2>_GRCh37.csv
``````````````````````````````````````````````````
If DNA is shared on both chromosomes in a pair, this file details the genes shared between the two
individuals on both chromosomes; these genes are located in the shared DNA segments specified in
`shared_dna_two_chroms_<name1>_<name2>.csv`_.

The file has the same columns as `shared_genes_one_chrom_<name1>_<name2>.csv`_.
