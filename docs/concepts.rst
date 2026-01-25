Understanding Genetic Relationships
====================================

This guide explains the genetic concepts behind ``lineage`` and how to interpret your results.

How DNA is Inherited
--------------------
Humans have 23 pairs of chromosomes. For each pair, you inherit one chromosome from your mother
and one from your father. This means:

- **You share exactly 50% of your DNA with each parent** — one complete chromosome from each pair
- **Siblings share variable amounts of DNA** — on average ~50%, but it varies due to the
  randomness of inheritance

Identity By Descent (IBD)
`````````````````````````
When ``lineage`` analyzes shared DNA, it uses **Identity By Descent (IBD)** to classify sharing
patterns at each genomic segment:

- **IBD0**: No chromosomes shared — the individuals inherited different chromosomes from both parents
- **IBD1**: One chromosome shared (half-identical region) — the individuals share DNA on one chromosome of a pair (shown in light blue in plots)
- **IBD2**: Both chromosomes shared (fully identical region) — the individuals share DNA on both chromosomes of a pair (shown in dark purple in plots)

Parent-Child Relationships
``````````````````````````
A child inherits exactly one chromosome from each parent at every position. This means:

- All positions between parent and child are **IBD1** (one chromosome shared)
- There is no IBD2 sharing (a child cannot inherit both chromosomes from one parent)
- The total shared DNA is approximately **3400–3700 centiMorgans**

.. image:: https://raw.githubusercontent.com/apriha/lineage/main/docs/images/shared_dna_Parent_Child_0p75cM_1100snps_GRCh37_HapMap2.png

Sibling Relationships
`````````````````````
Siblings inherit chromosomes from the same two parents, but which specific segments they inherit
is random due to a process called **meiotic recombination**. During the formation of egg and
sperm cells, chromosomes exchange segments, creating unique combinations.

For full siblings, at any given position in the genome:

- ~25% chance: **IBD0** (no sharing) — each sibling inherited different chromosomes from both parents
- ~50% chance: **IBD1** (one chromosome shared) — siblings inherited the same chromosome from one parent
- ~25% chance: **IBD2** (both chromosomes shared) — siblings inherited the same chromosomes from both parents

At any genomic position, these states are mutually exclusive: IBD1 means *exactly one* chromosome
copy is shared, while IBD2 means *both* copies are shared. Total shared DNA equals IBD1 + IBD2.

.. image:: https://raw.githubusercontent.com/apriha/lineage/main/docs/images/shared_dna_Sibling1_Sibling2_0p75cM_1100snps_GRCh37_CEU.png

What are CentiMorgans?
----------------------
A **centiMorgan (cM)** is a unit that measures genetic distance based on recombination frequency,
not physical distance. One centiMorgan corresponds to a 1% chance that a segment will be
separated by recombination in a single generation.

Why use centiMorgans instead of base pairs?

- Recombination rates vary across the genome — some regions recombine frequently, others rarely
- CentiMorgans account for these differences, giving a more accurate measure of genetic relatedness
- The total human genome is approximately **3400–3700 cM** (depending on the genetic map used)

Typical shared DNA amounts:

===============  =================  ==================  ==================
Relationship     Shared DNA (cM)    IBD1 (cM)           IBD2 (cM)
===============  =================  ==================  ==================
Parent/Child     ~3400–3700         ~3400–3700          0
Full Siblings    ~2300–2900         ~1700–2000          ~600–900
Half Siblings    ~1300–2300         ~1300–2300          0
First Cousins    ~550–1200          ~550–1200           0
===============  =================  ==================  ==================

.. note:: Total shared DNA = IBD1 + IBD2. Only relationships through both parents (full siblings,
          double first cousins, etc.) can have IBD2 sharing, since IBD2 requires inheriting the
          same DNA from both the mother and father.

Genetic Maps
------------
A **genetic map** provides recombination rates across the genome, allowing ``lineage`` to
convert physical distances (base pairs) into genetic distances (centiMorgans).

``lineage`` supports two sources of genetic maps:

HapMap Phase II
```````````````
The `International HapMap Project <https://www.genome.gov/10001688/international-hapmap-project/>`_
mapped genetic variation across human populations. This is the default genetic map in ``lineage``.

- Covers all autosomes (chromosomes 1–22) and the X chromosome
- Good general-purpose choice for most analyses

1000 Genomes Project
````````````````````
The `1000 Genomes Project <https://www.internationalgenome.org>`_ provides population-specific
genetic maps with more detailed recombination data.

Available populations include:

- **CEU**: Utah residents with Northern/Western European ancestry
- **YRI**: Yoruba in Ibadan, Nigeria
- **CHB**: Han Chinese in Beijing, China
- And `many others <https://www.internationalgenome.org/faq/which-populations-are-part-your-study/>`_

.. note:: The 1000 Genomes maps do not include the X chromosome, so shared DNA on chromosome X
          will not be computed when using these maps.

Using a population-specific map can provide more accurate centiMorgan calculations if you know
the ancestral background of the individuals being compared:

>>> results = l.find_shared_dna([ind1, ind2], genetic_map="CEU")  # doctest: +SKIP

Understanding Thresholds
------------------------
``lineage`` uses two thresholds to filter out noise and reduce false positives when detecting
shared DNA segments:

- **cM_threshold** (default: 0.75 cM): Minimum genetic length of a shared segment
- **snp_threshold** (default: 1100 SNPs): Minimum number of SNPs in a shared segment

These conservative defaults help ensure that detected segments represent true shared DNA rather
than matching by chance. You can adjust these thresholds:

>>> # More permissive — reveals smaller segments but may include false positives
>>> results = l.find_shared_dna([ind1, ind2], cM_threshold=0.5, snp_threshold=500)  # doctest: +SKIP

>>> # More conservative — only shows larger, high-confidence segments
>>> results = l.find_shared_dna([ind1, ind2], cM_threshold=1.0, snp_threshold=1500)  # doctest: +SKIP

Shared Genes
------------
Beyond identifying shared DNA segments, ``lineage`` can identify **shared genes** — genes
located within those shared segments where individuals have the same genetic variations.

The Science
```````````
The `Central Dogma of Molecular Biology <https://en.wikipedia.org/wiki/Central_dogma_of_molecular_biology>`_
describes how genetic information flows:

1. **DNA** is transcribed into **mRNA**
2. **mRNA** is translated into **proteins**

When two individuals share a DNA segment containing a gene, and they have identical genetic
variations in that region, they are likely producing the **same protein** from that gene.

Using Shared Genes
``````````````````
Enable shared gene analysis with the ``shared_genes`` parameter:

>>> results = l.find_shared_dna([ind1, ind2], shared_genes=True)  # doctest: +SKIP

This outputs additional CSV files listing genes in IBD1 and IBD2 regions, including:

- Gene names and symbols
- Chromosomal locations
- RefSeq and protein IDs
- Gene descriptions

.. note:: Shared DNA segments should theoretically produce the same proteins, but biological
          complexities like copy number variation (CNV) and gene expression differences mean
          this is an approximation.

Discordant SNPs
---------------
**Discordant SNPs** are genetic variations that don't follow expected inheritance patterns
(Mendelian inheritance). Finding discordant SNPs between a parent and child can indicate:

- Genotyping errors in the raw data
- Mutations that occurred during DNA replication
- Data quality issues from the testing company

A small number of discordant SNPs is normal. A large number may indicate data quality problems
or that the individuals may not have the expected biological relationship.

>>> discordant = l.find_discordant_snps(parent, child)  # doctest: +SKIP
>>> len(discordant)  # doctest: +SKIP
42  # A small number is typical

Interpreting Results
--------------------
The ``find_shared_dna`` method returns a dictionary with these keys:

================================  ===========
Key                               Description
================================  ===========
``one_chrom_shared_dna``          IBD1 segments (DNA shared on one chromosome)
``two_chrom_shared_dna``          IBD2 segments (DNA shared on both chromosomes)
``one_chrom_shared_genes``        Genes in IBD1 segments (if requested)
``two_chrom_shared_genes``        Genes in IBD2 segments (if requested)
``one_chrom_discrepant_snps``     SNPs with unexpected patterns in IBD1 regions
``two_chrom_discrepant_snps``     SNPs with unexpected patterns in IBD2 regions
================================  ===========

Tips for Analysis
`````````````````
1. **Start with defaults** — The default thresholds work well for most analyses

2. **Compare total centiMorgans** — Sum the ``cMs`` column to get total shared DNA and compare
   with expected values for the suspected relationship

3. **Look at the plot first** — The visualization quickly shows the pattern of sharing across
   all chromosomes

4. **Consider the relationship type**:

   - Parent-child: Expect IBD1 sharing across the entire genome
   - Full siblings: Expect a mix of IBD1 and IBD2 sharing
   - Half siblings: Expect IBD1 sharing only (no IBD2), with gaps showing IBD0 regions
   - More distant relatives: Expect smaller, scattered IBD1 segments

5. **Check for data quality** — Large gaps in the plot may indicate regions where SNPs weren't
   tested rather than regions of no sharing
