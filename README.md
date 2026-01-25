![lineage](https://raw.githubusercontent.com/apriha/lineage/main/docs/images/lineage_banner.png)

[![ci](https://github.com/apriha/lineage/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/apriha/lineage/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/apriha/lineage/branch/main/graph/badge.svg)](https://codecov.io/gh/apriha/lineage)
[![docs](https://readthedocs.org/projects/lineage/badge/?version=stable)](https://lineage.readthedocs.io/)
[![pypi](https://img.shields.io/pypi/v/lineage.svg)](https://pypi.python.org/pypi/lineage)
[![python](https://img.shields.io/pypi/pyversions/lineage.svg)](https://www.python.org)
[![downloads](https://pepy.tech/badge/lineage)](https://pepy.tech/project/lineage)
[![license](https://img.shields.io/pypi/l/lineage.svg)](https://github.com/apriha/lineage/blob/main/LICENSE.txt)
[![ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

# lineage
tools for analyzing and exploring genetic relationships 🧬

*`lineage` strives to be an easy-to-use and accessible open-source library for genetic genealogy*

## Features

### Shared DNA Analysis
- Find shared DNA between individuals using genetic maps from the
  [HapMap Project](https://www.genome.gov/10001688/international-hapmap-project/) and
  [1000 Genomes Project](https://www.internationalgenome.org)
- Compute centiMorgans (cMs) of shared DNA with configurable thresholds
- Detect IBD1 (half-identical) and IBD2 (fully identical) regions
- Visualize shared DNA segments across all chromosomes
- Find discordant SNPs inconsistent with Mendelian inheritance patterns

### Shared Gene Analysis
- Identify genes shared between individuals with the same genetic variations
- Determine which genes produce the same proteins across related individuals

### Synthetic Data Generation
- Generate realistic synthetic genotype data for related individuals
- Create parent-child pairs with proper single-allele inheritance
- Create sibling pairs with realistic meiotic recombination patterns

## Supported Genotype Files
`lineage` supports all genotype files supported by [snps](https://github.com/apriha/snps).

## Installation
`lineage` is [available](https://pypi.org/project/lineage/) on the
[Python Package Index](https://pypi.org). Install `lineage` (and its required
Python dependencies) via `pip`:

```bash
$ pip install lineage
```

## Examples
The examples below demonstrate the core features of `lineage`. For detailed explanations of
genetic concepts like IBD (Identity By Descent), centiMorgans, and how to interpret results,
see the [Concepts Guide](https://lineage.readthedocs.io/en/stable/concepts.html).

**Optional:** To see file save notifications, configure logging before running examples:

```python
import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
```

To try these examples, first generate some sample data:

```python
>>> from lineage import Lineage
>>> l = Lineage()
>>> paths = l.create_example_datasets()
```

### Load Individuals
Load genotype files and create `Individual` objects:

```python
>>> parent = l.create_individual('Parent', paths['parent'])
>>> child = l.create_individual('Child', paths['child'])
```

Each `Individual` inherits from [snps.SNPs](https://snps.readthedocs.io/en/stable/snps.html),
providing access to all SNPs properties and methods:

```python
>>> parent.build
37
>>> parent.assembly
'GRCh37'
>>> parent.count
899992
```

### Find Discordant SNPs
Identify SNPs inconsistent with Mendelian inheritance between parent and child:

```python
>>> discordant_snps = l.find_discordant_snps(parent, child, save_output=True)  # doctest: +SKIP
```

The example datasets include a small number of simulated genotyping errors (~0.01%) to
demonstrate discordant SNP detection.

### Find Shared DNA

#### Parent-Child Example
Compute shared DNA segments and generate a visualization:

```python
>>> results = l.find_shared_dna([parent, child])  # doctest: +SKIP
```

For parent-child relationships, all shared DNA appears on one chromosome only (IBD1),
representing the chromosome inherited from that parent (~3400-3700 cM total):

![Parent-Child shared DNA](https://raw.githubusercontent.com/apriha/lineage/main/docs/images/shared_dna_Parent_Child_0p75cM_1100snps_GRCh37_HapMap2.png)

#### Sibling Example
Analyze siblings using a population-specific genetic map:

```python
>>> sibling1 = l.create_individual('Sibling1', paths['sibling1'])  # doctest: +SKIP
>>> sibling2 = l.create_individual('Sibling2', paths['sibling2'])  # doctest: +SKIP
>>> results = l.find_shared_dna([sibling1, sibling2], shared_genes=True, genetic_map="CEU")  # doctest: +SKIP
```

Siblings share DNA on one chromosome (IBD1) and both chromosomes (IBD2), reflecting
segments where they inherited the same DNA from one or both parents:

![Sibling shared DNA](https://raw.githubusercontent.com/apriha/lineage/main/docs/images/shared_dna_Sibling1_Sibling2_0p75cM_1100snps_GRCh37_CEU.png)

## Documentation
Documentation is available [here](https://lineage.readthedocs.io/).

## Acknowledgements
Thanks to Whit Athey, Ryan Dale, Binh Bui, Jeff Gill, Gopal Vashishtha,
[CS50](https://cs50.harvard.edu). This project was historically validated using
data from [openSNP](https://opensnp.org).

`lineage` incorporates code and concepts generated with the assistance of various
generative AI tools (including but not limited to [ChatGPT](https://chatgpt.com),
[Grok](https://grok.com), and [Claude](https://claude.ai)). ✨

## License
`lineage` is licensed under the [MIT License](https://github.com/apriha/lineage/blob/main/LICENSE.txt).
