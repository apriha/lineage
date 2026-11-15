from unittest import TestCase

import numpy as np
import pandas as pd
from snps.testing import SNPsTestMixin, create_simulated_snp_df

from lineage import Lineage


def synthetic_genetic_map() -> dict[str, pd.DataFrame]:
    """Make a synthetic genetic map of chromosomes 1-22 and X.

    Each map runs from 1 Mb to 251 Mb, with a recombination rate that alternates between
    0.5 and 1.5 cM/Mb every 10 Mb. The rate starts at 0.5 on chromosomes 1, 3, ..., 21 and
    X, and at 1.5 on chromosomes 2, 4, ..., 22.

    Returns
    -------
    dict of pandas.DataFrame
        map of each chromosome, with the columns of
        ``Resources.get_genetic_map_HapMapII_GRCh37``: ``pos``, ``rate`` (cM/Mb), and
        ``map`` (cM)
    """
    pos = np.arange(1_000_000, 252_000_000, 10_000_000, dtype=np.int64)
    chroms = [str(chrom) for chrom in range(1, 23)] + ["X"]
    genetic_map = {}
    for i, chrom in enumerate(chroms):
        rate = np.resize([0.5, 1.5] if i % 2 == 0 else [1.5, 0.5], len(pos))
        cM = np.r_[0, np.cumsum(rate[:-1] * np.diff(pos) / 1e6)]
        genetic_map[chrom] = pd.DataFrame({"pos": pos, "rate": rate, "map": cM})
    return genetic_map


class BaseLineageTestCase(SNPsTestMixin, TestCase):
    def setUp(self):
        self.l = Lineage()

    def simulate_snps(
        self,
        ind,
        chrom="1",
        pos_start=1,
        pos_max=111700002,
        pos_step=10000,
        genotype="AA",
        insert_nulls=True,
        null_snp_step=101,
        complement_genotype_one_chrom=False,
        complement_genotype_two_chroms=False,
        complement_snp_step=50,
    ):
        """Simulate SNP data for an individual."""
        ind._build = 37
        ind._snps = create_simulated_snp_df(
            chrom=chrom,
            pos_start=pos_start,
            pos_max=pos_max,
            pos_step=pos_step,
            pos_dtype=np.int64,
            genotype=genotype,
            insert_nulls=insert_nulls,
            null_snp_step=null_snp_step,
            complement_genotype_one_allele=complement_genotype_one_chrom,
            complement_genotype_two_alleles=complement_genotype_two_chroms,
            complement_snp_step=complement_snp_step,
        )
        return ind
