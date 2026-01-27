from unittest import TestCase

import numpy as np
from snps.testing import SNPsTestMixin, create_simulated_snp_df

from lineage import Lineage


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
