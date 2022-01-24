"""
MIT License

Copyright (c) 2019 Andrew Riha

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import os
import shutil
from unittest import TestCase

import numpy as np
import pandas as pd
from snps.io.reader import NORMALIZED_DTYPES

from lineage import Lineage


class BaseLineageTestCase(TestCase):
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
        ind._build = 37

        positions = np.arange(pos_start, pos_max, pos_step, dtype=np.int64)
        snps = pd.DataFrame(
            {"chrom": chrom},
            index=pd.Index(
                ["rs" + str(x + 1) for x in range(len(positions))], name="rsid"
            ),
        )
        snps["pos"] = positions
        snps["genotype"] = genotype

        if insert_nulls:
            snps.loc[snps.iloc[0::null_snp_step, :].index, "genotype"] = np.nan

        indices = snps.iloc[0::complement_snp_step, :].index
        if complement_genotype_two_chroms:
            snps.loc[indices, "genotype"] = snps.loc[indices, "genotype"].apply(
                self.complement_two_chroms
            )
        elif complement_genotype_one_chrom:
            snps.loc[indices, "genotype"] = snps.loc[indices, "genotype"].apply(
                self.complement_one_chrom
            )

        ind._snps = snps

        return ind

    @staticmethod
    def get_complement(base):
        if base == "A":
            return "T"
        elif base == "G":
            return "C"
        elif base == "C":
            return "G"
        elif base == "T":
            return "A"
        else:
            return base

    def complement_one_chrom(self, genotype):
        if pd.isnull(genotype):
            return np.nan

        complement = ""

        for base in list(genotype):
            complement += self.get_complement(base)
            complement += genotype[1]
            return complement

    def complement_two_chroms(self, genotype):
        if pd.isnull(genotype):
            return np.nan

        complement = ""

        for base in list(genotype):
            complement += self.get_complement(base)

        return complement

    @staticmethod
    def create_snp_df(rsid, chrom, pos, genotype):
        df = pd.DataFrame(
            {"rsid": rsid, "chrom": chrom, "pos": pos, "genotype": genotype},
            columns=["rsid", "chrom", "pos", "genotype"],
        )
        df = df.astype(NORMALIZED_DTYPES)
        df = df.set_index("rsid")
        return df

    def generic_snps(self):
        return self.create_snp_df(
            rsid=["rs" + str(i) for i in range(1, 9)],
            chrom=["1"] * 8,
            pos=list(range(101, 109)),
            genotype=["AA", "CC", "GG", "TT", np.nan, "GC", "TC", "AT"],
        )

    @property
    def downloads_enabled(self):
        """Property indicating if downloads are enabled.

        Only download from external resources when an environment variable named
        "DOWNLOADS_ENABLED" is set to "true".

        Returns
        -------
        bool
        """
        return True if os.getenv("DOWNLOADS_ENABLED") == "true" else False
