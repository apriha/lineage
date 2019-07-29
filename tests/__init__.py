"""
Copyright (C) 2019 Andrew Riha

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

"""

import os
import shutil
from unittest import TestCase

import numpy as np
import pandas as pd

from lineage import Lineage


class BaseLineageTestCase(TestCase):
    def setUp(self):
        self.l = Lineage()
        self.del_output_dir_helper()

    def tearDown(self):
        self.del_output_dir_helper()

    @staticmethod
    def del_output_dir_helper():
        if os.path.exists("output"):
            shutil.rmtree("output")

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
        df = df.set_index("rsid")
        return df
