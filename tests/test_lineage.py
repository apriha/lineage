"""
Copyright (C) 2018 Andrew Riha

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
import warnings

import numpy as np
import pandas as pd

from tests import BaseLineageTestCase


class TestLineage(BaseLineageTestCase):
    def get_discordant_snps(self, ind, df):
        ind._build = 37

        snps = df.loc[:, ["chrom", "pos", ind.name]]
        snps = snps.rename(columns={ind.name: "genotype"})

        ind._snps = snps

        return ind

    def test_download_example_datasets(self):
        paths = self.l.download_example_datasets()

        for path in paths:
            if path is None or not os.path.exists(path):
                warnings.warn("Example dataset(s) not currently available")
                return

        assert True

    def test_find_discordant_snps(self):
        df = pd.read_csv(
            "tests/input/discordant_snps.csv",
            skiprows=1,
            na_values="--",
            names=["rsid", "chrom", "pos", "ind1", "ind2", "ind3"],
            index_col=0,
            dtype={"chrom": object, "pos": np.int64},
        )

        ind1 = self.get_discordant_snps(self.l.create_individual("ind1"), df)
        ind2 = self.get_discordant_snps(self.l.create_individual("ind2"), df)
        ind3 = self.get_discordant_snps(self.l.create_individual("ind3"), df)

        df_ind1_ind2 = self.l.find_discordant_snps(ind1, ind2, save_output=True)
        df_ind1_ind2_ind3 = self.l.find_discordant_snps(
            ind1, ind2, ind3, save_output=True
        )

        pd.testing.assert_index_equal(
            pd.Index(
                [
                    "rs34",
                    "rs47",
                    "rs48",
                    "rs49",
                    "rs50",
                    "rs106",
                    "rs110",
                    "rs111",
                    "rs112",
                    "rs113",
                    "rs137",
                    "rs140",
                    "rs141",
                    "rs144",
                    "rs146",
                    "rs147",
                ],
                name="rsid",
            ),
            df_ind1_ind2.index,
        )

        pd.testing.assert_index_equal(
            pd.Index(
                [
                    "rs30",
                    "rs34",
                    "rs38",
                    "rs42",
                    "rs46",
                    "rs47",
                    "rs48",
                    "rs49",
                    "rs50",
                    "rs60",
                    "rs75",
                    "rs85",
                    "rs100",
                    "rs102",
                    "rs106",
                    "rs110",
                    "rs111",
                    "rs112",
                    "rs113",
                    "rs114",
                    "rs118",
                    "rs122",
                    "rs135",
                    "rs137",
                    "rs139",
                    "rs140",
                    "rs141",
                    "rs142",
                    "rs144",
                    "rs146",
                    "rs147",
                    "rs148",
                ],
                name="rsid",
            ),
            df_ind1_ind2_ind3.index,
        )

        assert os.path.exists("output/discordant_snps_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/discordant_snps_ind1_ind2_ind3_GRCh37.csv")

    def test_find_shared_dna_one_ind(self):
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))

        d = self.l.find_shared_dna([ind1], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 0
        assert len(d["two_chrom_shared_dna"]) == 0
        assert len(d["one_chrom_shared_genes"]) == 0
        assert len(d["two_chrom_shared_genes"]) == 0
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        assert not os.path.exists("output/shared_dna_one_chrom_ind1_GRCh37.csv")
        assert not os.path.exists("output/shared_dna_two_chroms_ind1_GRCh37.csv")
        assert not os.path.exists("output/shared_genes_one_chrom_ind1_GRCh37.csv")
        assert not os.path.exists("output/shared_genes_two_chroms_ind1_GRCh37.csv")
        assert not os.path.exists("output/shared_dna_ind1.png")

    def test_find_shared_dna_two_chrom_shared(self):
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))
        ind2 = self.simulate_snps(self.l.create_individual("ind2"))

        d = self.l.find_shared_dna([ind1, ind2], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 1
        assert len(d["two_chrom_shared_dna"]) == 1
        assert len(d["one_chrom_shared_genes"]) == 3811
        assert len(d["two_chrom_shared_genes"]) == 3811
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        np.testing.assert_allclose(d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        np.testing.assert_allclose(d["two_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        assert os.path.exists("output/shared_dna_one_chrom_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_one_chrom_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_ind1_ind2.png")

    def test_find_shared_dna_two_chrom_shared_three_ind(self):
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))
        ind2 = self.simulate_snps(self.l.create_individual("ind2"))
        ind3 = self.simulate_snps(self.l.create_individual("ind3"))

        d = self.l.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 1
        assert len(d["two_chrom_shared_dna"]) == 1
        assert len(d["one_chrom_shared_genes"]) == 3811
        assert len(d["two_chrom_shared_genes"]) == 3811
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        np.testing.assert_allclose(d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        np.testing.assert_allclose(d["two_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        assert os.path.exists("output/shared_dna_one_chrom_ind1_ind2_ind3_GRCh37.csv")
        assert os.path.exists("output/shared_dna_two_chroms_ind1_ind2_ind3_GRCh37.csv")
        assert os.path.exists("output/shared_genes_one_chrom_ind1_ind2_ind3_GRCh37.csv")
        assert os.path.exists(
            "output/shared_genes_two_chroms_ind1_ind2_ind3_GRCh37.csv"
        )
        assert os.path.exists("output/shared_dna_ind1_ind2_ind3.png")

    def test_find_shared_dna_two_chrom_shared_no_output(self):
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))
        ind2 = self.simulate_snps(self.l.create_individual("ind2"))

        d = self.l.find_shared_dna([ind1, ind2], shared_genes=True, save_output=False)

        assert len(d["one_chrom_shared_dna"]) == 1
        assert len(d["two_chrom_shared_dna"]) == 1
        assert len(d["one_chrom_shared_genes"]) == 3811
        assert len(d["two_chrom_shared_genes"]) == 3811
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        np.testing.assert_allclose(d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        np.testing.assert_allclose(d["two_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        assert not os.path.exists("output/shared_dna_one_chrom_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_dna_two_chroms_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_genes_one_chrom_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_genes_two_chroms_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_dna_ind1_ind2.png")

    def test_find_shared_dna_one_chrom_shared(self):
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))
        ind2 = self.simulate_snps(
            self.l.create_individual("ind2"), complement_genotype_one_chrom=True
        )

        d = self.l.find_shared_dna([ind1, ind2], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 1
        assert len(d["two_chrom_shared_dna"]) == 0
        assert len(d["one_chrom_shared_genes"]) == 3811
        assert len(d["two_chrom_shared_genes"]) == 0
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        np.testing.assert_allclose(d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        assert os.path.exists("output/shared_dna_one_chrom_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_dna_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_one_chrom_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_genes_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_ind1_ind2.png")

    def test_find_shared_dna_one_chrom_shared_three_ind(self):
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))
        ind2 = self.simulate_snps(
            self.l.create_individual("ind2"), complement_genotype_one_chrom=True
        )
        ind3 = self.simulate_snps(self.l.create_individual("ind3"))

        d = self.l.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 1
        assert len(d["two_chrom_shared_dna"]) == 0
        assert len(d["one_chrom_shared_genes"]) == 3811
        assert len(d["two_chrom_shared_genes"]) == 0
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        np.testing.assert_allclose(d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        assert os.path.exists("output/shared_dna_one_chrom_ind1_ind2_ind3_GRCh37.csv")
        assert not os.path.exists(
            "output/shared_dna_two_chroms_ind1_ind2_ind3_GRCh37.csv"
        )
        assert os.path.exists("output/shared_genes_one_chrom_ind1_ind2_ind3_GRCh37.csv")
        assert not os.path.exists(
            "output/shared_genes_two_chroms_ind1_ind2_ind3_GRCh37.csv"
        )
        assert os.path.exists("output/shared_dna_ind1_ind2_ind3.png")

    def test_find_shared_dna_X_chrom_two_individuals_male(self):
        ind1 = self.simulate_snps(
            self.l.create_individual("ind1"),
            chrom="X",
            pos_max=155270560,
            pos_step=1000,
            genotype="AA",
        )
        ind2 = self.simulate_snps(
            self.l.create_individual("ind2"),
            chrom="X",
            pos_max=155270560,
            pos_step=1000,
            genotype="AA",
        )

        d = self.l.find_shared_dna([ind1, ind2], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
        assert len(d["two_chrom_shared_dna"]) == 1  # PAR1
        assert len(d["one_chrom_shared_genes"]) == 3022
        assert len(d["two_chrom_shared_genes"]) == 54
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        np.testing.assert_allclose(d["one_chrom_shared_dna"].loc[1]["cMs"], 202.022891)
        np.testing.assert_allclose(d["two_chrom_shared_dna"].loc[1]["cMs"], 20.837792)
        assert os.path.exists("output/shared_dna_one_chrom_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_one_chrom_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_ind1_ind2.png")

    def test_find_shared_dna_X_chrom_two_individuals_female(self):
        ind1 = self.simulate_snps(
            self.l.create_individual("ind1"),
            chrom="X",
            pos_max=155270560,
            genotype="AC",
        )
        ind2 = self.simulate_snps(
            self.l.create_individual("ind2"),
            chrom="X",
            pos_max=155270560,
            genotype="AC",
        )

        d = self.l.find_shared_dna([ind1, ind2], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
        assert len(d["two_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
        assert len(d["one_chrom_shared_genes"]) == 3022
        assert len(d["two_chrom_shared_genes"]) == 3022
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        np.testing.assert_allclose(d["one_chrom_shared_dna"].loc[1]["cMs"], 202.022891)
        np.testing.assert_allclose(d["two_chrom_shared_dna"].loc[1]["cMs"], 202.022891)
        assert os.path.exists("output/shared_dna_one_chrom_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_one_chrom_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_ind1_ind2.png")

    def test_find_shared_dna_two_chrom_shared_discrepant_snps(self):
        # simulate discrepant SNPs so that stitching of adjacent shared DNA segments is performed
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))
        ind2 = self.simulate_snps(
            self.l.create_individual("ind2"),
            complement_genotype_one_chrom=True,
            complement_snp_step=5000,
        )

        d = self.l.find_shared_dna([ind1, ind2], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 1
        assert len(d["two_chrom_shared_dna"]) == 1
        assert len(d["one_chrom_shared_genes"]) == 3811
        assert len(d["two_chrom_shared_genes"]) == 3811
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 2
        np.testing.assert_allclose(d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        np.testing.assert_allclose(d["two_chrom_shared_dna"].loc[1]["cMs"], 140.443968)
        assert os.path.exists("output/shared_dna_one_chrom_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_one_chrom_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_genes_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_ind1_ind2.png")

    def test_find_shared_dna_no_shared_dna(self):
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))
        ind2 = self.simulate_snps(
            self.l.create_individual("ind2"), complement_genotype_two_chroms=True
        )

        d = self.l.find_shared_dna([ind1, ind2], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 0
        assert len(d["two_chrom_shared_dna"]) == 0
        assert len(d["one_chrom_shared_genes"]) == 0
        assert len(d["two_chrom_shared_genes"]) == 0
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        assert not os.path.exists("output/shared_dna_one_chrom_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_dna_two_chroms_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_genes_one_chrom_ind1_ind2_GRCh37.csv")
        assert not os.path.exists("output/shared_genes_two_chroms_ind1_ind2_GRCh37.csv")
        assert os.path.exists("output/shared_dna_ind1_ind2.png")

    def test_find_shared_dna_no_shared_dna_three_ind(self):
        ind1 = self.simulate_snps(self.l.create_individual("ind1"))
        ind2 = self.simulate_snps(
            self.l.create_individual("ind2"), complement_genotype_two_chroms=True
        )
        ind3 = self.simulate_snps(self.l.create_individual("ind3"))

        d = self.l.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

        assert len(d["one_chrom_shared_dna"]) == 0
        assert len(d["two_chrom_shared_dna"]) == 0
        assert len(d["one_chrom_shared_genes"]) == 0
        assert len(d["two_chrom_shared_genes"]) == 0
        assert len(d["one_chrom_discrepant_snps"]) == 0
        assert len(d["two_chrom_discrepant_snps"]) == 0
        assert not os.path.exists(
            "output/shared_dna_one_chrom_ind1_ind2_ind3_GRCh37.csv"
        )
        assert not os.path.exists(
            "output/shared_dna_two_chroms_ind1_ind2_ind3_GRCh37.csv"
        )
        assert not os.path.exists(
            "output/shared_genes_one_chrom_ind1_ind2_ind3_GRCh37.csv"
        )
        assert not os.path.exists(
            "output/shared_genes_two_chroms_ind1_ind2_ind3_GRCh37.csv"
        )
        assert os.path.exists("output/shared_dna_ind1_ind2_ind3.png")
