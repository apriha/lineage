"""
MIT License

Copyright (c) 2018 Andrew Riha

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
import tempfile
from unittest.mock import Mock, mock_open, patch
import warnings

import numpy as np
import pandas as pd

from lineage import Lineage
from tests import BaseLineageTestCase


class TestLineage(BaseLineageTestCase):
    def get_discordant_snps(self, ind, df):
        ind._build = 37

        snps = df.loc[:, ["chrom", "pos", ind.name]]
        snps = snps.rename(columns={ind.name: "genotype"})

        ind._snps = snps

        return ind

    def _generate_test_genetic_map(
        self,
        chrom="1",
        pos=(1, 111700001),
        rate=(140.443968 / (111700001 / 1e6), 0),
        map_cMs=(0.000000, 140.443968),
        **kwargs,
    ):
        return {
            chrom: pd.DataFrame(
                {
                    "pos": pos,
                    "rate": rate,
                    "map": map_cMs,
                }
            ),
        }

    def _generate_test_cytoBand_hg19(self):
        return pd.DataFrame(
            {
                "chrom": ["1"],
                "start": [0],
                "end": [111800001],
                "name": ["test"],
                "gie_stain": ["gneg"],
            }
        )

    def _generate_test_gene_dfs(
        self,
        chrom="1",
        len1=3811,
        len2=4000,
        txStart1=1000000,
        txEnd1=2000000,
        txStart2=111600000,
        txEnd2=111800000,
        **kwargs,
    ):
        diff = len2 - len1
        txStart = [txStart1] * len1
        txStart.extend([txStart2] * diff)
        txEnd = [txEnd1] * len1
        txEnd.extend([txEnd2] * diff)

        kg = pd.DataFrame(
            {
                "name": [f"g{i}" for i in range(len2)],
                "chrom": [chrom] * len2,
                "strand": ["+"] * len2,
                "txStart": txStart,
                "txEnd": txEnd,
                "proteinID": ["g"] * len2,
            }
        )
        kg.set_index("name", inplace=True)

        kgXref = pd.DataFrame(
            {
                "kgID": [f"g{i}" for i in range(len2)],
                "geneSymbol": ["s"] * len2,
                "refseq": ["s"] * len2,
                "description": ["s"] * len2,
            }
        )
        kgXref.set_index("kgID", inplace=True)

        return kg, kgXref

    def run_find_shared_dna_test_X(self, f):
        self.run_find_shared_dna_test(
            f,
            chrom="X",
            pos=(1, 2695340, 154929412, 155270560),
            rate=(
                20.837792 / (2695340 / 1e6),
                180.837755 / ((154929412 - 2695340) / 1e6),
                0.347344 / ((155270560 - 154929412) / 1e6),
                0,
            ),
            map_cMs=(
                0.000000,
                20.837792,
                20.837792 + 180.837755,
                20.837792 + 180.837755 + 0.347344,
            ),
            len1=54,
            len2=3022,
            txStart1=2400000,
            txEnd1=2600000,
            txStart2=150000000,
            txEnd2=155000000,
        )

    def run_find_shared_dna_test_1000G(self, f):
        self.run_find_shared_dna_test(
            f,
            HapMap2=False,
            pos=(1, 43800001),
            rate=(63.0402663602 / (43800001 / 1e6), 0),
            map_cMs=(0.0, 63.0402663602),
            len1=2188,
        )

    def run_find_shared_dna_test(self, f, HapMap2=True, **kwargs):
        if self.downloads_enabled:
            f()
        else:
            genetic_map_patch = (
                "lineage.resources.Resources.get_genetic_map_HapMapII_GRCh37"
                if HapMap2
                else "lineage.resources.Resources.get_genetic_map_1000G_GRCh37"
            )
            genetic_map = self._generate_test_genetic_map(**kwargs)
            cytoband = self._generate_test_cytoBand_hg19()
            kg, kgXref = self._generate_test_gene_dfs(**kwargs)

            with patch(
                genetic_map_patch,
                Mock(return_value=genetic_map),
            ):
                with patch(
                    "lineage.resources.Resources.get_cytoBand_hg19",
                    Mock(return_value=cytoband),
                ):
                    with patch(
                        "lineage.resources.Resources.get_knownGene_hg19",
                        Mock(return_value=kg),
                    ):
                        with patch(
                            "lineage.resources.Resources.get_kgXref_hg19",
                            Mock(return_value=kgXref),
                        ):
                            f()

    def _assert_exists(self, files, idx):
        for i, file in enumerate(files):
            if i in idx:
                self.assertTrue(os.path.exists(file))

    def _assert_does_not_exist(self, files, idx):
        for i, file in enumerate(files):
            if i in idx:
                self.assertFalse(os.path.exists(file))

    def _make_file_exist_assertions(
        self, inds, exist="all", genetic_map="HapMap2", output_dir="output"
    ):
        filename_details = f"{inds}_0p75cM_1100snps_GRCh37_{genetic_map}"
        files = [
            os.path.join(output_dir, f"shared_dna_one_chrom_{filename_details}.csv"),
            os.path.join(output_dir, f"shared_dna_two_chroms_{filename_details}.csv"),
            os.path.join(output_dir, f"shared_genes_one_chrom_{filename_details}.csv"),
            os.path.join(output_dir, f"shared_genes_two_chroms_{filename_details}.csv"),
            os.path.join(output_dir, f"shared_dna_{filename_details}.png"),
        ]

        if exist == "all":
            self._assert_exists(files, list(range(5)))
        elif exist == "none":
            self._assert_does_not_exist(files, list(range(5)))
        elif exist == "one_chrom":
            self._assert_exists(files, [0, 2, 4])
            self._assert_does_not_exist(files, [1, 3])
        elif exist == "plots":
            self._assert_exists(files, [4])
            self._assert_does_not_exist(files, list(range(4)))

    def _check_example_paths(self, paths):
        for path in paths:
            if path is None or not os.path.exists(path):
                warnings.warn("Example dataset(s) not currently available")
                return

    def test_download_example_datasets(self):
        if self.downloads_enabled:
            paths = self.l.download_example_datasets()
            self._check_example_paths(paths)
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                self.l._resources._resources_dir = tmpdir

                # use a temporary directory for test resource data
                with patch("urllib.request.urlopen", mock_open(read_data=b"")):
                    paths = self.l.download_example_datasets()

                self._check_example_paths(paths)

                self.l._resources._resources_dir = "resources"

    def test_find_discordant_snps(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)
            df = pd.read_csv(
                "tests/input/discordant_snps.csv",
                skiprows=1,
                na_values="--",
                names=["rsid", "chrom", "pos", "ind1", "ind2", "ind3"],
                index_col=0,
                dtype={"chrom": object, "pos": np.int64},
            )

            ind1 = self.get_discordant_snps(l.create_individual("ind1"), df)
            ind2 = self.get_discordant_snps(l.create_individual("ind2"), df)
            ind3 = self.get_discordant_snps(l.create_individual("ind3"), df)

            df_ind1_ind2 = l.find_discordant_snps(ind1, ind2, save_output=True)
            df_ind1_ind2_ind3 = l.find_discordant_snps(
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

            assert os.path.exists(
                os.path.join(tmpdir, "discordant_snps_ind1_ind2_GRCh37.csv")
            )
            assert os.path.exists(
                os.path.join(tmpdir, "discordant_snps_ind1_ind2_ind3_GRCh37.csv")
            )

    def test_find_shared_dna_one_ind(self):
        self.run_find_shared_dna_test(self._test_find_shared_dna_one_ind)

    def _test_find_shared_dna_one_ind(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))

            d = l.find_shared_dna([ind1], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 0
            assert len(d["two_chrom_shared_dna"]) == 0
            assert len(d["one_chrom_shared_genes"]) == 0
            assert len(d["two_chrom_shared_genes"]) == 0
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            self._make_file_exist_assertions("ind1", exist="none", output_dir=tmpdir)

    def test_find_shared_dna_invalid_genetic_map(self):
        self.run_find_shared_dna_test(self._test_find_shared_dna_invalid_genetic_map)

    def _test_find_shared_dna_invalid_genetic_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(l.create_individual("ind2"))

            d = l.find_shared_dna([ind1, ind2], genetic_map="test")

            assert len(d["one_chrom_shared_dna"]) == 0
            assert len(d["two_chrom_shared_dna"]) == 0
            assert len(d["one_chrom_shared_genes"]) == 0
            assert len(d["two_chrom_shared_genes"]) == 0
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            self._make_file_exist_assertions(
                "ind1_ind2", exist="none", output_dir=tmpdir
            )

    def test_find_shared_dna_two_chrom_shared(self):
        self.run_find_shared_dna_test(self._test_find_shared_dna_two_chrom_shared)

    def _test_find_shared_dna_two_chrom_shared(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(l.create_individual("ind2"))

            d = l.find_shared_dna(
                [ind1, ind2], shared_genes=True, genetic_map="HapMap2"
            )

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            assert len(d["one_chrom_shared_genes"]) == 3811
            assert len(d["two_chrom_shared_genes"]) == 3811
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            self._make_file_exist_assertions("ind1_ind2", output_dir=tmpdir)

    def test_find_shared_dna_two_chrom_shared_1000G(self):
        self.run_find_shared_dna_test_1000G(
            self._test_find_shared_dna_two_chrom_shared_1000G
        )

    def _test_find_shared_dna_two_chrom_shared_1000G(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"), pos_max=43800002)
            ind2 = self.simulate_snps(l.create_individual("ind2"), pos_max=43800002)

            d = l.find_shared_dna([ind1, ind2], shared_genes=True, genetic_map="CEU")

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            assert len(d["one_chrom_shared_genes"]) == 2188
            assert len(d["two_chrom_shared_genes"]) == 2188
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 63.0402663602
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"], 63.0402663602
            )
            self._make_file_exist_assertions(
                "ind1_ind2", genetic_map="CEU", output_dir=tmpdir
            )

    def test_find_shared_dna_two_chrom_shared_three_ind(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_two_chrom_shared_three_ind
        )

    def _test_find_shared_dna_two_chrom_shared_three_ind(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(l.create_individual("ind2"))
            ind3 = self.simulate_snps(l.create_individual("ind3"))

            d = l.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            assert len(d["one_chrom_shared_genes"]) == 3811
            assert len(d["two_chrom_shared_genes"]) == 3811
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            self._make_file_exist_assertions("ind1_ind2_ind3", output_dir=tmpdir)

    def test_find_shared_dna_two_chrom_shared_no_output(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_two_chrom_shared_no_output
        )

    def _test_find_shared_dna_two_chrom_shared_no_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(l.create_individual("ind2"))

            d = l.find_shared_dna([ind1, ind2], shared_genes=True, save_output=False)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            assert len(d["one_chrom_shared_genes"]) == 3811
            assert len(d["two_chrom_shared_genes"]) == 3811
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            self._make_file_exist_assertions(
                "ind1_ind2", exist="none", output_dir=tmpdir
            )

    def test_find_shared_dna_one_chrom_shared(self):
        self.run_find_shared_dna_test(self._test_find_shared_dna_one_chrom_shared)

    def _test_find_shared_dna_one_chrom_shared(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(
                l.create_individual("ind2"), complement_genotype_one_chrom=True
            )

            d = l.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 0
            assert len(d["one_chrom_shared_genes"]) == 3811
            assert len(d["two_chrom_shared_genes"]) == 0
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            self._make_file_exist_assertions(
                "ind1_ind2", exist="one_chrom", output_dir=tmpdir
            )

    def test_find_shared_dna_one_chrom_shared_three_ind(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_one_chrom_shared_three_ind
        )

    def _test_find_shared_dna_one_chrom_shared_three_ind(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(
                l.create_individual("ind2"), complement_genotype_one_chrom=True
            )
            ind3 = self.simulate_snps(l.create_individual("ind3"))

            d = l.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 0
            assert len(d["one_chrom_shared_genes"]) == 3811
            assert len(d["two_chrom_shared_genes"]) == 0
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            self._make_file_exist_assertions(
                "ind1_ind2_ind3", exist="one_chrom", output_dir=tmpdir
            )

    def test_find_shared_dna_X_chrom_two_individuals_male(self):
        self.run_find_shared_dna_test_X(
            self._test_find_shared_dna_X_chrom_two_individuals_male
        )

    def _test_find_shared_dna_X_chrom_two_individuals_male(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(
                l.create_individual("ind1"),
                chrom="X",
                pos_max=155270560,
                pos_step=1000,
                genotype="AA",
            )
            ind2 = self.simulate_snps(
                l.create_individual("ind2"),
                chrom="X",
                pos_max=155270560,
                pos_step=1000,
                genotype="AA",
            )

            d = l.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
            assert len(d["two_chrom_shared_dna"]) == 1  # PAR1
            assert len(d["one_chrom_shared_genes"]) == 3022
            assert len(d["two_chrom_shared_genes"]) == 54
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"],
                202.022891,
                rtol=1e-5,
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"],
                20.837792,
                rtol=1e-3,
            )
            self._make_file_exist_assertions("ind1_ind2", output_dir=tmpdir)

    def test_find_shared_dna_X_chrom_two_individuals_female(self):
        self.run_find_shared_dna_test_X(
            self._test_find_shared_dna_X_chrom_two_individuals_female
        )

    def _test_find_shared_dna_X_chrom_two_individuals_female(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(
                l.create_individual("ind1"),
                chrom="X",
                pos_max=155270560,
                genotype="AC",
            )
            ind2 = self.simulate_snps(
                l.create_individual("ind2"),
                chrom="X",
                pos_max=155270560,
                genotype="AC",
            )

            d = l.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
            assert len(d["two_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
            assert len(d["one_chrom_shared_genes"]) == 3022
            assert len(d["two_chrom_shared_genes"]) == 3022
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"],
                202.022891,
                rtol=1e-5,
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"],
                202.022891,
                rtol=1e-5,
            )
            self._make_file_exist_assertions("ind1_ind2", output_dir=tmpdir)

    def test_find_shared_dna_two_chrom_shared_discrepant_snps(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_two_chrom_shared_discrepant_snps
        )

    def _test_find_shared_dna_two_chrom_shared_discrepant_snps(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            # simulate discrepant SNPs so that stitching of adjacent shared DNA segments is performed
            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(
                l.create_individual("ind2"),
                complement_genotype_one_chrom=True,
                complement_snp_step=5000,
            )

            d = l.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            assert len(d["one_chrom_shared_genes"]) == 3811
            assert len(d["two_chrom_shared_genes"]) == 3811
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 2
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"], 140.443968
            )
            self._make_file_exist_assertions("ind1_ind2", output_dir=tmpdir)

    def test_find_shared_dna_no_shared_dna(self):
        self.run_find_shared_dna_test(self._test_find_shared_dna_no_shared_dna)

    def _test_find_shared_dna_no_shared_dna(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(
                l.create_individual("ind2"), complement_genotype_two_chroms=True
            )

            d = l.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 0
            assert len(d["two_chrom_shared_dna"]) == 0
            assert len(d["one_chrom_shared_genes"]) == 0
            assert len(d["two_chrom_shared_genes"]) == 0
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            self._make_file_exist_assertions(
                "ind1_ind2", exist="plots", output_dir=tmpdir
            )

    def test_find_shared_dna_no_shared_dna_three_ind(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_no_shared_dna_three_ind
        )

    def _test_find_shared_dna_no_shared_dna_three_ind(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            l = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(l.create_individual("ind1"))
            ind2 = self.simulate_snps(
                l.create_individual("ind2"), complement_genotype_two_chroms=True
            )
            ind3 = self.simulate_snps(l.create_individual("ind3"))

            d = l.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 0
            assert len(d["two_chrom_shared_dna"]) == 0
            assert len(d["one_chrom_shared_genes"]) == 0
            assert len(d["two_chrom_shared_genes"]) == 0
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            self._make_file_exist_assertions(
                "ind1_ind2_ind3", exist="plots", output_dir=tmpdir
            )
