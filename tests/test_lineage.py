import os
import tempfile
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

from lineage import PILEUP_REGIONS, Lineage
from tests import BaseLineageTestCase
from tests.test_config import SHARED_GENE_COUNTS


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
        test_case="default",
        txStart1=1000000,
        txEnd1=2000000,
        txStart2=111600000,
        txEnd2=111800000,
        **kwargs,
    ):
        len1 = SHARED_GENE_COUNTS[test_case]["len1"]
        len2 = SHARED_GENE_COUNTS[test_case].get("len2", len1)
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
            test_case="X_chrom_male",
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
            test_case="1000G",
        )

    def run_find_shared_dna_test(self, f, HapMap2=True, test_case="default", **kwargs):
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
            kg, kgXref = self._generate_test_gene_dfs(test_case=test_case, **kwargs)

            with (
                patch(
                    genetic_map_patch,
                    Mock(return_value=genetic_map),
                ),
                patch(
                    "lineage.resources.Resources.get_cytoBand_hg19",
                    Mock(return_value=cytoband),
                ),
                patch(
                    "lineage.resources.Resources.get_knownGene_hg19",
                    Mock(return_value=kg),
                ),
                patch(
                    "lineage.resources.Resources.get_kgXref_hg19",
                    Mock(return_value=kgXref),
                ),
            ):
                f()

    def _assert_shared_genes(self, shared_genes, test_case, key="len1"):
        if self.downloads_enabled:
            # downloaded gene resources are updated upstream, so only verify genes were found
            self.assertGreater(len(shared_genes), 0)
        else:
            self.assertEqual(len(shared_genes), SHARED_GENE_COUNTS[test_case][key])

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

    def test_find_discordant_snps(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)
            df = pd.read_csv(
                "tests/input/discordant_snps.csv",
                skiprows=1,
                na_values="--",
                names=["rsid", "chrom", "pos", "ind1", "ind2", "ind3"],
                index_col=0,
                dtype={"chrom": object, "pos": np.int64},
            )

            ind1 = self.get_discordant_snps(ln.create_individual("ind1"), df)
            ind2 = self.get_discordant_snps(ln.create_individual("ind2"), df)
            ind3 = self.get_discordant_snps(ln.create_individual("ind3"), df)

            df_ind1_ind2 = ln.find_discordant_snps(ind1, ind2, save_output=True)
            df_ind1_ind2_ind3 = ln.find_discordant_snps(
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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))

            d = ln.find_shared_dna([ind1], shared_genes=True)

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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(ln.create_individual("ind2"))

            d = ln.find_shared_dna([ind1, ind2], genetic_map="test")

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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(ln.create_individual("ind2"))

            d = ln.find_shared_dna(
                [ind1, ind2], shared_genes=True, genetic_map="HapMap2"
            )

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            self._assert_shared_genes(d["one_chrom_shared_genes"], "default")
            self._assert_shared_genes(d["two_chrom_shared_genes"], "default")
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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"), pos_max=43800002)
            ind2 = self.simulate_snps(ln.create_individual("ind2"), pos_max=43800002)

            d = ln.find_shared_dna([ind1, ind2], shared_genes=True, genetic_map="CEU")

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            self._assert_shared_genes(d["one_chrom_shared_genes"], "1000G")
            self._assert_shared_genes(d["two_chrom_shared_genes"], "1000G")
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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(ln.create_individual("ind2"))
            ind3 = self.simulate_snps(ln.create_individual("ind3"))

            d = ln.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            self._assert_shared_genes(d["one_chrom_shared_genes"], "default")
            self._assert_shared_genes(d["two_chrom_shared_genes"], "default")
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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(ln.create_individual("ind2"))

            d = ln.find_shared_dna([ind1, ind2], shared_genes=True, save_output=False)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            self._assert_shared_genes(d["one_chrom_shared_genes"], "default")
            self._assert_shared_genes(d["two_chrom_shared_genes"], "default")
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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"), complement_genotype_one_chrom=True
            )

            d = ln.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 0
            self._assert_shared_genes(d["one_chrom_shared_genes"], "default")
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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"), complement_genotype_one_chrom=True
            )
            ind3 = self.simulate_snps(ln.create_individual("ind3"))

            d = ln.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 0
            self._assert_shared_genes(d["one_chrom_shared_genes"], "default")
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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(
                ln.create_individual("ind1"),
                chrom="X",
                pos_max=155270560,
                pos_step=1000,
                genotype="AA",
            )
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"),
                chrom="X",
                pos_max=155270560,
                pos_step=1000,
                genotype="AA",
            )

            d = ln.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
            assert len(d["two_chrom_shared_dna"]) == 1  # PAR1
            self._assert_shared_genes(
                d["one_chrom_shared_genes"], "X_chrom_male", key="len2"
            )
            self._assert_shared_genes(d["two_chrom_shared_genes"], "X_chrom_male")
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 202.022891, rtol=1e-5
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"], 20.837792, rtol=1e-3
            )
            self._make_file_exist_assertions("ind1_ind2", output_dir=tmpdir)

    def test_find_shared_dna_X_chrom_two_individuals_female(self):
        self.run_find_shared_dna_test_X(
            self._test_find_shared_dna_X_chrom_two_individuals_female
        )

    def _test_find_shared_dna_X_chrom_two_individuals_female(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(
                ln.create_individual("ind1"),
                chrom="X",
                pos_max=155270560,
                genotype="AC",
            )
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"),
                chrom="X",
                pos_max=155270560,
                genotype="AC",
            )

            d = ln.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
            assert len(d["two_chrom_shared_dna"]) == 1  # PAR1, non-PAR, PAR2
            self._assert_shared_genes(d["one_chrom_shared_genes"], "X_chrom_female")
            self._assert_shared_genes(
                d["two_chrom_shared_genes"], "X_chrom_female", key="len2"
            )
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            np.testing.assert_allclose(
                d["one_chrom_shared_dna"].loc[1]["cMs"], 202.022891, rtol=1e-5
            )
            np.testing.assert_allclose(
                d["two_chrom_shared_dna"].loc[1]["cMs"], 202.022891, rtol=1e-5
            )
            self._make_file_exist_assertions("ind1_ind2", output_dir=tmpdir)

    def test_find_shared_dna_two_chrom_shared_discrepant_snps(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_two_chrom_shared_discrepant_snps
        )

    def _test_find_shared_dna_two_chrom_shared_discrepant_snps(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)

            # simulate discrepant SNPs so that stitching of adjacent shared DNA segments is performed
            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"),
                complement_genotype_one_chrom=True,
                complement_snp_step=5000,
            )

            d = ln.find_shared_dna([ind1, ind2], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            self._assert_shared_genes(d["one_chrom_shared_genes"], "default")
            self._assert_shared_genes(d["two_chrom_shared_genes"], "default")
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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"), complement_genotype_two_chroms=True
            )

            d = ln.find_shared_dna([ind1, ind2], shared_genes=True)

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
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"), complement_genotype_two_chroms=True
            )
            ind3 = self.simulate_snps(ln.create_individual("ind3"))

            d = ln.find_shared_dna([ind1, ind2, ind3], shared_genes=True)

            assert len(d["one_chrom_shared_dna"]) == 0
            assert len(d["two_chrom_shared_dna"]) == 0
            assert len(d["one_chrom_shared_genes"]) == 0
            assert len(d["two_chrom_shared_genes"]) == 0
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0
            self._make_file_exist_assertions(
                "ind1_ind2_ind3", exist="plots", output_dir=tmpdir
            )

    def test_allele_comparisons(self):
        genotypes1 = pd.Series(["AG", "AG", "AA", "A", "A", "A", "ID", None])
        genotypes2 = pd.Series(["GA", "CT", "AG", "A", "G", "AG", "DI", "AA"])
        alleles1 = Lineage._get_alleles(genotypes1)
        alleles2 = Lineage._get_alleles(genotypes2)

        np.testing.assert_array_equal(
            Lineage._no_shared_allele(alleles1, alleles2),
            [False, True, False, False, True, False, False, True],
        )
        np.testing.assert_array_equal(
            Lineage._same_genotype(alleles1, alleles2),
            [True, False, False, True, False, False, True, False],
        )
        np.testing.assert_array_equal(
            Lineage._is_discordant(genotypes1, genotypes2),
            [False, True, False, False, True, False, False, False],
        )

    def _generate_test_cytoBand_hg19_bands(self):
        return pd.DataFrame(
            {
                "chrom": ["1", "1", "1"],
                "start": [0, 50000000, 52000000],
                "end": [50000000, 52000000, 111800001],
                "name": ["p36.33", "p36.32", "q11"],
                "gie_stain": ["gneg", "gpos25", "gneg"],
            }
        )

    def _find_shared_dna_exclude_regions(self, exclude_regions, cytobands=None):
        if cytobands is None:
            cytobands = self._generate_test_cytoBand_hg19_bands()

        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(ln.create_individual("ind2"))

            with patch(
                "lineage.resources.Resources.get_cytoBand_hg19",
                Mock(return_value=cytobands),
            ):
                return ln.find_shared_dna(
                    [ind1, ind2], save_output=False, exclude_regions=exclude_regions
                )

    def test_find_shared_dna_exclude_regions(self):
        self.run_find_shared_dna_test(self._test_find_shared_dna_exclude_regions)

    def _test_find_shared_dna_exclude_regions(self):
        for exclude_regions in [
            [("1", 50000000, 52000000)],
            [("chr1", 50000000, 52000000)],
            (region for region in [("1", 50000000, 52000000)]),
            ["1p36.32"],
        ]:
            d = self._find_shared_dna_exclude_regions(exclude_regions)

            # the shared DNA segment is split at the excluded region
            for key in ["one_chrom_shared_dna", "two_chrom_shared_dna"]:
                assert len(d[key]) == 2
                assert (d[key]["end"] < 50000000).sum() == 1
                assert (d[key]["start"] > 52000000).sum() == 1
                assert d[key]["cMs"].sum() < 140.443968
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0

    def test_find_shared_dna_exclude_regions_sub_bands(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_exclude_regions_sub_bands
        )

    def _test_find_shared_dna_exclude_regions_sub_bands(self):
        # "1p36.3" includes sub-bands 1p36.33 and 1p36.32
        d = self._find_shared_dna_exclude_regions(["1p36.3"])

        for key in ["one_chrom_shared_dna", "two_chrom_shared_dna"]:
            assert len(d[key]) == 1
            assert d[key].loc[1]["start"] > 52000000

    def test_find_shared_dna_exclude_regions_invalid(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_exclude_regions_invalid
        )

    def _test_find_shared_dna_exclude_regions_invalid(self):
        for region in [
            "chr1p36",
            "1p12",
            ("1", 50000000),
            ("Z", 50000000, 52000000),
            ("1", "50000000", "52000000"),
            ("1", 50000000.0, 52000000.0),
            ("1", 0, 52000000),
            ("1", 52000000, 50000000),
            ("1", True, 52000000),
        ]:
            with self.assertRaises(ValueError) as cm:
                self._find_shared_dna_exclude_regions([region])
            assert str(region) in str(cm.exception)

        # a single region instead of an iterable of regions
        for exclude_regions in [
            "1p36.32",
            ("1", 50000000, 52000000),
            ["1", 50000000, 52000000],
            ("1", "50000000", "52000000"),
            ("1", 50000000.0, 52000000.0),
            {"chrom": "1", "start": 50000000, "end": 52000000},
        ]:
            with self.assertRaises(ValueError) as cm:
                self._find_shared_dna_exclude_regions(exclude_regions)
            assert "single region" in str(cm.exception)
            assert str(exclude_regions) in str(cm.exception)

        with self.assertRaises(TypeError) as cm:
            self._find_shared_dna_exclude_regions(None)
        assert "None" in str(cm.exception)

    def test_find_shared_dna_exclude_regions_no_cytobands(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(ln.create_individual("ind2"))

            def assert_cytobands_not_available():
                with self.assertRaises(ValueError) as cm:
                    ln.find_shared_dna(
                        [ind1, ind2], save_output=False, exclude_regions=["8p23.1"]
                    )
                assert "8p23.1" in str(cm.exception)
                assert "(chrom, start, end)" in str(cm.exception)

            # the cytoBand table is empty
            with patch(
                "lineage.resources.Resources.get_cytoBand_hg19",
                Mock(return_value=pd.DataFrame()),
            ):
                assert_cytobands_not_available()

            # the cytoBand table couldn't be downloaded, so there's no file to load
            with (
                patch(
                    "lineage.resources.Resources._get_path_cytoBand_hg19",
                    Mock(return_value=""),
                ),
                patch.object(ln._resources, "_cytoBand_hg19", pd.DataFrame()),
            ):
                assert_cytobands_not_available()

    def test_find_shared_dna_exclude_regions_one_snp(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_exclude_regions_one_snp
        )

    def _test_find_shared_dna_exclude_regions_one_snp(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(ln.create_individual("ind2"))

            # a discordant SNP at position 20000001, so adjacent segments are stitched
            ind2._snps.loc["rs2001", "genotype"] = "CC"

            # the region contains only the SNP at position 50000001
            d = ln.find_shared_dna(
                [ind1, ind2],
                save_output=False,
                exclude_regions=[("1", 50000001, 50000001)],
            )

            # the shared DNA segment is stitched across the discordant SNP, but it's split at
            # the SNP in the region, which isn't reported as discrepant
            for key in ["one_chrom", "two_chrom"]:
                shared_dna = d[f"{key}_shared_dna"]
                assert len(shared_dna) == 2
                assert shared_dna.loc[1]["start"] == 1
                assert shared_dna.loc[1]["end"] == 49990001
                assert shared_dna.loc[2]["start"] == 50010001
                assert d[f"{key}_discrepant_snps"].tolist() == ["rs2001"]

    def test_find_shared_dna_exclude_regions_band_start(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_exclude_regions_band_start
        )

    def _test_find_shared_dna_exclude_regions_band_start(self):
        # cytoBand starts are 0-based, so a band with start 50000001 begins after the SNP at
        # position 50000001, and a band with start 50000000 begins at that SNP
        for band_start, segment_end in [(50000001, 50000001), (50000000, 49990001)]:
            cytobands = self._generate_test_cytoBand_hg19_bands()
            cytobands.loc[1, "start"] = band_start

            d = self._find_shared_dna_exclude_regions(["1p36.32"], cytobands)

            for key in ["one_chrom_shared_dna", "two_chrom_shared_dna"]:
                assert len(d[key]) == 2
                assert d[key].loc[1]["end"] == segment_end
                assert d[key].loc[2]["start"] == 52000001

    def test_find_shared_dna_exclude_regions_X_chrom(self):
        self.run_find_shared_dna_test_X(
            self._test_find_shared_dna_exclude_regions_X_chrom
        )

    def _test_find_shared_dna_exclude_regions_X_chrom(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.simulate_snps(
                ln.create_individual("ind1"),
                chrom="X",
                pos_max=155270560,
                genotype="AC",
            )
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"),
                chrom="X",
                pos_max=155270560,
                genotype="AC",
            )

            d = ln.find_shared_dna(
                [ind1, ind2],
                save_output=False,
                exclude_regions=[("chrX", 50000000, 52000000)],
            )

            for key in ["one_chrom_shared_dna", "two_chrom_shared_dna"]:
                assert len(d[key]) == 2
                assert d[key].loc[1]["end"] == 49990001
                assert d[key].loc[2]["start"] == 52000001
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0

    def test_find_shared_dna_exclude_regions_pileup_regions(self):
        chroms = ["6", "8", "17"]
        pos = np.arange(1, 60000001, 2000, dtype=np.int64)
        df = pd.DataFrame(
            {
                "chrom": np.repeat(chroms, len(pos)),
                "pos": np.tile(pos, len(chroms)),
                "ind1": "AA",
                "ind2": "AA",
            },
            index=pd.Index(
                [f"rs{i + 1}" for i in range(len(chroms) * len(pos))], name="rsid"
            ),
        )

        genetic_map = {}
        for chrom in chroms:
            genetic_map.update(
                self._generate_test_genetic_map(
                    chrom=chrom, pos=(1, 60000001), rate=(1, 0), map_cMs=(0, 60)
                )
            )

        # synthetic cytobands (not the UCSC coordinates of these bands)
        cytobands = pd.DataFrame(
            {
                "chrom": ["8", "8", "8", "17", "17", "17"],
                "start": [2000000, 6000000, 12000000, 38000000, 40000000, 45000000],
                "end": [6000000, 12000000, 16000000, 40000000, 45000000, 46000000],
                "name": ["p23.2", "p23.1", "p22", "q21.2", "q21.31", "q21.32"],
                "gie_stain": ["gpos75", "gneg", "gpos100", "gneg", "gpos25", "gneg"],
            }
        )

        with (
            tempfile.TemporaryDirectory() as tmpdir,
            patch(
                "lineage.resources.Resources.get_genetic_map_HapMapII_GRCh37",
                Mock(return_value=genetic_map),
            ),
            patch(
                "lineage.resources.Resources.get_cytoBand_hg19",
                Mock(return_value=cytobands),
            ),
        ):
            ln = Lineage(output_dir=tmpdir)

            ind1 = self.get_discordant_snps(ln.create_individual("ind1"), df)
            ind2 = self.get_discordant_snps(ln.create_individual("ind2"), df)

            # by default, shared DNA spans each chromosome
            d = ln.find_shared_dna([ind1, ind2], save_output=False)

            for key in ["one_chrom_shared_dna", "two_chrom_shared_dna"]:
                assert d[key]["chrom"].tolist() == chroms
                assert d[key]["start"].tolist() == [1, 1, 1]
                assert d[key]["end"].tolist() == [59998001, 59998001, 59998001]

            # shared DNA is split at the MHC region (chr6:28477797-33448354) and at the
            # synthetic bands 8p23.1 (chr8:6000001-12000000) and 17q21.31
            # (chr17:40000001-45000000)
            d = ln.find_shared_dna(
                [ind1, ind2], save_output=False, exclude_regions=PILEUP_REGIONS
            )

            for key in ["one_chrom_shared_dna", "two_chrom_shared_dna"]:
                assert d[key]["chrom"].tolist() == ["6", "6", "8", "8", "17", "17"]
                assert d[key]["start"].tolist() == [
                    1,
                    33450001,
                    1,
                    12000001,
                    1,
                    45000001,
                ]
                assert d[key]["end"].tolist() == [
                    28476001,
                    59998001,
                    5998001,
                    59998001,
                    39998001,
                    59998001,
                ]
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0

    def test_find_shared_dna_exclude_strand_ambiguous(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_exclude_strand_ambiguous
        )

    def _test_find_shared_dna_exclude_strand_ambiguous(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)

            # complemented genotypes are A/T SNPs, which are found as discrepant SNPs unless
            # strand-ambiguous SNPs are excluded
            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(
                ln.create_individual("ind2"),
                complement_genotype_one_chrom=True,
                complement_snp_step=5000,
            )

            d = ln.find_shared_dna(
                [ind1, ind2], save_output=False, exclude_strand_ambiguous=True
            )

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0

    def test_find_shared_dna_exclude_strand_ambiguous_homozygous(self):
        self.run_find_shared_dna_test(
            self._test_find_shared_dna_exclude_strand_ambiguous_homozygous
        )

    def _test_find_shared_dna_exclude_strand_ambiguous_homozygous(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir)

            # AA and TT can't be distinguished from a strand flip, so the SNP is
            # strand-ambiguous even though the genotypes differ
            ind1 = self.simulate_snps(ln.create_individual("ind1"))
            ind2 = self.simulate_snps(ln.create_individual("ind2"))
            ind2._snps.loc["rs5001", "genotype"] = "TT"

            # by default, segments are stitched across the SNP, which is discrepant
            d = ln.find_shared_dna([ind1, ind2], save_output=False)

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            assert d["one_chrom_discrepant_snps"].tolist() == ["rs5001"]
            assert d["two_chrom_discrepant_snps"].tolist() == ["rs5001"]

            d = ln.find_shared_dna(
                [ind1, ind2], save_output=False, exclude_strand_ambiguous=True
            )

            assert len(d["one_chrom_shared_dna"]) == 1
            assert len(d["two_chrom_shared_dna"]) == 1
            assert len(d["one_chrom_discrepant_snps"]) == 0
            assert len(d["two_chrom_discrepant_snps"]) == 0

    def test_find_discordant_snps_exclude_strand_ambiguous(self):
        ln = Lineage()
        df = pd.DataFrame(
            {
                "chrom": "1",
                "pos": [1, 2, 3, 4],
                # rs1 (A/T) and rs3 (C/G) are strand-ambiguous; rs2 (A/G) is not
                "ind1": ["AA", "AA", "CC", "AG"],
                "ind2": ["TT", "GG", "GG", "AG"],
                "ind3": ["TT", "GG", "GG", "AG"],
            },
            index=pd.Index(["rs1", "rs2", "rs3", "rs4"], name="rsid"),
        )

        ind1 = self.get_discordant_snps(ln.create_individual("ind1"), df)
        ind2 = self.get_discordant_snps(ln.create_individual("ind2"), df)
        ind3 = self.get_discordant_snps(ln.create_individual("ind3"), df)

        pd.testing.assert_index_equal(
            ln.find_discordant_snps(ind1, ind2).index,
            pd.Index(["rs1", "rs2", "rs3"], name="rsid"),
        )
        for result in [
            ln.find_discordant_snps(ind1, ind2, exclude_strand_ambiguous=True),
            ln.find_discordant_snps(ind1, ind2, ind3, exclude_strand_ambiguous=True),
        ]:
            pd.testing.assert_index_equal(result.index, pd.Index(["rs2"], name="rsid"))


class TestCreateExampleDatasets(BaseLineageTestCase):
    """Tests for the create_example_datasets method."""

    def test_create_example_datasets(self):
        """Test that create_example_datasets generates valid files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ln = Lineage(output_dir=tmpdir, resources_dir=tmpdir)
            paths = ln.create_example_datasets()

            # Verify all files created
            assert "parent" in paths
            assert "child" in paths
            assert "sibling1" in paths
            assert "sibling2" in paths

            for key, path in paths.items():
                assert os.path.exists(path), f"File not created for {key}: {path}"
