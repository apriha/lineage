"""
Copyright (C) 2017 Andrew Riha

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

import gzip
import os
import shutil
import zipfile

import numpy as np
import pandas as pd
import pytest

from lineage import sort_snps
from tests import BaseLineageTestCase


class TestIndividual(BaseLineageTestCase):
    def generic_snps(self):
        return self.create_snp_df(
            rsid=["rs1", "rs2", "rs3", "rs4", "rs5"],
            chrom=["1", "1", "1", "1", "1"],
            pos=[1, 2, 3, 4, 5],
            genotype=["AA", "CC", "GG", "TT", np.nan],
        )

    def snps_NCBI36(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[742429, 143649677, 143649678, 50908372],
            genotype=["AA", np.nan, "ID", "AG"],
        )

    def snps_NCBI36_discrepant_snps(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[742429, 143649677, 143649678, 50908372],
            genotype=["AA", np.nan, "ID", np.nan],
        )

    def snps_GRCh37(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[752566, 144938320, 144938321, 50927009],
            genotype=["AA", np.nan, "ID", "TC"],
        )

    def snps_GRCh38(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rsIndelTest", "rs2500347", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[817186, 148946168, 148946169, 50889578],
            genotype=["AA", "ID", np.nan, "TC"],
        )

    def snps_GRCh38_PAR(self):
        return self.create_snp_df(
            rsid=["rs28736870", "rs113313554"],
            chrom=["X", "Y"],
            pos=[304103, 624523],
            genotype=["AA", "AA"],
        )

    def test_name(self):
        ind = self.l.create_individual("test")
        assert ind.name == "test"

    def test_snps_23andme(self):
        # https://www.23andme.com
        ind = self.l.create_individual("", "tests/input/23andme.txt")
        assert ind.source == "23andMe"
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps())

    def test_snps_23andme_zip(self):
        with zipfile.ZipFile("tests/input/23andme.txt.zip", "w") as f:
            # https://stackoverflow.com/a/16104667
            f.write("tests/input/23andme.txt", arcname="23andme.txt")
        ind = self.l.create_individual("", "tests/input/23andme.txt.zip")
        assert ind.source == "23andMe"
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps())

    def test_snps_ftdna(self):
        # https://www.familytreedna.com
        ind = self.l.create_individual("", "tests/input/ftdna.csv")
        assert ind.source == "FTDNA"
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps())

    def test_snps_ftdna_gzip(self):
        with open("tests/input/ftdna.csv", "rb") as f_in:
            with gzip.open("tests/input/ftdna.csv.gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        ind = self.l.create_individual("", "tests/input/ftdna.csv.gz")
        assert ind.source == "FTDNA"
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps())

    def test_snps_ftdna_famfinder(self):
        # https://www.familytreedna.com
        ind = self.l.create_individual("", "tests/input/ftdna_famfinder.csv")
        assert ind.source == "FTDNA"
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps())

    def test_snps_ancestry(self):
        # https://www.ancestry.com
        ind = self.l.create_individual("", "tests/input/ancestry.txt")
        assert ind.source == "AncestryDNA"
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps())

    def test_source_lineage_file(self):
        ind = self.l.create_individual("", "tests/input/GRCh37.csv")
        assert ind.source == "generic"
        ind.load_snps("tests/input/23andme.txt")
        assert ind.source == "generic, 23andMe"
        file = ind.save_snps()
        ind_saved_snps = self.l.create_individual("", file)
        assert ind_saved_snps.source == "generic, 23andMe"
        pd.testing.assert_frame_equal(ind.snps, ind_saved_snps.snps)

    def test_source_lineage_file_gzip(self):
        ind = self.l.create_individual("", "tests/input/GRCh37.csv")
        assert ind.source == "generic"
        ind.load_snps("tests/input/23andme.txt")
        assert ind.source == "generic, 23andMe"
        file = ind.save_snps()
        with open(file, "rb") as f_in:
            with gzip.open(file + ".gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        ind_saved_snps = self.l.create_individual("", file + ".gz")
        assert ind_saved_snps.source == "generic, 23andMe"
        pd.testing.assert_frame_equal(ind.snps, ind_saved_snps.snps)

    def test_source_generic(self):
        ind = self.l.create_individual("", "tests/input/NCBI36.csv")
        assert ind.source == "generic"

    def test_snps_None(self):
        ind = self.l.create_individual("")
        assert ind.snps is None

    def test_snp_count(self):
        ind = self.l.create_individual("", "tests/input/NCBI36.csv")
        assert ind.snp_count == 4

    def test_snp_count_None(self):
        ind = self.l.create_individual("")
        assert ind.snp_count == 0

    def test_chromosomes(self):
        ind = self.l.create_individual("", "tests/input/chromosomes.csv")
        assert ind.chromosomes == ["1", "2", "3", "5", "PAR", "MT"]

    def test_chromosomes_None(self):
        ind = self.l.create_individual("")
        assert ind.chromosomes == []

    def test_chromosomes_summary(self):
        ind = self.l.create_individual("", "tests/input/chromosomes.csv")
        assert ind.chromosomes_summary == "1-3, 5, PAR, MT"

    def test_chromosomes_summary_None(self):
        ind = self.l.create_individual("")
        assert ind.chromosomes_summary == ""

    def test_build(self):
        ind = self.l.create_individual("", "tests/input/NCBI36.csv")
        assert ind.build == 36
        assert ind.assembly == "NCBI36"

    def test_load_snps_list(self):
        ind = self.l.create_individual("")
        ind.load_snps(["tests/input/GRCh37.csv", "tests/input/GRCh37.csv"])
        pd.testing.assert_frame_equal(ind.snps, self.snps_GRCh37())
        assert ind.source == "generic, generic"

    def test_load_snps_None(self):
        ind = self.l.create_individual("")
        with pytest.raises(TypeError):
            ind.load_snps(None)

    def test_sex_Male_Y_chrom(self):
        ind = self.simulate_snps(
            self.l.create_individual(""),
            chrom="Y",
            pos_start=1,
            pos_max=59373566,
            pos_step=10000,
        )
        assert ind.sex == "Male"

    def test_sex_Female_Y_chrom(self):
        ind = self.simulate_snps(
            self.l.create_individual(""),
            chrom="Y",
            pos_start=1,
            pos_max=59373566,
            pos_step=10000,
            null_snp_step=1,
        )
        assert ind.sex == "Female"

    def test_sex_Female_X_chrom(self):
        ind = self.simulate_snps(
            self.l.create_individual(""),
            chrom="X",
            pos_start=1,
            pos_max=155270560,
            pos_step=10000,
            genotype="AC",
        )
        assert ind.sex == "Female"

    def test_sex_Male_X_chrom(self):
        ind = self.simulate_snps(
            self.l.create_individual(""),
            chrom="X",
            pos_start=1,
            pos_max=155270560,
            pos_step=10000,
            genotype="AA",
        )
        assert ind.sex == "Male"

    def test_sex_not_determined(self):
        ind = self.simulate_snps(
            self.l.create_individual(""),
            chrom="1",
            pos_start=1,
            pos_max=249250621,
            pos_step=10000,
        )
        assert ind.sex == ""

    def test_discrepant_positions(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_positions) == 4

    def test_discrepant_genotypes(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_genotypes) == 1

    def test_discrepant_snps(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_snps) == 4

    def test_load_snps_non_existent_file(self):
        ind = self.l.create_individual("")
        ind.load_snps(["tests/input/GRCh37.csv", "tests/input/non_existent_file.csv"])
        pd.testing.assert_frame_equal(ind.snps, self.snps_GRCh37())

    def test_load_snps_invalid_file(self):
        ind = self.l.create_individual("")
        with open("tests/input/empty.txt", "w"):
            pass
        ind.load_snps(["tests/input/GRCh37.csv", "tests/input/empty.txt"])
        pd.testing.assert_frame_equal(ind.snps, self.snps_GRCh37())

    def test_load_snps_assembly_mismatch(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert not os.path.exists("output/ind_discrepant_positions_1.csv")
        assert not os.path.exists("output/ind_discrepant_genotypes_1.csv")
        assert len(ind.discrepant_positions) == 4
        assert len(ind.discrepant_genotypes) == 1
        pd.testing.assert_frame_equal(ind.snps, self.snps_NCBI36_discrepant_snps())

    def test_load_snps_assembly_mismatch_save_output(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"], save_output=True
        )
        assert os.path.exists("output/ind_discrepant_positions_1.csv")
        assert os.path.exists("output/ind_discrepant_genotypes_1.csv")
        assert len(ind.discrepant_positions) == 4
        assert len(ind.discrepant_genotypes) == 1
        pd.testing.assert_frame_equal(ind.snps, self.snps_NCBI36_discrepant_snps())

    def test_load_snps_assembly_mismatch_exceed_discrepant_positions_threshold(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"],
            discrepant_snp_positions_threshold=0,
        )
        assert not os.path.exists("output/ind_discrepant_positions_1.csv")
        assert not os.path.exists("output/ind_discrepant_genotypes_1.csv")
        assert len(ind.discrepant_positions) == 4
        assert len(ind.discrepant_genotypes) == 0
        pd.testing.assert_frame_equal(ind.snps, self.snps_NCBI36())

    def test_load_snps_assembly_mismatch_exceed_discrepant_genotypes_threshold(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"],
            discrepant_genotypes_threshold=0,
        )
        assert not os.path.exists("output/ind_discrepant_positions_1.csv")
        assert not os.path.exists("output/ind_discrepant_genotypes_1.csv")
        assert len(ind.discrepant_positions) == 4
        assert len(ind.discrepant_genotypes) == 1
        pd.testing.assert_frame_equal(ind.snps, self.snps_NCBI36())

    def test_merging_files_discrepant_snps(self):
        df = pd.read_csv(
            "tests/input/discrepant_snps.csv",
            skiprows=1,
            na_values="--",
            names=[
                "rsid",
                "chrom",
                "pos_file1",
                "pos_file2",
                "genotype_file1",
                "genotype_file2",
                "discrepant_position",
                "discrepant_genotype",
                "expected_position",
                "expected_genotype",
            ],
            index_col=0,
            dtype={
                "chrom": object,
                "pos_file1": np.int64,
                "pos_file2": np.int64,
                "discrepant_position": bool,
                "discrepant_genotype": bool,
            },
        )

        df1 = df[["chrom", "pos_file1", "genotype_file1"]]
        df2 = df[["chrom", "pos_file2", "genotype_file2"]]

        df1.to_csv(
            "tests/input/discrepant_snps1.csv",
            na_rep="--",
            header=["chromosome", "position", "genotype"],
        )

        df2.to_csv(
            "tests/input/discrepant_snps2.csv",
            na_rep="--",
            header=["chromosome", "position", "genotype"],
        )

        ind = self.l.create_individual(
            "", ["tests/input/discrepant_snps1.csv", "tests/input/discrepant_snps2.csv"]
        )

        expected = df[
            [
                "chrom",
                "discrepant_position",
                "discrepant_genotype",
                "expected_position",
                "expected_genotype",
            ]
        ]
        expected = expected.rename(
            columns={"expected_position": "pos", "expected_genotype": "genotype"}
        )
        expected = sort_snps(expected)

        pd.testing.assert_index_equal(
            ind.discrepant_positions.index,
            expected.loc[expected["discrepant_position"] == True].index,
        )

        pd.testing.assert_index_equal(
            ind.discrepant_genotypes.index,
            expected.loc[expected["discrepant_genotype"] == True].index,
        )

        pd.testing.assert_series_equal(ind.snps["pos"], expected["pos"])
        pd.testing.assert_series_equal(ind.snps["genotype"], expected["genotype"])

    def test_save_snps(self):
        ind = self.l.create_individual("test save snps", "tests/input/GRCh37.csv")
        assert (
            os.path.relpath(ind.save_snps())
            == "output/test_save_snps_lineage_GRCh37.csv"
        )
        ind_saved_snps = self.l.create_individual(
            "", "output/test_save_snps_lineage_GRCh37.csv"
        )
        pd.testing.assert_frame_equal(ind_saved_snps.snps, self.snps_GRCh37())

    def test_save_snps_specify_file(self):
        ind = self.l.create_individual("test save snps", "tests/input/GRCh37.csv")
        assert os.path.relpath(ind.save_snps("snps.csv")) == "output/snps.csv"
        ind_saved_snps = self.l.create_individual("", "output/snps.csv")
        pd.testing.assert_frame_equal(ind_saved_snps.snps, self.snps_GRCh37())

    def test_save_snps_no_snps(self):
        ind = self.l.create_individual("")
        assert not ind.save_snps()

    def test_save_snps_invalid_output_dir(self):
        ind = self.l.create_individual("", "tests/input/GRCh37.csv")
        ind._output_dir = None
        assert not ind.save_snps()

    def test_save_snps_exception(self):
        ind = self.l.create_individual("")
        ind._snps = "invalid"
        assert not ind.save_snps()

    def test_save_discrepant_positions(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_positions) == 4
        discrepant_positions_file = ind.save_discrepant_positions()
        assert (
            os.path.relpath(discrepant_positions_file)
            == "output/ind_discrepant_positions.csv"
        )
        assert os.path.exists(discrepant_positions_file)

    def test_save_discrepant_positions_specify_file(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_positions) == 4
        discrepant_positions_file = ind.save_discrepant_positions(
            "discrepant_positions.csv"
        )
        assert (
            os.path.relpath(discrepant_positions_file)
            == "output/discrepant_positions.csv"
        )
        assert os.path.exists(discrepant_positions_file)

    def test_save_discrepant_positions_no_discrepant_snps(self):
        ind = self.l.create_individual("ind")
        assert len(ind.discrepant_positions) == 0
        assert not ind.save_discrepant_positions()

    def test_save_discrepant_positions_invalid_output_dir(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_positions) == 4
        ind._output_dir = None
        assert not ind.save_discrepant_positions()

    def test_save_discrepant_positions_exception(self):
        ind = self.l.create_individual("")
        ind._discrepant_positions = "invalid"
        assert not ind.save_discrepant_positions()

    def test_save_discrepant_genotypes(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_genotypes) == 1
        discrepant_genotypes_file = ind.save_discrepant_genotypes()
        assert (
            os.path.relpath(discrepant_genotypes_file)
            == "output/ind_discrepant_genotypes.csv"
        )
        assert os.path.exists(discrepant_genotypes_file)

    def test_save_discrepant_genotypes_specify_file(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_genotypes) == 1
        discrepant_genotypes_file = ind.save_discrepant_genotypes(
            "discrepant_genotypes.csv"
        )
        assert (
            os.path.relpath(discrepant_genotypes_file)
            == "output/discrepant_genotypes.csv"
        )
        assert os.path.exists(discrepant_genotypes_file)

    def test_save_discrepant_genotypes_no_discrepant_snps(self):
        ind = self.l.create_individual("ind")
        assert len(ind.discrepant_genotypes) == 0
        assert not ind.save_discrepant_genotypes()

    def test_save_discrepant_genotypes_invalid_output_dir(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_genotypes) == 1
        ind._output_dir = None
        assert not ind.save_discrepant_genotypes()

    def test_save_discrepant_genotypes_exception(self):
        ind = self.l.create_individual("")
        ind._discrepant_genotypes = "invalid"
        assert not ind.save_discrepant_genotypes()

    def test_save_discrepant_snps(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_snps) == 4
        discrepant_snps_file = ind.save_discrepant_snps()
        assert os.path.relpath(discrepant_snps_file) == "output/ind_discrepant_snps.csv"
        assert os.path.exists(discrepant_snps_file)

    def test_save_discrepant_snps_specify_file(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_snps) == 4
        discrepant_snps_file = ind.save_discrepant_snps("discrepant_snps.csv")
        assert os.path.relpath(discrepant_snps_file) == "output/discrepant_snps.csv"
        assert os.path.exists(discrepant_snps_file)

    def test_save_discrepant_snps_no_discrepant_snps(self):
        ind = self.l.create_individual("ind")
        assert len(ind.discrepant_snps) == 0
        assert not ind.save_discrepant_snps()

    def test_save_discrepant_snps_invalid_output_dir(self):
        ind = self.l.create_individual("ind")
        ind.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(ind.discrepant_snps) == 4
        ind._output_dir = None
        assert not ind.save_discrepant_snps()

    def test_save_discrepant_snps_exception(self):
        ind = self.l.create_individual("")
        ind._discrepant_snps = "invalid"
        assert not ind.save_discrepant_snps()

    def test_get_var_name(self):
        ind = self.l.create_individual("test a. name")
        assert ind.get_var_name() == "test_a__name"

    def test_remap_snps_36_to_37(self):
        ind = self.l.create_individual("", "tests/input/NCBI36.csv")
        chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(37)
        assert ind.build == 37
        assert ind.assembly == "GRCh37"
        if len(chromosomes_remapped) == 2:
            assert len(chromosomes_not_remapped) == 0
            pd.testing.assert_frame_equal(ind.snps, self.snps_GRCh37())

    def test_remap_snps_37_to_36(self):
        ind = self.l.create_individual("", "tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(36)
        assert ind.build == 36
        assert ind.assembly == "NCBI36"
        if len(chromosomes_remapped) == 2:
            assert len(chromosomes_not_remapped) == 0
            pd.testing.assert_frame_equal(ind.snps, self.snps_NCBI36())

    def test_remap_snps_37_to_38(self):
        ind = self.l.create_individual("", "tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(38)
        assert ind.build == 38
        assert ind.assembly == "GRCh38"
        if len(chromosomes_remapped) == 2:
            assert len(chromosomes_not_remapped) == 0
            pd.testing.assert_frame_equal(ind.snps, self.snps_GRCh38())

    def test_remap_snps_37_to_38_via_lineage(self):
        ind = self.l.create_individual("", "tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = self.l.remap_snps(ind, 38)
        assert ind.build == 38
        assert ind.assembly == "GRCh38"
        if len(chromosomes_remapped) == 2:
            assert len(chromosomes_not_remapped) == 0
            pd.testing.assert_frame_equal(ind.snps, self.snps_GRCh38())

    def test_remap_snps_37_to_38_with_PAR_SNP(self):
        ind = self.l.create_individual("", "tests/input/GRCh37_PAR.csv")
        assert ind.snp_count == 3
        chromosomes_remapped, chromosomes_not_remapped = self.l.remap_snps(ind, 38)
        assert ind.build == 38
        assert ind.assembly == "GRCh38"
        if len(chromosomes_remapped) == 2:
            assert len(chromosomes_not_remapped) == 1
            assert ind.snp_count == 2
            pd.testing.assert_frame_equal(ind.snps, self.snps_GRCh38_PAR())

    def test_remap_snps_37_to_37(self):
        ind = self.l.create_individual("", "tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(37)
        assert ind.build == 37
        assert ind.assembly == "GRCh37"
        assert len(chromosomes_remapped) == 0
        assert len(chromosomes_not_remapped) == 2
        pd.testing.assert_frame_equal(ind.snps, self.snps_GRCh37())

    def test_remap_snps_no_snps(self):
        ind = self.l.create_individual("")
        chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(38)
        assert ind.build is None
        assert len(chromosomes_remapped) == 0
        assert len(chromosomes_not_remapped) == 0

    def test_remap_snps_invalid_assembly(self):
        ind = self.l.create_individual("", "tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(-1)
        assert ind.build == 37
        assert ind.assembly == "GRCh37"
        assert len(chromosomes_remapped) == 0
        assert len(chromosomes_not_remapped) == 2

    def test___repr__(self):
        ind = self.l.create_individual("test")
        assert "Individual('test')" == ind.__repr__()
