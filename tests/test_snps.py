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

from lineage import Lineage, SNPs
from tests import BaseLineageTestCase


class TestSnps(BaseLineageTestCase):
    def setUp(self):
        self.l = Lineage()
        self.snps_GRCh38 = SNPs("tests/input/GRCh38.csv")
        self.snps = SNPs("tests/input/chromosomes.csv")
        self.snps_none = SNPs(None)
        self.del_output_dir_helper()

    def snps_discrepant_pos(self):
        return self.create_snp_df(
            rsid=["rs3094315"], chrom=["1"], pos=[1], genotype=["AA"]
        )

    def test_assembly(self):
        assert self.snps_GRCh38.assembly == "GRCh38"

    def test_assembly_no_snps(self):
        assert self.snps_none.assembly == ""

    def test_snp_count(self):
        assert self.snps.snp_count == 6

    def test_snp_count_no_snps(self):
        assert self.snps_none.snp_count == 0

    def test_chromosomes(self):
        assert self.snps.chromosomes == ["1", "2", "3", "5", "PAR", "MT"]

    def test_chromosomes_no_snps(self):
        assert self.snps_none.chromosomes == []

    def test_chromosomes_summary(self):
        assert self.snps.chromosomes_summary == "1-3, 5, PAR, MT"

    def test_chromosomes_summary_no_snps(self):
        assert self.snps_none.chromosomes_summary == ""

    def test_build_no_snps(self):
        assert self.snps_none.build is None

    def test_build_detected_no_snps(self):
        assert not self.snps_none.build_detected

    def test_build_detected_PAR_snps(self):
        snps = SNPs("tests/input/GRCh37_PAR.csv")
        assert snps.build == 37
        assert snps.build_detected

    def test_sex_no_snps(self):
        assert self.snps_none.sex == ""

    def test_sex_Male_Y_chrom(self):
        ind = self.simulate_snps(
            self.l.create_individual("test_snps_sex_Male_Y_chrom"),
            chrom="Y",
            pos_start=1,
            pos_max=59373566,
            pos_step=10000,
        )
        file = ind.save_snps()
        from lineage.snps import SNPs

        snps = SNPs(file)
        assert snps.sex == "Male"

    def test_get_summary(self):
        assert self.snps_GRCh38.get_summary() == {
            "source": "generic",
            "assembly": "GRCh38",
            "build": 38,
            "build_detected": True,
            "snp_count": 4,
            "chromosomes": "1, 3",
            "sex": "",
        }

    def test_get_summary_no_snps(self):
        assert self.snps_none.get_summary() is None

    def test_is_valid_True(self):
        assert self.snps_GRCh38.is_valid()

    def test_is_valid_False(self):
        assert not self.snps_none.is_valid()

    def test__read_raw_data(self):
        assert self.snps_none.snps is None
        assert self.snps_none.source == ""

    def test__lookup_build_with_snp_pos_None(self):
        from lineage.snps import detect_build

        assert detect_build(self.snps_discrepant_pos()) is None

    def test_get_assembly_None(self):
        from lineage.snps import get_assembly

        assert get_assembly(None) is ""
