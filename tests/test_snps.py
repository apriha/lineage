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

import pandas as pd
import pytest

from tests.test_lineage import simulate_snps

def create_snp_df(rsid, chrom, pos, genotype):
    df = pd.DataFrame({'rsid': rsid, 'chrom': chrom, 'pos': pos, 'genotype': genotype},
                      columns=['rsid', 'chrom', 'pos', 'genotype'])
    df = df.set_index('rsid')
    return df


@pytest.fixture(scope='module')
def snps_discrepant_pos():
    return create_snp_df(rsid=['rs3094315'], chrom=['1'], pos=[1], genotype=['AA'])


@pytest.fixture(scope='module')
def snps_GRCh38():
    from lineage.snps import SNPs
    return SNPs('tests/input/GRCh38.csv')


@pytest.fixture(scope='module')
def snps():
    from lineage.snps import SNPs
    return SNPs('tests/input/chromosomes.csv')


@pytest.fixture(scope='module')
def snps_none():
    from lineage.snps import SNPs
    return SNPs(None)


def test_assembly(snps_GRCh38):
    assert snps_GRCh38.assembly == 'GRCh38'


def test_assembly_no_snps(snps_none):
    assert snps_none.assembly == ''


def test_snp_count(snps):
    assert snps.snp_count == 6


def test_snp_count_no_snps(snps_none):
    assert snps_none.snp_count == 0


def test_chromosomes(snps):
    assert snps.chromosomes == ['1', '2', '3', '5', 'PAR', 'MT']


def test_chromosomes_no_snps(snps_none):
    assert snps_none.chromosomes == []


def test_chromosomes_summary(snps):
    assert snps.chromosomes_summary == '1-3, 5, PAR, MT'


def test_chromosomes_summary_no_snps(snps_none):
    assert snps_none.chromosomes_summary == ''


def test_build_no_snps(snps_none):
    assert snps_none.build is None


def test_build_detected_no_snps(snps_none):
    assert snps_none.build_detected == False


def test_build_detected_PAR_snps():
    from lineage.snps import SNPs
    snps = SNPs('tests/input/GRCh37_PAR.csv')
    assert snps.build == 37
    assert snps.build_detected


def test_sex_no_snps(snps_none):
    assert snps_none.sex == ''


def test_sex_Male_Y_chrom(l):
    ind = simulate_snps(l.create_individual('test_snps_sex_Male_Y_chrom'), chrom='Y', pos_start=1,
                        pos_max=59373566, pos_step=10000)
    ind.save_snps()
    from lineage.snps import SNPs
    snps = SNPs('output/test_snps_sex_Male_Y_chrom.csv')
    assert snps.sex == 'Male'


def test__read_raw_data(snps_none):
    assert snps_none.snps is None
    assert snps_none.source == ''


def test__lookup_build_with_snp_pos_None(snps_discrepant_pos):
    from lineage.snps import detect_build
    assert detect_build(snps_discrepant_pos) is None


def test_get_assembly_None():
    from lineage.snps import get_assembly
    assert get_assembly(None) is ''
