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

from lineage.snps import sort_snps
from tests.test_lineage import simulate_snps


def create_snp_df(rsid, chrom, pos, genotype):
    df = pd.DataFrame({'rsid': rsid, 'chrom': chrom, 'pos': pos, 'genotype': genotype},
                      columns=['rsid', 'chrom', 'pos', 'genotype'])
    df = df.set_index('rsid')
    return df


@pytest.fixture(scope='module')
def generic_snps():
    return create_snp_df(rsid=['rs1', 'rs2', 'rs3', 'rs4', 'rs5'], chrom=['1', '1', '1', '1', '1'],
                pos=[1, 2, 3, 4, 5], genotype=['AA', 'CC', 'GG', 'TT', np.nan])


@pytest.fixture(scope='module')
def snps_NCBI36():
    return create_snp_df(rsid=['rs3094315', 'rs2500347', 'rsIndelTest', 'rs11928389'],
                         chrom=['1', '1', '1', '3'], pos=[742429, 143649677, 143649678, 50908372],
                         genotype=['AA', np.nan, 'ID', 'AG'])


@pytest.fixture(scope='module')
def snps_NCBI36_discrepant_snps():
    return create_snp_df(rsid=['rs3094315', 'rs2500347', 'rsIndelTest', 'rs11928389'],
                         chrom=['1', '1', '1', '3'], pos=[742429, 143649677, 143649678, 50908372],
                         genotype=['AA', np.nan, 'ID', np.nan])


@pytest.fixture(scope='module')
def snps_GRCh37():
    return create_snp_df(rsid=['rs3094315', 'rs2500347', 'rsIndelTest', 'rs11928389'],
                         chrom=['1', '1', '1', '3'], pos=[752566, 144938320, 144938321, 50927009],
                         genotype=['AA', np.nan, 'ID', 'TC'])


@pytest.fixture(scope='module')
def snps_GRCh38():
    return create_snp_df(rsid=['rs3094315', 'rsIndelTest', 'rs2500347', 'rs11928389'],
                         chrom=['1', '1', '1', '3'], pos=[817186, 148946168, 148946169, 50889578],
                         genotype=['AA', 'ID', np.nan, 'TC'])


@pytest.fixture(scope='module')
def snps_GRCh38_PAR():
    return create_snp_df(rsid=['rs28736870', 'rs113313554'], chrom=['X', 'Y'],
                         pos=[304103, 624523], genotype=['AA', 'AA'])


def test_name(l):
    ind = l.create_individual('test')
    assert ind.name == 'test'


def test_snps_23andme(l, generic_snps):
    # https://www.23andme.com
    ind = l.create_individual('', 'tests/input/23andme.txt')
    assert ind.source == '23andMe'
    pd.testing.assert_frame_equal(ind.snps, generic_snps)


def test_snps_23andme_zip(l, generic_snps):
    with zipfile.ZipFile('tests/input/23andme.txt.zip', 'w') as f:
        # https://stackoverflow.com/a/16104667
        f.write('tests/input/23andme.txt', arcname='23andme.txt')
    ind = l.create_individual('', 'tests/input/23andme.txt.zip')
    assert ind.source == '23andMe'
    pd.testing.assert_frame_equal(ind.snps, generic_snps)


def test_snps_ftdna(l, generic_snps):
    # https://www.familytreedna.com
    ind = l.create_individual('', 'tests/input/ftdna.csv')
    assert ind.source == 'FTDNA'
    pd.testing.assert_frame_equal(ind.snps, generic_snps)


def test_snps_ftdna_gzip(l, generic_snps):
    with open('tests/input/ftdna.csv', 'rb') as f_in:
        with gzip.open('tests/input/ftdna.csv.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    ind = l.create_individual('', 'tests/input/ftdna.csv.gz')
    assert ind.source == 'FTDNA'
    pd.testing.assert_frame_equal(ind.snps, generic_snps)


def test_snps_ancestry(l, generic_snps):
    # https://www.ancestry.com
    ind = l.create_individual('', 'tests/input/ancestry.txt')
    assert ind.source == 'AncestryDNA'
    pd.testing.assert_frame_equal(ind.snps, generic_snps)


def test_source_lineage(l):
    ind = l.create_individual('', 'tests/input/chromosomes.csv')
    assert ind.source == 'generic'
    file = ind.save_snps()
    ind_saved_snps = l.create_individual('', file)
    assert ind_saved_snps.source == 'lineage'
    pd.testing.assert_frame_equal(ind.snps, ind_saved_snps.snps)


def test_source_generic(l):
    ind = l.create_individual('', 'tests/input/NCBI36.csv')
    assert ind.source == 'generic'


def test_snps_None(l):
    ind = l.create_individual('')
    assert ind.snps is None


def test_snp_count(l):
    ind = l.create_individual('', 'tests/input/NCBI36.csv')
    assert ind.snp_count == 4


def test_snp_count_None(l):
    ind = l.create_individual('')
    assert ind.snp_count == 0


def test_chromosomes(l):
    ind = l.create_individual('', 'tests/input/chromosomes.csv')
    assert ind.chromosomes == ['1', '2', '3', '5', 'PAR', 'MT']


def test_chromosomes_None(l):
    ind = l.create_individual('')
    assert ind.chromosomes == []


def test_chromosomes_summary(l):
    ind = l.create_individual('', 'tests/input/chromosomes.csv')
    assert ind.chromosomes_summary == '1-3, 5, PAR, MT'


def test_chromosomes_summary_None(l):
    ind = l.create_individual('')
    assert ind.chromosomes_summary == ''


def test_build(l):
    ind = l.create_individual('', 'tests/input/NCBI36.csv')
    assert ind.build == 36
    assert ind.assembly_name == 'NCBI36'


def test_load_snps_list(l, snps_GRCh37):
    ind = l.create_individual('')
    ind.load_snps(['tests/input/GRCh37.csv', 'tests/input/GRCh37.csv'])
    pd.testing.assert_frame_equal(ind.snps, snps_GRCh37)
    assert ind.source == 'generic, generic'


def test_load_snps_None(l):
    ind = l.create_individual('')
    with pytest.raises(TypeError):
        ind.load_snps(None)


def test_sex_Male_Y_chrom(l):
    ind = simulate_snps(l.create_individual(''), chrom='Y', pos_start=1, pos_max=59373566,
                        pos_step=10000)
    assert ind.sex == 'Male'


def test_sex_Female_Y_chrom(l):
    ind = simulate_snps(l.create_individual(''), chrom='Y', pos_start=1, pos_max=59373566,
                        pos_step=10000, null_snp_step=1)
    assert ind.sex == 'Female'


def test_sex_Female_X_chrom(l):
    ind = simulate_snps(l.create_individual(''), chrom='X', pos_start=1, pos_max=155270560,
                        pos_step=10000, genotype='AC')
    assert ind.sex == 'Female'


def test_sex_Male_X_chrom(l):
    ind = simulate_snps(l.create_individual(''), chrom='X', pos_start=1, pos_max=155270560,
                        pos_step=10000, genotype='AA')
    assert ind.sex == 'Male'


def test_sex_not_determined(l):
    ind = simulate_snps(l.create_individual(''), chrom='1', pos_start=1, pos_max=249250621,
                        pos_step=10000)
    assert ind.sex == ''


def test_load_snps_non_existent_file(l, snps_GRCh37):
    ind = l.create_individual('')
    ind.load_snps(['tests/input/GRCh37.csv', 'tests/input/non_existent_file.csv'])
    pd.testing.assert_frame_equal(ind.snps, snps_GRCh37)


def test_load_snps_invalid_file(l, snps_GRCh37):
    ind = l.create_individual('')
    with open('tests/input/empty.txt', 'w'):
        pass
    ind.load_snps(['tests/input/GRCh37.csv', 'tests/input/empty.txt'])
    pd.testing.assert_frame_equal(ind.snps, snps_GRCh37)


def test_load_snps_assembly_mismatch(l, snps_NCBI36_discrepant_snps):
    ind = l.create_individual('ind')
    ind.load_snps(['tests/input/NCBI36.csv', 'tests/input/GRCh37.csv'])
    assert not os.path.exists('output/ind_discrepant_positions_1.csv')
    assert not os.path.exists('output/ind_discrepant_genotypes_1.csv')
    assert len(ind.discrepant_positions) == 4
    assert len(ind.discrepant_genotypes) == 1
    pd.testing.assert_frame_equal(ind.snps, snps_NCBI36_discrepant_snps)


def test_load_snps_assembly_mismatch_save_output(l, snps_NCBI36_discrepant_snps):
    ind = l.create_individual('ind')
    ind.load_snps(['tests/input/NCBI36.csv', 'tests/input/GRCh37.csv'], save_output=True)
    assert os.path.exists('output/ind_discrepant_positions_1.csv')
    assert os.path.exists('output/ind_discrepant_genotypes_1.csv')
    assert len(ind.discrepant_positions) == 4
    assert len(ind.discrepant_genotypes) == 1
    pd.testing.assert_frame_equal(ind.snps, snps_NCBI36_discrepant_snps)


def test_load_snps_assembly_mismatch_exceed_discrepant_positions_threshold(l, snps_NCBI36):
    ind = l.create_individual('ind')
    ind.load_snps(['tests/input/NCBI36.csv', 'tests/input/GRCh37.csv'],
                  discrepant_snp_positions_threshold=0)
    assert not os.path.exists('output/ind_discrepant_positions_1.csv')
    assert not os.path.exists('output/ind_discrepant_genotypes_1.csv')
    assert len(ind.discrepant_positions) == 4
    assert len(ind.discrepant_genotypes) == 0
    pd.testing.assert_frame_equal(ind.snps, snps_NCBI36)


def test_load_snps_assembly_mismatch_exceed_discrepant_genotypes_threshold(l, snps_NCBI36):
    ind = l.create_individual('ind')
    ind.load_snps(['tests/input/NCBI36.csv', 'tests/input/GRCh37.csv'],
                  discrepant_genotypes_threshold=0)
    assert not os.path.exists('output/ind_discrepant_positions_1.csv')
    assert not os.path.exists('output/ind_discrepant_genotypes_1.csv')
    assert len(ind.discrepant_positions) == 4
    assert len(ind.discrepant_genotypes) == 1
    pd.testing.assert_frame_equal(ind.snps, snps_NCBI36)


def test_merging_files_discrepant_snps(l):
    df = pd.read_csv('tests/input/discrepant_snps.csv', skiprows=1, na_values='--',
                     names=['rsid', 'chrom', 'pos_file1', 'pos_file2',
                            'genotype_file1', 'genotype_file2', 'discrepant_position',
                            'discrepant_genotype', 'expected_position', 'expected_genotype'],
                     index_col=0,
                     dtype={'chrom': object, 'pos_file1': np.int64, 'pos_file2': np.int64,
                            'discrepant_position': bool, 'discrepant_genotype': bool})

    df1 = df[['chrom', 'pos_file1', 'genotype_file1']]
    df2 = df[['chrom', 'pos_file2', 'genotype_file2']]

    df1.to_csv('tests/input/discrepant_snps1.csv', na_rep='--',
               header=['chromosome', 'position', 'genotype'])

    df2.to_csv('tests/input/discrepant_snps2.csv', na_rep='--',
               header=['chromosome', 'position', 'genotype'])

    ind = l.create_individual('', ['tests/input/discrepant_snps1.csv',
                                   'tests/input/discrepant_snps2.csv'])

    expected = df[['chrom', 'discrepant_position', 'discrepant_genotype', 'expected_position',
              'expected_genotype']]
    expected = expected.rename(columns={'expected_position': 'pos',
                                        'expected_genotype': 'genotype'})
    expected = sort_snps(expected)

    pd.testing.assert_index_equal(ind.discrepant_positions.index,
                                  expected.loc[expected['discrepant_position'] == True].index)

    pd.testing.assert_index_equal(ind.discrepant_genotypes.index,
                                  expected.loc[expected['discrepant_genotype'] == True].index)

    pd.testing.assert_series_equal(ind.snps['pos'], expected['pos'])
    pd.testing.assert_series_equal(ind.snps['genotype'], expected['genotype'])

def test_save_snps(l, snps_GRCh37):
    ind = l.create_individual('test save snps', 'tests/input/GRCh37.csv')
    assert os.path.relpath(ind.save_snps()) == 'output/test_save_snps.csv'
    ind_saved_snps = l.create_individual('', 'output/test_save_snps.csv')
    pd.testing.assert_frame_equal(ind_saved_snps.snps, snps_GRCh37)


def test_save_snps_no_snps(l):
    ind = l.create_individual('')
    assert not ind.save_snps()


def test_save_snps_invalid_output_dir(l):
    ind = l.create_individual('', 'tests/input/GRCh37.csv')
    ind._output_dir = None
    assert not ind.save_snps()


def test_save_snps_exception(l):
    ind = l.create_individual('')
    ind._snps = 'invalid'
    assert not ind.save_snps()


def test_get_var_name(l):
    ind = l.create_individual('test a. name')
    assert ind.get_var_name() == 'test_a__name'


def test_remap_snps_36_to_37(l, snps_GRCh37):
    ind = l.create_individual('', 'tests/input/NCBI36.csv')
    chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(37)
    assert ind.build == 37
    assert ind.assembly_name == 'GRCh37'
    if len(chromosomes_remapped) == 2:
        assert len(chromosomes_not_remapped) == 0
        pd.testing.assert_frame_equal(ind.snps, snps_GRCh37)


def test_remap_snps_37_to_36(l, snps_NCBI36):
    ind = l.create_individual('', 'tests/input/GRCh37.csv')
    chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(36)
    assert ind.build == 36
    assert ind.assembly_name == 'NCBI36'
    if len(chromosomes_remapped) == 2:
        assert len(chromosomes_not_remapped) == 0
        pd.testing.assert_frame_equal(ind.snps, snps_NCBI36)


def test_remap_snps_37_to_38(l, snps_GRCh38):
    ind = l.create_individual('', 'tests/input/GRCh37.csv')
    chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(38)
    assert ind.build == 38
    assert ind.assembly_name == 'GRCh38'
    if len(chromosomes_remapped) == 2:
        assert len(chromosomes_not_remapped) == 0
        pd.testing.assert_frame_equal(ind.snps, snps_GRCh38)


def test_remap_snps_37_to_38_via_lineage(l, snps_GRCh38):
    ind = l.create_individual('', 'tests/input/GRCh37.csv')
    chromosomes_remapped, chromosomes_not_remapped = l.remap_snps(ind, 38)
    assert ind.build == 38
    assert ind.assembly_name == 'GRCh38'
    if len(chromosomes_remapped) == 2:
        assert len(chromosomes_not_remapped) == 0
        pd.testing.assert_frame_equal(ind.snps, snps_GRCh38)


def test_remap_snps_37_to_38_with_PAR_SNP(l, snps_GRCh38_PAR):
    ind = l.create_individual('', 'tests/input/GRCh37_PAR.csv')
    assert ind.snp_count == 3
    chromosomes_remapped, chromosomes_not_remapped = l.remap_snps(ind, 38)
    assert ind.build == 38
    assert ind.assembly_name == 'GRCh38'
    if len(chromosomes_remapped) == 2:
        assert len(chromosomes_not_remapped) == 1
        assert ind.snp_count == 2
        pd.testing.assert_frame_equal(ind.snps, snps_GRCh38_PAR)


def test_remap_snps_37_to_37(l, snps_GRCh37):
    ind = l.create_individual('', 'tests/input/GRCh37.csv')
    chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(37)
    assert ind.build == 37
    assert ind.assembly_name == 'GRCh37'
    assert len(chromosomes_remapped) == 0
    assert len(chromosomes_not_remapped) == 2
    pd.testing.assert_frame_equal(ind.snps, snps_GRCh37)


def test_remap_snps_no_snps(l):
    ind = l.create_individual('')
    chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(38)
    assert ind.build is None
    assert len(chromosomes_remapped) == 0
    assert len(chromosomes_not_remapped) == 0


def test_remap_snps_invalid_assembly(l):
    ind = l.create_individual('', 'tests/input/GRCh37.csv')
    chromosomes_remapped, chromosomes_not_remapped = ind.remap_snps(-1)
    assert ind.build == 37
    assert ind.assembly_name == 'GRCh37'
    assert len(chromosomes_remapped) == 0
    assert len(chromosomes_not_remapped) == 2


def test___repr__(l):
    ind = l.create_individual('test')
    assert "Individual('test')" == ind.__repr__()
