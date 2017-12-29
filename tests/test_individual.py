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

import numpy as np
import pandas as pd
import pytest

@pytest.fixture(scope='module')
def normalized_snps():
    rsid = ['rs1', 'rs2', 'rs3', 'rs4', 'rs5']
    chrom = ['1', '1', '1', '1', '1']
    pos = [1, 2, 3, 4, 5]
    genotype = ['AA', 'CC', 'GG', 'TT', np.nan]

    df = pd.DataFrame({'rsid': rsid, 'chrom': chrom, 'pos': pos, 'genotype': genotype},
                      columns=['rsid', 'chrom', 'pos', 'genotype'])

    df = df.set_index('rsid')
    return df

def test__read_23andme(l, normalized_snps):
    # https://www.23andme.com
    test_user = l.create_individual('test user', 'tests/input/23andme.txt')
    pd.testing.assert_frame_equal(test_user.snps, normalized_snps)

def test__read_ftdna(l, normalized_snps):
    # https://www.familytreedna.com
    test_user = l.create_individual('test user', 'tests/input/ftdna.csv')
    pd.testing.assert_frame_equal(test_user.snps, normalized_snps)

def test__read_ancestry(l, normalized_snps):
    # http://www.ancestry.com
    test_user = l.create_individual('test user', 'tests/input/ancestry.txt')
    pd.testing.assert_frame_equal(test_user.snps, normalized_snps)
