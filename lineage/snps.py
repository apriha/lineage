""" Utilities for reading and analyzing SNPs within the `lineage` framework. """

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

import gzip
import os
import zipfile

import numpy as np
import pandas as pd


class SNPs(object):

    def __init__(self, file):
        self.snps = self._read_raw_data(file)

    def _read_raw_data(self, file):
        if not os.path.exists(file):
            print(file + ' does not exist; skipping')
            return None

        # peek into files to determine the data format
        if '.zip' in file:
            with zipfile.ZipFile(file) as z:
                with z.open(z.namelist()[0], 'r') as f:
                    # https://stackoverflow.com/a/606199
                    line = f.readline().decode('utf-8')
        elif '.gz' in file:
            with gzip.open(file, 'rt') as f:
                line = f.readline()
        else:
            with open(file, 'r') as f:
                line = f.readline()

        if '23andMe' in line:
            return self._read_23andme(file)
        elif 'Ancestry' in line:
            return self._read_ancestry(file)
        elif line[:4] == 'RSID':
            return self._read_ftdna(file)
        elif line[:4] == 'rsid':
            return self._read_generic_csv(file)
        else:
            return None

    @staticmethod
    def _read_23andme(file):
        """ Read and parse 23andMe file.

        https://www.23andme.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            individual's genetic data normalized for use with `lineage`
        """
        try:
            return pd.read_csv(file, comment='#', sep='\t', na_values='--',
                               names=['rsid', 'chrom', 'pos', 'genotype'],
                               index_col=0, dtype={'chrom': object})
        except Exception as err:
            print(err)
            return None

    @staticmethod
    def _read_ftdna(file):
        """ Read and parse Family Tree DNA (FTDNA) file.

        https://www.familytreedna.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            individual's genetic data normalized for use with `lineage`
        """
        try:
            df = pd.read_csv(file, skiprows=1, na_values='--',
                             names=['rsid', 'chrom', 'pos', 'genotype'],
                             index_col=0, dtype={'chrom': object})

            # remove incongruous data
            df = df.drop(df.loc[df['chrom'] == '0'].index)
            df = df.drop(df.loc[df.index == 'RSID'].index)  # second header for concatenated data

            # if second header existed, pos dtype will be object (should be np.int64)
            df['pos'] = df['pos'].astype(np.int64)

            return df
        except Exception as err:
            print(err)
            return None

    @staticmethod
    def _read_ancestry(file):
        """ Read and parse Ancestry.com file.

        http://www.ancestry.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            individual's genetic data normalized for use with `lineage`
        """
        try:
            df = pd.read_csv(file, comment='#', header=0, sep='\t', na_values=0,
                             names=['rsid', 'chrom', 'pos', 'allele1', 'allele2'],
                             index_col=0, dtype={'chrom': object})

            # create genotype column from allele columns
            df['genotype'] = df['allele1'] + df['allele2']

            # delete allele columns
            # http://stackoverflow.com/a/13485766
            del df['allele1']
            del df['allele2']

            # https://redd.it/5y90un
            df.ix[np.where(df['chrom'] == '23')[0], 'chrom'] = 'X'
            df.ix[np.where(df['chrom'] == '24')[0], 'chrom'] = 'Y'
            df.ix[np.where(df['chrom'] == '25')[0], 'chrom'] = 'PAR'
            df.ix[np.where(df['chrom'] == '26')[0], 'chrom'] = 'MT'

            return df
        except Exception as err:
            print(err)
            return None

    @staticmethod
    def _read_generic_csv(file):
        """ Read and parse generic CSV file.

        Notes
        -----
        Assumes columns are 'rsid', 'chrom' / 'chromosome', 'pos' / 'position', and 'genotype';
        values are comma separated; unreported genotypes are indicated by '--'; and one header row
        precedes data. For example:

            rsid,chromosome,position,genotype
            rs1,1,1,AA
            rs2,1,2,CC
            rs3,1,3,--

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            individual's genetic data normalized for use with `lineage`
        """
        try:
            df = pd.read_csv(file, skiprows=1, na_values='--',
                             names=['rsid', 'chrom', 'pos', 'genotype'],
                             index_col=0, dtype={'chrom': object, 'pos': np.int64})

            return df
        except Exception as err:
            print(err)
            return None
