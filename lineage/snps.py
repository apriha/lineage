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

from itertools import groupby, count
import gzip
import os
import zipfile

import numpy as np
import pandas as pd


class SNPs(object):

    def __init__(self, file):
        self.snps = self._read_raw_data(file)
        self.assembly = None
        self.assembly_detected = False

        if self.snps is not None:
            self.assembly = detect_assembly(self.snps)

            if self.assembly is None:
                self.assembly = 37  # assume GRCh37 if not detected
            else:
                self.assembly_detected = True

    @property
    def assembly_name(self):
        """ Get the name of the assembly of ``SNPs``.

        Returns
        -------
        str
        """
        return get_assembly_name(self.assembly)

    @property
    def snp_count(self):
        """ Count of SNPs.

        Returns
        -------
        int
        """
        return get_snp_count(self.snps)

    @property
    def chromosomes(self):
        """ Chromosomes of ``SNPs``.

        Returns
        -------
        list
            list of str chromosomes (e.g., ['1', '2', '3', 'MT'], empty list if no chromosomes
        """
        return get_chromosomes(self.snps)

    @property
    def chromosomes_summary(self):
        """ Summary of the chromosomes of ``SNPs``.

        Returns
        -------
        str
            human-readable listing of chromosomes (e.g., '1-3, MT'), empty str if no chromosomes
        """
        return get_chromosomes_summary(self.snps)

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


def detect_assembly(snps):
    """ Detect assembly of SNPs.

    Use the coordinates of common SNPs to identify the assembly / build of a genotype file
    that is being loaded.

    Notes
    -----
    rs3094315 : plus strand in 36, 37, and 38
    rs11928389 : plus strand in 36, minus strand in 37 and 38
    rs2500347 : plus strand in 36 and 37, minus strand in 38
    rs964481 : plus strand in 36, 37, and 38

    Parameters
    ----------
    snps : pandas.DataFrame
        SNPs to add

    Returns
    -------
    int
        detected assembly of SNPs, else None

    References
    ----------
    ..[1] Yates et. al. (doi:10.1093/bioinformatics/btu613),
      http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613
    ..[2] Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098
    ..[3] Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K.
      dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;29(1):308-11.
    ..[4] Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center
      for Biotechnology Information, National Library of Medicine. dbSNP accession: rs3094315,
      rs11928389, rs2500347, and rs964481 (dbSNP Build ID: 151). Available from:
      http://www.ncbi.nlm.nih.gov/SNP/
    """
    def lookup_assembly_with_snp_pos(pos, s):
        try:
            return s.loc[s == pos].index[0]
        except:
            return None

    assembly = None

    rsids = ['rs3094315', 'rs11928389', 'rs2500347', 'rs964481']
    df = pd.DataFrame({36: [742429, 50908372, 143649677, 27566744],
                       37: [752566, 50927009, 144938320, 27656823],
                       38: [817186, 50889578, 148946169, 27638706]}, index=rsids)

    for rsid in rsids:
        if rsid in snps.index:
            assembly = lookup_assembly_with_snp_pos(snps.loc[rsid].pos, df.loc[rsid])

        if assembly is not None:
            break

    return assembly


def get_assembly_name(assembly):
    """ Get the name of an assembly.

    Parameters
    ----------
    assembly : int {36, 37, 38}

    Returns
    -------
    str
        empty str if `assembly` is None
    """

    if assembly is None:
        return ''
    elif assembly == 36:
        return 'NCBI36'
    elif assembly == 38:
        return 'GRCh38'
    else:
        return 'GRCh37'


def get_snp_count(snps):
    """ Count of SNPs.

    Parameters
    ----------
    snps : pandas.DataFrame

    Returns
    -------
    int
    """

    if snps is not None:
        return len(snps)
    else:
        return 0


def get_chromosomes(snps):
    """ Get the chromosomes of SNPs.

    Parameters
    ----------
    snps : pandas.DataFrame

    Returns
    -------
    list
        list of str chromosomes (e.g., ['1', '2', '3', 'MT'], empty list if no chromosomes
    """

    if snps is not None:
        return list(pd.unique(snps['chrom']))
    else:
        return []


def get_chromosomes_summary(snps):
    """ Summary of the chromosomes of SNPs.

    Parameters
    ----------
    snps : pandas.DataFrame

    Returns
    -------
    str
        human-readable listing of chromosomes (e.g., '1-3, MT'), empty str if no chromosomes
    """

    if snps is not None:
        chroms = list(pd.unique(snps['chrom']))

        int_chroms = [int(chrom) for chrom in chroms if chrom.isdigit()]
        str_chroms = [chrom for chrom in chroms if not chrom.isdigit()]

        # https://codereview.stackexchange.com/a/5202
        def as_range(iterable):
            l = list(iterable)
            if len(l) > 1:
                return '{0}-{1}'.format(l[0], l[-1])
            else:
                return '{0}'.format(l[0])

        # create str representations
        int_chroms = ', '.join(as_range(g) for _, g in
                               groupby(int_chroms, key=lambda n, c=count(): n - next(c)))
        str_chroms = ', '.join(str_chroms)

        if str_chroms != '':
            int_chroms += ', '

        return int_chroms + str_chroms
    else:
        return ''
