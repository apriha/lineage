""" Class for representing individuals within the `lineage` framework. """

"""
Copyright (C) 2016 Andrew Riha

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
import re
import zipfile

import numpy as np
import pandas as pd

import lineage

class Individual(object):
    """ Object used to represent and interact with an individual.

    The `Individual` object maintains information about an individual. The object provides
    methods for loading an individual's genetic data (SNPs) and normalizing it for use with the
    `lineage` framework.

    """

    def __init__(self, name, raw_data=None, output_dir='', ensembl_rest_client=None):
        """ Initialize an individual.

        Parameters
        ----------
        name : str
            name of the individual
        raw_data : list or str
            path(s) to file(s) with raw genotype data
        output_dir : str
            path to output directory
        ensembl_rest_client : EnsemblRestClient
            client for making requests to the Ensembl REST API

        """

        self._name = name
        self._output_dir = output_dir
        self._ensembl_rest_client = ensembl_rest_client
        self._snps = None
        self._discrepant_positions_file_count = 0
        self._discrepant_genotypes_file_count = 0

        if raw_data is not None:
            self.load_snps(raw_data)

    def __repr__(self):
        return 'Individual({!r})'.format(self._name)

    @property
    def name(self):
        return self._name

    @property
    def snps(self):
        if self._snps is not None:
            return self._snps.copy()
        else:
            return None

    def load_snps(self, raw_data, discrepant_snp_positions_threshold=100,
                  discrepant_genotypes_threshold=10000):
        """ Load raw genotype data.

        Parameters
        ----------
        raw_data : list or str
            path(s) to file(s) with raw genotype data
        discrepant_snp_positions_threshold : int
            threshold for discrepant SNP positions between existing data and data to be loaded,
            a large value could indicate mismatched genome assemblies
        discrepant_genotypes_threshold : int
            threshold for discrepant genotype data between existing data and data to be loaded,
            a large value could indicated mismatched individuals

        """
        if type(raw_data) is list:
            for file in raw_data:
                print('Loading ' + os.path.relpath(file))
                self._add_snps(self._read_raw_data(file), discrepant_snp_positions_threshold,
                               discrepant_genotypes_threshold)
        elif type(raw_data) is str:
            print('Loading ' + os.path.relpath(raw_data))
            self._add_snps(self._read_raw_data(raw_data), discrepant_snp_positions_threshold,
                           discrepant_genotypes_threshold)
        else:
            raise TypeError('invalid filetype')

    def save_snps(self):
        """ Save SNPs to file. """
        if self._snps is not None:
            try:
                if lineage.dir_exists(self._output_dir):
                    output_dir = self._output_dir
                else:
                    output_dir = ''

                file = os.path.join(output_dir, self.get_var_name() + '.csv')
                print('Saving ' + os.path.relpath(file))
                self._snps.to_csv(file, na_rep="--", header=['chromosome', 'position', 'genotype'])
            except Exception as err:
                print(err)
        else:
            print('no SNPs to save...')

    def get_var_name(self):
        """ Clean a string so that it can be a valid Python variable
        name.

        Returns
        -------
        str
            cleaned string that can be used as a variable name

        """

        return re.sub('\W|^(?=\d)', '_', self.name)  # [SO-02]

    def remap_snps(self, source_assembly, target_assembly, complement_bases=True):
        """ Remap the SNP coordinates of this ``Individual`` from one assembly to another.

        This method uses the assembly map endpoint of the Ensembl REST API service (via this
        ``Individual``'s ``EnsemblRestClient``) to convert SNP coordinates / positions from one
        assembly to another. After remapping, the coordinates / positions for the ``Individual``'s
        SNPs will be that of the target assembly.

        Parameters
        ----------
        source_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            starting assembly of an Individual's SNPs
        target_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap to
        complement_bases : bool
            complement bases when remapping SNPs to the minus strand

        Notes
        -----
        An assembly is also know as a "build." For example:

        Assembly NCBI36 = Build 36
        Assembly GRCh37 = Build 37
        Assembly GRCh38 = Build 38

        See https://www.ncbi.nlm.nih.gov/assembly for more information about assemblies and
        remapping.

        References
        ----------
        ..[1] Ensembl, Assembly Map Endpoint,
          http://rest.ensembl.org/documentation/info/assembly_map

        """

        if self._ensembl_rest_client is None:
            print('Need an ``EnsemblRestClient`` to remap SNPs')
            return

        valid_assemblies = ['NCBI36', 'GRCh37', 'GRCh38']

        if source_assembly not in valid_assemblies:
            print('Invalid source assembly')
            return
        elif target_assembly not in valid_assemblies:
            print('Invalid target assembly')
            return

        for chrom in self._snps['chrom'].unique():
            print('Remapping chromosome ' + chrom + '...')

            # extract SNPs for this chrom for faster remapping
            temp = pd.DataFrame(self._snps.loc[self._snps['chrom'] == chrom])

            temp['remapped'] = False

            pos_start = str(int(temp['pos'].describe()['min']))
            pos_end = str(int(temp['pos'].describe()['max']))

            endpoint = '/map/human/' + source_assembly + '/' + chrom + ':' + \
                       pos_start + '..' + pos_end + '/' + target_assembly + '?'

            # get remapping information
            response = self._ensembl_rest_client.perform_rest_action(endpoint)

            if response is None:
                continue

            for mapping in response['mappings']:
                orig_range_len = mapping['original']['end'] - mapping['original']['start']
                mapped_range_len = mapping['mapped']['end'] - mapping['mapped']['start']

                orig_region = mapping['original']['seq_region_name']
                mapped_region = mapping['mapped']['seq_region_name']

                if orig_region != mapped_region:
                    print('discrepant chroms')
                    continue

                if orig_range_len != mapped_range_len:
                    print('discrepant coords')  # observed when mapping NCBI36 -> GRCh38
                    continue

                # find the SNPs that are being remapped for this mapping
                snp_indices = temp.loc[~temp['remapped'] &
                                   (temp['pos'] >= mapping['original']['start']) &
                                   (temp['pos'] <= mapping['original']['end'])].index

                if len(snp_indices) > 0:
                    # remap the SNPs
                    if mapping['mapped']['strand'] == -1:
                        # flip and (optionally) complement since we're mapping to minus strand
                        diff_from_start = temp.loc[snp_indices, 'pos'] - mapping['original']['start']
                        temp.loc[snp_indices, 'pos'] = mapping['mapped']['end'] - diff_from_start

                        if complement_bases:
                            self._snps.loc[snp_indices, 'genotype'] = \
                                temp.loc[snp_indices, 'genotype'].apply(self._complement_bases)
                    else:
                        # mapping is on same (plus) strand, so just remap based on offset
                        offset = mapping['mapped']['start'] - mapping['original']['start']
                        temp.loc[snp_indices, 'pos'] = temp['pos'] + offset

                    # mark these SNPs as remapped
                    temp.loc[snp_indices, 'remapped'] = True

            # update SNP positions for this chrom
            self._snps.loc[temp.index, 'pos'] = temp['pos']

        self._sort_snps()

    def _complement_bases(self, genotype):
        if pd.isnull(genotype):
            return np.nan

        complement = ''

        for base in list(genotype):
            if base == 'A':
                complement += 'T'
            elif base == 'G':
                complement += 'C'
            elif base == 'C':
                complement += 'G'
            elif base == 'T':
                complement += 'A'

        return complement

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

            # delete allele columns [SO-01]
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

    def _add_snps(self, snps, discrepant_snp_positions_threshold, discrepant_genotypes_threshold):
        """ Add SNPs to this Individual.

        Parameters
        ----------
        snps : pandas.DataFrame
            SNPs to add
        discrepant_snp_positions_threshold : int
            see above
        discrepant_genotypes_threshold : int
            see above

        """

        if snps is None:
            return

        # ensure there area always two X alleles
        snps = self._double_single_alleles(snps, 'X')

        if self._snps is None:
            self._snps = snps
        else:
            common_snps = self._snps.join(snps, how='inner', rsuffix='_added')

            discrepant_positions = common_snps.loc[
                (common_snps['chrom'] != common_snps['chrom_added']) |
                (common_snps['pos'] != common_snps['pos_added'])]

            if 0 < len(discrepant_positions) < discrepant_snp_positions_threshold:
                print(str(len(discrepant_positions)) + ' SNP positions being added differ; '
                      'keeping original positions')

                self._discrepant_positions_file_count += 1
                lineage.save_df_as_csv(discrepant_positions, self._output_dir,
                                       self.get_var_name() + '_discrepant_positions_' + str(
                                           self._discrepant_positions_file_count) + '.csv')
            elif len(discrepant_positions) >= discrepant_snp_positions_threshold:
                print('too many SNPs differ in position; ensure same genome build is being used')
                return

            # remove null genotypes
            common_snps = common_snps.loc[~common_snps['genotype'].isnull() &
                                          ~common_snps['genotype_added'].isnull()]

            # discrepant genotypes are where alleles are not equivalent (i.e., alleles are not the
            # same and not swapped)
            discrepant_genotypes = common_snps.loc[
                ((common_snps['genotype'].str.len() == 1) &
                 (common_snps['genotype_added'].str.len() == 1) &
                 ~(common_snps['genotype'].str[0] == common_snps['genotype_added'].str[0])) |
                ((common_snps['genotype'].str.len() == 2) &
                 (common_snps['genotype_added'].str.len() == 2) &
                 ~((common_snps['genotype'].str[0] == common_snps['genotype_added'].str[0]) &
                   (common_snps['genotype'].str[1] == common_snps['genotype_added'].str[1])) &
                 ~((common_snps['genotype'].str[0] == common_snps['genotype_added'].str[1]) &
                   (common_snps['genotype'].str[1] == common_snps['genotype_added'].str[0])))]

            if 0 < len(discrepant_genotypes) < discrepant_genotypes_threshold:
                print(str(len(discrepant_genotypes)) + ' genotypes were discrepant; '
                      'marking those as null')

                self._discrepant_genotypes_file_count += 1
                lineage.save_df_as_csv(discrepant_genotypes, self._output_dir,
                                       self.get_var_name() + '_discrepant_genotypes_' + str(
                                           self._discrepant_genotypes_file_count) + '.csv')
            elif len(discrepant_genotypes) >= discrepant_genotypes_threshold:
                print('too many SNPs differ in their genotype; ensure file is for same '
                      'individual')
                return

            # add new SNPs
            self._snps = self._snps.combine_first(snps)
            self._snps.loc[discrepant_genotypes.index, 'genotype'] = np.nan

            # combine_first converts position to float64, so convert it back to int64
            self._snps['pos'] = self._snps['pos'].astype(np.int64)

        self._sort_snps()

    @staticmethod
    def _double_single_alleles(df, chrom):
        """ Double any single alleles in the specified chromosome.

        Parameters
        ----------
        df : pandas.DataFrame
            SNPs
        chrom : str
            chromosome of alleles to double

        Returns
        -------
        df : pandas.DataFrame
            SNPs with specified chromosome's single alleles doubled

        """
        # find all single alleles of the specified chromosome
        single_alleles = np.where((df['chrom'] == chrom) &
                                  (df['genotype'].str.len() == 1))[0]

        # double those alleles
        df.ix[single_alleles, 'genotype'] = df.ix[single_alleles, 'genotype'] * 2

        return df

    def _sort_snps(self):
        """ Sort this individual's SNPs. """

        sorted_list = sorted(self._snps['chrom'].unique(), key=self._natural_sort_key)

        # move PAR and MT to the end of the dataframe
        if 'PAR' in sorted_list:
            sorted_list.remove('PAR')
            sorted_list.append('PAR')

        if 'MT' in sorted_list:
            sorted_list.remove('MT')
            sorted_list.append('MT')

        # https://stackoverflow.com/a/26707444
        # convert chrom column to category for sorting
        self._snps['chrom'] = \
            self._snps['chrom'].astype('category', categories=sorted_list, ordered=True)

        # sort based on ordered chromosome list and position
        self._snps = self._snps.sort_values(['chrom', 'pos'])

        # convert chromosome back to object
        self._snps['chrom'] = self._snps['chrom'].astype(object)

    # https://stackoverflow.com/a/16090640
    @staticmethod
    def _natural_sort_key(s, natural_sort_re=re.compile('([0-9]+)')):
        return [int(text) if text.isdigit() else text.lower()
                for text in re.split(natural_sort_re, s)]


"""
Stack Overflow Attributions
---------------------------

[SO-##] references throughout the code above correspond to the following
content from Stack Overflow (http://stackoverflow.com):

[SO-01] "Delete column from pandas DataFrame"
        http://stackoverflow.com/q/13411544
        John : http://stackoverflow.com/users/390388/john
        http://stackoverflow.com/a/13485766
        Wes McKinney : http://stackoverflow.com/users/776560/wes-mckinney

[SO-02] "how do I convert a string to a valid variable name in python?"
        http://stackoverflow.com/q/3303312
        George Profenza : http://stackoverflow.com/users/89766/george-profenza
        http://stackoverflow.com/a/3305731
        Nas Banov : http://stackoverflow.com/users/226086/nas-banov

"""
