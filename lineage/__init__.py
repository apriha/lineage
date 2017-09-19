""" lineage

tools for genetic genealogy and the analysis of consumer DNA test results

"""

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

import os
import numpy as np
import pandas as pd

# http://mikegrouchy.com/blog/2012/05/be-pythonic-__init__py.html
from lineage.ensembl import EnsemblRestClient
from lineage.individual import Individual
from lineage.resources import Resources
from lineage.visualization import plot_chromosomes


class Lineage(object):
    def __init__(self, output_dir='output', resources_dir='resources'):
        self._output_dir = os.path.abspath(output_dir)
        self._resources = Resources(resources_dir=resources_dir)
        self._ensembl_rest_client = EnsemblRestClient()

    def create_individual(self, name, raw_data=None):
        """ Initialize an individual in the context of the `lineage` framework.

        Parameters
        ----------
        name : str
            name of the individual
        raw_data : list or str
            path(s) to file(s) with raw genotype data

        Returns
        -------
        Individual
            ``Individual`` initialized in the context of the `lineage` framework

        """
        return Individual(name, raw_data, self._output_dir, self._ensembl_rest_client)

    def find_discordant_snps(self, individual1, individual2, individual3=None, save_output=False):
        """ Find discordant SNPs between two or three individuals.

        Parameters
        ----------
        individual1 : Individual
            reference individual (child if `individual2` and `individual3` are parents)
        individual2 : Individual
            comparison individual
        individual3 : Individual
            other parent if `individual1` is child and `individual2` is a parent
        save_output : bool
            specifies whether to save output to a CSV file in the output directory

        Returns
        -------
        pandas.DataFrame
            discordant SNPs and associated genetic data

        References
        ----------
        ..[1] David Pike, "Search for Discordant SNPs in Parent-Child
          Raw Data Files," David Pike's Utilities,
          http://www.math.mun.ca/~dapike/FF23utils/pair-discord.php
        ..[2] David Pike, "Search for Discordant SNPs when given data
          for child and both parents," David Pike's Utilities,
          http://www.math.mun.ca/~dapike/FF23utils/trio-discord.php

        """

        df = individual1.snps

        # remove nulls for reference individual
        df = df.loc[df['genotype'].notnull()]

        # add SNPs shared with `individual2`
        df = df.join(individual2.snps['genotype'], rsuffix='2')

        genotype1 = 'genotype_' + individual1.get_var_name()
        genotype2 = 'genotype_' + individual2.get_var_name()

        if individual3 is None:
            df = df.rename(columns={'genotype': genotype1, 'genotype2': genotype2})

            # find discordant SNPs between reference and comparison individuals
            df = df.loc[df[genotype2].notnull() &
                        ((df[genotype1].str.len() == 1) &
                         (df[genotype2].str.len() == 1) &
                         (df[genotype1] != df[genotype2])) |
                        ((df[genotype1].str.len() == 2) &
                         (df[genotype2].str.len() == 2) &
                         (df[genotype1].str[0] != df[genotype2].str[1]) &
                         (df[genotype1].str[0] != df[genotype2].str[1]) &
                         (df[genotype1].str[1] != df[genotype2].str[0]) &
                         (df[genotype1].str[1] != df[genotype2].str[1]))]
            if save_output:
                save_df_as_csv(df, self._output_dir, 'discordant_snps_' +
                               individual1.get_var_name() + '_' +
                               individual2.get_var_name() + '.csv')
        else:
            # add SNPs shared with `individual3`
            df = df.join(individual3.snps['genotype'], rsuffix='3')

            genotype3 = 'genotype_' + individual3.get_var_name()

            df = df.rename(columns={'genotype': genotype1, 'genotype2': genotype2,
                                    'genotype3': genotype3})

            # find discordant SNPs between child and two parents
            df = df.loc[(df[genotype2].notnull() &
                         ((df[genotype1].str.len() == 1) &
                          (df[genotype2].str.len() == 1) &
                          (df[genotype1] != df[genotype2])) |
                         ((df[genotype1].str.len() == 2) &
                          (df[genotype2].str.len() == 2) &
                          (df[genotype1].str[0] != df[genotype2].str[0]) &
                          (df[genotype1].str[0] != df[genotype2].str[1]) &
                          (df[genotype1].str[1] != df[genotype2].str[0]) &
                          (df[genotype1].str[1] != df[genotype2].str[1]))) |
                        (df[genotype3].notnull() &
                         ((df[genotype1].str.len() == 1) &
                          (df[genotype3].str.len() == 1) &
                          (df[genotype1] != df[genotype3])) |
                         ((df[genotype1].str.len() == 2) &
                          (df[genotype3].str.len() == 2) &
                          (df[genotype1].str[0] != df[genotype3].str[0]) &
                          (df[genotype1].str[0] != df[genotype3].str[1]) &
                          (df[genotype1].str[1] != df[genotype3].str[0]) &
                          (df[genotype1].str[1] != df[genotype3].str[1]))) |
                        (df[genotype2].notnull() & df[genotype3].notnull() &
                         (df[genotype2].str.len() == 2) &
                         (df[genotype2].str[0] == df[genotype2].str[1]) &
                         (df[genotype2] == df[genotype3]) &
                         (df[genotype1] != df[genotype2]))]

            if save_output:
                save_df_as_csv(df, self._output_dir, 'discordant_snps_' +
                               individual1.get_var_name() + '_' +
                               individual2.get_var_name() + '_' +
                               individual3.get_var_name() + '.csv')

        return df

    def find_shared_dna(self, individual1, individual2, build=37, cM_threshold=0.75,
                        snp_threshold=1100, shared_genes=False):
        """ Find the shared DNA between two individuals.

        Computes the genetic distance in centimorgans (cMs) between SNPs using the specified
        assembly's HapMap. Applies thresholds to determine the shared DNA. Plots shared DNA.
        Optionally determines shared genes (i.e., genes that are transcribed from the shared DNA).

        All output is saved to the output directory as `CSV` or `PNG` files.

        Parameters
        ----------
        individual1 : Individual
        individual2 : Individual
        build : {36, 37}
            human genome assembly
        cM_threshold : float
            minimum centimorgans for each shared DNA segment
        snp_threshold : int
            minimum SNPs for each shared DNA segment
        shared_genes : bool
            determine shared genes

        """
        df = individual1.snps

        df = df.join(individual2.snps['genotype'], rsuffix='2', how='inner')

        genotype1 = 'genotype_' + individual1.get_var_name()
        genotype2 = 'genotype_' + individual2.get_var_name()

        df = df.rename(columns={'genotype': genotype1, 'genotype2': genotype2})

        one_x_chrom = self._is_one_individual_male(df, genotype1, genotype2)

        # determine the genetic distance between each SNP using HapMap tables
        hapmap, df = self._compute_snp_distances(df, build)

        # determine where individuals share an allele on one chromosome
        df['one_chrom_match'] = np.where(
            df[genotype1].isnull() | df[genotype2].isnull() |
            (df[genotype1].str[0] == df[genotype2].str[0]) |
            (df[genotype1].str[0] == df[genotype2].str[1]) |
            (df[genotype1].str[1] == df[genotype2].str[0]) |
            (df[genotype1].str[1] == df[genotype2].str[1]), True, False)

        # determine where individuals share alleles on both chromosomes
        df['two_chrom_match'] = np.where(
            df[genotype1].isnull() | df[genotype2].isnull() |
            ((df[genotype1].str.len() == 2) &
             (df[genotype2].str.len() == 2) &
             ((df[genotype1] == df[genotype2]) |
              ((df[genotype1].str[0] == df[genotype2].str[1]) &
               (df[genotype1].str[1] == df[genotype2].str[0])))), True, False)

        # compute shared DNA between individuals
        one_chrom_shared_dna = self._compute_shared_dna(df, hapmap, 'one_chrom_match',
                                                        cM_threshold, snp_threshold, one_x_chrom)

        two_chrom_shared_dna = self._compute_shared_dna(df, hapmap, 'two_chrom_match',
                                                        cM_threshold, snp_threshold, one_x_chrom)

        if build == 36:
            cytobands = self._resources.get_cytoband_h36()
        else:
            cytobands = self._resources.get_cytoband_h37()

        # plot data
        if dir_exists(self._output_dir):
            plot_chromosomes(one_chrom_shared_dna, two_chrom_shared_dna, cytobands,
                             os.path.join(self._output_dir, 'shared_dna_' +
                                          individual1.get_var_name() + '_' +
                                          individual2.get_var_name() + '.png'),
                             individual1.name + ' / ' + individual2.name + ' shared DNA', build)

        if shared_genes and build == 36:
            print('Remap SNPs to Build 37 to determine shared genes')
            shared_genes = False

        # save results in CSV format
        if len(one_chrom_shared_dna) > 0:
            self._save_shared_dna_csv_format(one_chrom_shared_dna, 'one',
                                             individual1.get_var_name(),
                                             individual2.get_var_name())
            if shared_genes:
                self._compute_shared_genes(one_chrom_shared_dna, 'one',
                                           individual1.get_var_name(),
                                           individual2.get_var_name())

        if len(two_chrom_shared_dna) > 0:
            self._save_shared_dna_csv_format(two_chrom_shared_dna, 'two',
                                             individual1.get_var_name(),
                                             individual2.get_var_name())
            if shared_genes:
                self._compute_shared_genes(two_chrom_shared_dna, 'two',
                                           individual1.get_var_name(),
                                           individual2.get_var_name())

    def _compute_shared_genes(self, shared_dna, type, individual1_name, individual2_name):
        knownGenes = self._resources.get_knownGene_h37()
        kgXref = self._resources.get_kgXref_h37()

        # http://seqanswers.com/forums/showthread.php?t=22336
        df = knownGenes.join(kgXref)

        shared_genes = []

        for shared_segment in shared_dna:
            # determine genes transcribed from this shared DNA segment
            temp = df.loc[(df['chrom'] == shared_segment['chrom']) &
                          (df['txStart'] >= shared_segment['start']) &
                          (df['txEnd'] <= shared_segment['end'])].copy()

            # select subset of columns
            temp = temp[['geneSymbol', 'chrom', 'strand', 'txStart', 'txEnd', 'refseq',
                         'proteinID', 'description']]

            if len(temp) > 0:
                shared_genes.append(temp)

        if len(shared_genes) > 0:
            if type == 'one':
                chroms = 'one_chrom'
            else:
                chroms = 'two_chroms'

            file = os.path.join(self._output_dir, 'shared_genes_' + chroms + '_' +
                                individual1_name + '_' + individual2_name + '.csv')

            save_df_as_csv(pd.concat(shared_genes), self._output_dir, file)

    def _save_shared_dna_csv_format(self, shared_dna, type, individual1_name, individual2_name):
        if type == 'one':
            chroms = 'one_chrom'
        else:
            chroms = 'two_chroms'

        file = os.path.join(self._output_dir, 'shared_dna_' + chroms + '_' + individual1_name +
                            '_' + individual2_name + '.csv')

        print('Saving ' + os.path.relpath(file))

        with open(file, 'w') as f:
            f.write('chrom,start,stop,cMs,snps\n')

            for shared_segment in shared_dna:
                f.write('{chrom},{start:d},{stop:d},{cMs:.2f},{snps:d}\n'.format(
                            chrom=shared_segment['chrom'], start=shared_segment['start'],
                            stop=shared_segment['end'], cMs=shared_segment['cMs'],
                            snps=shared_segment['snps']))

    def _is_one_individual_male(self, df, genotype1, genotype2, heterozygous_x_snp_threshold=100):
        # determine if at least one individual is male by counting heterozygous X SNPs

        if len(df.loc[(df['chrom'] == 'X') &
                      (df[genotype1].notnull()) &
                      (df[genotype1].str[0] != df[genotype1].str[1])]) < heterozygous_x_snp_threshold:
            return True
        else:
            if len(df.loc[(df['chrom'] == 'X') &
                          (df[genotype2].notnull()) &
                          (df[genotype2].str[0] != df[genotype2].str[1])]) < heterozygous_x_snp_threshold:
                return True
            else:
                return False


    def _compute_snp_distances(self, df, build):
        if build == 36:
            hapmap = self._resources.get_hapmap_h36()
        else:
            hapmap = self._resources.get_hapmap_h37()

        for chrom in df['chrom'].unique():
            if chrom not in hapmap.keys():
                continue

            # create a new dataframe from the positions for the current chromosome
            temp = pd.DataFrame(df.loc[(df['chrom'] == chrom)]['pos'].values, columns=['pos'])

            # merge HapMap for this chrom
            temp = temp.append(hapmap[chrom], ignore_index=True)

            # sort based on pos
            temp = temp.sort_values('pos')

            # fill cM rates forward and backward
            temp['rate'] = temp['rate'].fillna(method='ffill')
            temp['rate'] = temp['rate'].fillna(method='bfill')

            # get difference between positions
            pos_diffs = np.ediff1d(temp['pos'])

            # compute cMs between each pos based on probablistic recombination rate
            # https://www.biostars.org/p/123539/
            cMs_match_segment = (temp['rate'] * np.r_[pos_diffs, 0] / 1e6).values

            # add back into temp
            temp['cMs'] = np.r_[0, cMs_match_segment][:-1]

            temp = temp.reset_index()
            del temp['index']

            # use null `map` values to find locations of SNPs
            snp_indices = temp.loc[temp['map'].isnull()].index

            # use SNP indices to determine boundaries over which to sum cMs
            start_snp_ix = snp_indices + 1
            end_snp_ix = np.r_[snp_indices, snp_indices[-1]][1:] + 1
            snp_boundaries = np.c_[start_snp_ix, end_snp_ix]

            # sum cMs between SNPs to get total cM distance between SNPs
            # http://stackoverflow.com/a/7471967
            c = np.r_[0, temp['cMs'].cumsum()][snp_boundaries]
            cM_from_prev_snp = c[:, 1] - c[:, 0]

            # debug
            # temp.loc[snp_indices, 'cM_from_prev_snp'] = np.r_[0, cM_from_prev_snp][:-1]
            # temp.to_csv('debug.csv')

            # add back into df
            df.loc[(df['chrom'] == chrom), 'cM_from_prev_snp'] = np.r_[0, cM_from_prev_snp][:-1]

        return hapmap, df

    def _compute_shared_dna(self, df, hapmap, col, cM_threshold, snp_threshold, one_x_chrom):
        shared_dna = []

        for chrom in df['chrom'].unique():
            if chrom not in hapmap.keys():
                continue

            # skip calculating two chrom shared on the X chromosome if an individual is male
            if chrom == 'X' and col == 'two_chrom_match' and one_x_chrom:
                continue

            # get consecutive strings of trues
            # http://stackoverflow.com/a/17151327
            a = df.loc[(df['chrom'] == chrom)][col].values
            a = np.r_[a, False]
            a_rshifted = np.roll(a, 1)
            starts = a & ~a_rshifted
            ends = ~a & a_rshifted
            a_starts = np.nonzero(starts)[0]
            a_starts = np.reshape(a_starts, (len(a_starts), 1))
            a_ends = np.nonzero(ends)[0]
            a_ends = np.reshape(a_ends, (len(a_ends), 1))

            matches = np.hstack((a_starts, a_ends))

            c = np.r_[0, df.loc[(df['chrom'] == chrom)]['cM_from_prev_snp'].cumsum()][matches]
            cMs_match_segment = c[:, 1] - c[:, 0]

            # get matching segments where total cMs is greater than the threshold
            matches_passed = matches[np.where(cMs_match_segment > cM_threshold)]

            # get indices where a_end boundary is adjacent to a_start, perhaps indicating a
            # discrepant SNP
            adjacent_ix = np.where(np.roll(matches_passed[:, 0], -1) - matches_passed[:, 1] == 1)

            # if there are adjacent segments
            if len(adjacent_ix[0]) != 0:
                matches_stitched = np.array([0, 0])
                prev = -1
                counter = 0

                # stitch together adjacent segments
                for x in matches_passed:
                    if prev != -1 and prev + 1 == x[0]:
                        matches_stitched[counter, 1] = x[1]
                    else:
                        matches_stitched = np.vstack((matches_stitched, x))
                        counter += 1

                    prev = x[1]

                matches_passed = matches_stitched[1:]

            snp_counts = matches_passed[:, 1] - matches_passed[:, 0]

            # apply SNP count threshold for each matching segment; we now have the shared DNA
            # segments that pass the centimorgan and SNP thresholds
            matches_passed = matches_passed[np.where(snp_counts > snp_threshold)]

            # compute total cMs for each match segment
            c = np.r_[0, df.loc[(df['chrom'] == chrom)]['cM_from_prev_snp'].cumsum()][matches_passed]
            cMs_match_segment = c[:, 1] - c[:, 0]

            counter = 0
            # save matches for this chromosome
            if col == 'one_chrom_match':
                chrom_stain = 'one_chrom'
            else:
                chrom_stain = 'two_chrom'

            for x in matches_passed:
                shared_dna.append({"chrom": chrom,
                                   "start": df.loc[(df['chrom'] == chrom)].ix[x[0]].pos,
                                   "end": df.loc[(df['chrom'] == chrom)].ix[x[1] - 1].pos,
                                   "cMs": cMs_match_segment[counter],
                                   "snps": x[1] - x[0],
                                   "gie_stain": chrom_stain})
                counter += 1
        return shared_dna

def dir_exists(path):
    # https://stackoverflow.com/a/5032238
    try:
        os.makedirs(path, exist_ok=True)
    except Exception as err:
        print(err)
        return False

    if os.path.exists(path):
        return True
    else:
        return False


def save_df_as_csv(df, path, filename):
    if dir_exists(path):
        destination = os.path.join(path, filename)
        print('Saving ' + os.path.relpath(destination))
        df.to_csv(destination, na_rep='--')
