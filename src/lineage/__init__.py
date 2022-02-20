""" lineage

tools for genetic genealogy and the analysis of consumer DNA test results

"""

"""
MIT License

Copyright (c) 2016 Andrew Riha

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import datetime
from itertools import chain, combinations
import logging
import os

import numpy as np
import pandas as pd
from snps.utils import Parallelizer, create_dir, save_df_as_csv

from lineage.individual import Individual
from lineage.resources import Resources
from lineage.visualization import plot_chromosomes

# set version string with Versioneer
from lineage._version import get_versions

__version__ = get_versions()["version"]
del get_versions

logger = logging.getLogger(__name__)


class Lineage:
    """Object used to interact with the `lineage` framework."""

    def __init__(
        self,
        output_dir="output",
        resources_dir="resources",
        parallelize=False,
        processes=os.cpu_count(),
    ):
        """Initialize a ``Lineage`` object.

        Parameters
        ----------
        output_dir : str
            name / path of output directory
        resources_dir : str
            name / path of resources directory
        parallelize : bool
            utilize multiprocessing to speedup calculations
        processes : int
            processes to launch if multiprocessing
        """
        self._output_dir = output_dir
        self._resources_dir = resources_dir
        self._resources = Resources(resources_dir=resources_dir)
        self._parallelizer = Parallelizer(parallelize=parallelize, processes=processes)

    def create_individual(self, name, raw_data=(), **kwargs):
        """Initialize an individual in the context of the `lineage` framework.

        Parameters
        ----------
        name : str
            name of the individual
        raw_data : str, bytes, ``SNPs`` (or list or tuple thereof)
            path(s) to file(s), bytes, or ``SNPs`` object(s) with raw genotype data
        **kwargs
            parameters to ``snps.SNPs`` and/or ``snps.SNPs.merge``

        Returns
        -------
        Individual
            ``Individual`` initialized in the context of the `lineage` framework
        """
        if "output_dir" not in kwargs:
            kwargs["output_dir"] = self._output_dir
        if "resources_dir" not in kwargs:
            kwargs["resources_dir"] = self._resources_dir

        return Individual(name, raw_data, **kwargs)

    def download_example_datasets(self):
        """Download example datasets from `openSNP <https://opensnp.org>`_.

        Per openSNP, "the data is donated into the public domain using `CC0 1.0
        <http://creativecommons.org/publicdomain/zero/1.0/>`_."

        Returns
        -------
        paths : list of str or empty str
            paths to example datasets

        References
        ----------
        1. Greshake B, Bayer PE, Rausch H, Reda J (2014), "openSNP-A Crowdsourced Web Resource
           for Personal Genomics," PLOS ONE, 9(3): e89204,
           https://doi.org/10.1371/journal.pone.0089204
        """
        paths = self._resources.download_example_datasets()

        if "" in paths:
            logger.warning("Example dataset(s) not currently available")

        return paths

    def find_discordant_snps(
        self, individual1, individual2, individual3=None, save_output=False
    ):
        """Find discordant SNPs between two or three individuals.

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
        1. David Pike, "Search for Discordant SNPs in Parent-Child
           Raw Data Files," David Pike's Utilities,
           http://www.math.mun.ca/~dapike/FF23utils/pair-discord.php
        2. David Pike, "Search for Discordant SNPs when given data
           for child and both parents," David Pike's Utilities,
           http://www.math.mun.ca/~dapike/FF23utils/trio-discord.php
        """
        self._remap_snps_to_GRCh37([individual1, individual2, individual3])

        df = individual1.snps

        # remove nulls for reference individual
        df = df.loc[df["genotype"].notnull()]

        # add SNPs shared with `individual2`
        df = df.join(individual2.snps["genotype"], rsuffix="2")

        genotype1 = "genotype_" + individual1.get_var_name()
        genotype2 = "genotype_" + individual2.get_var_name()

        if individual3 is None:
            df = df.rename(columns={"genotype": genotype1, "genotype2": genotype2})

            # find discordant SNPs between reference and comparison individuals
            df = df.loc[
                df[genotype2].notnull()
                & (
                    (df[genotype1].str.len() == 1)
                    & (df[genotype2].str.len() == 1)
                    & (df[genotype1] != df[genotype2])
                )
                | (
                    (df[genotype1].str.len() == 2)
                    & (df[genotype2].str.len() == 2)
                    & (df[genotype1].str[0] != df[genotype2].str[0])
                    & (df[genotype1].str[0] != df[genotype2].str[1])
                    & (df[genotype1].str[1] != df[genotype2].str[0])
                    & (df[genotype1].str[1] != df[genotype2].str[1])
                )
            ]
            if save_output:
                save_df_as_csv(
                    df,
                    self._output_dir,
                    "discordant_snps_{}_{}_GRCh37.csv".format(
                        individual1.get_var_name(), individual2.get_var_name()
                    ),
                    comment=self._get_csv_header(),
                    prepend_info=False,
                )
        else:
            # add SNPs shared with `individual3`
            df = df.join(individual3.snps["genotype"], rsuffix="3")

            genotype3 = "genotype_" + individual3.get_var_name()

            df = df.rename(
                columns={
                    "genotype": genotype1,
                    "genotype2": genotype2,
                    "genotype3": genotype3,
                }
            )

            # find discordant SNPs between child and two parents
            df = df.loc[
                (
                    df[genotype2].notnull()
                    & (
                        (df[genotype1].str.len() == 1)
                        & (df[genotype2].str.len() == 1)
                        & (df[genotype1] != df[genotype2])
                    )
                    | (
                        (df[genotype1].str.len() == 2)
                        & (df[genotype2].str.len() == 2)
                        & (df[genotype1].str[0] != df[genotype2].str[0])
                        & (df[genotype1].str[0] != df[genotype2].str[1])
                        & (df[genotype1].str[1] != df[genotype2].str[0])
                        & (df[genotype1].str[1] != df[genotype2].str[1])
                    )
                )
                | (
                    df[genotype3].notnull()
                    & (
                        (df[genotype1].str.len() == 1)
                        & (df[genotype3].str.len() == 1)
                        & (df[genotype1] != df[genotype3])
                    )
                    | (
                        (df[genotype1].str.len() == 2)
                        & (df[genotype3].str.len() == 2)
                        & (df[genotype1].str[0] != df[genotype3].str[0])
                        & (df[genotype1].str[0] != df[genotype3].str[1])
                        & (df[genotype1].str[1] != df[genotype3].str[0])
                        & (df[genotype1].str[1] != df[genotype3].str[1])
                    )
                )
                | (
                    df[genotype2].notnull()
                    & df[genotype3].notnull()
                    & (df[genotype2].str.len() == 2)
                    & (df[genotype2].str[0] == df[genotype2].str[1])
                    & (df[genotype2] == df[genotype3])
                    & (df[genotype1] != df[genotype2])
                )
            ]

            if save_output:
                save_df_as_csv(
                    df,
                    self._output_dir,
                    "discordant_snps_{}_{}_{}_GRCh37.csv".format(
                        individual1.get_var_name(),
                        individual2.get_var_name(),
                        individual3.get_var_name(),
                    ),
                    comment=self._get_csv_header(),
                    prepend_info=False,
                )

        return df

    def find_shared_dna(
        self,
        individuals=(),
        cM_threshold=0.75,
        snp_threshold=1100,
        shared_genes=False,
        save_output=True,
        genetic_map="HapMap2",
    ):
        """Find the shared DNA between individuals.

        Computes the genetic distance in centiMorgans (cMs) between SNPs using the specified genetic
        map. Applies thresholds to determine the shared DNA. Plots shared DNA. Optionally determines
        shared genes (i.e., genes transcribed from the shared DNA).

        All output is saved to the output directory as `CSV` or `PNG` files.

        Notes
        -----
        The code is commented throughout to help describe the algorithm and its operation.

        To summarize, the algorithm first computes the genetic distance in cMs between SNPs
        common to all individuals using the specified genetic map.

        Then, individuals are compared for whether they share one or two alleles for each SNP in
        common; in this manner, where all individuals share one chromosome, for example, there
        will be several SNPs in a row where at least one allele is shared between individuals for
        each SNP. The ``cM_threshold`` is then applied to each of these "matching segments" to
        determine whether the segment could be a potential shared DNA segment (i.e., whether each
        segment has a cM value greater than the threshold).

        The matching segments that passed the ``cM_threshold`` are then checked to see if they
        are adjacent to another matching segment, and if so, the segments are stitched together,
        and the single SNP separating the segments is flagged as potentially discrepant. (This
        means that multiple smaller matching segments passing the ``cM_threshold`` could be
        stitched, identifying the SNP between each segment as discrepant.)

        Next, the ``snp_threshold`` is applied to each segment to ensure there are enough SNPs in
        the segment and the segment is not only a few SNPs in a region with a high recombination
        rate; for each segment that passes this test, we have a segment of shared DNA, and the
        total cMs for this segment are computed.

        Finally, discrepant SNPs are checked to ensure that only SNPs internal to a shared DNA
        segment are reported as discrepant (i.e., don't report a discrepant SNP if it was part of a
        segment that didn't pass the ``snp_threshold``). Currently, no action other than reporting
        is taken on discrepant SNPs.

        Parameters
        ----------
        individuals : iterable of Individuals
        cM_threshold : float
            minimum centiMorgans for each shared DNA segment
        snp_threshold : int
            minimum SNPs for each shared DNA segment
        shared_genes : bool
            determine shared genes
        save_output : bool
            specifies whether to save output files in the output directory
        genetic_map : {'HapMap2', 'ACB', 'ASW', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'FIN', 'GBR', 'GIH', 'IBS', 'JPT', 'KHV', 'LWK', 'MKK', 'MXL', 'PEL', 'PUR', 'TSI', 'YRI'}
            genetic map to use for computation of shared DNA; `HapMap2` corresponds to the HapMap
            Phase II genetic map from the
            `International HapMap Project <https://www.genome.gov/10001688/international-hapmap-project/>`_
            and all others correspond to the
            `population-specific <https://www.internationalgenome.org/faq/which-populations-are-part-your-study/>`_
            genetic maps generated from the
            `1000 Genomes Project <https://www.internationalgenome.org>`_ phased OMNI data.
            Note that shared DNA is not computed on the X chromosome with the 1000 Genomes
            Project genetic maps since the X chromosome is not included in these genetic maps.

        Returns
        -------
        dict
            dict with the following items:

            one_chrom_shared_dna (pandas.DataFrame)
                segments of shared DNA on one chromosome
            two_chrom_shared_dna (pandas.DataFrame)
                segments of shared DNA on two chromosomes
            one_chrom_shared_genes (pandas.DataFrame)
                shared genes on one chromosome
            two_chrom_shared_genes (pandas.DataFrame)
                shared genes on two chromosomes
            one_chrom_discrepant_snps (pandas.Index)
                discrepant SNPs discovered while finding shared DNA on one chromosome
            two_chrom_discrepant_snps (pandas.Index)
                discrepant SNPs discovered while finding shared DNA on two chromosomes
        """
        # initialize all objects to be returned to be empty to start
        one_chrom_shared_dna = pd.DataFrame()
        two_chrom_shared_dna = pd.DataFrame()
        one_chrom_shared_genes = pd.DataFrame()
        two_chrom_shared_genes = pd.DataFrame()
        one_chrom_discrepant_snps = pd.Index([])
        two_chrom_discrepant_snps = pd.Index([])

        # ensure that all individuals have SNPs that are mapped relative to Build 37
        self._remap_snps_to_GRCh37(individuals)

        # return if there aren't enough individuals to compare
        if len(individuals) < 2:
            logger.warning("find_shared_dna requires two or more individuals...")
            return self._find_shared_dna_return_helper(
                one_chrom_shared_dna,
                two_chrom_shared_dna,
                one_chrom_shared_genes,
                two_chrom_shared_genes,
                one_chrom_discrepant_snps,
                two_chrom_discrepant_snps,
            )

        # load the specified genetic map (one genetic map for each chromosome)
        genetic_map_dfs = self._resources.get_genetic_map(genetic_map)

        if len(genetic_map_dfs) == 0:
            return self._find_shared_dna_return_helper(
                one_chrom_shared_dna,
                two_chrom_shared_dna,
                one_chrom_shared_genes,
                two_chrom_shared_genes,
                one_chrom_discrepant_snps,
                two_chrom_discrepant_snps,
            )

        # generate a list of dynamically named columns for each individual's genotype
        # (e.g., genotype0, genotype1, etc).
        cols = [f"genotype{str(i)}" for i in range(len(individuals))]

        # set the reference SNPs to compare to be that of the first individual
        df = individuals[0].snps
        df = df.rename(columns={"genotype": cols[0]})

        # build-up a dataframe of SNPs that are common to all individuals
        for i, ind in enumerate(individuals[1:]):
            # join SNPs for all individuals
            df = df.join(ind.snps["genotype"], how="inner")
            df = df.rename(columns={"genotype": cols[i + 1]})

        # set a flag for if one individual is male (i.e., only one chromosome match on the X
        # chromosome is possible in the non-PAR region)
        one_x_chrom = self._is_one_individual_male(individuals)

        # create tasks to compute the genetic distances (cMs) between each SNP on each chromosome
        tasks = []
        chroms_to_drop = []
        for chrom in df["chrom"].unique():
            if chrom not in genetic_map_dfs.keys():
                chroms_to_drop.append(chrom)
                continue

            # each task requires the genetic map for the chromosome and the positions of all SNPs
            # in common on that chromosome
            tasks.append(
                {
                    "genetic_map": genetic_map_dfs[chrom],
                    # get positions for the current chromosome
                    "snps": pd.DataFrame(df.loc[(df["chrom"] == chrom)]["pos"]),
                }
            )

        # drop chromosomes without genetic distance data (e.g., chroms MT, PAR, etc.)
        for chrom in chroms_to_drop:
            df = df.drop(df.loc[df["chrom"] == chrom].index)

        # determine the genetic distance between each SNP using the specified genetic map
        snp_distances = map(self._compute_snp_distances, tasks)
        snp_distances = pd.concat(snp_distances)

        # extract the column "cM_from_prev_snp" from the result and add that to the dataframe
        # of SNPs common to all individuals; now we have the genetic distance between each SNP
        df["cM_from_prev_snp"] = snp_distances["cM_from_prev_snp"]

        # now we apply a mask for whether all individuals match on one or two chromosomes...
        # first, set all rows for these columns to True
        df["one_chrom_match"] = True
        df["two_chrom_match"] = True
        # determine where individuals share an allele on one chromosome (i.e., set to False when
        # at least one allele doesn't match for all individuals)
        for genotype1, genotype2 in combinations(cols, 2):
            df.loc[
                ~df[genotype1].isnull()
                & ~df[genotype2].isnull()
                & (df[genotype1].str[0] != df[genotype2].str[0])
                & (df[genotype1].str[0] != df[genotype2].str[1])
                & (df[genotype1].str[1] != df[genotype2].str[0])
                & (df[genotype1].str[1] != df[genotype2].str[1]),
                "one_chrom_match",
            ] = False

        # determine where individuals share alleles on two chromosomes (i.e., set to False when
        # two alleles don't match for all individuals)
        for genotype1, genotype2 in combinations(cols, 2):
            df.loc[
                ~df[genotype1].isnull()
                & ~df[genotype2].isnull()
                & (df[genotype1] != df[genotype2])
                & ~(
                    (df[genotype1].str[0] == df[genotype2].str[1])
                    & (df[genotype1].str[1] == df[genotype2].str[0])
                ),
                "two_chrom_match",
            ] = False

        # genotype columns are no longer required for calculation
        df = df.drop(cols, axis=1)

        # find shared DNA on one chrom
        one_chrom_shared_dna, one_chrom_discrepant_snps = self._find_shared_dna_helper(
            df[["chrom", "pos", "cM_from_prev_snp", "one_chrom_match"]],
            cM_threshold,
            snp_threshold,
            one_x_chrom,
        )
        # find shared DNA on two chroms
        two_chrom_shared_dna, two_chrom_discrepant_snps = self._find_shared_dna_helper(
            df[["chrom", "pos", "cM_from_prev_snp", "two_chrom_match"]],
            cM_threshold,
            snp_threshold,
            one_x_chrom,
        )

        if shared_genes:
            one_chrom_shared_genes = self._compute_shared_genes(one_chrom_shared_dna)
            two_chrom_shared_genes = self._compute_shared_genes(two_chrom_shared_dna)

        if save_output:
            self._find_shared_dna_output_helper(
                individuals,
                cM_threshold,
                snp_threshold,
                one_chrom_shared_dna,
                two_chrom_shared_dna,
                one_chrom_shared_genes,
                two_chrom_shared_genes,
                genetic_map,
            )

        return self._find_shared_dna_return_helper(
            one_chrom_shared_dna,
            two_chrom_shared_dna,
            one_chrom_shared_genes,
            two_chrom_shared_genes,
            one_chrom_discrepant_snps,
            two_chrom_discrepant_snps,
        )

    def _find_shared_dna_helper(self, df, cM_threshold, snp_threshold, one_x_chrom):
        tasks = []

        for chrom in df["chrom"].unique():
            tasks.append(
                {
                    "df": df.loc[df["chrom"] == chrom],
                    "chrom": chrom,
                    "cM_threshold": cM_threshold,
                    "snp_threshold": snp_threshold,
                    "one_x_chrom": one_x_chrom,
                }
            )

        # compute shared DNA between individuals
        results = map(self._compute_shared_dna, tasks)

        shared_dna = []
        discrepant_snps = pd.Index([], name="rsid")
        for result in list(results):
            shared_dna.append(result["shared_dna"])
            discrepant_snps = discrepant_snps.append(result["discrepant_snps"])

        # https://stackoverflow.com/a/952952
        shared_dna = list(chain.from_iterable(shared_dna))

        return (self._convert_shared_dna_list_to_df(shared_dna), discrepant_snps)

    def _find_shared_dna_output_helper(
        self,
        individuals,
        cM_threshold,
        snp_threshold,
        one_chrom_shared_dna,
        two_chrom_shared_dna,
        one_chrom_shared_genes,
        two_chrom_shared_genes,
        genetic_map,
    ):
        def output_csv(df, file, float_format="%.2f"):
            save_df_as_csv(
                df,
                self._output_dir,
                file,
                comment=self._get_csv_header(),
                prepend_info=False,
                float_format=float_format,
            )

        cytobands = self._resources.get_cytoBand_hg19()

        individuals_filename = ""
        individuals_plot_title = ""

        for individual in individuals:
            individuals_filename += individual.get_var_name() + "_"
            individuals_plot_title += individual.name + " / "

        individuals_filename = individuals_filename[:-1]
        individuals_plot_title = individuals_plot_title[:-3]

        cM = "{:.2f}".format(cM_threshold).replace(".", "p")
        filename_details = (
            f"{individuals_filename}_{cM}cM_{snp_threshold}snps_GRCh37_{genetic_map}"
        )

        if create_dir(self._output_dir):
            plot_chromosomes(
                one_chrom_shared_dna,
                two_chrom_shared_dna,
                cytobands,
                os.path.join(
                    self._output_dir,
                    f"shared_dna_{filename_details}.png",
                ),
                f"{individuals_plot_title} shared DNA",
                37,
            )

        if len(one_chrom_shared_dna) > 0:
            output_csv(
                one_chrom_shared_dna,
                f"shared_dna_one_chrom_{filename_details}.csv",
            )

        if len(two_chrom_shared_dna) > 0:
            output_csv(
                two_chrom_shared_dna,
                f"shared_dna_two_chroms_{filename_details}.csv",
            )

        if len(one_chrom_shared_genes) > 0:
            output_csv(
                one_chrom_shared_genes,
                f"shared_genes_one_chrom_{filename_details}.csv",
                None,
            )

        if len(two_chrom_shared_genes) > 0:
            output_csv(
                two_chrom_shared_genes,
                f"shared_genes_two_chroms_{filename_details}.csv",
                None,
            )

    def _find_shared_dna_return_helper(
        self,
        one_chrom_shared_dna,
        two_chrom_shared_dna,
        one_chrom_shared_genes,
        two_chrom_shared_genes,
        one_chrom_discrepant_snps,
        two_chrom_discrepant_snps,
    ):

        return {
            "one_chrom_shared_dna": one_chrom_shared_dna,
            "two_chrom_shared_dna": two_chrom_shared_dna,
            "one_chrom_shared_genes": one_chrom_shared_genes,
            "two_chrom_shared_genes": two_chrom_shared_genes,
            "one_chrom_discrepant_snps": one_chrom_discrepant_snps,
            "two_chrom_discrepant_snps": two_chrom_discrepant_snps,
        }

    def _convert_shared_dna_list_to_df(self, shared_dna):
        df = pd.DataFrame(shared_dna, columns=["chrom", "start", "end", "cMs", "snps"])
        df.index.name = "segment"
        df.index = df.index + 1
        return df

    def _compute_shared_genes(self, shared_dna):
        knownGenes = self._resources.get_knownGene_hg19()
        kgXref = self._resources.get_kgXref_hg19()

        # http://seqanswers.com/forums/showthread.php?t=22336
        df = knownGenes.join(kgXref)

        shared_genes_dfs = []
        shared_genes = pd.DataFrame()

        for shared_segment in shared_dna.itertuples():
            # determine genes transcribed from this shared DNA segment
            temp = df.loc[
                (df["chrom"] == shared_segment.chrom)
                & (df["txStart"] >= shared_segment.start)
                & (df["txEnd"] <= shared_segment.end)
            ].copy()

            # select subset of columns
            temp = temp[
                [
                    "geneSymbol",
                    "chrom",
                    "strand",
                    "txStart",
                    "txEnd",
                    "refseq",
                    "proteinID",
                    "description",
                ]
            ]

            if len(temp) > 0:
                shared_genes_dfs.append(temp)

        if len(shared_genes_dfs) > 0:
            shared_genes = pd.concat(shared_genes_dfs, sort=True)

        return shared_genes

    def _is_one_individual_male(self, individuals):
        for ind in individuals:
            if ind.sex == "Male":
                return True
        return False

    def _compute_snp_distances(self, task):
        """Compute genetic distance in cMs between SNPs.

        Parameters
        ----------
        task : dict
            dict with `snps` to compute distance between using `genetic_map`

        Returns
        -------
        pandas.DataFrame
            genetic distances between SNPs
        """
        genetic_map = task["genetic_map"]
        temp = task["snps"]

        # merge genetic map for this chrom
        temp = pd.concat([temp, genetic_map], ignore_index=False, sort=True)

        # sort based on pos
        temp = temp.sort_values("pos")

        # fill recombination rates forward
        temp["rate"] = temp["rate"].fillna(method="ffill")

        # assume recombination rate of 0 for SNPs upstream of first defined rate
        temp["rate"] = temp["rate"].fillna(0)

        # get difference between positions
        pos_diffs = np.ediff1d(temp["pos"])

        # compute cMs between each pos based on probabilistic recombination rate
        # https://www.biostars.org/p/123539/
        cMs_match_segment = (temp["rate"] * np.r_[pos_diffs, 0] / 1e6).values

        # add back into temp
        temp["cMs"] = np.r_[0, cMs_match_segment][:-1]

        temp = temp.reset_index()

        # use null `map` values to find locations of SNPs
        snp_indices = temp.loc[temp["map"].isnull()].index

        # use SNP indices to determine boundaries over which to sum cMs
        start_snp_ix = snp_indices + 1
        end_snp_ix = np.r_[snp_indices, snp_indices[-1]][1:] + 1
        snp_boundaries = np.c_[start_snp_ix, end_snp_ix]

        # sum cMs between SNPs to get total cM distance between SNPs
        # http://stackoverflow.com/a/7471967
        c = np.r_[0, temp["cMs"].cumsum()][snp_boundaries]
        cM_from_prev_snp = c[:, 1] - c[:, 0]

        temp = temp.loc[temp["map"].isna()]

        # add back into temp
        temp["cM_from_prev_snp"] = np.r_[0, cM_from_prev_snp][:-1]

        # restore index
        temp = temp.set_index("index")

        return pd.DataFrame(temp["cM_from_prev_snp"])

    def _compute_shared_dna(self, task):
        df = task["df"]
        chrom = task["chrom"]
        cM_threshold = task["cM_threshold"]
        snp_threshold = task["snp_threshold"]
        one_x_chrom = task["one_x_chrom"]

        if "one_chrom_match" in df.keys():
            match_col = "one_chrom_match"
        else:
            match_col = "two_chrom_match"

        shared_dna = []
        discrepant_snps = pd.Index([], name="rsid")

        # set two_chrom_match in non-PAR region to False if an individual is male
        if chrom == "X" and match_col == "two_chrom_match" and one_x_chrom:
            df = df.copy()
            # https://www.ncbi.nlm.nih.gov/grc/human
            df.loc[
                (df["chrom"] == "X") & (df["pos"] > 2699520) & (df["pos"] < 154931044),
                "two_chrom_match",
            ] = False

        # get consecutive strings of Trues, for where there's a one or two chrom match between
        # individuals, depending on the task; http://stackoverflow.com/a/17151327
        a = df.loc[(df["chrom"] == chrom)][match_col].values
        a = np.r_[a, False]
        a_rshifted = np.roll(a, 1)
        starts = a & ~a_rshifted
        ends = ~a & a_rshifted
        a_starts = np.nonzero(starts)[0]
        a_starts = np.reshape(a_starts, (len(a_starts), 1))
        a_ends = np.nonzero(ends)[0]
        a_ends = np.reshape(a_ends, (len(a_ends), 1))

        # get the matching segments
        matches = np.hstack((a_starts, a_ends))

        # compute total cMs for each matching segment
        c = np.r_[0, df.loc[(df["chrom"] == chrom)]["cM_from_prev_snp"].cumsum()][
            matches
        ]
        cMs_match_segment = c[:, 1] - c[:, 0]

        # get matching segments where total cMs is greater than the threshold
        matches_passed = matches[np.where(cMs_match_segment > cM_threshold)]

        # get indices where a_end boundary is adjacent to a_start, perhaps indicating a
        # discrepant SNP
        adjacent_ix = np.where(
            np.roll(matches_passed[:, 0], -1) - matches_passed[:, 1] == 1
        )

        # if there are adjacent segments
        if len(adjacent_ix[0]) != 0:
            discrepant_snps = df.iloc[matches_passed[adjacent_ix[0], 1]].index
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
        # segments that pass the centiMorgan and SNP thresholds
        matches_passed = matches_passed[np.where(snp_counts > snp_threshold)]

        # compute total cMs for each match segment
        c = np.r_[0, df.loc[(df["chrom"] == chrom)]["cM_from_prev_snp"].cumsum()][
            matches_passed
        ]
        cMs_match_segment = c[:, 1] - c[:, 0]

        discrepant_snps_passed = pd.Index([], name="rsid")
        counter = 0
        # save matches for this chromosome
        for x in matches_passed:
            d = {
                "chrom": chrom,
                "start": df.loc[(df["chrom"] == chrom)].iloc[x[0]].pos,
                "end": df.loc[(df["chrom"] == chrom)].iloc[x[1] - 1].pos,
                "cMs": cMs_match_segment[counter],
                "snps": x[1] - x[0],
            }

            # ensure discrepant SNPs are in shared DNA segments
            for discrepant_snp in discrepant_snps:
                if d["start"] <= df.loc[discrepant_snp].pos <= d["end"]:
                    discrepant_snps_passed = discrepant_snps_passed.append(
                        df.loc[[discrepant_snp]].index
                    )

            # remove found discrepant SNPs from search on next iteration
            discrepant_snps = discrepant_snps.drop(
                discrepant_snps_passed, errors="ignore"
            )

            shared_dna.append(d)
            counter += 1
        return {"shared_dna": shared_dna, "discrepant_snps": discrepant_snps_passed}

    def _remap_snps_to_GRCh37(self, individuals):
        for ind in individuals:
            if ind is None:
                continue

            ind.remap(37)

    def _get_csv_header(self):
        return (
            f"# Generated by lineage v{__version__}; https://pypi.org/project/lineage/{os.linesep}"
            f"# Generated at {datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC{os.linesep}"
        )
