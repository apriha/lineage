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

import os
import re

import numpy as np
import pandas as pd

import lineage
from lineage.snps import (
    SNPs,
    get_assembly,
    get_chromosomes,
    get_chromosomes_summary,
    get_snp_count,
    sort_snps,
    determine_sex,
)


class Individual(object):
    """ Object used to represent and interact with an individual.

    The ``Individual`` object maintains information about an individual. The object provides
    methods for loading an individual's genetic data (SNPs) and normalizing it for use with the
    `lineage` framework.

    """

    def __init__(self, name, raw_data=None, output_dir="output"):
        """ Initialize an ``Individual`` object.

        Parameters
        ----------
        name : str
            name of the individual
        raw_data : list or str
            path(s) to file(s) with raw genotype data
        output_dir : str
            path to output directory
        """
        self._name = name
        self._output_dir = output_dir
        self._snps = None
        self._build = None
        self._source = []
        self._discrepant_positions_file_count = 0
        self._discrepant_genotypes_file_count = 0
        self._discrepant_positions = pd.DataFrame()
        self._discrepant_genotypes = pd.DataFrame()

        if raw_data is not None:
            self.load_snps(raw_data)

    def __repr__(self):
        return "Individual({!r})".format(self._name)

    @property
    def name(self):
        """ Get this ``Individual``'s name.

        Returns
        -------
        str
        """
        return self._name

    @property
    def snps(self):
        """ Get a copy of this ``Individual``'s SNPs.

        Returns
        -------
        pandas.DataFrame
        """
        if self._snps is not None:
            return self._snps.copy()
        else:
            return None

    @property
    def snp_count(self):
        """ Count of SNPs loaded for this ``Individual``.

        Returns
        -------
        int
        """
        return get_snp_count(self._snps)

    @property
    def chromosomes(self):
        """ Get the chromosomes of this ``Individual``'s SNPs.

        Returns
        -------
        list
            list of str chromosomes (e.g., ['1', '2', '3', 'MT'], empty list if no chromosomes
        """
        return get_chromosomes(self._snps)

    @property
    def chromosomes_summary(self):
        """ Summary of the chromosomes represented by this ``Individual``'s SNPs.

        Returns
        -------
        str
            human-readable listing of chromosomes (e.g., '1-3, MT'), empty str if no chromosomes
        """
        return get_chromosomes_summary(self._snps)

    @property
    def build(self):
        """ Get the build of this ``Individual``'s SNPs.

        Returns
        -------
        int
        """
        return self._build

    @property
    def assembly(self):
        """ Get the assembly of this ``Individual``'s SNPs.

        Returns
        -------
        str
        """
        return get_assembly(self._build)

    @property
    def source(self):
        """ Summary of the SNP data source(s) for this ``Individual``'s SNPs.

        Returns
        -------
        str
        """
        return ", ".join(self._source)

    @property
    def sex(self):
        """ Sex of this ``Individual`` derived from this ``Individual``'s SNPs.

        Returns
        -------
        str
            'Male' or 'Female' if detected, else empty str
        """
        return determine_sex(self._snps)

    @property
    def discrepant_positions(self):
        """ SNPs with discrepant positions discovered while loading SNPs.

        Returns
        -------
        pandas.DataFrame
        """
        return self._discrepant_positions

    @property
    def discrepant_genotypes(self):
        """ SNPs with discrepant genotypes discovered while loading SNPs.

        Returns
        -------
        pandas.DataFrame
        """
        return self._discrepant_genotypes

    @property
    def discrepant_snps(self):
        """ SNPs with discrepant positions and / or genotypes discovered while loading SNPs.

        Returns
        -------
        pandas.DataFrame
        """
        df = self._discrepant_positions.append(self._discrepant_genotypes)
        if len(df) > 1:
            df = df.drop_duplicates()
        return df

    def load_snps(
        self,
        raw_data,
        discrepant_snp_positions_threshold=100,
        discrepant_genotypes_threshold=500,
        save_output=False,
    ):
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
        save_output : bool
            specifies whether to save discrepant SNP output to CSV files in the output directory
        """
        if type(raw_data) is list:
            for file in raw_data:
                self._load_snps_helper(
                    file,
                    discrepant_snp_positions_threshold,
                    discrepant_genotypes_threshold,
                    save_output,
                )
        elif type(raw_data) is str:
            self._load_snps_helper(
                raw_data,
                discrepant_snp_positions_threshold,
                discrepant_genotypes_threshold,
                save_output,
            )
        else:
            raise TypeError("invalid filetype")

    def _load_snps_helper(
        self,
        file,
        discrepant_snp_positions_threshold,
        discrepant_genotypes_threshold,
        save_output,
    ):
        print("Loading " + os.path.relpath(file))
        discrepant_positions, discrepant_genotypes = self._add_snps(
            SNPs(file),
            discrepant_snp_positions_threshold,
            discrepant_genotypes_threshold,
            save_output,
        )

        self._discrepant_positions = self._discrepant_positions.append(
            discrepant_positions, sort=True
        )
        self._discrepant_genotypes = self._discrepant_genotypes.append(
            discrepant_genotypes, sort=True
        )

    def save_snps(self, filename=None):
        """ Save SNPs to file.

        Parameters
        ----------
        filename : str
            filename for file to save

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        comment = (
            "# Source(s): {}\n"
            "# Assembly: {}\n"
            "# SNPs: {}\n"
            "# Chromosomes: {}\n".format(
                self.source, self.assembly, self.snp_count, self.chromosomes_summary
            )
        )

        if filename is None:
            filename = self.get_var_name() + "_lineage_" + self.assembly + ".csv"

        return lineage.save_df_as_csv(
            self._snps,
            self._output_dir,
            filename,
            comment=comment,
            header=["chromosome", "position", "genotype"],
        )

    def save_discrepant_positions(self, filename=None):
        """ Save SNPs with discrepant positions to file.

        Parameters
        ----------
        filename : str
            filename for file to save

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        return self._save_discrepant_snps_file(
            self.discrepant_positions, "discrepant_positions", filename
        )

    def save_discrepant_genotypes(self, filename=None):
        """ Save SNPs with discrepant genotypes to file.

        Parameters
        ----------
        filename : str
            filename for file to save

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        return self._save_discrepant_snps_file(
            self.discrepant_genotypes, "discrepant_genotypes", filename
        )

    def save_discrepant_snps(self, filename=None):
        """ Save SNPs with discrepant positions and / or genotypes to file.

        Parameters
        ----------
        filename : str
            filename for file to save

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        return self._save_discrepant_snps_file(
            self.discrepant_snps, "discrepant_snps", filename
        )

    def _save_discrepant_snps_file(self, df, name, filename):
        if filename is None:
            filename = self.get_var_name() + "_" + name + ".csv"

        return lineage.save_df_as_csv(
            df,
            self._output_dir,
            filename,
            comment="# Source(s): {}\n".format(self.source),
        )

    def get_var_name(self):
        return self._clean_string(self.name)

    def _clean_string(self, s):
        """ Clean a string so that it can be a valid Python variable
        name.

        Parameters
        ----------
        s : str
            string to clean

        Returns
        -------
        str
            cleaned string that can be used as a variable name
        """
        # http://stackoverflow.com/a/3305731
        return re.sub("\W|^(?=\d)", "_", s)

    def remap_snps(self, target_assembly, complement_bases=True):
        """ Remap the SNP coordinates of this ``Individual`` from one assembly to another.

        This method is a wrapper for `remap_snps` in the ``Lineage`` class.

        This method uses the assembly map endpoint of the Ensembl REST API service to convert SNP
        coordinates / positions from one assembly to another. After remapping, the coordinates /
        positions for the ``Individual``'s SNPs will be that of the target assembly.

        If the SNPs are already mapped relative to the target assembly, remapping will not be
        performed.

        Parameters
        ----------
        target_assembly : {'NCBI36', 'GRCh37', 'GRCh38', 36, 37, 38}
            assembly to remap to
        complement_bases : bool
            complement bases when remapping SNPs to the minus strand

        Returns
        -------
        chromosomes_remapped : list of str
            chromosomes remapped; empty if None
        chromosomes_not_remapped : list of str
            chromosomes not remapped; empty if None

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
        from lineage import Lineage

        l = Lineage()
        return l.remap_snps(self, target_assembly, complement_bases)

    def _set_snps(self, snps, build=37):
        """ Set `_snps` and `_build` properties of this ``Individual``.

        Notes
        -----
        Intended to be used internally to `lineage`.

        Parameters
        ----------
        snps : pandas.DataFrame
            individual's genetic data normalized for use with `lineage`
        build : int
            build of this ``Individual``'s SNPs
        """
        self._snps = snps
        self._build = build

    def _add_snps(
        self,
        snps,
        discrepant_snp_positions_threshold,
        discrepant_genotypes_threshold,
        save_output,
    ):
        """ Add SNPs to this Individual.

        Parameters
        ----------
        snps : SNPs
            SNPs to add
        discrepant_snp_positions_threshold : int
            see above
        discrepant_genotypes_threshold : int
            see above
        save_output
            see above

        Returns
        -------
        discrepant_positions : pandas.DataFrame
        discrepant_genotypes : pandas.DataFrame
        """
        discrepant_positions = pd.DataFrame()
        discrepant_genotypes = pd.DataFrame()

        if snps.snps is None:
            return discrepant_positions, discrepant_genotypes

        build = snps.build
        source = [s.strip() for s in snps.source.split(",")]

        if not snps.build_detected:
            print("build not detected, assuming build {}".format(snps.build))

        if self._build is None:
            self._build = build
        elif self._build != build:
            print(
                "build / assembly mismatch between current build of SNPs and SNPs being loaded"
            )

        # ensure there area always two X alleles
        snps = self._double_single_alleles(snps.snps, "X")

        if self._snps is None:
            self._source.extend(source)
            self._snps = snps
        else:
            common_snps = self._snps.join(snps, how="inner", rsuffix="_added")

            discrepant_positions = common_snps.loc[
                (common_snps["chrom"] != common_snps["chrom_added"])
                | (common_snps["pos"] != common_snps["pos_added"])
            ]

            if 0 < len(discrepant_positions) < discrepant_snp_positions_threshold:
                print(
                    str(len(discrepant_positions)) + " SNP positions were discrepant; "
                    "keeping original positions"
                )

                if save_output:
                    self._discrepant_positions_file_count += 1
                    lineage.save_df_as_csv(
                        discrepant_positions,
                        self._output_dir,
                        self.get_var_name()
                        + "_discrepant_positions_"
                        + str(self._discrepant_positions_file_count)
                        + ".csv",
                    )
            elif len(discrepant_positions) >= discrepant_snp_positions_threshold:
                print(
                    "too many SNPs differ in position; ensure same genome build is being used"
                )
                return discrepant_positions, discrepant_genotypes

            # remove null genotypes
            common_snps = common_snps.loc[
                ~common_snps["genotype"].isnull()
                & ~common_snps["genotype_added"].isnull()
            ]

            # discrepant genotypes are where alleles are not equivalent (i.e., alleles are not the
            # same and not swapped)
            discrepant_genotypes = common_snps.loc[
                (
                    (common_snps["genotype"].str.len() == 1)
                    & (common_snps["genotype_added"].str.len() == 1)
                    & ~(
                        common_snps["genotype"].str[0]
                        == common_snps["genotype_added"].str[0]
                    )
                )
                | (
                    (common_snps["genotype"].str.len() == 2)
                    & (common_snps["genotype_added"].str.len() == 2)
                    & ~(
                        (
                            common_snps["genotype"].str[0]
                            == common_snps["genotype_added"].str[0]
                        )
                        & (
                            common_snps["genotype"].str[1]
                            == common_snps["genotype_added"].str[1]
                        )
                    )
                    & ~(
                        (
                            common_snps["genotype"].str[0]
                            == common_snps["genotype_added"].str[1]
                        )
                        & (
                            common_snps["genotype"].str[1]
                            == common_snps["genotype_added"].str[0]
                        )
                    )
                )
            ]

            if 0 < len(discrepant_genotypes) < discrepant_genotypes_threshold:
                print(
                    str(len(discrepant_genotypes)) + " SNP genotypes were discrepant; "
                    "marking those as null"
                )

                if save_output:
                    self._discrepant_genotypes_file_count += 1
                    lineage.save_df_as_csv(
                        discrepant_genotypes,
                        self._output_dir,
                        self.get_var_name()
                        + "_discrepant_genotypes_"
                        + str(self._discrepant_genotypes_file_count)
                        + ".csv",
                    )
            elif len(discrepant_genotypes) >= discrepant_genotypes_threshold:
                print(
                    "too many SNPs differ in their genotype; ensure file is for same "
                    "individual"
                )
                return discrepant_positions, discrepant_genotypes

            # add new SNPs
            self._source.extend(source)
            self._snps = self._snps.combine_first(snps)
            self._snps.loc[discrepant_genotypes.index, "genotype"] = np.nan

            # combine_first converts position to float64, so convert it back to int64
            self._snps["pos"] = self._snps["pos"].astype(np.int64)

        self._snps = sort_snps(self._snps)

        return discrepant_positions, discrepant_genotypes

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
        single_alleles = np.where(
            (df["chrom"] == chrom) & (df["genotype"].str.len() == 1)
        )[0]

        # double those alleles
        df.ix[single_alleles, "genotype"] = df.ix[single_alleles, "genotype"] * 2

        return df
