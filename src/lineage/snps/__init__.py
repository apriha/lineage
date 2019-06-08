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
import os
import re

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

from lineage.ensembl import EnsemblRestClient
from lineage.resources import Resources
from lineage.snps.io import Reader, Writer
from lineage.utils import save_df_as_csv, Parallelizer, clean_str


class SNPs:
    def __init__(
        self,
        file="",
        assign_par_snps=True,
        output_dir="output",
        resources_dir="resources",
        parallelize=False,
        processes=os.cpu_count(),
    ):
        """ Object used to read and parse genotype / raw data files.

        Parameters
        ----------
        file : str
            path to file to load
        assign_par_snps : bool
            assign PAR SNPs to the X and Y chromosomes
        output_dir : str
            path to output directory
        resources_dir : str
            name / path of resources directory
        parallelize : bool
            utilize multiprocessing to speedup calculations
        processes : int
            processes to launch if multiprocessing
        """
        self._snps = pd.DataFrame()
        self._source = ""
        self._build = 0
        self._build_detected = False
        self._output_dir = output_dir
        self._resources = Resources(resources_dir=resources_dir)
        self._parallelizer = Parallelizer(parallelize=parallelize, processes=processes)

        if file:
            self._snps, self._source = self._read_raw_data(file)

            if not self._snps.empty:
                self.sort_snps()

                self._build = self.detect_build()

                if not self._build:
                    self._build = 37  # assume Build 37 / GRCh37 if not detected
                else:
                    self._build_detected = True

                if assign_par_snps:
                    self._assign_par_snps()

    @property
    def source(self):
        """ Summary of the SNP data source for ``SNPs``.

        Returns
        -------
        str
        """
        return self._source

    @property
    def snps(self):
        """ Get a copy of SNPs.

        Returns
        -------
        pandas.DataFrame
        """
        return self._snps.copy()

    @property
    def build(self):
        """ Get the build of ``SNPs``.

        Returns
        -------
        int
        """
        return self._build

    @property
    def build_detected(self):
        """ Get status indicating if build of ``SNPs`` was detected.

        Returns
        -------
        bool
        """
        return self._build_detected

    @property
    def assembly(self):
        """ Get the assembly of ``SNPs``.

        Returns
        -------
        str
        """
        return self.get_assembly()

    @property
    def snp_count(self):
        """ Count of SNPs.

        Returns
        -------
        int
        """
        return self.get_snp_count()

    @property
    def chromosomes(self):
        """ Chromosomes of ``SNPs``.

        Returns
        -------
        list
            list of str chromosomes (e.g., ['1', '2', '3', 'MT'], empty list if no chromosomes
        """
        return self.get_chromosomes()

    @property
    def chromosomes_summary(self):
        """ Summary of the chromosomes of ``SNPs``.

        Returns
        -------
        str
            human-readable listing of chromosomes (e.g., '1-3, MT'), empty str if no chromosomes
        """
        return self.get_chromosomes_summary()

    @property
    def sex(self):
        """ Sex derived from ``SNPs``.

        Returns
        -------
        str
            'Male' or 'Female' if detected, else empty str
        """
        return self.determine_sex()

    def get_summary(self):
        """ Get summary of ``SNPs``.

        Returns
        -------
        dict
            summary info if ``SNPs`` is valid, else {}
        """
        if not self.is_valid():
            return {}
        else:
            return {
                "source": self.source,
                "assembly": self.assembly,
                "build": self.build,
                "build_detected": self.build_detected,
                "snp_count": self.snp_count,
                "chromosomes": self.chromosomes_summary,
                "sex": self.sex,
            }

    def is_valid(self):
        """ Determine if ``SNPs`` is valid.

        ``SNPs`` is valid when the input file has been successfully parsed.

        Returns
        -------
        bool
            True if ``SNPs`` is valid
        """
        if self._snps.empty:
            return False
        else:
            return True

    def save_snps(self, filename=""):
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
        return Writer.write_file(snps=self, filename=filename)

    def _read_raw_data(self, file):
        return Reader.read_file(file)

    def _assign_par_snps(self):
        """ Assign PAR SNPs to the X or Y chromosome using SNP position.

        References
        -----
        ..[1] National Center for Biotechnology Information, Variation Services, RefSNP,
          https://api.ncbi.nlm.nih.gov/variation/v0/
        ..[2] Yates et. al. (doi:10.1093/bioinformatics/btu613),
          http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613
        ..[3] Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098
        ..[4] Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K.
          dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;
          29(1):308-11.
        ..[5] Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center
          for Biotechnology Information, National Library of Medicine. dbSNP accession: rs28736870,
          rs113313554, and rs758419898 (dbSNP Build ID: 151). Available from:
          http://www.ncbi.nlm.nih.gov/SNP/
        """
        rest_client = EnsemblRestClient(
            server="https://api.ncbi.nlm.nih.gov", reqs_per_sec=1
        )
        for rsid in self._snps.loc[self._snps["chrom"] == "PAR"].index.values:
            if "rs" in rsid:
                try:
                    id = rsid.split("rs")[1]
                    response = rest_client.perform_rest_action(
                        "/variation/v0/refsnp/" + id
                    )

                    if response is not None:
                        for item in response["primary_snapshot_data"][
                            "placements_with_allele"
                        ]:
                            if "NC_000023" in item["seq_id"]:
                                assigned = self._assign_snp(rsid, item["alleles"], "X")
                            elif "NC_000024" in item["seq_id"]:
                                assigned = self._assign_snp(rsid, item["alleles"], "Y")
                            else:
                                assigned = False

                            if assigned:
                                if not self._build_detected:
                                    self._build = self._extract_build(item)
                                    self._build_detected = True
                                break

                except Exception as err:
                    print(err)

    def _assign_snp(self, rsid, alleles, chrom):
        # only assign SNP if positions match (i.e., same build)
        for allele in alleles:
            allele_pos = allele["allele"]["spdi"]["position"]
            # ref SNP positions seem to be 0-based...
            if allele_pos == self._snps.loc[rsid].pos - 1:
                self._snps.loc[rsid, "chrom"] = chrom
                return True
        return False

    def _extract_build(self, item):
        assembly_name = item["placement_annot"]["seq_id_traits_by_assembly"][0][
            "assembly_name"
        ]
        assembly_name = assembly_name.split(".")[0]
        return int(assembly_name[-2:])

    def detect_build(self):
        """ Detect build of SNPs.

        Use the coordinates of common SNPs to identify the build / assembly of a genotype file
        that is being loaded.

        Notes
        -----
        rs3094315 : plus strand in 36, 37, and 38
        rs11928389 : plus strand in 36, minus strand in 37 and 38
        rs2500347 : plus strand in 36 and 37, minus strand in 38
        rs964481 : plus strand in 36, 37, and 38
        rs2341354 : plus strand in 36, 37, and 38

        Returns
        -------
        int
            detected build of SNPs, else 0

        References
        ----------
        ..[1] Yates et. al. (doi:10.1093/bioinformatics/btu613),
          http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613
        ..[2] Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098
        ..[3] Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K.
          dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;29(1):308-11.
        ..[4] Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center
          for Biotechnology Information, National Library of Medicine. dbSNP accession: rs3094315,
          rs11928389, rs2500347, rs964481, and rs2341354 (dbSNP Build ID: 151). Available from:
          http://www.ncbi.nlm.nih.gov/SNP/
        """

        def lookup_build_with_snp_pos(pos, s):
            try:
                return s.loc[s == pos].index[0]
            except:
                return 0

        build = 0

        rsids = ["rs3094315", "rs11928389", "rs2500347", "rs964481", "rs2341354"]
        df = pd.DataFrame(
            {
                36: [742429, 50908372, 143649677, 27566744, 908436],
                37: [752566, 50927009, 144938320, 27656823, 918573],
                38: [817186, 50889578, 148946169, 27638706, 983193],
            },
            index=rsids,
        )

        for rsid in rsids:
            if rsid in self._snps.index:
                build = lookup_build_with_snp_pos(
                    self._snps.loc[rsid].pos, df.loc[rsid]
                )

            if build:
                break

        return build

    def get_assembly(self):
        """ Get the assembly of a build.

        Returns
        -------
        str
        """

        if self._build == 37:
            return "GRCh37"
        elif self._build == 36:
            return "NCBI36"
        elif self._build == 38:
            return "GRCh38"
        else:
            return ""

    def get_snp_count(self):
        """ Count of SNPs.

        Returns
        -------
        int
        """

        return len(self._snps)

    def get_chromosomes(self):
        """ Get the chromosomes of SNPs.

        Returns
        -------
        list
            list of str chromosomes (e.g., ['1', '2', '3', 'MT'], empty list if no chromosomes
        """

        if not self._snps.empty:
            return list(pd.unique(self._snps["chrom"]))
        else:
            return []

    def get_chromosomes_summary(self):
        """ Summary of the chromosomes of SNPs.

        Returns
        -------
        str
            human-readable listing of chromosomes (e.g., '1-3, MT'), empty str if no chromosomes
        """

        if not self._snps.empty:
            chroms = list(pd.unique(self._snps["chrom"]))

            int_chroms = [int(chrom) for chrom in chroms if chrom.isdigit()]
            str_chroms = [chrom for chrom in chroms if not chrom.isdigit()]

            # https://codereview.stackexchange.com/a/5202
            def as_range(iterable):
                l = list(iterable)
                if len(l) > 1:
                    return "{0}-{1}".format(l[0], l[-1])
                else:
                    return "{0}".format(l[0])

            # create str representations
            int_chroms = ", ".join(
                as_range(g)
                for _, g in groupby(int_chroms, key=lambda n, c=count(): n - next(c))
            )
            str_chroms = ", ".join(str_chroms)

            if int_chroms != "" and str_chroms != "":
                int_chroms += ", "

            return int_chroms + str_chroms
        else:
            return ""

    def determine_sex(
        self, y_snps_not_null_threshold=0.1, heterozygous_x_snps_threshold=0.01
    ):
        """ Determine sex from SNPs using thresholds.

        Parameters
        ----------
        y_snps_not_null_threshold : float
            percentage Y SNPs that are not null; above this threshold, Male is determined
        heterozygous_x_snps_threshold : float
            percentage heterozygous X SNPs; above this threshold, Female is determined

        Returns
        -------
        str
            'Male' or 'Female' if detected, else empty str
        """

        if not self._snps.empty:
            y_snps = len(self._snps.loc[(self._snps["chrom"] == "Y")])

            if y_snps > 0:
                y_snps_not_null = len(
                    self._snps.loc[
                        (self._snps["chrom"] == "Y")
                        & (self._snps["genotype"].notnull())
                    ]
                )

                if y_snps_not_null / y_snps > y_snps_not_null_threshold:
                    return "Male"
                else:
                    return "Female"

            x_snps = len(self._snps.loc[self._snps["chrom"] == "X"])

            if x_snps == 0:
                return ""

            heterozygous_x_snps = len(
                self._snps.loc[
                    (self._snps["chrom"] == "X")
                    & (self._snps["genotype"].notnull())
                    & (self._snps["genotype"].str[0] != self._snps["genotype"].str[1])
                ]
            )

            if heterozygous_x_snps / x_snps > heterozygous_x_snps_threshold:
                return "Female"
            else:
                return "Male"
        else:
            return ""

    def sort_snps(self):
        """ Sort SNPs based on ordered chromosome list and position. """

        sorted_list = sorted(self._snps["chrom"].unique(), key=self._natural_sort_key)

        # move PAR and MT to the end of the dataframe
        if "PAR" in sorted_list:
            sorted_list.remove("PAR")
            sorted_list.append("PAR")

        if "MT" in sorted_list:
            sorted_list.remove("MT")
            sorted_list.append("MT")

        # convert chrom column to category for sorting
        # https://stackoverflow.com/a/26707444
        self._snps["chrom"] = self._snps["chrom"].astype(
            CategoricalDtype(categories=sorted_list, ordered=True)
        )

        # sort based on ordered chromosome list and position
        snps = self._snps.sort_values(["chrom", "pos"])

        # convert chromosome back to object
        snps["chrom"] = snps["chrom"].astype(object)

        self._snps = snps

    def remap_snps(self, target_assembly, complement_bases=True):
        """ Remap SNP coordinates from one assembly to another.

        This method uses the assembly map endpoint of the Ensembl REST API service (via
        ``Resources``'s ``EnsemblRestClient``) to convert SNP coordinates / positions from one
        assembly to another. After remapping, the coordinates / positions for the
        SNPs will be that of the target assembly.

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
            chromosomes remapped
        chromosomes_not_remapped : list of str
            chromosomes not remapped

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
        chromosomes_remapped = []
        chromosomes_not_remapped = []

        snps = self.snps

        if snps.empty:
            print("No SNPs to remap")
            return chromosomes_remapped, chromosomes_not_remapped
        else:
            chromosomes = snps["chrom"].unique()
            chromosomes_not_remapped = list(chromosomes)

        valid_assemblies = ["NCBI36", "GRCh37", "GRCh38", 36, 37, 38]

        if target_assembly not in valid_assemblies:
            print("Invalid target assembly")
            return chromosomes_remapped, chromosomes_not_remapped

        if isinstance(target_assembly, int):
            if target_assembly == 36:
                target_assembly = "NCBI36"
            else:
                target_assembly = "GRCh" + str(target_assembly)

        if self.build == 36:
            source_assembly = "NCBI36"
        else:
            source_assembly = "GRCh" + str(self.build)

        if source_assembly == target_assembly:
            return chromosomes_remapped, chromosomes_not_remapped

        assembly_mapping_data = self._resources.get_assembly_mapping_data(
            source_assembly, target_assembly
        )

        if not assembly_mapping_data:
            return chromosomes_remapped, chromosomes_not_remapped

        tasks = []

        for chrom in chromosomes:
            if chrom in assembly_mapping_data:
                chromosomes_remapped.append(chrom)
                chromosomes_not_remapped.remove(chrom)
                mappings = assembly_mapping_data[chrom]
                tasks.append(
                    {
                        "snps": snps.loc[snps["chrom"] == chrom],
                        "mappings": mappings,
                        "complement_bases": complement_bases,
                    }
                )
            else:
                print(
                    "Chromosome {} not remapped; "
                    "removing chromosome from SNPs for consistency".format(chrom)
                )
                snps = snps.drop(snps.loc[snps["chrom"] == chrom].index)

        # remap SNPs
        remapped_snps = self._parallelizer(self._remapper, tasks)
        remapped_snps = pd.concat(remapped_snps)

        # update SNP positions and genotypes
        snps.loc[remapped_snps.index, "pos"] = remapped_snps["pos"]
        snps.loc[remapped_snps.index, "genotype"] = remapped_snps["genotype"]

        self._snps = snps
        self.sort_snps()
        self._build = int(target_assembly[-2:])

        return chromosomes_remapped, chromosomes_not_remapped

    def _remapper(self, task):
        """ Remap SNPs for a chromosome.

        Parameters
        ----------
        task : dict
            dict with `snps` to remap per `mappings`, optionally `complement_bases`

        Returns
        -------
        pandas.DataFrame
            remapped SNPs
        """
        temp = task["snps"].copy()
        mappings = task["mappings"]
        complement_bases = task["complement_bases"]

        temp["remapped"] = False

        pos_start = int(temp["pos"].describe()["min"])
        pos_end = int(temp["pos"].describe()["max"])

        for mapping in mappings["mappings"]:
            # skip if mapping is outside of range of SNP positions
            if (
                mapping["original"]["end"] <= pos_start
                or mapping["original"]["start"] >= pos_end
            ):
                continue

            orig_range_len = mapping["original"]["end"] - mapping["original"]["start"]
            mapped_range_len = mapping["mapped"]["end"] - mapping["mapped"]["start"]

            orig_region = mapping["original"]["seq_region_name"]
            mapped_region = mapping["mapped"]["seq_region_name"]

            if orig_region != mapped_region:
                print("discrepant chroms")
                continue

            if orig_range_len != mapped_range_len:
                print("discrepant coords")  # observed when mapping NCBI36 -> GRCh38
                continue

            # find the SNPs that are being remapped for this mapping
            snp_indices = temp.loc[
                ~temp["remapped"]
                & (temp["pos"] >= mapping["original"]["start"])
                & (temp["pos"] <= mapping["original"]["end"])
            ].index

            if len(snp_indices) > 0:
                # remap the SNPs
                if mapping["mapped"]["strand"] == -1:
                    # flip and (optionally) complement since we're mapping to minus strand
                    diff_from_start = (
                        temp.loc[snp_indices, "pos"] - mapping["original"]["start"]
                    )
                    temp.loc[snp_indices, "pos"] = (
                        mapping["mapped"]["end"] - diff_from_start
                    )

                    if complement_bases:
                        temp.loc[snp_indices, "genotype"] = temp.loc[
                            snp_indices, "genotype"
                        ].apply(self._complement_bases)
                else:
                    # mapping is on same (plus) strand, so just remap based on offset
                    offset = mapping["mapped"]["start"] - mapping["original"]["start"]
                    temp.loc[snp_indices, "pos"] = temp["pos"] + offset

                # mark these SNPs as remapped
                temp.loc[snp_indices, "remapped"] = True

        return temp

    def _complement_bases(self, genotype):
        if pd.isnull(genotype):
            return np.nan

        complement = ""

        for base in list(genotype):
            if base == "A":
                complement += "T"
            elif base == "G":
                complement += "C"
            elif base == "C":
                complement += "G"
            elif base == "T":
                complement += "A"
            else:
                complement += base

        return complement

    # https://stackoverflow.com/a/16090640
    @staticmethod
    def _natural_sort_key(s, natural_sort_re=re.compile("([0-9]+)")):
        return [
            int(text) if text.isdigit() else text.lower()
            for text in re.split(natural_sort_re, s)
        ]


class SNPsCollection(SNPs):
    def __init__(self, raw_data=None, output_dir="output", name="", **kwargs):
        """

        Parameters
        ----------
        raw_data : list or str
            path(s) to file(s) with raw genotype data
        output_dir : str
            path to output directory
        name : str
            name for this ``SNPsCollection``
        """
        super().__init__(file="", output_dir=output_dir, **kwargs)

        self._source = []
        self._discrepant_positions_file_count = 0
        self._discrepant_genotypes_file_count = 0
        self._discrepant_positions = pd.DataFrame()
        self._discrepant_genotypes = pd.DataFrame()
        self._name = name

        if raw_data is not None:
            self.load_snps(raw_data)

    @property
    def source(self):
        """ Summary of the SNP data source for ``SNPs``.

        Returns
        -------
        str
        """
        return ", ".join(self._source)

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

    def save_snps(self, filename=""):
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
        if not filename:
            filename = "{}_lineage_{}{}".format(
                clean_str(self._name), self.assembly, ".csv"
            )
        return super().save_snps(filename)

    def save_discrepant_positions(self, filename=""):
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

    def save_discrepant_genotypes(self, filename=""):
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

    def save_discrepant_snps(self, filename=""):
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
        if not filename:
            filename = clean_str(self._name) + "_" + name + ".csv"

        return save_df_as_csv(
            df,
            self._output_dir,
            filename,
            comment="# Source(s): {}\n".format(self.source),
        )

    def _add_snps(
        self,
        snps,
        discrepant_snp_positions_threshold,
        discrepant_genotypes_threshold,
        save_output,
    ):
        """ Add SNPs to this ``SNPsCollection``.

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

        if snps._snps.empty:
            return discrepant_positions, discrepant_genotypes

        build = snps._build
        source = [s.strip() for s in snps._source.split(",")]

        if not snps._build_detected:
            print("build not detected, assuming build {}".format(snps._build))

        if not self._build:
            self._build = build
        elif self._build != build:
            print(
                "build / assembly mismatch between current build of SNPs and SNPs being loaded"
            )

        # ensure there area always two X alleles
        snps = self._double_single_alleles(snps._snps, "X")

        if self._snps.empty:
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
                    save_df_as_csv(
                        discrepant_positions,
                        self._output_dir,
                        clean_str(self._name)
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
                    save_df_as_csv(
                        discrepant_genotypes,
                        self._output_dir,
                        clean_str(self._name)
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

        self.sort_snps()

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
        df.iloc[single_alleles, 2] = df.iloc[single_alleles, 2] * 2

        return df
