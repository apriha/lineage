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
import re
import zipfile

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

import lineage
from lineage.ensembl import EnsemblRestClient


class SNPs:
    def __init__(self, file="", assign_par_snps=True, output_dir="output"):
        """ Object used to read and parse genotype / raw data files.

        Parameters
        ----------
        file : str
            path to file to load
        assign_par_snps : bool
            assign PAR SNPs to the X and Y chromosomes
        output_dir : str
            path to output directory
        """
        self._snps = None
        self._source = ""
        self._build = None
        self._build_detected = False
        self._output_dir = output_dir

        if file:
            self._snps, self._source = self._read_raw_data(file)

            if self._snps is not None:
                self.sort_snps()

                self._build = self.detect_build()

                if self._build is None:
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
        if self._snps is not None:
            return self._snps.copy()
        else:
            return None

    @property
    def build(self):
        """ Get the build of ``SNPs``.

        Returns
        -------
        int
        """
        return self._build

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
        return self.determine_sex(self._snps)

    def get_summary(self):
        """ Get summary of ``SNPs``.

        Returns
        -------
        dict
            summary info, else None if ``SNPs`` is not valid
        """
        if not self.is_valid():
            return None
        else:
            return {
                "source": self._source,
                "assembly": self.assembly,
                "build": self._build,
                "build_detected": self._build_detected,
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
        if self._snps is None:
            return False
        else:
            return True

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
            filename = (
                self.get_var_repr(self._source) + "_lineage_" + self.assembly + ".csv"
            )

        return lineage.save_df_as_csv(
            self._snps,
            self._output_dir,
            filename,
            comment=comment,
            header=["chromosome", "position", "genotype"],
        )

    def get_var_repr(self, s):
        return self._clean_string(s)

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

    def _read_raw_data(self, file):
        try:
            if not os.path.exists(file):
                print(file + " does not exist; skipping")
                return None, ""

            # peek into files to determine the data format
            if ".zip" in file:
                with zipfile.ZipFile(file) as z:
                    with z.open(z.namelist()[0], "r") as f:
                        first_line, comments = self._extract_comments(f, True)
            elif ".gz" in file:
                with gzip.open(file, "rt") as f:
                    first_line, comments = self._extract_comments(f, False)
            else:
                with open(file, "r") as f:
                    first_line, comments = self._extract_comments(f, False)

            if "23andMe" in first_line:
                return self._read_23andme(file)
            elif "Ancestry" in first_line:
                return self._read_ancestry(file)
            elif first_line.startswith("RSID"):
                return self._read_ftdna(file)
            elif "famfinder" in first_line:
                return self._read_ftdna_famfinder(file)
            elif "lineage" in first_line:
                return self._read_lineage_csv(file, comments)
            elif first_line.startswith("rsid"):
                return self._read_generic_csv(file)
            else:
                return None, ""
        except Exception as err:
            print(err)
            return None, ""

    def _extract_comments(self, f, decode):
        line = self._read_line(f, decode)
        first_line = line
        comments = ""

        while line.startswith("#"):
            comments += line
            line = self._read_line(f, decode)

        return first_line, comments

    def _read_line(self, f, decode):
        if decode:
            # https://stackoverflow.com/a/606199
            return f.readline().decode("utf-8")
        else:
            return f.readline()

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
        str
            name of data source
        """
        df = pd.read_csv(
            file,
            comment="#",
            sep="\t",
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype={"chrom": object},
        )

        return df, "23andMe"

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
        str
            name of data source
        """
        df = pd.read_csv(
            file,
            skiprows=1,
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype={"chrom": object},
        )

        # remove incongruous data
        df = df.drop(df.loc[df["chrom"] == "0"].index)
        df = df.drop(
            df.loc[df.index == "RSID"].index
        )  # second header for concatenated data

        # if second header existed, pos dtype will be object (should be np.int64)
        df["pos"] = df["pos"].astype(np.int64)

        return df, "FTDNA"

    @staticmethod
    def _read_ftdna_famfinder(file):
        """ Read and parse Family Tree DNA (FTDNA) "famfinder" file.

        https://www.familytreedna.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            individual's genetic data normalized for use with `lineage`
        str
            name of data source
        """
        df = pd.read_csv(
            file,
            comment="#",
            na_values="-",
            names=["rsid", "chrom", "pos", "allele1", "allele2"],
            index_col=0,
            dtype={"chrom": object},
        )

        # create genotype column from allele columns
        df["genotype"] = df["allele1"] + df["allele2"]

        # delete allele columns
        # http://stackoverflow.com/a/13485766
        del df["allele1"]
        del df["allele2"]

        return df, "FTDNA"

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
        str
            name of data source
        """
        df = pd.read_csv(
            file,
            comment="#",
            header=0,
            sep="\t",
            na_values=0,
            names=["rsid", "chrom", "pos", "allele1", "allele2"],
            index_col=0,
            dtype={"chrom": object},
        )

        # create genotype column from allele columns
        df["genotype"] = df["allele1"] + df["allele2"]

        # delete allele columns
        # http://stackoverflow.com/a/13485766
        del df["allele1"]
        del df["allele2"]

        # https://redd.it/5y90un
        df.ix[np.where(df["chrom"] == "23")[0], "chrom"] = "X"
        df.ix[np.where(df["chrom"] == "24")[0], "chrom"] = "Y"
        df.ix[np.where(df["chrom"] == "25")[0], "chrom"] = "PAR"
        df.ix[np.where(df["chrom"] == "26")[0], "chrom"] = "MT"

        return df, "AncestryDNA"

    @staticmethod
    def _read_lineage_csv(file, comments):
        """ Read and parse CSV file generated by lineage.

        Parameters
        ----------
        file : str
            path to file
        comments : str
            comments at beginning of file

        Returns
        -------
        pandas.DataFrame
            individual's genetic data normalized for use with `lineage`
        str
            name of data source(s)
        """
        source = ""
        for comment in comments.split("\n"):
            if "Source(s):" in comment:
                source = comment.split("Source(s):")[1].strip()
                break

        df = pd.read_csv(
            file,
            comment="#",
            header=0,
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype={"chrom": object, "pos": np.int64},
        )

        return df, source

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
        str
            name of data source
        """
        df = pd.read_csv(
            file,
            skiprows=1,
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype={"chrom": object, "pos": np.int64},
        )

        return df, "generic"

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
            detected build of SNPs, else None

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
                return None

        build = None

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

            if build is not None:
                break

        return build

    def get_assembly(self):
        """ Get the assembly of a build.

        Returns
        -------
        str
            empty str if `build` is None
        """

        if self._build is None:
            return ""
        elif self._build == 36:
            return "NCBI36"
        elif self._build == 38:
            return "GRCh38"
        else:
            return "GRCh37"

    def get_snp_count(self):
        """ Count of SNPs.

        Returns
        -------
        int
        """

        if self._snps is not None:
            return len(self._snps)
        else:
            return 0

    def get_chromosomes(self):
        """ Get the chromosomes of SNPs.

        Returns
        -------
        list
            list of str chromosomes (e.g., ['1', '2', '3', 'MT'], empty list if no chromosomes
        """

        if isinstance(self._snps, pd.DataFrame):
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

        if isinstance(self._snps, pd.DataFrame):
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
        self, snps, y_snps_not_null_threshold=0.1, heterozygous_x_snps_threshold=0.01
    ):
        """ Determine sex from SNPs using thresholds.

        Parameters
        ----------
        snps : pandas.DataFrame
        y_snps_not_null_threshold : float
            percentage Y SNPs that are not null; above this threshold, Male is determined
        heterozygous_x_snps_threshold : float
            percentage heterozygous X SNPs; above this threshold, Female is determined

        Returns
        -------
        str
            'Male' or 'Female' if detected, else empty str
        """

        if isinstance(snps, pd.DataFrame):
            y_snps = len(snps.loc[(snps["chrom"] == "Y")])

            if y_snps > 0:
                y_snps_not_null = len(
                    snps.loc[(snps["chrom"] == "Y") & (snps["genotype"].notnull())]
                )

                if y_snps_not_null / y_snps > y_snps_not_null_threshold:
                    return "Male"
                else:
                    return "Female"

            x_snps = len(snps.loc[snps["chrom"] == "X"])

            if x_snps == 0:
                return ""

            heterozygous_x_snps = len(
                snps.loc[
                    (snps["chrom"] == "X")
                    & (snps["genotype"].notnull())
                    & (snps["genotype"].str[0] != snps["genotype"].str[1])
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

    # https://stackoverflow.com/a/16090640
    @staticmethod
    def _natural_sort_key(s, natural_sort_re=re.compile("([0-9]+)")):
        return [
            int(text) if text.isdigit() else text.lower()
            for text in re.split(natural_sort_re, s)
        ]


class SNPsCollection(SNPs):
    def __init__(self, raw_data=None, output_dir="output", name=""):
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
        super().__init__(file="", output_dir=output_dir)

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
        if filename is None:
            filename = (
                self.get_var_repr(self._name) + "_lineage_" + self.assembly + ".csv"
            )
        return super().save_snps(filename)

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
            filename = self.get_var_repr(self._name) + "_" + name + ".csv"

        return lineage.save_df_as_csv(
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

        if snps._snps is None:
            return discrepant_positions, discrepant_genotypes

        build = snps._build
        source = [s.strip() for s in snps._source.split(",")]

        if not snps._build_detected:
            print("build not detected, assuming build {}".format(snps._build))

        if self._build is None:
            self._build = build
        elif self._build != build:
            print(
                "build / assembly mismatch between current build of SNPs and SNPs being loaded"
            )

        # ensure there area always two X alleles
        snps = self._double_single_alleles(snps._snps, "X")

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
                        self.get_var_repr(self._name)
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
                        self.get_var_repr(self._name)
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
        df.ix[single_alleles, "genotype"] = df.ix[single_alleles, "genotype"] * 2

        return df
