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

from lineage.ensembl import EnsemblRestClient


class SNPs(object):
    def __init__(self, file, assign_par_snps=True):
        """ Object used to read and parse genotype / raw data files.

        Parameters
        ----------
        file : str
            path to file to load
        assign_par_snps : bool
            assign PAR SNPs to the X and Y chromosomes
        """
        self.snps, self.source = self._read_raw_data(file)
        self.build = None
        self.build_detected = False

        if self.snps is not None:
            self.build = detect_build(self.snps)

            if self.build is None:
                self.build = 37  # assume Build 37 / GRCh37 if not detected
            else:
                self.build_detected = True

            if assign_par_snps:
                self._assign_par_snps()

    @property
    def assembly(self):
        """ Get the assembly of ``SNPs``.

        Returns
        -------
        str
        """
        return get_assembly(self.build)

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

    @property
    def sex(self):
        """ Sex derived from ``SNPs``.

        Returns
        -------
        str
            'Male' or 'Female' if detected, else empty str
        """
        return determine_sex(self.snps)

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
        if self.snps is None:
            return False
        else:
            return True

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

        return sort_snps(df), "23andMe"

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

        return sort_snps(df), "FTDNA"

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

        return sort_snps(df), "FTDNA"

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

        return sort_snps(df), "AncestryDNA"

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

        return sort_snps(df), source

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

        return sort_snps(df), "generic"

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
        rest_client = EnsemblRestClient(server="https://api.ncbi.nlm.nih.gov")
        for rsid in self.snps.loc[self.snps["chrom"] == "PAR"].index.values:
            if "rs" in rsid:
                try:
                    id = rsid.split("rs")[1]
                    response = rest_client.perform_rest_action(
                        "/variation/v0/beta/refsnp/" + id
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
                                if not self.build_detected:
                                    self.build = self._extract_build(item)
                                    self.build_detected = True
                                continue

                except Exception as err:
                    print(err)

    def _assign_snp(self, rsid, alleles, chrom):
        for allele in alleles:
            allele_pos = allele["allele"]["spdi"]["position"]
            # ref SNP positions seem to be 0-based...
            if allele_pos == self.snps.loc[rsid].pos - 1:
                self.snps.loc[rsid, "chrom"] = chrom
                return True
        return False

    def _extract_build(self, item):
        assembly_name = item["placement_annot"]["seq_id_traits_by_assembly"][0][
            "assembly_name"
        ]
        assembly_name = assembly_name.split(".")[0]
        return int(assembly_name[-2:])


def detect_build(snps):
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

    Parameters
    ----------
    snps : pandas.DataFrame
        SNPs to add

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
        if rsid in snps.index:
            build = lookup_build_with_snp_pos(snps.loc[rsid].pos, df.loc[rsid])

        if build is not None:
            break

    return build


def get_assembly(build):
    """ Get the assembly of a build.

    Parameters
    ----------
    build : int {36, 37, 38}

    Returns
    -------
    str
        empty str if `build` is None
    """

    if build is None:
        return ""
    elif build == 36:
        return "NCBI36"
    elif build == 38:
        return "GRCh38"
    else:
        return "GRCh37"


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

    if isinstance(snps, pd.DataFrame):
        return list(pd.unique(snps["chrom"]))
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

    if isinstance(snps, pd.DataFrame):
        chroms = list(pd.unique(snps["chrom"]))

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
    snps, y_snps_not_null_threshold=0.1, heterozygous_x_snps_threshold=0.01
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


def sort_snps(snps):
    """ Sort SNPs based on ordered chromosome list and position. """

    sorted_list = sorted(snps["chrom"].unique(), key=_natural_sort_key)

    # move PAR and MT to the end of the dataframe
    if "PAR" in sorted_list:
        sorted_list.remove("PAR")
        sorted_list.append("PAR")

    if "MT" in sorted_list:
        sorted_list.remove("MT")
        sorted_list.append("MT")

    # convert chrom column to category for sorting
    # https://stackoverflow.com/a/26707444
    snps["chrom"] = snps["chrom"].astype(
        CategoricalDtype(categories=sorted_list, ordered=True)
    )

    # sort based on ordered chromosome list and position
    snps = snps.sort_values(["chrom", "pos"])

    # convert chromosome back to object
    snps["chrom"] = snps["chrom"].astype(object)

    return snps


# https://stackoverflow.com/a/16090640
def _natural_sort_key(s, natural_sort_re=re.compile("([0-9]+)")):
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(natural_sort_re, s)
    ]
