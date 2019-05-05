"""
Copyright (C) 2019 Andrew Riha

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
import gzip
import zipfile

import numpy as np
import pandas as pd


class Reader:
    """ Class for reading and parsing raw data / genotype files. """

    def __init__(self, file=""):
        """ Initialize a `Reader`.

        Parameters
        ----------
        file : str
            path to file to read
        """
        self.file = file

    def __call__(self):
        """ Read and parse a raw data / genotype file.

        Returns
        -------
        tuple : (pandas.DataFrame, str)
            dataframe of parsed SNPs, detected source of SNPs
        """
        file = self.file

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
                return self.read_23andme(file)
            elif "Ancestry" in first_line:
                return self.read_ancestry(file)
            elif first_line.startswith("RSID"):
                return self.read_ftdna(file)
            elif "famfinder" in first_line:
                return self.read_ftdna_famfinder(file)
            elif "MyHeritage" in first_line:
                return self.read_myheritage(file)
            elif "lineage" in first_line:
                return self.read_lineage_csv(file, comments)
            elif first_line.startswith("rsid"):
                return self.read_generic_csv(file)
            else:
                return None, ""
        except Exception as err:
            print(err)
            return None, ""

    @classmethod
    def read_file(cls, file):
        """ Read `file`.

        Parameters
        ----------
        file : str
            path to file to read

        Returns
        -------
        tuple : (pandas.DataFrame, str)
            dataframe of parsed SNPs, detected source of SNPs
        """
        p = cls(file)
        return p()

    def _extract_comments(self, f, decode):
        line = self._read_line(f, decode)
        first_line = line
        comments = ""

        while line.startswith("#"):
            comments += line
            line = self._read_line(f, decode)

        return first_line, comments

    @staticmethod
    def _read_line(f, decode):
        if decode:
            # https://stackoverflow.com/a/606199
            return f.readline().decode("utf-8")
        else:
            return f.readline()

    @staticmethod
    def read_23andme(file):
        """ Read and parse 23andMe file.

        https://www.23andme.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            genetic data normalized for use with `lineage`
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
    def read_ftdna(file):
        """ Read and parse Family Tree DNA (FTDNA) file.

        https://www.familytreedna.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            genetic data normalized for use with `lineage`
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
    def read_ftdna_famfinder(file):
        """ Read and parse Family Tree DNA (FTDNA) "famfinder" file.

        https://www.familytreedna.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            genetic data normalized for use with `lineage`
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
    def read_ancestry(file):
        """ Read and parse Ancestry.com file.

        http://www.ancestry.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            genetic data normalized for use with `lineage`
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
    def read_myheritage(file):
        """ Read and parse MyHeritage file.

        https://www.myheritage.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        pandas.DataFrame
            genetic data normalized for use with `lineage`
        str
            name of data source
        """
        df = pd.read_csv(
            file,
            comment="#",
            header=0,
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype={"chrom": object, "pos": np.int64},
        )

        return df, "MyHeritage"

    @staticmethod
    def read_lineage_csv(file, comments):
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
            genetic data normalized for use with `lineage`
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
    def read_generic_csv(file):
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
            genetic data normalized for use with `lineage`
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
