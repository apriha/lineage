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

import pandas as pd
import re


class Individual(object):
    """ Object used to represent and interact with an individual.

    The `Individual` object maintains information about an individual.
    Additionally, the object provides methods for loading and analyzing
    an individual's genetic data.

    """

    def __init__(self, name, filename, relation="self"):
        """ Initialize an individual.

        Parameters
        ----------
        name : str
            name of the individual
        filename : str
            path to CSV file with raw genetic data
        relation : {"self", "parent", "child"}, optional
            relation of the individual in the context of a family

        """

        self.name = name

        if relation not in ["self", "parent", "child"]:
            raise ValueError("Invalid relation for Individual")

        self.relation = relation
        self.genome = self._load_raw_data(filename)

    @staticmethod
    def _load_raw_data(filename):
        """ Load raw genetic data.

        Parameters
        ----------
        filename : str
            path to CSV file with raw genetic data

        Returns
        -------
        pandas.DataFrame
            individual's genetic data normalized for use with `lineage`

        """

        with open(filename, "r") as f:
            line = f.readline()

        try:
            if "23andMe" in line:
                # file is from 23andMe (https://www.23andme.com)
                return pd.read_csv(filename, comment="#", sep="\t", na_values="--",
                                   names=["rsid", "chromosome", "position", "genotype"],
                                   index_col=0, dtype={"chromosome": str})
            elif "Ancestry" in line:
                # file is from Ancestry (http://www.ancestry.com)

                # read file
                df = pd.read_csv(filename, comment="#", header=0, sep="\t", na_values=0,
                                 names=["rsid", "chromosome", "position", "allele1", "allele2"],
                                 index_col=0, dtype={"chromosome": str})

                # create genotype column from allele columns [SO-01]
                df["genotype"] = df["allele1"].str.cat(df["allele2"])

                # delete allele columns [SO-02]
                del df["allele1"]
                del df["allele2"]

                return df
            elif line[:4] == "RSID":
                # assume file is from FTDNA (https://www.familytreedna.com)
                return pd.read_csv(filename, skiprows=1, na_values="--",
                                   names=["rsid", "chromosome", "position", "genotype"],
                                   index_col=0, dtype={"chromosome": str})
            else:
                return None
        except Exception as e:
            print(e)
            return None

    @staticmethod
    def _clean(s):
        """ Clean a string so that it can be a valid Python variable
        name.

        Returns
        -------
        str
            cleaned string that can be used as a variable name

        """

        return re.sub('\W|^(?=\d)', '_', s)  # [SO-03]

    def find_discordant_snps(self, individual_to_compare):
        """ Find discordant SNPs between this individual and another
        individual.

        Parameters
        ----------
        individual_to_compare : Individual
            Individual to compare against

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

        reference = "genotype_" + self._clean(self.name)
        comparison = "genotype_" + self._clean(individual_to_compare.name)

        df = self.genome.rename(columns={"genotype": reference})
        df = df.join(individual_to_compare.genome["genotype"])
        df.rename(columns={"genotype": comparison}, inplace=True)

        return df.loc[df[reference].notnull() & df[comparison].notnull() &
                      ~(df[reference].str[0] == df[comparison].str[0]) &
                      ~(df[reference].str[0] == df[comparison].str[1]) &
                      ~(df[reference].str[1] == df[comparison].str[0]) &
                      ~(df[reference].str[1] == df[comparison].str[1])]

"""
Stack Overflow Attributions
---------------------------

[SO-##] references throughout the code above correspond to the following
content from Stack Overflow (http://stackoverflow.com):

[SO-01] "Combine two columns of text in dataframe in pandas/python"
        http://stackoverflow.com/q/19377969
        user2866103 : http://stackoverflow.com/users/2866103/user2866103
        http://stackoverflow.com/a/35850749
        LeoRochael : http://stackoverflow.com/users/1273938/leorochael

[SO-02] "Delete column from pandas DataFrame"
        http://stackoverflow.com/q/13411544
        John : http://stackoverflow.com/users/390388/john
        http://stackoverflow.com/a/13485766
        Wes McKinney : http://stackoverflow.com/users/776560/wes-mckinney

[SO-03] "how do I convert a string to a valid variable name in python?"
        http://stackoverflow.com/q/3303312
        George Profenza : http://stackoverflow.com/users/89766/george-profenza
        http://stackoverflow.com/a/3305731
        Nas Banov : http://stackoverflow.com/users/226086/nas-banov

"""
