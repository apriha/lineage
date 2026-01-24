"""
MIT License

Copyright (c) 2019 Andrew Riha

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

import os
from unittest import TestCase

import numpy as np
import pandas as pd
from pandas.api.types import is_object_dtype, is_string_dtype
from snps.io.reader import NORMALIZED_DTYPES

from lineage import Lineage


def get_complement(base):
    """Get the complement of a DNA base."""
    if base == "A":
        return "T"
    elif base == "G":
        return "C"
    elif base == "C":
        return "G"
    elif base == "T":
        return "A"
    else:
        return base


def complement_one_chrom(genotype):
    """Complement the genotype for one chromosome."""
    if pd.isnull(genotype):
        return np.nan

    complement = ""
    for base in list(genotype):
        complement += get_complement(base)
        complement += genotype[1]
        return complement


def complement_two_chroms(genotype):
    """Complement the genotype for both chromosomes."""
    if pd.isnull(genotype):
        return np.nan

    complement = ""
    for base in list(genotype):
        complement += get_complement(base)
    return complement


def simulate_snps(
    ind,
    chrom="1",
    pos_start=1,
    pos_max=111700002,
    pos_step=10000,
    genotype="AA",
    insert_nulls=True,
    null_snp_step=101,
    complement_genotype_one_chrom=False,
    complement_genotype_two_chroms=False,
    complement_snp_step=50,
):
    """Simulate SNP data for an individual."""
    ind._build = 37

    positions = np.arange(pos_start, pos_max, pos_step, dtype=np.int64)
    snps = pd.DataFrame(
        {"chrom": chrom},
        index=pd.Index(["rs" + str(x + 1) for x in range(len(positions))], name="rsid"),
    )
    snps["pos"] = positions
    snps["genotype"] = genotype

    if insert_nulls:
        snps.loc[snps.iloc[0::null_snp_step, :].index, "genotype"] = np.nan

    indices = snps.iloc[0::complement_snp_step, :].index
    if complement_genotype_two_chroms:
        snps.loc[indices, "genotype"] = snps.loc[indices, "genotype"].apply(
            complement_two_chroms
        )
    elif complement_genotype_one_chrom:
        snps.loc[indices, "genotype"] = snps.loc[indices, "genotype"].apply(
            complement_one_chrom
        )

    ind._snps = snps

    return ind


class BaseLineageTestCase(TestCase):
    def setUp(self):
        self.l = Lineage()

    def simulate_snps(
        self,
        ind,
        chrom="1",
        pos_start=1,
        pos_max=111700002,
        pos_step=10000,
        genotype="AA",
        insert_nulls=True,
        null_snp_step=101,
        complement_genotype_one_chrom=False,
        complement_genotype_two_chroms=False,
        complement_snp_step=50,
    ):
        """Simulate SNP data for an individual (wrapper for standalone simulate_snps)."""
        return simulate_snps(
            ind,
            chrom,
            pos_start,
            pos_max,
            pos_step,
            genotype,
            insert_nulls,
            null_snp_step,
            complement_genotype_one_chrom,
            complement_genotype_two_chroms,
            complement_snp_step,
        )

    @staticmethod
    def create_snp_df(rsid, chrom, pos, genotype):
        df = pd.DataFrame(
            {"rsid": rsid, "chrom": chrom, "pos": pos, "genotype": genotype},
            columns=["rsid", "chrom", "pos", "genotype"],
        )
        df = df.astype(NORMALIZED_DTYPES)
        df = df.set_index("rsid")
        return df

    def generic_snps(self):
        return self.create_snp_df(
            rsid=["rs" + str(i) for i in range(1, 9)],
            chrom=["1"] * 8,
            pos=list(range(101, 109)),
            genotype=["AA", "CC", "GG", "TT", np.nan, "GC", "TC", "AT"],
        )

    def assert_series_equal_with_string_dtype(self, left, right, **kwargs):
        """Assert Series are equal, accepting both object and StringDtype for string data.

        In Python 3.14+, pandas infers StringDtype for string data instead of object.
        This wrapper compares Series without strict dtype matching for string data.

        Parameters
        ----------
        left : pd.Series
            First Series to compare
        right : pd.Series
            Second Series to compare
        **kwargs : dict
            Additional arguments passed to pd.testing.assert_series_equal
        """
        # Verify string series have string or object dtypes
        if is_string_dtype(left.dtype) or is_object_dtype(left.dtype):
            self.assertTrue(
                is_string_dtype(right.dtype) or is_object_dtype(right.dtype),
                f"Right series dtype {right.dtype} should be string/object type",
            )
        # Compare Series without strict dtype matching
        pd.testing.assert_series_equal(left, right, check_dtype=False, **kwargs)

    def assert_frame_equal_with_string_index(self, left, right, **kwargs):
        """Assert DataFrames are equal, accepting both object and StringDtype for string columns.

        In Python 3.14+, pandas infers StringDtype for string columns/indices instead of object.
        This wrapper validates that string columns have string types, then compares the
        DataFrames without strict dtype matching for object/string columns.

        Parameters
        ----------
        left : pd.DataFrame
            First DataFrame to compare
        right : pd.DataFrame
            Second DataFrame to compare
        **kwargs : dict
            Additional arguments passed to pd.testing.assert_frame_equal
        """
        # Verify index dtypes are string types if they're named 'rsid'
        if left.index.name == "rsid":
            self.assertTrue(
                is_string_dtype(left.index.dtype),
                f"Left index dtype {left.index.dtype} is not a string type",
            )
        if right.index.name == "rsid":
            self.assertTrue(
                is_string_dtype(right.index.dtype),
                f"Right index dtype {right.index.dtype} is not a string type",
            )

        # Verify string columns (chrom, genotype) have string dtypes
        for col in ["chrom", "genotype"]:
            if col in left.columns:
                self.assertTrue(
                    is_string_dtype(left[col].dtype)
                    or is_object_dtype(left[col].dtype),
                    f"Left column '{col}' dtype {left[col].dtype} is not a string/object type",
                )
            if col in right.columns:
                self.assertTrue(
                    is_string_dtype(right[col].dtype)
                    or is_object_dtype(right[col].dtype),
                    f"Right column '{col}' dtype {right[col].dtype} is not a string/object type",
                )

        # Compare DataFrames without strict dtype matching for string columns
        pd.testing.assert_frame_equal(
            left, right, check_index_type=False, check_dtype=False, **kwargs
        )

    @property
    def downloads_enabled(self):
        """Property indicating if downloads are enabled.

        Only download from external resources when an environment variable named
        "DOWNLOADS_ENABLED" is set to "true".

        Returns
        -------
        bool
        """
        return True if os.getenv("DOWNLOADS_ENABLED") == "true" else False
