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

import datetime
import os
import gzip
import zipfile

import numpy as np
import pandas as pd
import vcf

import lineage
from lineage.utils import save_df_as_csv, clean_str


class Reader:
    """ Class for reading and parsing raw data / genotype files. """

    def __init__(self, file=""):
        """ Initialize a `Reader`.

        Parameters
        ----------
        file : str
            path to file to read
        """
        self._file = file

    def __call__(self):
        """ Read and parse a raw data / genotype file.

        Returns
        -------
        tuple : (pandas.DataFrame, str)
            dataframe of parsed SNPs, detected source of SNPs
        """
        file = self._file

        try:
            if not os.path.exists(file):
                print(file + " does not exist; skipping")
                return pd.DataFrame(), ""

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
            elif "vcf" in comments.lower():
                return self.read_vcf(file)
            else:
                return pd.DataFrame(), ""
        except Exception as err:
            print(err)
            return pd.DataFrame(), ""

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
        r = cls(file)
        return r()

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
        df.iloc[np.where(df["chrom"] == "23")[0], 0] = "X"
        df.iloc[np.where(df["chrom"] == "24")[0], 0] = "Y"
        df.iloc[np.where(df["chrom"] == "25")[0], 0] = "PAR"
        df.iloc[np.where(df["chrom"] == "26")[0], 0] = "MT"

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

    @staticmethod
    def read_vcf(file):
        """ Read and parse VCF file.

        Notes
        -----
        This function uses the PyVCF python module to parse the genotypes from VCF files:
        https://pyvcf.readthedocs.io/en/latest/index.html


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
        df = pd.DataFrame(columns=["rsid", "chrom", "pos", "genotype"])
        df = df.astype(
            {"rsid": object, "chrom": object, "pos": np.int64, "genotype": object}
        )

        vcf_reader = vcf.Reader(open(file, "r"))

        # lineage does not yet support multi-sample vcf.
        if len(vcf_reader.samples) > 1:
            print(
                "Multiple samples detected in the vcf file, please use a single sample vcf."
            )
            return df, "vcf"

        for i, record in enumerate(vcf_reader):
            # assign null genotypes if either allele is None
            # Could capture full genotype, if REF is None, but genotype is 1/1 or
            # if ALT is None, but genotype is 0/0
            if record.REF is None or record.ALT[0] is None:
                genotype = np.nan
            # skip SNPs with missing rsIDs.
            elif record.ID is None:
                continue
            # skip insertions and deletions
            elif len(record.REF) > 1 or len(record.ALT[0]) > 1:
                continue
            else:
                alleles = record.genotype(vcf_reader.samples[0]).gt_bases
                a1 = alleles[0]
                a2 = alleles[-1]
                genotype = "{}{}".format(a1, a2)

            record_info = {
                "rsid": record.ID,
                "chrom": "{}".format(record.CHROM).strip("chr"),
                "pos": record.POS,
                "genotype": genotype,
            }
            # append the record to the DataFrame
            df = df.append(pd.DataFrame([record_info]), ignore_index=True, sort=False)
        df.set_index("rsid", inplace=True, drop=True)

        return df, "vcf"


class Writer:
    """ Class for writing SNPs to files. """

    def __init__(self, snps=None, filename="", vcf=False):
        """ Initialize a `Writer`.

        Parameters
        ----------
        snps : SNPs
            SNPs to save to file
        filename : str
            filename for file to save
        vcf : bool
            flag to save file as VCF
        """
        self._snps = snps
        self._filename = filename
        self._vcf = vcf

    def __call__(self):
        if self._vcf:
            return self._write_vcf()
        else:
            return self._write_csv()

    @classmethod
    def write_file(cls, snps=None, filename="", vcf=False):
        """ Save SNPs to file.

        Parameters
        ----------
        snps : SNPs
            SNPs to save to file
        filename : str
            filename for file to save
        vcf : bool
            flag to save file as VCF

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        w = cls(snps=snps, filename=filename, vcf=vcf)
        return w()

    def _write_csv(self):
        """ Write SNPs to a CSV file.

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        filename = self._filename
        if not filename:
            filename = "{}_lineage_{}{}".format(
                clean_str(self._snps._source), self._snps.assembly, ".csv"
            )

        comment = (
            "# Source(s): {}\n"
            "# Assembly: {}\n"
            "# SNPs: {}\n"
            "# Chromosomes: {}\n".format(
                self._snps.source,
                self._snps.assembly,
                self._snps.snp_count,
                self._snps.chromosomes_summary,
            )
        )

        return save_df_as_csv(
            self._snps._snps,
            self._snps._output_dir,
            filename,
            comment=comment,
            header=["chromosome", "position", "genotype"],
        )

    def _write_vcf(self):
        """ Write SNPs to a VCF file.

        References
        ----------
        ..[1] The Variant Call Format (VCF) Version 4.2 Specification, 8 Mar 2019,
          https://samtools.github.io/hts-specs/VCFv4.2.pdf

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        filename = self._filename
        if not filename:
            filename = "{}_lineage_{}{}".format(
                clean_str(self._snps._source), self._snps.assembly, ".vcf"
            )

        comment = (
            "##fileformat=VCFv4.2\n"
            "##fileDate={}\n"
            '##source="lineage v{}; https://github.com/apriha/lineage"\n'.format(
                datetime.datetime.utcnow().strftime("%Y%m%d"), lineage.__version__
            )
        )

        reference_sequence_chroms = (
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "X",
            "Y",
            "MT",
        )

        df = self._snps.snps

        tasks = []

        # skip insertions and deletions
        df = df.drop(
            df.loc[
                df["genotype"].notnull()
                & (
                    (df["genotype"].str[0] == "I")
                    | (df["genotype"].str[0] == "D")
                    | (df["genotype"].str[1] == "I")
                    | (df["genotype"].str[1] == "D")
                )
            ].index
        )

        chroms_to_drop = []
        for chrom in df["chrom"].unique():
            if chrom not in reference_sequence_chroms:
                chroms_to_drop.append(chrom)
                continue

            tasks.append(
                {
                    "resources": self._snps._resources,
                    "assembly": self._snps.assembly,
                    "chrom": chrom,
                    "snps": pd.DataFrame(df.loc[(df["chrom"] == chrom)]),
                }
            )

        # drop chromosomes without reference sequence data (e.g., unassigned PAR)
        for chrom in chroms_to_drop:
            df = df.drop(df.loc[df["chrom"] == chrom].index)

        # create the VCF representation for SNPs
        results = map(self._create_vcf_representation, tasks)

        contigs = []
        vcf = []
        for result in list(results):
            contigs.append(result["contig"])
            vcf.append(result["vcf"])

        vcf = pd.concat(vcf)

        comment += "".join(contigs)
        comment += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        comment += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"

        return save_df_as_csv(
            vcf,
            self._snps._output_dir,
            filename,
            comment=comment,
            prepend_info=False,
            header=False,
            index=False,
            na_rep=".",
            sep="\t",
        )

    def _create_vcf_representation(self, task):
        resources = task["resources"]
        assembly = task["assembly"]
        chrom = task["chrom"]
        snps = task["snps"]

        seqs = resources.get_reference_sequences(assembly, [chrom])
        seq = seqs[chrom]

        contig = '##contig=<ID={},URL={},length={},assembly={},md5={},species="{}">\n'.format(
            seq.ID, seq.url, seq.length, seq.build, seq.md5, seq.species
        )

        snps = snps.reset_index()

        df = pd.DataFrame(
            columns=[
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                "SAMPLE",
            ]
        )
        df = df.astype(
            {
                "CHROM": object,
                "POS": np.int64,
                "ID": object,
                "REF": object,
                "ALT": object,
                "QUAL": np.int64,
                "FILTER": object,
                "INFO": object,
                "FORMAT": object,
                "SAMPLE": object,
            }
        )

        df["CHROM"] = snps["chrom"]
        df["POS"] = snps["pos"]
        df["ID"] = snps["rsid"]

        # https://stackoverflow.com/a/24838429
        df["REF"] = list(map(chr, seq.sequence[snps.pos - seq.start]))

        df["FORMAT"] = "GT"

        seq.clear()

        df["genotype"] = snps["genotype"]

        temp = df.loc[df["genotype"].notnull()]

        # https://stackoverflow.com/a/19976286
        df.loc[df["genotype"].notnull(), "ALT"] = np.vectorize(self._compute_alt)(
            temp["REF"], temp["genotype"]
        )

        temp = df.loc[df["genotype"].notnull()]

        df.loc[df["genotype"].notnull(), "SAMPLE"] = np.vectorize(
            self._compute_genotype
        )(temp["REF"], temp["ALT"], temp["genotype"])

        df.loc[df["SAMPLE"].isnull(), "SAMPLE"] = "./."

        del df["genotype"]

        return {"contig": contig, "vcf": df}

    def _compute_alt(self, ref, genotype):
        genotype_alleles = list(set(genotype))

        if ref in genotype_alleles:
            if len(genotype_alleles) == 1:
                return "N"
            else:
                genotype_alleles.remove(ref)
                return genotype_alleles.pop(0)
        else:
            return ",".join(genotype_alleles)

    def _compute_genotype(self, ref, alt, genotype):
        alleles = [ref]

        if pd.notna(alt):
            alleles.extend(alt.split(","))

        if len(genotype) == 2:
            return "{}/{}".format(
                alleles.index(genotype[0]), alleles.index(genotype[1])
            )
        else:
            return "{}".format(alleles.index(genotype[0]))
