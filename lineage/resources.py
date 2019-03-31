""" Class for downloading and loading required external resources.

`lineage` uses tables and data from UCSC's Genome Browser:

* http://genome.ucsc.edu/
* http://genome.ucsc.edu/cgi-bin/hgTables

References
----------
..[1] Karolchik D, Hinrichs AS, Furey TS, Roskin KM, Sugnet CW, Haussler D, Kent WJ.
  The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 2004 Jan
  1;32(Database issue):D493-6. PubMed PMID: 14681465; PubMed Central PMCID:
  PMC308837. https://www.ncbi.nlm.nih.gov/pubmed/14681465
..[2] Tyner C, Barber GP, Casper J, Clawson H, Diekhans M, Eisenhart C, Fischer CM,
  Gibson D, Gonzalez JN, Guruvadoo L, Haeussler M, Heitner S, Hinrichs AS,
  Karolchik D, Lee BT, Lee CM, Nejad P, Raney BJ, Rosenbloom KR, Speir ML,
  Villarreal C, Vivian J, Zweig AS, Haussler D, Kuhn RM, Kent WJ. The UCSC Genome
  Browser database: 2017 update. Nucleic Acids Res. 2017 Jan 4;45(D1):D626-D634.
  doi: 10.1093/nar/gkw1134. Epub 2016 Nov 29. PubMed PMID: 27899642; PubMed Central
  PMCID: PMC5210591. https://www.ncbi.nlm.nih.gov/pubmed/27899642
..[3] International Human Genome Sequencing Consortium. Initial sequencing and
  analysis of the human genome. Nature. 2001 Feb 15;409(6822):860-921.
  http://dx.doi.org/10.1038/35057062
..[4] hg19 (GRCh37): Hiram Clawson, Brooke Rhead, Pauline Fujita, Ann Zweig, Katrina
  Learned, Donna Karolchik and Robert Kuhn, https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19
..[5] Yates et. al. (doi:10.1093/bioinformatics/btu613),
  http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613
..[6] Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098

"""

"""
Copyright (C) 2017 Andrew Riha

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

import gzip
import itertools
import json
import os
import tarfile
import tempfile
import urllib.error
import urllib.request
import zlib

import pandas as pd

import lineage


class Resources(object):
    """ Object used to manage resources required by `lineage`. """

    def __init__(self, resources_dir="resources", ensembl_rest_client=None):
        """ Initialize a ``Resources`` object.

        Parameters
        ----------
        resources_dir : str
            name / path of resources directory
        """
        self._resources_dir = os.path.abspath(resources_dir)
        self._genetic_map_HapMapII_GRCh37 = None
        self._cytoBand_hg19 = None
        self._knownGene_hg19 = None
        self._kgXref_hg19 = None
        self._ensembl_rest_client = ensembl_rest_client

    def get_genetic_map_HapMapII_GRCh37(self):
        """ Get International HapMap Consortium HapMap Phase II genetic map for Build 37.

        Returns
        -------
        dict
            dict of pandas.DataFrame HapMapII genetic maps if loading was successful, else None
        """
        if self._genetic_map_HapMapII_GRCh37 is None:
            self._genetic_map_HapMapII_GRCh37 = self._load_genetic_map(
                self._get_path_genetic_map_HapMapII_GRCh37()
            )

        return self._genetic_map_HapMapII_GRCh37

    def get_cytoBand_hg19(self):
        """ Get UCSC cytoBand table for Build 37.

        Returns
        -------
        pandas.DataFrame
            cytoBand table if loading was successful, else None
        """
        if self._cytoBand_hg19 is None:
            self._cytoBand_hg19 = self._load_cytoBand(self._get_path_cytoBand_hg19())

        return self._cytoBand_hg19

    def get_knownGene_hg19(self):
        """ Get UCSC knownGene table for Build 37.

        Returns
        -------
        pandas.DataFrame
            knownGene table if loading was successful, else None
        """
        if self._knownGene_hg19 is None:
            self._knownGene_hg19 = self._load_knownGene(self._get_path_knownGene_hg19())

        return self._knownGene_hg19

    def get_kgXref_hg19(self):
        """ Get UCSC kgXref table for Build 37.

        Returns
        -------
        pandas.DataFrame
            kgXref table if loading was successful, else None
        """
        if self._kgXref_hg19 is None:
            self._kgXref_hg19 = self._load_kgXref(self._get_path_kgXref_hg19())

        return self._kgXref_hg19

    def get_assembly_mapping_data(self, source_assembly, target_assembly):
        """ Get assembly mapping data.

        Parameters
        ----------
        source_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap from
        target_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap to

        Returns
        -------
        dict
            dict of json assembly mapping data if loading was successful, else None
        """
        return self._load_assembly_mapping_data(
            self._get_path_assembly_mapping_data(source_assembly, target_assembly)
        )

    def download_example_datasets(self):
        """ Download example datasets from `openSNP <https://opensnp.org>`_.

        Per openSNP, "the data is donated into the public domain using `CC0 1.0
        <http://creativecommons.org/publicdomain/zero/1.0/>`_."

        Returns
        -------
        paths : list of str or None
            paths to example datasets

        References
        ----------
        ..[1] Greshake B, Bayer PE, Rausch H, Reda J (2014), "openSNP-A Crowdsourced Web Resource
          for Personal Genomics," PLOS ONE, 9(3): e89204,
          https://doi.org/10.1371/journal.pone.0089204
        """
        paths = []
        paths.append(
            self._download_file(
                "https://opensnp.org/data/662.23andme.304",
                "662.23andme.304.txt.gz",
                compress=True,
            )
        )
        paths.append(
            self._download_file(
                "https://opensnp.org/data/662.23andme.340",
                "662.23andme.340.txt.gz",
                compress=True,
            )
        )
        paths.append(
            self._download_file(
                "https://opensnp.org/data/662.ftdna-illumina.341",
                "662.ftdna-illumina.341.csv.gz",
                compress=True,
            )
        )
        paths.append(
            self._download_file(
                "https://opensnp.org/data/663.23andme.305",
                "663.23andme.305.txt.gz",
                compress=True,
            )
        )

        # these two files consist of concatenated gzip files and therefore need special handling
        paths.append(
            self._download_file(
                "https://opensnp.org/data/4583.ftdna-illumina.3482",
                "4583.ftdna-illumina.3482.csv.gz",
            )
        )
        paths.append(
            self._download_file(
                "https://opensnp.org/data/4584.ftdna-illumina.3483",
                "4584.ftdna-illumina.3483.csv.gz",
            )
        )

        try:
            for gzip_path in paths[-2:]:
                # https://stackoverflow.com/q/4928560
                # https://stackoverflow.com/a/37042747
                with open(gzip_path, "rb") as f:
                    decompressor = zlib.decompressobj(31)

                    # decompress data from first concatenated gzip file
                    data = decompressor.decompress(f.read())

                    if len(decompressor.unused_data) > 0:
                        # decompress data from second concatenated gzip file, if any
                        additional_data = zlib.decompress(decompressor.unused_data, 31)
                        data += additional_data[33:]  # skip over second header

                # recompress data
                with gzip.open(gzip_path, "wb") as f:
                    f.write(data)
        except Exception as err:
            print(err)

        return paths

    def get_all_resources(self):
        """ Get / download all resources used throughout `lineage`.

        Returns
        -------
        dict
            dict of resources
        """
        resources = {}
        resources[
            "genetic_map_HapMapII_GRCh37"
        ] = self.get_genetic_map_HapMapII_GRCh37()
        resources["cytoBand_hg19"] = self.get_cytoBand_hg19()
        resources["knownGene_hg19"] = self.get_knownGene_hg19()
        resources["kgXref_hg19"] = self.get_kgXref_hg19()
        for source, target in itertools.permutations(["NCBI36", "GRCh37", "GRCh38"], 2):
            resources[source + "_" + target] = self.get_assembly_mapping_data(
                source, target
            )
        return resources

    @staticmethod
    def _load_genetic_map(filename):
        """ Load genetic map (e.g. HapMapII).

        Parameters
        ----------
        filename : str
            path to compressed archive with genetic map data

        Returns
        -------
        genetic_map : dict
            dict of pandas.DataFrame genetic maps if loading was successful, else None

        Notes
        -----
        Keys of returned dict are chromosomes and values are the corresponding genetic map.
        """
        try:
            genetic_map = {}

            with tarfile.open(filename, "r") as tar:
                # http://stackoverflow.com/a/2018576
                for member in tar.getmembers():
                    if "genetic_map" in member.name:
                        df = pd.read_csv(tar.extractfile(member), sep="\t")
                        df = df.rename(
                            columns={
                                "Position(bp)": "pos",
                                "Rate(cM/Mb)": "rate",
                                "Map(cM)": "map",
                            }
                        )
                        del df["Chromosome"]
                        start_pos = member.name.index("chr") + 3
                        end_pos = member.name.index(".")
                        genetic_map[member.name[start_pos:end_pos]] = df

            # X chrom consists of X PAR regions and X non-PAR region
            genetic_map["X"] = pd.concat(
                [genetic_map["X_par1"], genetic_map["X"], genetic_map["X_par2"]]
            )
            del genetic_map["X_par1"]
            del genetic_map["X_par2"]

            return genetic_map
        except Exception as err:
            print(err)
            return None

    @staticmethod
    def _load_assembly_mapping_data(filename):
        """ Load assembly mapping data.

        Parameters
        ----------
        filename : str
            path to compressed archive with assembly mapping data

        Returns
        -------
        assembly_mapping_data : dict
            dict of assembly maps if loading was successful, else None

        Notes
        -----
        Keys of returned dict are chromosomes and values are the corresponding assembly map.
        """
        try:
            assembly_mapping_data = {}

            with tarfile.open(filename, "r") as tar:
                # http://stackoverflow.com/a/2018576
                for member in tar.getmembers():
                    if ".json" in member.name:
                        with tar.extractfile(member) as tar_file:
                            tar_bytes = tar_file.read()
                        # https://stackoverflow.com/a/42683509/4727627
                        assembly_mapping_data[member.name.split(".")[0]] = json.loads(
                            tar_bytes.decode("utf-8")
                        )

            return assembly_mapping_data
        except Exception as err:
            print(err)
            return None

    @staticmethod
    def _load_cytoBand(filename):
        """ Load UCSC cytoBand table.

        Parameters
        ----------
        filename : str
            path to cytoBand file

        Returns
        -------
        df : pandas.DataFrame
            cytoBand table if loading was successful, else None

        References
        ----------
        ..[1] Ryan Dale, GitHub Gist,
          https://gist.github.com/daler/c98fc410282d7570efc3#file-ideograms-py
        """
        try:
            # adapted from chromosome plotting code (see [1]_)
            df = pd.read_table(
                filename, names=["chrom", "start", "end", "name", "gie_stain"]
            )
            df["chrom"] = df["chrom"].str[3:]
            return df
        except Exception as err:
            print(err)
            return None

    @staticmethod
    def _load_knownGene(filename):
        """ Load UCSC knownGene table.

        Parameters
        ----------
        filename : str
            path to knownGene file

        Returns
        -------
        df : pandas.DataFrame
            knownGene table if loading was successful, else None
        """
        try:
            df = pd.read_table(
                filename,
                names=[
                    "name",
                    "chrom",
                    "strand",
                    "txStart",
                    "txEnd",
                    "cdsStart",
                    "cdsEnd",
                    "exonCount",
                    "exonStarts",
                    "exonEnds",
                    "proteinID",
                    "alignID",
                ],
                index_col=0,
            )
            df["chrom"] = df["chrom"].str[3:]
            return df
        except Exception as err:
            print(err)
            return None

    @staticmethod
    def _load_kgXref(filename):
        """ Load UCSC kgXref table.

        Parameters
        ----------
        filename : str
            path to kgXref file

        Returns
        -------
        df : pandas.DataFrame
            kgXref table if loading was successful, else None
        """
        try:
            df = pd.read_table(
                filename,
                names=[
                    "kgID",
                    "mRNA",
                    "spID",
                    "spDisplayID",
                    "geneSymbol",
                    "refseq",
                    "protAcc",
                    "description",
                    "rfamAcc",
                    "tRnaName",
                ],
                index_col=0,
                dtype=object,
            )
            return df
        except Exception as err:
            print(err)
            return None

    def _get_path_cytoBand_hg19(self):
        """ Get local path to cytoBand file for hg19 / GRCh37 from UCSC, downloading if necessary.

        Returns
        -------
        str
            path to cytoBand_hg19.txt.gz
        """
        return self._download_file(
            "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
            "cytoBand_hg19.txt.gz",
        )

    def _get_path_genetic_map_HapMapII_GRCh37(self):
        """ Get local path to HapMap Phase II genetic map for hg19 / GRCh37 (HapMapII),
        downloading if necessary.

        Returns
        -------
        str
            path to genetic_map_HapMapII_GRCh37.tar.gz

        References
        ----------
        ..[1] "The International HapMap Consortium (2007).  A second generation human haplotype
          map of over 3.1 million SNPs.  Nature 449: 851-861."
        ..[2] "The map was generated by lifting the HapMap Phase II genetic map from build 35 to
          GRCh37. The original map was generated using LDhat as described in the 2007 HapMap
          paper (Nature, 18th Sept 2007). The conversion from b35 to GRCh37 was achieved using
          the UCSC liftOver tool. Adam Auton, 08/12/2010"
        """
        return self._download_file(
            "ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/"
            "genetic_map_HapMapII_GRCh37.tar.gz",
            "genetic_map_HapMapII_GRCh37.tar.gz",
        )

    def _get_path_knownGene_hg19(self):
        """ Get local path to knownGene file for hg19 / GRCh37 from UCSC, downloading if necessary.

        Returns
        -------
        str
            path to knownGene_hg19.txt.gz
        """
        return self._download_file(
            "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz",
            "knownGene_hg19.txt.gz",
        )

    def _get_path_kgXref_hg19(self):
        """ Get local path to kgXref file for hg19 / GRCh37 from UCSC, downloading if necessary.

        Returns
        -------
        str
            path to kgXref_hg19.txt.gz
        """
        return self._download_file(
            "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz",
            "kgXref_hg19.txt.gz",
        )

    def _get_path_assembly_mapping_data(
        self, source_assembly, target_assembly, retries=10
    ):
        """ Get local path to assembly mapping data, downloading if necessary.

        Parameters
        ----------
        source_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap from
        target_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap to
        retries : int
            number of retries per chromosome to download assembly mapping data

        Returns
        -------
        str
            path to <source_assembly>_<target_assembly>.tar.gz

        References
        ----------
        ..[1] Ensembl, Assembly Information Endpoint,
          https://rest.ensembl.org/documentation/info/assembly_info
        ..[2] Ensembl, Assembly Map Endpoint,
          http://rest.ensembl.org/documentation/info/assembly_map

        """

        if not lineage.create_dir(self._resources_dir):
            return None

        chroms = [
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
        ]

        assembly_mapping_data = source_assembly + "_" + target_assembly
        destination = os.path.join(
            self._resources_dir, assembly_mapping_data + ".tar.gz"
        )

        if not os.path.exists(destination) or not self._all_chroms_in_tar(
            chroms, destination
        ):
            print("Downloading {}".format(os.path.relpath(destination)))

            try:
                with tarfile.open(destination, "w:gz") as out_tar:
                    for chrom in chroms:
                        file = chrom + ".json"

                        map_endpoint = (
                            "/map/human/"
                            + source_assembly
                            + "/"
                            + chrom
                            + "/"
                            + target_assembly
                            + "?"
                        )

                        # get assembly mapping data
                        response = None
                        retry = 0
                        while response is None and retry < retries:
                            response = self._ensembl_rest_client.perform_rest_action(
                                map_endpoint
                            )
                            retry += 1

                        if response is not None:
                            # open temp file, save json response to file, close temp file
                            with tempfile.NamedTemporaryFile(
                                delete=False, mode="w"
                            ) as f:
                                json.dump(response, f)

                            # add temp file to archive
                            out_tar.add(f.name, arcname=file)

                            # remove temp file
                            os.remove(f.name)
            except Exception as err:
                print(err)
                return None

        return destination

    def _all_chroms_in_tar(self, chroms, filename):
        try:
            with tarfile.open(filename, "r") as tar:
                members = tar.getnames()

            for chrom in chroms:
                if chrom + ".json" not in members:
                    return False
        except Exception as err:
            print(err)
            return False

        return True

    def _download_file(self, url, filename, compress=False, timeout=30):
        """ Download a file to the resources folder.

        Download data from `url`, save as `filename`, and optionally compress with gzip.

        Parameters
        ----------
        url : str
            URL to download data from
        filename : str
            name of file to save; if compress, ensure '.gz' is appended
        compress : bool
            compress with gzip
        timeout : int
            seconds for timeout of download request

        Returns
        -------
        str
            path to downloaded file, None if error
        """
        if not lineage.create_dir(self._resources_dir):
            return None

        if compress and filename[-3:] != ".gz":
            filename += ".gz"

        destination = os.path.join(self._resources_dir, filename)

        if not os.path.exists(destination):
            try:
                if compress:
                    open_func = gzip.open
                else:
                    open_func = open

                # get file if it hasn't already been downloaded
                # http://stackoverflow.com/a/7244263
                with urllib.request.urlopen(
                    url, timeout=timeout
                ) as response, open_func(destination, "wb") as f:
                    self._print_download_msg(destination)
                    data = response.read()  # a `bytes` object
                    f.write(data)
            except urllib.error.URLError as err:
                print(err)
                destination = None
                # try HTTP if an FTP error occurred
                if "ftp://" in url:
                    destination = self._download_file(
                        url.replace("ftp://", "http://"),
                        filename,
                        compress=compress,
                        timeout=timeout,
                    )
            except Exception as err:
                print(err)
                return None

        return destination

    @staticmethod
    def _print_download_msg(path):
        """ Print download message.

        Parameters
        ----------
        path : str
            path to file being downloaded
        """
        print("Downloading " + os.path.relpath(path))
