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
import hashlib
import itertools
import json
import os
import tarfile
import tempfile
import urllib.error
import urllib.request
import zlib

from atomicwrites import atomic_write
import numpy as np
import pandas as pd

from lineage.ensembl import EnsemblRestClient
from lineage.utils import create_dir, Singleton


class Resources(metaclass=Singleton):
    """ Object used to manage resources required by `lineage`. """

    def __init__(self, resources_dir="resources"):
        """ Initialize a ``Resources`` object.

        Parameters
        ----------
        resources_dir : str
            name / path of resources directory
        """
        self._resources_dir = os.path.abspath(resources_dir)
        self._genetic_map_HapMapII_GRCh37 = {}
        self._cytoBand_hg19 = pd.DataFrame()
        self._knownGene_hg19 = pd.DataFrame()
        self._kgXref_hg19 = pd.DataFrame()
        self._ensembl_rest_client = EnsemblRestClient()
        self._reference_sequences = {}

    def get_genetic_map_HapMapII_GRCh37(self):
        """ Get International HapMap Consortium HapMap Phase II genetic map for Build 37.

        Returns
        -------
        dict
            dict of pandas.DataFrame HapMapII genetic maps if loading was successful, else {}
        """
        if not self._genetic_map_HapMapII_GRCh37:
            self._genetic_map_HapMapII_GRCh37 = self._load_genetic_map(
                self._get_path_genetic_map_HapMapII_GRCh37()
            )

        return self._genetic_map_HapMapII_GRCh37

    def get_reference_sequences(
        self,
        assembly="GRCh37",
        chroms=(
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
        ),
    ):
        """ Get Homo sapiens reference sequences for `chroms` of `assembly`.

        Notes
        -----
        This function can download over 800MB of data for each assembly.

        Parameters
        ----------
        assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            reference sequence assembly
        chroms : list of str
            reference sequence chromosomes

        Returns
        -------
        dict
            dict of ReferenceSequence, else {}
        """
        valid_assemblies = ["NCBI36", "GRCh37", "GRCh38"]

        if assembly not in valid_assemblies:
            print("Invalid assembly")
            return {}

        if not self._reference_chroms_available(assembly, chroms):
            self._reference_sequences[assembly] = self._create_reference_sequences(
                *self._get_paths_reference_sequences(assembly=assembly, chroms=chroms)
            )

        return self._reference_sequences[assembly]

    def _reference_chroms_available(self, assembly, chroms):
        if assembly in self._reference_sequences:
            for chrom in chroms:
                if chrom not in self._reference_sequences[assembly]:
                    return False
            return True
        else:
            return False

    def get_cytoBand_hg19(self):
        """ Get UCSC cytoBand table for Build 37.

        Returns
        -------
        pandas.DataFrame
            cytoBand table if loading was successful, else empty DataFrame
        """
        if self._cytoBand_hg19.empty:
            self._cytoBand_hg19 = self._load_cytoBand(self._get_path_cytoBand_hg19())

        return self._cytoBand_hg19

    def get_knownGene_hg19(self):
        """ Get UCSC knownGene table for Build 37.

        Returns
        -------
        pandas.DataFrame
            knownGene table if loading was successful, else empty DataFrame
        """
        if self._knownGene_hg19.empty:
            self._knownGene_hg19 = self._load_knownGene(self._get_path_knownGene_hg19())

        return self._knownGene_hg19

    def get_kgXref_hg19(self):
        """ Get UCSC kgXref table for Build 37.

        Returns
        -------
        pandas.DataFrame
            kgXref table if loading was successful, else empty DataFrame
        """
        if self._kgXref_hg19.empty:
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
            dict of json assembly mapping data if loading was successful, else {}
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
        paths : list of str or empty str
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
                with atomic_write(gzip_path, mode="wb", overwrite=True) as f:
                    self._write_data_to_gzip(f, data)
        except Exception as err:
            print(err)

        return paths

    def get_all_resources(self):
        """ Get / download all resources (except reference sequences) used throughout `lineage`.

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

    def get_all_reference_sequences(self, **kwargs):
        """ Get Homo sapiens reference sequences for Builds 36, 37, and 38 from Ensembl.

        Notes
        -----
        This function can download over 2.5GB of data.

        Returns
        -------
        dict
            dict of ReferenceSequence, else {}
        """
        for assembly in ("NCBI36", "GRCh37", "GRCh38"):
            self.get_reference_sequences(assembly=assembly, **kwargs)
        return self._reference_sequences

    @staticmethod
    def _write_data_to_gzip(f, data):
        """ Write `data` to `f` in `gzip` format.

        Parameters
        ----------
        f : file object opened with `mode="wb"`
        data : `bytes` object
        """
        with gzip.open(f, "wb") as f_gzip:
            f_gzip.write(data)

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
            dict of pandas.DataFrame genetic maps if loading was successful, else {}

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
            return {}

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
            dict of assembly maps if loading was successful, else {}

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
            return {}

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
            cytoBand table if loading was successful, else empty DataFrame

        References
        ----------
        ..[1] Ryan Dale, GitHub Gist,
          https://gist.github.com/daler/c98fc410282d7570efc3#file-ideograms-py
        """
        try:
            # adapted from chromosome plotting code (see [1]_)
            df = pd.read_csv(
                filename, names=["chrom", "start", "end", "name", "gie_stain"], sep="\t"
            )
            df["chrom"] = df["chrom"].str[3:]
            return df
        except Exception as err:
            print(err)
            return pd.DataFrame()

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
            knownGene table if loading was successful, else empty DataFrame
        """
        try:
            df = pd.read_csv(
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
                sep="\t",
            )
            df["chrom"] = df["chrom"].str[3:]
            return df
        except Exception as err:
            print(err)
            return pd.DataFrame()

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
            kgXref table if loading was successful, else empty DataFrame
        """
        try:
            df = pd.read_csv(
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
                sep="\t",
                dtype=object,
            )
            return df
        except Exception as err:
            print(err)
            return pd.DataFrame()

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

    def _get_paths_reference_sequences(
        self, sub_dir="reference_sequences", assembly="GRCh37", chroms=()
    ):
        """ Get local paths to Homo sapiens reference sequences from Ensembl.

        Notes
        -----
        This function can download over 800MB of data for each assembly.

        Parameters
        ----------
        sub_dir : str
            directory under resources to store reference sequence data
        assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            reference sequence assembly
        chroms : list of str
            reference sequence chromosomes

        Returns
        -------
        assembly : str
            reference sequence assembly
        chroms : list of str
            reference sequence chromosomes
        urls : list of str
            urls to Ensembl reference sequences
        paths : list of str
            paths to local reference sequences

        References
        ----------
        ..[1] Daniel R. Zerbino, Premanand Achuthan, Wasiu Akanni, M. Ridwan Amode,
          Daniel Barrell, Jyothish Bhai, Konstantinos Billis, Carla Cummins, Astrid Gall,
          Carlos García Giro´n, Laurent Gil, Leo Gordon, Leanne Haggerty, Erin Haskell,
          Thibaut Hourlier, Osagie G. Izuogu, Sophie H. Janacek, Thomas Juettemann,
          Jimmy Kiang To, Matthew R. Laird, Ilias Lavidas, Zhicheng Liu, Jane E. Loveland,
          Thomas Maurel, William McLaren, Benjamin Moore, Jonathan Mudge, Daniel N. Murphy,
          Victoria Newman, Michael Nuhn, Denye Ogeh, Chuang Kee Ong, Anne Parker,
          Mateus Patricio, Harpreet Singh Riat, Helen Schuilenburg, Dan Sheppard,
          Helen Sparrow, Kieron Taylor, Anja Thormann, Alessandro Vullo, Brandon Walts,
          Amonida Zadissa, Adam Frankish, Sarah E. Hunt, Myrto Kostadima, Nicholas Langridge,
          Fergal J. Martin, Matthieu Muffato, Emily Perry, Magali Ruffier, Dan M. Staines,
          Stephen J. Trevanion, Bronwen L. Aken, Fiona Cunningham, Andrew Yates, Paul Flicek
          Ensembl 2018.
          PubMed PMID: 29155950.
          doi:10.1093/nar/gkx1098
        ..[2] NCBI 36, Oct 2005, Ensembl release 54, Database version: 54.36p
        ..[3] GRCh37.p13 (Genome Reference Consortium Human Reference 37),
          INSDC Assembly GCA_000001405.14, Feb 2009, Ensembl GRCh37 release 96, Database
          version: 96.37
        ..[4] GRCh38.p12 (Genome Reference Consortium Human Build 38),
          INSDC Assembly GCA_000001405.27, Dec 2013, Ensembl release 96, Database
          version: 96.38
        """
        release = ""

        # https://www.biostars.org/p/374149/#374219
        if assembly == "GRCh37":
            base = "ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/"
        elif assembly == "NCBI36":
            base = "ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/"
            release = "54."
        elif assembly == "GRCh38":
            base = "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/"
        else:
            return ("", [], [], [])

        filenames = [
            "Homo_sapiens.{}.{}dna.chromosome.{}.fa.gz".format(assembly, release, chrom)
            for chrom in chroms
        ]

        urls = ["{}{}".format(base, filename) for filename in filenames]

        local_filenames = [
            "{}{}{}{}{}".format(sub_dir, os.sep, assembly, os.sep, filename)
            for filename in filenames
        ]

        return (
            assembly,
            chroms,
            urls,
            list(map(self._download_file, urls, local_filenames)),
        )

    def _create_reference_sequences(self, assembly, chroms, urls, paths):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        seqs = {}

        for i, path in enumerate(paths):
            if not path:
                continue

            d = {}
            d["ID"] = chroms[i]
            d["url"] = urls[i]
            d["path"] = os.path.relpath(path)
            d["assembly"] = "B{}".format(assembly[-2:])
            d["species"] = "Homo sapiens"
            d["taxonomy"] = "x"
            seqs[chroms[i]] = ReferenceSequence(**d)

        return seqs

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

        if not create_dir(self._resources_dir):
            return ""

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
                self._download_assembly_mapping_data(
                    destination, chroms, source_assembly, target_assembly, retries
                )
            except Exception as err:
                print(err)
                return ""

        return destination

    def _download_assembly_mapping_data(
        self, destination, chroms, source_assembly, target_assembly, retries
    ):
        with atomic_write(destination, mode="wb", overwrite=True) as f:
            with tarfile.open(fileobj=f, mode="w:gz") as out_tar:
                for chrom in chroms:
                    file = chrom + ".json"

                    map_endpoint = "/map/human/{}/{}/{}?".format(
                        source_assembly, chrom, target_assembly
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
                        ) as f_tmp:
                            json.dump(response, f_tmp)

                        # add temp file to archive
                        out_tar.add(f_tmp.name, arcname=file)

                        # remove temp file
                        os.remove(f_tmp.name)

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
            path to downloaded file, empty str if error
        """
        if compress and filename[-3:] != ".gz":
            filename += ".gz"

        destination = os.path.join(self._resources_dir, filename)

        if not create_dir(os.path.relpath(os.path.dirname(destination))):
            return ""

        if not os.path.exists(destination):
            try:
                # get file if it hasn't already been downloaded
                # http://stackoverflow.com/a/7244263
                with urllib.request.urlopen(
                    url, timeout=timeout
                ) as response, atomic_write(destination, mode="wb") as f:
                    self._print_download_msg(destination)
                    data = response.read()  # a `bytes` object

                    if compress:
                        self._write_data_to_gzip(f, data)
                    else:
                        f.write(data)
            except urllib.error.URLError as err:
                print(err)
                destination = ""
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
                return ""

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


class ReferenceSequence:
    """ Object used to represent and interact with a reference sequence. """

    def __init__(self, ID="", url="", path="", assembly="", species="", taxonomy=""):
        """ Initialize a ``ReferenceSequence`` object.

        Parameters
        ----------
        ID : str
            reference sequence chromosome
        url : str
            url to Ensembl reference sequence
        path : str
            path to local reference sequence
        assembly : str
            reference sequence assembly / build (e.g., "B37")
        species : str
            reference sequence species
        taxonomy : str
            reference sequence taxonomy

        References
        ----------
        ..[1] The Variant Call Format (VCF) Version 4.2 Specification, 8 Mar 2019,
          https://samtools.github.io/hts-specs/VCFv4.2.pdf
        """
        self._ID = ID
        self._url = url
        self._path = path
        self._assembly = assembly
        self._species = species
        self._taxonomy = taxonomy
        self._sequence = np.array([], dtype=np.uint8)
        self._md5 = ""
        self._start = 0
        self._end = 0
        self._length = 0

    def __repr__(self):
        return "ReferenceSequence(assembly={!r}, ID={!r})".format(
            self._assembly, self._ID
        )

    @property
    def ID(self):
        """ Get reference sequence chromosome.

        Returns
        -------
        str
        """
        return self._ID

    @property
    def chrom(self):
        """ Get reference sequence chromosome.

        Returns
        -------
        str
        """
        return self._ID

    @property
    def url(self):
        """ Get URL to Ensembl reference sequence.

        Returns
        -------
        str
        """
        return self._url

    @property
    def path(self):
        """ Get path to local reference sequence.

        Returns
        -------
        str
        """
        return self._path

    @property
    def assembly(self):
        """ Get reference sequence assembly.

        Returns
        -------
        str
        """
        return self._assembly

    @property
    def species(self):
        """ Get reference sequence species.

        Returns
        -------
        str
        """
        return self._species

    @property
    def taxonomy(self):
        """ Get reference sequence taxonomy.

        Returns
        -------
        str
        """
        return self._taxonomy

    @property
    def sequence(self):
        """ Get reference sequence.

        Returns
        -------
        np.array(dtype=np.uint8)
        """
        self._load_sequence()
        return self._sequence

    @property
    def md5(self):
        """ Get reference sequence MD5 hash.

        Returns
        -------
        str
        """
        self._load_sequence()
        return self._md5

    @property
    def start(self):
        """ Get reference sequence start position (1-based).

        Returns
        -------
        int
        """
        self._load_sequence()
        return self._start

    @property
    def end(self):
        """ Get reference sequence end position (1-based).

        Returns
        -------
        int
        """
        self._load_sequence()
        return self._end

    @property
    def length(self):
        """ Get reference sequence length.

        Returns
        -------
        int
        """
        self._load_sequence()
        return self._sequence.size

    def clear(self):
        """ Clear reference sequence. """
        self._sequence = np.array([], dtype=np.uint8)
        self._md5 = ""
        self._start = 0
        self._end = 0
        self._length = 0

    def _load_sequence(self):
        if not self._sequence.size:
            # decompress and read file
            with gzip.open(self._path, "rb") as f:
                data = f.read()

            # convert bytes to str
            data = str(data, encoding="utf-8", errors="strict")

            data = data.split("\n")

            self._start, self._end = self._parse_first_line(data[0])

            # convert str (FASTA sequence) to bytes
            data = bytearray("".join(data[1:]), encoding="utf-8", errors="strict")

            # get MD5 of FASTA sequence
            self._md5 = hashlib.md5(data).hexdigest()

            # store FASTA sequence as `np.uint8` array
            self._sequence = np.array(data, dtype=np.uint8)

    def _parse_first_line(self, first_line):
        items = first_line.split(":")
        return (
            int(items[items.index(self._ID) + 1]),
            int(items[items.index(self._ID) + 2]),
        )
