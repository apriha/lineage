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
..[5] hg18 (NCBI36): Engineering effort led by Fan Hsu; QA effort led by Ann Zweig,
  https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg18

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

import ftplib
import os
import tarfile
import tempfile
import urllib.request

import pandas as pd

import lineage


class Resources(object):
    """ Object used to manage resources required by `lineage`. """

    def __init__(self, resources_dir):
        self._resources_dir = os.path.abspath(resources_dir)
        self._hapmap_h36 = None
        self._hapmap_h37 = None
        self._cytoband_h36 = None
        self._cytoband_h37 = None
        self._knownGene_h37 = None
        self._kgXref_h37 = None

    def get_hapmap_h36(self):
        if self._hapmap_h36 is None:
            self._hapmap_h36 = self._load_hapmap(self._get_path_hapmap_h36())

        return self._hapmap_h36

    def get_hapmap_h37(self):
        if self._hapmap_h37 is None:
            self._hapmap_h37 = self._load_hapmap(self._get_path_hapmap_h37())

        return self._hapmap_h37

    def get_cytoband_h36(self):
        if self._cytoband_h36 is None:
            self._cytoband_h36 = self._load_cytoband(self._get_path_cytoband_h36())

        return self._cytoband_h36

    def get_cytoband_h37(self):
        if self._cytoband_h37 is None:
            self._cytoband_h37 = self._load_cytoband(self._get_path_cytoband_h37())

        return self._cytoband_h37

    def get_knownGene_h37(self):
        if self._knownGene_h37 is None:
            self._knownGene_h37 = self._load_knownGene(self._get_path_knownGene_h37())

        return self._knownGene_h37

    def get_kgXref_h37(self):
        if self._kgXref_h37 is None:
            self._kgXref_h37 = self._load_kgXref(self._get_path_kgXref_h37())

        return self._kgXref_h37

    @staticmethod
    def _load_hapmap(filename):
        """ Load HapMap.

        Parameters
        ----------
        filename : str
            path to compressed archive with HapMap data

        Returns
        -------
        hapmap : dict
            dict of pandas.DataFrame HapMap tables if loading was successful, else None

        """
        if filename is None:
            return None

        try:
            hapmap = {}

            if '36' in filename:
                with tarfile.open(filename, 'r') as tar:
                    # http://stackoverflow.com/a/2018576
                    for member in tar.getmembers():
                        if 'genetic_map' in member.name:
                            df = pd.read_csv(tar.extractfile(member), sep=' ')
                            df = df.rename(columns={'position': 'pos',
                                                    'COMBINED_rate(cM/Mb)': 'rate',
                                                    'Genetic_Map(cM)': 'map'})
                            start_pos = member.name.index('chr') + 3
                            end_pos = member.name.index('_b36')
                            hapmap[member.name[start_pos:end_pos]] = df
            elif '37' in filename:
                with tarfile.open(filename, 'r') as tar:
                    for member in tar.getmembers():
                        if 'genetic_map' in member.name:
                            df = pd.read_csv(tar.extractfile(member), sep='\t')
                            df = df.rename(columns={'Position(bp)': 'pos',
                                                    'Rate(cM/Mb)': 'rate',
                                                    'Map(cM)': 'map'})
                            del df['Chromosome']
                            start_pos = member.name.index('chr') + 3
                            end_pos = member.name.index('.')
                            hapmap[member.name[start_pos:end_pos]] = df
            else:
                hapmap = None
        except Exception as err:
            print(err)
            hapmap = None

        return hapmap

    @staticmethod
    def _load_cytoband(filename):
        """ Load UCSC cytoBand table.

        Parameters
        ----------
        filename : str
            path to cytoband file

        Returns
        -------
        df : pandas.DataFrame
            cytoband data if loading was successful, else None

        References
        ----------
        ..[1] Ryan Dale, GitHub Gist,
          https://gist.github.com/daler/c98fc410282d7570efc3#file-ideograms-py

        """
        if filename is None:
            return None

        try:
            # adapted from chromosome plotting code (see [1]_)
            df = pd.read_table(filename, names=['chrom', 'start', 'end', 'name', 'gie_stain'])
            df['chrom'] = df['chrom'].str[3:]
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
            knownGene data if loading was successful, else None
        """
        if filename is None:
            return None

        try:
            df = pd.read_table(filename, names=['name', 'chrom', 'strand', 'txStart', 'txEnd',
                                                'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts',
                                                'exonEnds', 'proteinID', 'alignID'], index_col=0)
            df['chrom'] = df['chrom'].str[3:]
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
            kgXref data if loading was successful, else None
        """
        if filename is None:
            return None

        try:
            df = pd.read_table(filename, names=['kgID', 'mRNA', 'spID', 'spDisplayID',
                                                'geneSymbol', 'refseq', 'protAcc',
                                                'description', 'rfamAcc', 'tRnaName'], index_col=0,
                               dtype=object)
            return df
        except Exception as err:
            print(err)
            return None

    def _get_path_cytoband_h36(self):
        """ Get local path to cytoBand file for hg18 / NCBI36 from UCSC, downloading if necessary.

        Returns
        -------
        str
            path to cytoband_h36.txt.gz

        """

        return self._download_file(
            'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz',
            'cytoband_h36.txt.gz')

    def _get_path_cytoband_h37(self):
        """ Get local path to cytoBand file for hg19 / GRCh37 from UCSC, downloading if necessary.

        Returns
        -------
        str
            path to cytoband_h37.txt.gz

        """

        return self._download_file(
            'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz',
            'cytoband_h37.txt.gz')

    def _get_path_hapmap_h36(self):
        """ Get local path to HapMap for hg18 / NCBI36, downloading if necessary.

        Returns
        -------
        str
            path to hapmap_h36.tar.gz

        References
        ----------
        ..[1] "The International HapMap Consortium (2007).  A second generation human haplotype
          map of over 3.1 million SNPs.  Nature 449: 851-861."

        """
        if not lineage.dir_exists(self._resources_dir):
            return None

        hapmap = 'hapmap_h36'
        destination = os.path.join(self._resources_dir, hapmap + '.tar.gz')

        if not os.path.exists(destination):
            try:
                # make FTP connection to NCBI
                with ftplib.FTP('ftp.ncbi.nlm.nih.gov') as ftp:
                    ftp.login()
                    ftp.cwd('hapmap/recombination/2008-03_rel22_B36/rates')

                    # download each HapMap file and add to compressed tar
                    with tarfile.open(destination, "w:gz") as out_tar:
                        for filename in ftp.nlst():
                            if '.txt' in filename:
                                path = os.path.join(destination, hapmap, filename)
                                self._print_download_msg(path)

                                # open temp file, download HapMap file, close temp file
                                with tempfile.NamedTemporaryFile(delete=False) as fp:
                                    ftp.retrbinary('RETR ' + filename, fp.write)

                                # add temp file to archive
                                out_tar.add(fp.name, arcname=os.path.join(hapmap, filename))

                                # remove temp file
                                os.remove(fp.name)
                    ftp.quit()
            except Exception as err:
                print(err)
                return None

        return destination

    def _get_path_hapmap_h37(self):
        """ Get local path to HapMap for hg19 / GRCh37, downloading if necessary.

        Returns
        -------
        str
            path to hapmap_h37.tar.gz

        References
        ----------
        ..[1] "The map was generated by lifting the HapMap Phase II genetic map from build 35 to
          GRCh37. The original map was generated using LDhat as described in the 2007 HapMap
          paper (Nature, 18th Sept 2007). The conversion from b35 to GRCh37 was achieved using
          the UCSC liftOver tool. Adam Auton, 08/12/2010"

        """

        return self._download_file(
            'ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37'
            '/genetic_map_HapMapII_GRCh37.tar.gz', 'hapmap_h37.tar.gz')

    def _get_path_knownGene_h37(self):
        """ Get local path to knownGene file for hg19 / GRCh37 from UCSC, downloading if necessary.

        Returns
        -------
        str
            path to knownGene_h37.txt.gz

        """

        return self._download_file(
            'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz',
            'knownGene_h37.txt.gz')


    def _get_path_kgXref_h37(self):
        """ Get local path to kgXref file for hg19 / GRCh37 from UCSC, downloading if necessary.

        Returns
        -------
        str
            path to kgXref_h37.txt.gz

        """

        return self._download_file(
            'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz',
            'kgXref_h37.txt.gz')

    def _download_file(self, url, filename):
        """ Download a file.

        Download data from `url` and save as `filename`.

        Parameters
        ----------
        url : str
        filename : str

        Returns
        -------
        str
            path to downloaded file, None if error

        """
        if not lineage.dir_exists(self._resources_dir):
            return None

        destination = os.path.join(self._resources_dir, filename)

        if not os.path.exists(destination):
            try:
                self._print_download_msg(destination)
                # get file if it hasn't already been downloaded [SO-01]
                with urllib.request.urlopen(url) as response, open(destination, 'wb') as out_file:
                    data = response.read()  # a `bytes` object
                    out_file.write(data)
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
        print('Downloading ' + os.path.relpath(path))


"""
Stack Overflow Attributions
---------------------------

[SO-01] "Download file from web in Python 3"
        http://stackoverflow.com/q/7243750
        Bo Milanovich : http://stackoverflow.com/users/647897/bo-milanovich
        http://stackoverflow.com/a/7244263
        Oleh Prypin : http://stackoverflow.com/users/241039/oleh-prypin

"""
