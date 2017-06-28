""" Class for downloading and loading required external resources. """

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
        """ Load UCSC cytoBandIdeo table.

        http://genome.ucsc.edu/cgi-bin/hgTables

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
            path to hapmap_36.tar.gz

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
                                self._print_download_msg(filename)

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
            path to hapmap_36.tar.gz

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
                self._print_download_msg(filename)
                # get file if it hasn't already been downloaded [SO-01]
                with urllib.request.urlopen(url) as response, open(destination, 'wb') as out_file:
                    data = response.read()  # a `bytes` object
                    out_file.write(data)
            except Exception as err:
                print(err)
                return None

        return destination

    @staticmethod
    def _print_download_msg(filename):
        """ Print download message.

        Parameters
        ----------
        filename : str
            filename to print
        """
        print('Downloading ' + filename + '...')


"""
Stack Overflow Attributions
---------------------------

[SO-01] "Download file from web in Python 3"
        http://stackoverflow.com/q/7243750
        Bo Milanovich : http://stackoverflow.com/users/647897/bo-milanovich
        http://stackoverflow.com/a/7244263
        Oleh Prypin : http://stackoverflow.com/users/241039/oleh-prypin

"""
