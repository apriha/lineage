"""
MIT License

Copyright (c) 2017 Andrew Riha

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

import gzip
import io
import os
import tarfile
import tempfile
from unittest.mock import mock_open, patch
import warnings

import pandas as pd

from lineage.resources import Resources
from tests import BaseLineageTestCase


class TestResources(BaseLineageTestCase):
    def _reset_resource(self):
        self.resource._genetic_map = {}
        self.resource._genetic_map_name = ""
        self.resource._cytoBand_hg19 = pd.DataFrame()
        self.resource._knownGene_hg19 = pd.DataFrame()
        self.resource._kgXref_hg19 = pd.DataFrame()

    def run(self, result=None):
        # set resources directory based on if downloads are being performed
        # https://stackoverflow.com/a/11180583

        self.resource = Resources()
        self._reset_resource()
        if self.downloads_enabled:
            self.resource._resources_dir = "resources"
            super().run(result)
        else:
            # use a temporary directory for test resource data
            with tempfile.TemporaryDirectory() as tmpdir:
                self.resource._resources_dir = tmpdir
                super().run(result)
                self.resource._resources_dir = "resources"

    def _generate_test_genetic_map_HapMapII_GRCh37_resource(self):
        filenames = [f"genetic_map_GRCh37_chr{chrom}.txt" for chrom in range(1, 23)]
        filenames.extend(
            [
                "genetic_map_GRCh37_chrX.txt",
                "genetic_map_GRCh37_chrX_par1.txt",
                "genetic_map_GRCh37_chrX_par2.txt",
                "README.txt",
            ]
        )

        # create compressed tar in memory
        tar_file = io.BytesIO()
        with tarfile.open(fileobj=tar_file, mode="w:gz") as out_tar:
            for filename in filenames:
                if filename != "README.txt":
                    chrom = filename[filename.find("chr") :].split(".")[0]
                    s = "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n"
                    s += f"{chrom}\t0\t0.0\t0.0\n"
                else:
                    s = "test"

                # add file to tar; https://stackoverflow.com/a/40392022
                data = s.encode()
                file = io.BytesIO(data)
                tar_info = tarfile.TarInfo(name=filename)
                tar_info.size = len(data)
                out_tar.addfile(tar_info, fileobj=file)

        mock = mock_open(read_data=tar_file.getvalue())
        with patch("urllib.request.urlopen", mock):
            self.resource._get_path_genetic_map_HapMapII_GRCh37()

    def test_get_genetic_map_HapMapII_GRCh37(self):
        def f():
            # mock download of test data
            self._generate_test_genetic_map_HapMapII_GRCh37_resource()
            return self.resource.get_genetic_map_HapMapII_GRCh37()

        genetic_map_HapMapII_GRCh37 = (
            self.resource.get_genetic_map_HapMapII_GRCh37()
            if self.downloads_enabled
            else f()
        )

        assert len(genetic_map_HapMapII_GRCh37) == 23

        # get already loaded resource
        genetic_map_HapMapII_GRCh37 = self.resource.get_genetic_map_HapMapII_GRCh37()
        assert len(genetic_map_HapMapII_GRCh37) == 23

        # get already loaded resource
        genetic_map = self.resource.get_genetic_map("HapMap2")
        assert len(genetic_map) == 23

    def _generate_test_genetic_map_1000G_GRCh37_resource(self):
        filenames = [f"CEU-{chrom}-final.txt.gz" for chrom in range(1, 23)]

        # create tar in memory
        tar_file = io.BytesIO()
        with tarfile.open(fileobj=tar_file, mode="w") as out_tar:
            for filename in filenames:
                s = "Position(bp)\tRate(cM/Mb)\tMap(cM)\tFiltered\n"
                s += f"   0\t0.0\t0.0\t0\n"

                # add file to tar; https://stackoverflow.com/a/40392022
                data = gzip.compress(s.encode())
                file = io.BytesIO(data)
                tar_info = tarfile.TarInfo(name=filename)
                tar_info.size = len(data)
                out_tar.addfile(tar_info, fileobj=file)

        mock = mock_open(read_data=tar_file.getvalue())
        with patch("urllib.request.urlopen", mock):
            self.resource._get_path_genetic_map_1000G_GRCh37("CEU")

    def test_get_genetic_map_1000G_GRCh37(self):
        def f():
            # mock download of test data
            self._generate_test_genetic_map_1000G_GRCh37_resource()
            return self.resource.get_genetic_map_1000G_GRCh37("CEU")

        genetic_map = (
            self.resource.get_genetic_map_1000G_GRCh37("CEU")
            if self.downloads_enabled
            else f()
        )

        assert len(genetic_map) == 22

        # get already loaded resource
        genetic_map = self.resource.get_genetic_map_1000G_GRCh37("CEU")
        assert len(genetic_map) == 22

        # get already loaded resource
        genetic_map = self.resource.get_genetic_map("CEU")
        assert len(genetic_map) == 22

    def test_invalid_genetic_map(self):
        # https://stackoverflow.com/a/46767037
        with self.assertLogs() as log:
            genetic_map = self.resource.get_genetic_map("test")
            self.assertEqual(len(log.output), 1)
            self.assertEqual(len(log.records), 1)
            self.assertIn("Invalid genetic map", log.output[0])
            self.assertEqual(len(genetic_map), 0)

    def _generate_test_cytoBand_hg19_resource(self):
        s = f"s\t0\t0\ts\ts\n" * 862
        mock = mock_open(read_data=gzip.compress(s.encode()))
        with patch("urllib.request.urlopen", mock):
            self.resource._get_path_cytoBand_hg19()

    def test_get_cytoBand_hg19(self):
        def f():
            # mock download of test data
            self._generate_test_cytoBand_hg19_resource()
            return self.resource.get_cytoBand_hg19()

        cytoBand_hg19 = (
            self.resource.get_cytoBand_hg19() if self.downloads_enabled else f()
        )

        assert len(cytoBand_hg19) == 862

        # get already loaded resource
        cytoBand_hg19 = self.resource.get_cytoBand_hg19()
        assert len(cytoBand_hg19) == 862

    def _generate_test_knownGene_hg19_resource(self):
        s = "s\ts\ts\t0\t0\t0\t0\t0\ts\ts\ts\ts\n" * 82960

        mock = mock_open(read_data=gzip.compress(s.encode()))
        with patch("urllib.request.urlopen", mock):
            self.resource._get_path_knownGene_hg19()

    def test_get_knownGene_hg19(self):
        def f():
            # mock download of test data
            self._generate_test_knownGene_hg19_resource()
            return self.resource.get_knownGene_hg19()

        knownGene_hg19 = (
            self.resource.get_knownGene_hg19() if self.downloads_enabled else f()
        )

        assert len(knownGene_hg19) == 82960

        # get already loaded resource
        knownGene_hg19 = self.resource.get_knownGene_hg19()
        assert len(knownGene_hg19) == 82960

    def _generate_test_kgXref_hg19_resource(self):
        s = "s\ts\ts\ts\ts\ts\ts\ts\n" * 82960

        mock = mock_open(read_data=gzip.compress(s.encode()))
        with patch("urllib.request.urlopen", mock):
            self.resource._get_path_kgXref_hg19()

    def test_get_kgXref_hg19(self):
        def f():
            # mock download of test data
            self._generate_test_kgXref_hg19_resource()
            return self.resource.get_kgXref_hg19()

        kgXref_hg19 = self.resource.get_kgXref_hg19() if self.downloads_enabled else f()

        assert len(kgXref_hg19) == 82960

        # get already loaded resource
        kgXref_hg19 = self.resource.get_kgXref_hg19()
        assert len(kgXref_hg19) == 82960

    def test_get_all_resources(self):
        def f():
            # mock download of test data for each resource
            self._generate_test_genetic_map_HapMapII_GRCh37_resource()
            self._generate_test_cytoBand_hg19_resource()
            self._generate_test_knownGene_hg19_resource()
            self._generate_test_kgXref_hg19_resource()

            return self.resource.get_all_resources()

        resources = self.resource.get_all_resources() if self.downloads_enabled else f()

        for k, v in resources.items():
            self.assertGreater(len(v), 0)

    def test_download_example_datasets(self):
        def f():
            with patch("urllib.request.urlopen", mock_open(read_data=b"")):
                return self.resource.download_example_datasets()

        paths = (
            self.resource.download_example_datasets() if self.downloads_enabled else f()
        )

        for path in paths:
            if path is None or not os.path.exists(path):
                warnings.warn("Example dataset(s) not currently available")
                return
