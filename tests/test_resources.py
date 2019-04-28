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

import os
import warnings

from lineage import Resources, EnsemblRestClient
from tests import BaseLineageTestCase


class TestResources(BaseLineageTestCase):
    def setUp(self):
        self.resource = Resources(resources_dir="resources")
        self.resource_assembly_mapping = Resources(
            resources_dir="resources", ensembl_rest_client=EnsemblRestClient()
        )
        self.del_output_dir_helper()

    def test_get_genetic_map_HapMapII_GRCh37(self):
        genetic_map_HapMapII_GRCh37 = self.resource.get_genetic_map_HapMapII_GRCh37()
        assert len(genetic_map_HapMapII_GRCh37) == 23

    def test_get_cytoBand_hg19(self):
        cytoBand_hg19 = self.resource.get_cytoBand_hg19()
        assert len(cytoBand_hg19) == 862

    def test_get_knownGene_hg19(self):
        knownGene_hg19 = self.resource.get_knownGene_hg19()
        assert len(knownGene_hg19) == 82960

    def test_get_kgXref_hg19(self):
        kgXref_hg19 = self.resource.get_kgXref_hg19()
        assert len(kgXref_hg19) == 82960

    def test_get_assembly_mapping_data_no_EnsemblRestClient(self):
        if os.path.exists("resources/NCBI36_GRCh37.tar.gz"):
            os.remove("resources/NCBI36_GRCh37.tar.gz")
        assembly_mapping_data = self.resource.get_assembly_mapping_data(
            "NCBI36", "GRCh37"
        )
        assert assembly_mapping_data is None

    def test_get_assembly_mapping_data_bad_tar(self):
        with open("resources/NCBI36_GRCh37.tar.gz", "w"):
            pass
        assembly_mapping_data = self.resource_assembly_mapping.get_assembly_mapping_data(
            "NCBI36", "GRCh37"
        )
        assert len(assembly_mapping_data) == 25

    def test_get_assembly_mapping_data(self):
        assembly_mapping_data = self.resource_assembly_mapping.get_assembly_mapping_data(
            "NCBI36", "GRCh37"
        )
        assert len(assembly_mapping_data) == 25

    def test_get_all_resources(self):
        resources = self.resource_assembly_mapping.get_all_resources()
        for k, v in resources.items():
            if v is None:
                assert False
        assert True

    def test__all_chroms_in_tar(self):
        assert not self.resource_assembly_mapping._all_chroms_in_tar(
            ["PAR"], "resources/NCBI36_GRCh37.tar.gz"
        )

    def test_get_assembly_mapping_data_invalid_dir(self):
        self.resource_assembly_mapping._resources_dir = None
        assembly_mapping_data = self.resource_assembly_mapping.get_assembly_mapping_data(
            "NCBI36", "GRCh37"
        )
        assert assembly_mapping_data is None

    def test_download_example_datasets(self):
        paths = self.resource.download_example_datasets()

        for path in paths:
            if path is None or not os.path.exists(path):
                warnings.warn("Example dataset(s) not currently available")
                return

        assert True

    def test__load_genetic_map_None(self):
        result = self.resource._load_genetic_map(None)
        assert result is None

    def test__load_cytoBand_None(self):
        result = self.resource._load_cytoBand(None)
        assert result is None

    def test__load_knownGene_None(self):
        result = self.resource._load_knownGene(None)
        assert result is None

    def test__load_kgXref_None(self):
        result = self.resource._load_kgXref(None)
        assert result is None

    def test__load_assembly_mapping_data_None(self):
        result = self.resource._load_assembly_mapping_data(None)
        assert result is None

    def test__download_file_compress(self):
        result = self.resource._download_file("", "", compress=True)
        assert result is None

    def test__download_file_invalid_dir(self):
        self.resource._resources_dir = None
        result = self.resource._download_file("", "")
        assert result is None
