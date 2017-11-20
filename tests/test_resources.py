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

import urllib.request


class TestResources(object):

    @staticmethod
    def url_exists(url):
        # test existence of external resource without downloading

        # https://stackoverflow.com/a/16778473
        try:
            with urllib.request.urlopen(url):
                return True
        except:
            return False

    def test__get_url_cytoband_h37(self, resource):
        assert self.url_exists(resource._get_url_cytoband_h37())

    def test__get_url_hapmap_h37(self, resource):
        assert self.url_exists(resource._get_url_hapmap_h37())

    def test__get_knownGene_h37(self, resource):
        assert self.url_exists(resource._get_url_knownGene_h37())

    def test__kgXref_h37(self, resource):
        assert self.url_exists(resource._get_url_kgXref_h37())

    def test__get_url_dataset_user662_304(self, resource):
        assert self.url_exists(resource._get_url_dataset_user662_304())

    def test__get_url_dataset_user662_340(self, resource):
        assert self.url_exists(resource._get_url_dataset_user662_340())

    def test__get_url_dataset_user662_341(self, resource):
        assert self.url_exists(resource._get_url_dataset_user662_341())

    def test__get_url_dataset_user663_305(self, resource):
        assert self.url_exists(resource._get_url_dataset_user663_305())

    def test__get_url_dataset_user4583_3482(self, resource):
        assert self.url_exists(resource._get_url_dataset_user4583_3482())

    def test__get_url_dataset_user4584_3483(self, resource):
        assert self.url_exists(resource._get_url_dataset_user4584_3483())
