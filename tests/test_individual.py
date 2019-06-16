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

from tests import BaseLineageTestCase


class TestIndividual(BaseLineageTestCase):
    def test_name(self):
        ind = self.l.create_individual("test")
        assert ind.name == "test"

    def test_get_var_name(self):
        ind = self.l.create_individual("test a. name")
        assert ind.get_var_name() == "test_a__name"

    def test___repr__(self):
        ind = self.l.create_individual("test")
        assert "Individual('test')" == ind.__repr__()
