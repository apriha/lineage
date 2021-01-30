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
import pandas as pd
from snps import SNPs

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

    def test_load_path(self):
        ind = self.l.create_individual("test", "tests/input/generic.csv")
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps(), check_exact=True)

    def test_load_SNPs(self):
        s = SNPs("tests/input/generic.csv")
        ind = self.l.create_individual("test", s)
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps(), check_exact=True)

    def test_load_list_bytes(self):
        with open("tests/input/generic.csv", "rb") as f:
            data = f.read()
        ind = self.l.create_individual("test", [SNPs(), data])
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps(), check_exact=True)

    def test_load_resource_output_dirs(self):
        ind = self.l.create_individual(
            "test",
            "tests/input/generic.csv",
            output_dir="output1",
            resources_dir="resources1",
        )
        self.assertEqual(self.l._output_dir, "output")
        self.assertEqual(self.l._resources_dir, "resources")
        self.assertEqual(ind._output_dir, "output1")
        pd.testing.assert_frame_equal(ind.snps, self.generic_snps(), check_exact=True)
