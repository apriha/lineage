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
        self.assert_frame_equal_with_string_index(
            ind.snps, self.generic_snps(), check_exact=True
        )

    def test_load_SNPs(self):
        s = SNPs("tests/input/generic.csv")
        ind = self.l.create_individual("test", s)
        self.assert_frame_equal_with_string_index(
            ind.snps, self.generic_snps(), check_exact=True
        )

    def test_load_list_bytes(self):
        with open("tests/input/generic.csv", "rb") as f:
            data = f.read()
        ind = self.l.create_individual("test", [SNPs(), data])
        self.assert_frame_equal_with_string_index(
            ind.snps, self.generic_snps(), check_exact=True
        )

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
        self.assert_frame_equal_with_string_index(
            ind.snps, self.generic_snps(), check_exact=True
        )
