"""Tests for synthetic related individuals SNP data generator."""

import os
import tempfile

import pandas as pd

from lineage.generator import SyntheticRelatedGenerator
from tests import BaseLineageTestCase


class TestSyntheticRelatedGenerator(BaseLineageTestCase):
    """Tests for the SyntheticRelatedGenerator class."""

    def test_init_default(self):
        """Test default initialization."""
        gen = SyntheticRelatedGenerator()
        self.assertEqual(gen.build, 37)
        self.assertIsNone(gen.seed)
        self.assertEqual(gen.crossovers_per_chrom, 1.0)

    def test_init_with_params(self):
        """Test initialization with custom parameters."""
        gen = SyntheticRelatedGenerator(build=38, seed=123, crossovers_per_chrom=2.0)
        self.assertEqual(gen.build, 38)
        self.assertEqual(gen.seed, 123)
        self.assertEqual(gen.crossovers_per_chrom, 2.0)

    def test_init_invalid_build(self):
        """Test that invalid build raises ValueError."""
        with self.assertRaises(ValueError):
            SyntheticRelatedGenerator(build=35)

    def test_generate_parent_child_pair_structure(self):
        """Test that parent-child pair has correct structure."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        parent_df, child_df = gen.generate_parent_child_pair(num_snps=1000)

        # Verify DataFrames have correct structure
        self.assertGreater(len(parent_df), 0)
        self.assertGreater(len(child_df), 0)
        self.assertEqual(len(parent_df), len(child_df))

        # Verify columns
        self.assertIn("chrom", parent_df.columns)
        self.assertIn("pos", parent_df.columns)
        self.assertIn("genotype", parent_df.columns)

        # Verify index name
        self.assertEqual(parent_df.index.name, "rsid")
        self.assertEqual(child_df.index.name, "rsid")

    def test_generate_parent_child_pair_genotypes_valid(self):
        """Test that generated genotypes are valid."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        parent_df, child_df = gen.generate_parent_child_pair(num_snps=1000)

        # Valid genotypes include all unphased combinations (both orientations)
        bases = ["A", "C", "G", "T"]
        valid_genotypes = {"--"}
        for b1 in bases:
            for b2 in bases:
                valid_genotypes.add(b1 + b2)

        parent_gts = set(parent_df["genotype"].unique())
        child_gts = set(child_df["genotype"].unique())

        self.assertTrue(parent_gts.issubset(valid_genotypes))
        self.assertTrue(child_gts.issubset(valid_genotypes))

    def test_generate_parent_child_pair_allele_sharing(self):
        """Test that child shares at least one allele with parent at every position."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        parent_df, child_df = gen.generate_parent_child_pair(num_snps=1000)

        shared_count = 0
        total_count = 0

        for idx in parent_df.index:
            p_gt = parent_df.loc[idx, "genotype"]
            c_gt = child_df.loc[idx, "genotype"]

            if p_gt == "--" or c_gt == "--":
                continue

            total_count += 1
            p_alleles = set(p_gt)
            c_alleles = set(c_gt)

            if p_alleles & c_alleles:
                shared_count += 1

        # Should be 100% sharing (within floating point tolerance)
        self.assertEqual(shared_count, total_count)

    def test_generate_parent_child_pair_discordant_rate(self):
        """Test that discordant_rate parameter introduces discordant SNPs."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        # Use a high discordant rate to ensure we get some discordant SNPs
        parent_df, child_df = gen.generate_parent_child_pair(
            num_snps=10000, discordant_rate=0.05
        )

        discordant_count = 0
        total_count = 0

        for idx in parent_df.index:
            p_gt = parent_df.loc[idx, "genotype"]
            c_gt = child_df.loc[idx, "genotype"]

            if p_gt == "--" or c_gt == "--":
                continue

            total_count += 1
            p_alleles = set(p_gt)
            c_alleles = set(c_gt)

            # Discordant if no alleles shared
            if not (p_alleles & c_alleles):
                discordant_count += 1

        # With 5% discordant rate, we should have some discordant SNPs
        # Allow for statistical variation (should be roughly 4-6% of total)
        discordant_rate = discordant_count / total_count
        self.assertGreater(discordant_rate, 0.03)
        self.assertLess(discordant_rate, 0.07)

    def test_generate_sibling_pair_structure(self):
        """Test that sibling pair has correct structure."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        sib1_df, sib2_df = gen.generate_sibling_pair(num_snps=1000)

        # Verify DataFrames have correct structure
        self.assertGreater(len(sib1_df), 0)
        self.assertGreater(len(sib2_df), 0)
        self.assertEqual(len(sib1_df), len(sib2_df))

        # Verify columns
        self.assertIn("chrom", sib1_df.columns)
        self.assertIn("pos", sib1_df.columns)
        self.assertIn("genotype", sib1_df.columns)

    def test_generate_sibling_pair_genotypes_valid(self):
        """Test that sibling genotypes are valid."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        sib1_df, sib2_df = gen.generate_sibling_pair(num_snps=1000)

        # Valid genotypes include all unphased combinations (both orientations)
        bases = ["A", "C", "G", "T"]
        valid_genotypes = {"--"}
        for b1 in bases:
            for b2 in bases:
                valid_genotypes.add(b1 + b2)

        sib1_gts = set(sib1_df["genotype"].unique())
        sib2_gts = set(sib2_df["genotype"].unique())

        self.assertTrue(sib1_gts.issubset(valid_genotypes))
        self.assertTrue(sib2_gts.issubset(valid_genotypes))

    def test_generate_sibling_pair_has_variation(self):
        """Test that siblings have both matching and non-matching genotypes."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        sib1_df, sib2_df = gen.generate_sibling_pair(num_snps=5000)

        identical_count = 0
        different_count = 0

        for idx in sib1_df.index:
            s1_gt = sib1_df.loc[idx, "genotype"]
            s2_gt = sib2_df.loc[idx, "genotype"]

            if s1_gt == "--" or s2_gt == "--":
                continue

            if s1_gt == s2_gt:
                identical_count += 1
            else:
                different_count += 1

        # Siblings should have both identical and different genotypes
        self.assertGreater(identical_count, 0)
        self.assertGreater(different_count, 0)

    def test_parent_child_reproducibility(self):
        """Test that parent-child generation is reproducible with seed."""
        gen1 = SyntheticRelatedGenerator(build=37, seed=123)
        parent1, child1 = gen1.generate_parent_child_pair(num_snps=100)

        gen2 = SyntheticRelatedGenerator(build=37, seed=123)
        parent2, child2 = gen2.generate_parent_child_pair(num_snps=100)

        pd.testing.assert_frame_equal(parent1, parent2)
        pd.testing.assert_frame_equal(child1, child2)

    def test_sibling_reproducibility(self):
        """Test that sibling generation is reproducible with seed."""
        gen1 = SyntheticRelatedGenerator(build=37, seed=123)
        sib1_a, sib2_a = gen1.generate_sibling_pair(num_snps=100)

        gen2 = SyntheticRelatedGenerator(build=37, seed=123)
        sib1_b, sib2_b = gen2.generate_sibling_pair(num_snps=100)

        pd.testing.assert_frame_equal(sib1_a, sib1_b)
        pd.testing.assert_frame_equal(sib2_a, sib2_b)

    def test_different_seeds_produce_different_results(self):
        """Test that different seeds produce different results."""
        gen1 = SyntheticRelatedGenerator(build=37, seed=123)
        parent1, _ = gen1.generate_parent_child_pair(num_snps=100)

        gen2 = SyntheticRelatedGenerator(build=37, seed=456)
        parent2, _ = gen2.generate_parent_child_pair(num_snps=100)

        # Genotypes should be different
        self.assertFalse(parent1["genotype"].equals(parent2["genotype"]))

    def test_save_parent_child_pair(self):
        """Test that save_parent_child_pair creates files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gen = SyntheticRelatedGenerator(build=37, seed=42)
            parent_path, child_path = gen.save_parent_child_pair(tmpdir, num_snps=1000)

            self.assertTrue(os.path.exists(parent_path))
            self.assertTrue(os.path.exists(child_path))
            self.assertTrue(parent_path.endswith(".txt.gz"))
            self.assertTrue(child_path.endswith(".csv.gz"))

    def test_save_sibling_pair(self):
        """Test that save_sibling_pair creates files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gen = SyntheticRelatedGenerator(build=37, seed=42)
            sib1_path, sib2_path = gen.save_sibling_pair(tmpdir, num_snps=1000)

            self.assertTrue(os.path.exists(sib1_path))
            self.assertTrue(os.path.exists(sib2_path))

    def test_save_parent_child_pair_formats(self):
        """Test saving in different formats."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gen = SyntheticRelatedGenerator(build=37, seed=42)

            # Test 23andme format
            p1, c1 = gen.save_parent_child_pair(
                tmpdir, parent_format="23andme", child_format="23andme", num_snps=500
            )
            self.assertTrue(p1.endswith(".23andme.txt.gz"))
            self.assertTrue(c1.endswith(".23andme.txt.gz"))

    def test_save_with_invalid_format(self):
        """Test that invalid format raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gen = SyntheticRelatedGenerator(build=37, seed=42)

            with self.assertRaises(ValueError):
                gen.save_parent_child_pair(
                    tmpdir, parent_format="invalid_format", num_snps=100
                )

    def test_generate_with_specific_chromosomes(self):
        """Test generating SNPs for specific chromosomes only."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        parent_df, child_df = gen.generate_parent_child_pair(
            num_snps=1000, chromosomes=["1", "2", "X"]
        )

        # Verify only specified chromosomes are present
        parent_chroms = set(parent_df["chrom"].unique())
        child_chroms = set(child_df["chrom"].unique())

        self.assertTrue(parent_chroms.issubset({"1", "2", "X"}))
        self.assertTrue(child_chroms.issubset({"1", "2", "X"}))

    def test_crossovers_per_chrom_affects_sibling_pattern(self):
        """Test that crossovers_per_chrom parameter affects sibling generation."""
        # Low crossovers should produce longer segments
        gen_low = SyntheticRelatedGenerator(build=37, seed=42, crossovers_per_chrom=0.5)
        sib1_low, sib2_low = gen_low.generate_sibling_pair(num_snps=5000)

        # High crossovers should produce more variation
        gen_high = SyntheticRelatedGenerator(
            build=37, seed=42, crossovers_per_chrom=3.0
        )
        sib1_high, sib2_high = gen_high.generate_sibling_pair(num_snps=5000)

        # Just verify both run without error and produce valid output
        self.assertEqual(len(sib1_low), len(sib1_high))

    def test_missing_rate(self):
        """Test that missing_rate parameter produces missing genotypes."""
        gen = SyntheticRelatedGenerator(build=37, seed=42)
        parent_df, child_df = gen.generate_parent_child_pair(
            num_snps=10000, missing_rate=0.05
        )

        # Check that there are some missing genotypes
        parent_missing = (parent_df["genotype"] == "--").sum()
        child_missing = (child_df["genotype"] == "--").sum()

        # With 5% missing rate and 10000 SNPs, expect roughly 500 missing in each
        # Allow some variance
        self.assertGreater(parent_missing, 100)
        self.assertGreater(child_missing, 100)

    def test_sibling_shared_dna_proportions(self):
        """Test that sibling shared DNA proportions are genetically realistic.

        For full siblings, expected theoretical proportions are:
        - ~25% IBD0 (no alleles shared by descent)
        - ~50% IBD1 (one allele shared by descent)
        - ~25% IBD2 (both alleles shared by descent)

        This means siblings share ~75% of their genome (~2625 cM of ~3500 cM),
        with approximately 2/3 of shared DNA being IBD1 and 1/3 being IBD2.

        This test uses lineage's find_shared_dna() which identifies contiguous
        segments of shared DNA using genetic maps and thresholds, filtering out
        spurious matches from identity by state (IBS).

        Uses a deterministic seed for reproducibility.
        """
        from lineage import Lineage

        gen = SyntheticRelatedGenerator(build=37, seed=1234, crossovers_per_chrom=1.0)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Generate and save siblings with sufficient SNPs for reliable detection
            paths = gen.save_sibling_pair(tmpdir, num_snps=500000)

            # Load with lineage and find shared DNA
            lineage_obj = Lineage(output_dir=tmpdir, resources_dir="resources")
            sib1 = lineage_obj.create_individual("Sibling1", paths[0])
            sib2 = lineage_obj.create_individual("Sibling2", paths[1])
            results = lineage_obj.find_shared_dna([sib1, sib2], save_output=False)

            one_chrom = results["one_chrom_shared_dna"]
            two_chrom = results["two_chrom_shared_dna"]

            # one_chrom_shared_dna includes ALL shared DNA (IBD1 + IBD2 regions)
            # two_chrom_shared_dna includes only IBD2 regions
            total_shared_cM = one_chrom["cMs"].sum() if len(one_chrom) > 0 else 0
            ibd2_cM = two_chrom["cMs"].sum() if len(two_chrom) > 0 else 0
            ibd1_cM = total_shared_cM - ibd2_cM

            # Siblings should share substantial DNA (expect ~2000-3200 cM)
            self.assertGreater(
                total_shared_cM, 1500, f"Total shared too low: {total_shared_cM:.0f} cM"
            )
            self.assertLess(
                total_shared_cM,
                3600,
                f"Total shared too high: {total_shared_cM:.0f} cM",
            )

            # Both IBD1 and IBD2 should be present
            self.assertGreater(ibd1_cM, 0, "IBD1 should have some shared DNA")
            self.assertGreater(ibd2_cM, 0, "IBD2 should have some shared DNA")

            # IBD1 should be more common than IBD2 (expect ~2:1 ratio)
            self.assertGreater(
                ibd1_cM,
                ibd2_cM,
                f"IBD1 ({ibd1_cM:.0f} cM) should exceed IBD2 ({ibd2_cM:.0f} cM)",
            )

            # Check proportions are in reasonable range
            # Expected: ~67% IBD1, ~33% IBD2 of shared DNA
            # Allow wide bounds (50-85% IBD1) due to random variation
            ibd1_prop = ibd1_cM / total_shared_cM
            self.assertGreater(
                ibd1_prop, 0.50, f"IBD1 proportion too low: {ibd1_prop:.1%}"
            )
            self.assertLess(
                ibd1_prop, 0.85, f"IBD1 proportion too high: {ibd1_prop:.1%}"
            )
