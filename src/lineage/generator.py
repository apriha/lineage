"""Generate synthetic genotype data for related individuals.

This module extends the SyntheticSNPGenerator from snps to create pairs of
genetically related individuals (parent-child and sibling pairs) suitable
for testing and demonstrating lineage's shared DNA analysis capabilities.
"""

from __future__ import annotations

import logging
import os

import numpy as np
import pandas as pd
from snps.io.generator import SyntheticSNPGenerator

logger = logging.getLogger(__name__)

# Valid nucleotide bases for genotype generation
BASES = np.array(["A", "C", "G", "T"])


class SyntheticRelatedGenerator(SyntheticSNPGenerator):
    """Generate synthetic genotype data for related individuals.

    This class extends SyntheticSNPGenerator to create pairs of genetically
    related individuals with realistic inheritance patterns.

    Genetic Relationships:

    - Parent-Child: Share exactly ONE allele at every position (IBD1 = 100%)
    - Siblings: Share DNA in segments due to recombination

      - ~25% IBD2 (both chromosomes identical)
      - ~50% IBD1 (one chromosome identical)
      - ~25% IBD0 (no sharing)

    Parameters
    ----------
    build : int
        Genome build (36, 37, or 38), default is 37
    seed : int, optional
        Random seed for reproducibility
    crossovers_per_chrom : float
        Average number of crossovers per chromosome per meiosis (default: 1.0)

    Examples
    --------
    >>> gen = SyntheticRelatedGenerator(build=37, seed=42)
    >>> parent_df, child_df = gen.generate_parent_child_pair(num_snps=10000)
    >>> sib1_df, sib2_df = gen.generate_sibling_pair(num_snps=10000)
    """

    def __init__(
        self,
        build: int = 37,
        seed: int | None = None,
        crossovers_per_chrom: float = 1.0,
    ) -> None:
        super().__init__(build=build, seed=seed)
        self.crossovers_per_chrom = crossovers_per_chrom

    def generate_parent_child_pair(
        self,
        num_snps: int = 100000,
        chromosomes: list[str] | None = None,
        missing_rate: float = 0.01,
        discordant_rate: float = 0.0,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Generate a parent-child pair with realistic inheritance.

        The child inherits exactly one allele from the parent at every position,
        resulting in 100% IBD1 (one chromosome shared) when analyzed.

        Parameters
        ----------
        num_snps : int
            Approximate number of SNPs to generate (default: 100000)
        chromosomes : list of str, optional
            Chromosomes to include (default: all autosomes plus X, Y, MT)
        missing_rate : float
            Proportion of SNPs with missing genotypes (default: 0.01)
        discordant_rate : float
            Proportion of SNPs that are discordant (child doesn't inherit
            from parent), simulating genotyping errors (default: 0.0)

        Returns
        -------
        tuple of (pd.DataFrame, pd.DataFrame)
            (parent_df, child_df) - DataFrames with rsid index, chrom, pos, genotype
        """
        # Generate parent genotypes
        parent_df = self.generate_snps(
            num_snps=num_snps,
            chromosomes=chromosomes,
            missing_rate=missing_rate,
        )

        # Generate child by inheriting one allele from parent
        child_df = self._generate_child_from_parent(
            parent_df, missing_rate, discordant_rate
        )

        return parent_df, child_df

    def _generate_child_from_parent(
        self,
        parent_df: pd.DataFrame,
        missing_rate: float,
        discordant_rate: float = 0.0,
    ) -> pd.DataFrame:
        """Generate child genotypes from parent.

        For each SNP:
        - Child inherits one allele from parent
        - Child gets one random allele from "other parent" (simulated)
        - With discordant_rate probability, child gets a genotype that doesn't
          match parent (simulating genotyping errors)
        """
        child_df = parent_df.copy()
        n = len(parent_df)

        parent_genotypes = parent_df["genotype"].values
        child_genotypes = np.empty(n, dtype=object)

        # Generate random alleles for "other parent"
        other_parent_alleles = BASES[self.rng.integers(0, 4, size=n)]

        # Determine which SNPs will be discordant
        discordant_mask = self.rng.random(n) < discordant_rate

        for i, parent_gt in enumerate(parent_genotypes):
            if parent_gt == "--":
                # If parent is missing, child is also missing
                child_genotypes[i] = "--"
            elif discordant_mask[i]:
                # Generate a discordant genotype (no allele matches parent)
                # Get alleles that are NOT in the parent genotype
                parent_alleles = set(parent_gt)
                non_parent_alleles = [b for b in BASES if b not in parent_alleles]
                if len(non_parent_alleles) >= 2:
                    # Pick two random non-parent alleles
                    chosen = self.rng.choice(non_parent_alleles, size=2, replace=True)
                    child_genotypes[i] = chosen[0] + chosen[1]
                elif len(non_parent_alleles) == 1:
                    # Homozygous for the non-parent allele
                    child_genotypes[i] = non_parent_alleles[0] + non_parent_alleles[0]
                else:
                    # Parent is heterozygous with all 4 bases (shouldn't happen)
                    # Fall back to normal inheritance
                    inherited_allele = parent_gt[self.rng.integers(0, 2)]
                    child_genotypes[i] = inherited_allele + other_parent_alleles[i]
            else:
                # Child inherits one allele from parent (randomly choose which)
                inherited_allele = parent_gt[self.rng.integers(0, 2)]
                other_allele = other_parent_alleles[i]
                # Combine alleles without sorting to mimic unphased genotype files
                child_genotypes[i] = inherited_allele + other_allele

        # Apply missing rate to child (additional missing data)
        missing_mask = self.rng.random(n) < missing_rate
        child_genotypes[missing_mask] = "--"

        child_df["genotype"] = child_genotypes
        return child_df

    def generate_sibling_pair(
        self,
        num_snps: int = 100000,
        chromosomes: list[str] | None = None,
        missing_rate: float = 0.01,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Generate a sibling pair with realistic recombination patterns.

        Siblings share DNA in segments due to recombination during meiosis:
        - ~25% IBD2 (both chromosomes identical from both parents)
        - ~50% IBD1 (one chromosome identical)
        - ~25% IBD0 (no sharing)

        Parameters
        ----------
        num_snps : int
            Approximate number of SNPs to generate (default: 100000)
        chromosomes : list of str, optional
            Chromosomes to include (default: all autosomes plus X, Y, MT)
        missing_rate : float
            Proportion of SNPs with missing genotypes (default: 0.01)

        Returns
        -------
        tuple of (pd.DataFrame, pd.DataFrame)
            (sibling1_df, sibling2_df) - DataFrames with rsid index, chrom, pos, genotype
        """
        from snps.constants import REFERENCE_SEQUENCE_CHROMS

        if chromosomes is None:
            chromosomes = list(REFERENCE_SEQUENCE_CHROMS)

        # Generate four parental haplotypes (two per parent)
        # Parent 1: haplotypes A and B
        # Parent 2: haplotypes C and D
        base_snps = self.generate_snps(
            num_snps=num_snps,
            chromosomes=chromosomes,
            missing_rate=0,  # Handle missing rate separately
            inject_build_markers=True,
        )

        # Generate haplotypes for all four parental chromosomes
        hap_a = self._generate_haplotype(len(base_snps))
        hap_b = self._generate_haplotype(len(base_snps))
        hap_c = self._generate_haplotype(len(base_snps))
        hap_d = self._generate_haplotype(len(base_snps))

        # For each sibling, simulate meiosis from each parent
        sib1_genotypes = self._simulate_sibling_genotypes(
            base_snps, hap_a, hap_b, hap_c, hap_d, chromosomes
        )
        sib2_genotypes = self._simulate_sibling_genotypes(
            base_snps, hap_a, hap_b, hap_c, hap_d, chromosomes
        )

        # Create sibling DataFrames
        sib1_df = base_snps.copy()
        sib2_df = base_snps.copy()
        sib1_df["genotype"] = sib1_genotypes
        sib2_df["genotype"] = sib2_genotypes

        # Apply missing rate
        n = len(base_snps)
        missing_mask_1 = self.rng.random(n) < missing_rate
        missing_mask_2 = self.rng.random(n) < missing_rate
        sib1_df.loc[sib1_df.index[missing_mask_1], "genotype"] = "--"
        sib2_df.loc[sib2_df.index[missing_mask_2], "genotype"] = "--"

        return sib1_df, sib2_df

    def _generate_haplotype(self, n: int) -> np.ndarray:
        """Generate a random haplotype (single alleles, not genotypes)."""
        return BASES[self.rng.integers(0, 4, size=n)]

    def _simulate_sibling_genotypes(
        self,
        base_snps: pd.DataFrame,
        hap_a: np.ndarray,
        hap_b: np.ndarray,
        hap_c: np.ndarray,
        hap_d: np.ndarray,
        chromosomes: list[str],
    ) -> np.ndarray:
        """Simulate genotypes for one sibling via meiosis from both parents.

        For each chromosome:
        1. Simulate recombination in parent 1 (between hap_a and hap_b)
        2. Simulate recombination in parent 2 (between hap_c and hap_d)
        3. Combine the resulting gametes to form the sibling's genotype
        """
        n = len(base_snps)
        genotypes = np.empty(n, dtype=object)

        # Get chromosome boundaries
        chrom_col = base_snps["chrom"].values

        for chrom in chromosomes:
            chrom_mask = chrom_col == chrom
            chrom_indices = np.where(chrom_mask)[0]

            if len(chrom_indices) == 0:
                continue

            # Simulate meiosis for this chromosome from each parent
            gamete_from_p1 = self._simulate_meiosis(
                hap_a[chrom_indices], hap_b[chrom_indices]
            )
            gamete_from_p2 = self._simulate_meiosis(
                hap_c[chrom_indices], hap_d[chrom_indices]
            )

            # Combine gametes to form genotypes (without sorting to mimic unphased data)
            for i, idx in enumerate(chrom_indices):
                genotypes[idx] = gamete_from_p1[i] + gamete_from_p2[i]

        return genotypes

    def _simulate_meiosis(self, hap1: np.ndarray, hap2: np.ndarray) -> np.ndarray:
        """Simulate meiosis with recombination between two haplotypes.

        Returns a gamete (single haplotype) that is a recombinant of the two inputs.
        """
        n = len(hap1)
        if n == 0:
            return np.array([], dtype=object)

        # Generate crossover points using Poisson distribution
        num_crossovers = self.rng.poisson(self.crossovers_per_chrom)
        if num_crossovers > 0 and n > 1:
            crossover_positions = np.sort(self.rng.integers(1, n, size=num_crossovers))
        else:
            crossover_positions = np.array([], dtype=int)

        # Build the gamete by alternating between haplotypes at crossover points
        gamete = np.empty(n, dtype=object)
        current_hap = hap1 if self.rng.random() < 0.5 else hap2
        other_hap = hap2 if current_hap is hap1 else hap1

        prev_pos = 0
        for pos in crossover_positions:
            gamete[prev_pos:pos] = current_hap[prev_pos:pos]
            # Swap haplotypes
            current_hap, other_hap = other_hap, current_hap
            prev_pos = pos

        # Fill in the rest
        gamete[prev_pos:] = current_hap[prev_pos:]

        return gamete

    def save_parent_child_pair(
        self,
        output_dir: str,
        parent_format: str = "23andme",
        child_format: str = "ftdna",
        num_snps: int = 100000,
        discordant_rate: float = 0.0,
        **kwargs,
    ) -> tuple[str, str]:
        """Generate and save a parent-child pair to files.

        Parameters
        ----------
        output_dir : str
            Directory for output files
        parent_format : str
            Format for parent file ('23andme', 'ftdna', 'ancestry')
        child_format : str
            Format for child file ('23andme', 'ftdna', 'ancestry')
        num_snps : int
            Approximate number of SNPs to generate
        discordant_rate : float
            Proportion of SNPs that are discordant (default: 0.0)
        **kwargs
            Additional arguments passed to generate_parent_child_pair

        Returns
        -------
        tuple of (str, str)
            Paths to (parent_file, child_file)
        """
        os.makedirs(output_dir, exist_ok=True)

        parent_df, child_df = self.generate_parent_child_pair(
            num_snps=num_snps, discordant_rate=discordant_rate, **kwargs
        )

        parent_path = self._save_dataframe(
            parent_df, output_dir, "parent", parent_format
        )
        child_path = self._save_dataframe(child_df, output_dir, "child", child_format)

        return parent_path, child_path

    def save_sibling_pair(
        self,
        output_dir: str,
        sib1_format: str = "23andme",
        sib2_format: str = "ftdna",
        num_snps: int = 100000,
        **kwargs,
    ) -> tuple[str, str]:
        """Generate and save a sibling pair to files.

        Parameters
        ----------
        output_dir : str
            Directory for output files
        sib1_format : str
            Format for sibling 1 file ('23andme', 'ftdna', 'ancestry')
        sib2_format : str
            Format for sibling 2 file ('23andme', 'ftdna', 'ancestry')
        num_snps : int
            Approximate number of SNPs to generate
        **kwargs
            Additional arguments passed to generate_sibling_pair

        Returns
        -------
        tuple of (str, str)
            Paths to (sibling1_file, sibling2_file)
        """
        os.makedirs(output_dir, exist_ok=True)

        sib1_df, sib2_df = self.generate_sibling_pair(num_snps=num_snps, **kwargs)

        sib1_path = self._save_dataframe(sib1_df, output_dir, "sibling1", sib1_format)
        sib2_path = self._save_dataframe(sib2_df, output_dir, "sibling2", sib2_format)

        return sib1_path, sib2_path

    def _save_dataframe(
        self, df: pd.DataFrame, output_dir: str, name: str, format: str
    ) -> str:
        """Save a DataFrame in the specified format."""
        format_info = {
            "23andme": ("txt.gz", self._write_snps_as_23andme),
            "ftdna": ("csv.gz", self._write_snps_as_ftdna),
            "ancestry": ("txt.gz", self._write_snps_as_ancestry),
        }

        if format not in format_info:
            raise ValueError(
                f"Unknown format: {format}. Use: {list(format_info.keys())}"
            )

        ext, writer = format_info[format]
        path = os.path.join(output_dir, f"{name}.{format}.{ext}")
        logger.info(f"Creating {os.path.relpath(path)}")
        writer(df, path)
        return path
