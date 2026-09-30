"""
Configuration file for test resource sizes and shared gene counts.
These values size the mocked resource data used when downloads are disabled.
"""

# Resource sizes for hg19 data
RESOURCE_SIZES = {"knownGene_hg19": 503413, "kgXref_hg19": 503369, "cytoBand_hg19": 862}

# Shared gene counts for specific test cases
SHARED_GENE_COUNTS = {
    "default": {"len1": 23148, "len2": 23148},
    "1000G": {"len1": 12812},
    "X_chrom_male": {"len1": 383, "len2": 15547},
    "X_chrom_female": {"len1": 15547, "len2": 15547},
}
