import os
from pathlib import Path

class DefaultValues():
    """Default values for filenames and common constants."""

    COMBINED_ASSEMBLY_NAME = 'combined_assemblies.fasta'
    FIRST_PASS_MAPPING_NAME = 'ava_initial_alignments.paf'
    SECOND_PASS_NAME = 'initial_subfilter'
    CIRCULAR_ALIGNMENTS_NAME = 'circular_alignments'

    DEFAULT_FASTA_HEADER_SEPARATOR = '@'

    FIRST_PASS_AVA_LR_CUTOFF = 0.9
    FIRST_PASS_AVA_AR_CUTOFF = 0.9
    FIRST_PASS_AVA_ANI_CUTOFF = 0.9

    MIN_CLUSTER_SIZE = 3

    EDGE_PRUNE_THRESHOLD = 6

    FILTER_CHUNKSIZE = 25000




