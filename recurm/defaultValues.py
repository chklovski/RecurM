import os
from pathlib import Path

class DefaultValues():
    """Default values for filenames and common constants."""

    COMBINED_ASSEMBLY_NAME = 'combined_assemblies.fasta'
    FIRST_PASS_MAPPING_NAME = 'ava_initial_alignments'
    SECOND_PASS_NAME = 'initial_subfilter.passed'
    CIRCULAR_ALIGNMENTS_NAME = 'circular_alignments.passed'
    HASH_FILE_NAME = 'combined_hash_data.hash'

    COMBINED_CONTIGS_FILE = 'combined_contigs_ref.fasta'
    COMBINED_CONTIGS_MAPPING_FILE = 'combined_ref_map.paf'

    COLLAPSED_LOG_FILE = 'collapsed_related_clusters.tsv'

    MAPPING_DIR = 'mapping_data'
    ALIGNMENT_DIR = 'alignment_data'
    HASHING_DIR = 'hashing_data'
    RESULTS_DIR = 'results'

    GRAPH_PICKLE_NAME = 'single_linkage_graph_data.graph'

    DEFAULT_FASTA_HEADER_SEPARATOR = 'Ω'

    FIRST_PASS_AVA_LR_CUTOFF = 0.9
    FIRST_PASS_AVA_AR_CUTOFF = 0.9
    FIRST_PASS_AVA_ANI_CUTOFF = 0.9

    SECOND_PASS_CONTIGS_AR_SHORT_CUTOFF = 0.70

    MIN_CLUSTER_SIZE = 3

    EDGE_PRUNE_THRESHOLD = 6

    FILTER_CHUNKSIZE = 25000

    PAF_CHUNK_SIZE = '250MB'




