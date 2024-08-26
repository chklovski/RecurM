import os
from pathlib import Path

class DefaultValues():
    """Default values for filenames and common constants."""


    COMBINED_ASSEMBLY_NAME = 'combined_assemblies.fasta'
    FIRST_PASS_MAPPING_NAME = 'ava_initial_alignments'
    SECOND_PASS_NAME = 'initial_subfilter.passed'
    CIRCULAR_ALIGNMENTS_NAME = 'circular_alignments.passed'
    HASH_FILE_NAME = 'combined_hash_data.hash'

    COLLAPSE_MAPPING_NAME = 'ava_collapse_alignments'

    COMBINED_CONTIGS_FILE = 'combined_contigs_ref.fasta'
    COMBINED_CONTIGS_MAPPING_FILE = 'combined_ref_map.paf'

    ASSEMBLY_DEREP_FILE = 'cluster_to_assembly_map.paf'

    COLLAPSED_LOG_FILE = 'collapsed_related_clusters.tsv'

    MAPPING_DIR = 'mapping_data'
    ALIGNMENT_DIR = 'alignment_data'
    HASHING_DIR = 'hashing_data'
    MULTIALIGN_DIR = 'multialign_data'
    RESULTS_DIR = 'results'
    COLLAPSE_DIR = 'collapse_mapping'



    CONTIG_SIZE_BINS = 10
    DEFAULT_MIN_BIN_SIZE = 8000

    SIZE_BIN_NAME = 'RECURM_SIZEBIN'
    MASTER_ASSEMBLY_FILE = 'MASTER_ASSEMBLY_FILE'

    EXTRACTED_CHUNK_CONTIGS_FILE = 'contigs_to_cluster.tsv'

    EXTRACTED_CHUNK_ALIGNMENT_FILE = 'all_passed_alignments.paf'

    GRAPH_PICKLE_NAME = 'single_linkage_graph_data.graph'


    DEFAULT_FASTA_HEADER_SEPARATOR = 'Ω'

    FIRST_PASS_AVA_LR_CUTOFF = 0.9
    FIRST_PASS_AVA_AR_CUTOFF = 0.9
    FIRST_PASS_AVA_ANI_CUTOFF = 0.9

    ASSEMBLY_DEREP_CONTIG_PERCENT_LONGER = 0.25

    # aligns to x% of all other contigs in cluster
    TIGHTNESS_CLUSTER_CUTOFF = 0.5

    # percent representative sequence size is allowed to deviate from mean
    # intended to filter out smeared assembly length type clusters
    ALLOWED_SIZE_DEVIATION_WITHIN_CLUSTER = 0.5

    SECOND_PASS_CONTIGS_AR_SHORT_CUTOFF = 0.70

    CIRCULAR_ANI_INCLUDING_CUTOFF = 0.956

    MIN_CLUSTER_SIZE = 3

    EDGE_PRUNE_THRESHOLD = 6

    FILTER_CHUNKSIZE = 25000

    PAF_CHUNK_SIZE = '50MB'




