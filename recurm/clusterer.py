import os
import errno
import sys
import logging
import shutil
import subprocess
import tempfile


from recurm import FASTA_manager
from recurm import fileManager
from recurm import mapper
from recurm.defaultValues import DefaultValues
from recurm import grapher


class Clusterer():
    def __init__(self, mincontiglen, threads, assemblies_list, out, filext, overwrite=False):
        self.nthreads = threads
        self.input_list = assemblies_list
        fileManager.check_empty_dir(os.path.abspath(out), overwrite)
        self.outdir = os.path.abspath(out)
        self.extension = filext
        self.minlen = mincontiglen


    def cluster_all_vs_all(self):
        logging.info('Concatenating all assemblies into master assembly file. Only contigs of lengh {} or more bp considered.'.format(self.minlen))
        combined_assemblies = FASTA_manager.setup_and_format_assemblies(self.minlen, self.input_list, self.outdir)
        AVA_mapper = mapper.AllVsAllMapper(self.nthreads)
        logging.info('Mapping all contigs against each other to identify putative clusters (warning: this can take a long time).')
        paf_first_out = AVA_mapper.mapAVA(combined_assemblies, self.outdir)
        logging.info('Hashing and filtering alignments based on following cut-offs:')
        logging.info('\t\tLength Ratio: More than {}'.format(DefaultValues.FIRST_PASS_AVA_LR_CUTOFF))
        logging.info('\t\tAlignment Ratio: More than {}'.format(DefaultValues.FIRST_PASS_AVA_AR_CUTOFF))
        logging.info('\t\tANI Ratio: More than {}'.format(DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF))
        hashfile, filteredfile = AVA_mapper.subfilter_paf_mapping_results(paf_first_out,
                                                 self.outdir,
                                                 DefaultValues.SECOND_PASS_NAME,
                                                 DefaultValues.FILTER_CHUNKSIZE,
                                                 DefaultValues.FIRST_PASS_AVA_LR_CUTOFF,
                                                 DefaultValues.FIRST_PASS_AVA_AR_CUTOFF,
                                                 DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF)

        logging.info('Sorting extracted hashes with {} threads'.format(self.nthreads))
        try:
            sortedfile = '{}/TMP.hashes'.format(self.outdir)
            cmd = "sort -t$'\t' -k2 --parallel {} -nr {} > {} && mv {} {}"\
                .format(self.nthreads, hashfile, sortedfile, sortedfile, hashfile)
            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished sorting')
        except Exception as e:
            logging.error('An error occured while sorting hashes: {}'.format(e))
            sys.exit(1)
        logging.info('Identifying multiple alignments from hashes and extracting mappings.')
        hashes_extracted_file = AVA_mapper.extract_hashes_with_multiple_hits(hashfile)
        circular_success = AVA_mapper.read_hashes_and_process_alignments(hashfile,
                                                      paf_first_out,
                                                      DefaultValues.FIRST_PASS_AVA_LR_CUTOFF,
                                                      DefaultValues.FIRST_PASS_AVA_AR_CUTOFF,
                                                      DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF,
                                                      '{}/{}.passed'.format(self.outdir, DefaultValues.CIRCULAR_ALIGNMENTS_NAME))

        logging.info('Drawing master graph from contig alignments.')
        graph_process = grapher.ClusterGraph(self.nthreads)
        master_graph, alignments_file = graph_process.create_graph_from_alignments('{}/{}.passed'.format(self.outdir, DefaultValues.CIRCULAR_ALIGNMENTS_NAME),
                                                                  '{}/{}.passed'.format(self.outdir, DefaultValues.SECOND_PASS_NAME), circular_success)
        logging.info('Extracting disconnected subgraphs.')
        cluster_information, sub_graphs = graph_process.retrieve_disconnected_subgraphs(master_graph, DefaultValues.MIN_CLUSTER_SIZE, alignments_file)

        logging.info('Writing summary information to file: {}/cluster_information.tsv'.format(self.outdir))

        cluster_information.to_csv('{}/cluster_information.tsv'.format(self.outdir), sep='\t', index=False)

        logging.info('Extracting representative contigs.')
        FASTA_manager.write_clusters_to_file(combined_assemblies, self.outdir, sub_graphs,cluster_information)

        logging.info('Cleaning up intermediate files.')
        fileManager.remove_intermediates([hashfile,
                                          combined_assemblies,
                                          filteredfile,
                                          paf_first_out,
                                          '{}/{}.passed'.format(self.outdir, DefaultValues.CIRCULAR_ALIGNMENTS_NAME),
                                          hashes_extracted_file])








