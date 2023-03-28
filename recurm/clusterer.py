import os
import errno
import sys
import logging
import shutil
import subprocess
import tempfile
from os import listdir


from recurm import FASTA_manager
from recurm import fileManager
from recurm import mapper
from recurm.defaultValues import DefaultValues
from recurm import grapher


class Clusterer():
    def __init__(self, mincontiglen, minclustersize, threads, out, filext, overwrite=False, resume=False,
                 keep_related=False, collapse_against_assembly=False):
        self.nthreads = threads
        #self.input_list = assemblies_list
        self.keep_related = keep_related
        self.collapse_against_assembly = collapse_against_assembly
        self.use_short_circular = False #this needs to be fixed as currently longer contigs are used to circularise shorter contigs which is not right
        if not resume:
            #fileManager.check_empty_dir(os.path.abspath(out), overwrite)
            self.outdir = os.path.abspath(out)
            self.extension = filext
            self.minlen = mincontiglen
            self.clustersize = int(minclustersize)
            self.mapping_dir = os.path.join(self.outdir, DefaultValues.MAPPING_DIR)
            fileManager.make_sure_path_exists(self.mapping_dir)
            self.alignment_dir = os.path.join(self.outdir, DefaultValues.ALIGNMENT_DIR)
            fileManager.make_sure_path_exists(self.alignment_dir)
            self.hashing_dir = os.path.join(self.outdir, DefaultValues.HASHING_DIR)
            fileManager.make_sure_path_exists(self.hashing_dir)
            self.results_dir = os.path.join(self.outdir, DefaultValues.RESULTS_DIR)
            fileManager.make_sure_path_exists(self.results_dir)

            self.collapse_dir = os.path.join(self.outdir, DefaultValues.COLLAPSE_DIR)
            fileManager.make_sure_path_exists(self.collapse_dir)

            self.combined_assemblies = None
            self.existing_mapping_directory = None
        else:
            logging.info('Resuming previous run of RecurM. Checking...')
            self.outdir = os.path.abspath(out)
            fileManager.make_sure_path_exists(self.outdir)

            self.combined_assemblies_file = os.path.join(self.outdir, DefaultValues.COMBINED_ASSEMBLY_NAME)
            if not os.path.isfile(self.combined_assemblies_file):
                self.combined_assemblies_file = None
                print('MADE ASSEMBLY NONE')

            self.minlen = mincontiglen
            self.clustersize = int(minclustersize)
            self.mapping_dir = os.path.join(self.outdir, DefaultValues.MAPPING_DIR)
            fileManager.make_sure_path_exists(self.mapping_dir)
            self.alignment_dir = os.path.join(self.outdir, DefaultValues.ALIGNMENT_DIR)
            fileManager.check_empty_dir(self.alignment_dir, overwrite=True)
            self.hashing_dir = os.path.join(self.outdir, DefaultValues.HASHING_DIR)
            fileManager.check_empty_dir(self.hashing_dir, overwrite=True)
            self.results_dir = os.path.join(self.outdir, DefaultValues.RESULTS_DIR)
            fileManager.check_empty_dir(self.results_dir, overwrite=True)

            for filename in os.listdir(self.outdir):
                if str(filename) != str(DefaultValues.COMBINED_ASSEMBLY_NAME):
                    file_path = os.path.join(self.outdir, filename)
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)

            cluster_folder = os.path.join(self.outdir, 'clusters')
            try:
                shutil.rmtree(cluster_folder)
            except Exception as e:
                pass
            

            if len(os.listdir(self.mapping_dir) ) == 0:
                self.existing_mapping_directory = None
            else:
                self.existing_mapping_directory = self.mapping_dir





    def cluster_all_vs_all(self, combined_assemblies, master_assembly, is_chunk):
        if self.existing_mapping_directory is None or self.combined_assemblies_file is None:
            #logging.info('Concatenating all assemblies into master assembly file. Only contigs of length {} or more bp considered.'.format(self.minlen))
            #combined_assemblies = FASTA_manager.setup_and_format_assemblies(self.minlen, self.input_list, self.outdir)

            AVA_mapper = mapper.AllVsAllMapper(self.nthreads)
            logging.info('Mapping all contigs against each other to identify putative clusters (warning: this can take a long time).')
            paf_alignment_list = AVA_mapper.mapAVA(combined_assemblies, self.mapping_dir, self.minlen)
        else:
            AVA_mapper = mapper.AllVsAllMapper(self.nthreads)
            combined_assemblies = os.path.abspath(combined_assemblies)
            paf_alignment_list = [os.path.join(self.existing_mapping_directory, mapping_file) for mapping_file in listdir(self.existing_mapping_directory) if mapping_file.endswith('.paf')]

        if not is_chunk:
            logging.info('Hashing and filtering alignments based on following cut-offs:')
            logging.info('')
            logging.info('\tLength Ratio: \t\t\t\tMore than {}'.format(DefaultValues.FIRST_PASS_AVA_LR_CUTOFF))
            logging.info('\tAlignment Ratio: \t\t\tMore than {}'.format(DefaultValues.FIRST_PASS_AVA_AR_CUTOFF))
            logging.info('\tANI Ratio: \t\t\t\tMore than {}'.format(DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF))
            logging.info('\tAllowed cluster size deviation: \tLess than {}'.format(DefaultValues.ALLOWED_SIZE_DEVIATION_WITHIN_CLUSTER))
            logging.info('\tCluster alignment tightness: \t\tMore than {}'.format(DefaultValues.TIGHTNESS_CLUSTER_CUTOFF))
            logging.info('')


        if len(paf_alignment_list) < 1:
            if is_chunk:
                return
            else:
                logging.error('ERROR: No alignment files were generated. Consider adding more metagenomes to RecurM input. Exiting.')
                sys.exit(1)


        hashfile, filteredfile, hash_success = AVA_mapper.subfilter_paf_mapping_results(self.nthreads,
                                                 paf_alignment_list,
                                                 self.alignment_dir,
                                                 self.hashing_dir,
                                                 DefaultValues.FILTER_CHUNKSIZE,
                                                 DefaultValues.FIRST_PASS_AVA_LR_CUTOFF,
                                                 DefaultValues.FIRST_PASS_AVA_AR_CUTOFF,
                                                 DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF)

        if not hash_success:
            if is_chunk:
                return
            else:
                logging.error('ERROR: No valid hash files were generated. Consider adding more metagenomes to RecurM input. Exiting.')
                sys.exit(1)

        if not is_chunk:
            logging.info('Sorting extracted hashes with {} threads'.format(self.nthreads))

        fileManager.bash_sort_file(hashfile, 2, self.outdir, self.nthreads, numerical=True)

        if not is_chunk:
            logging.info('Identifying multiple alignments from hashes and extracting mappings.')
        hashes_extracted_file = AVA_mapper.extract_hashes_with_multiple_hits(hashfile)

        #sort to pop off by sample
        fileManager.bash_sort_file(hashes_extracted_file, 3, self.outdir, self.nthreads, numerical=False)

        circular_success = AVA_mapper.read_hashes_and_process_alignments(self.mapping_dir,
                                                      hashes_extracted_file,
                                                      DefaultValues.FIRST_PASS_AVA_LR_CUTOFF,
                                                      DefaultValues.FIRST_PASS_AVA_AR_CUTOFF,
                                                      DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF,
                                                      '{}/{}'.format(self.outdir, DefaultValues.CIRCULAR_ALIGNMENTS_NAME),
                                                      self.use_short_circular)

        if is_chunk:
            return

        logging.info('Drawing master graph from contig alignments.')
        graph_process = grapher.ClusterGraph(self.nthreads)
        master_graph, alignments_file = \
            graph_process.create_graph_from_alignments('{}/{}'.format(self.outdir, DefaultValues.CIRCULAR_ALIGNMENTS_NAME),
                                                                  '{}/{}'.format(self.alignment_dir, DefaultValues.SECOND_PASS_NAME), circular_success)
        logging.info('Extracting disconnected subgraphs.')
        cluster_information, sub_graphs = \
            graph_process.retrieve_disconnected_subgraphs(master_graph, self.clustersize, alignments_file, self.outdir)

        logging.info('Extracting representative contigs.')
        cluster_folder, cluster_contigs_info, leftover_contigs_info, filtered_out, x1, x2, x3, x4 = \
            FASTA_manager.write_clusters_to_file(combined_assemblies, self.outdir, sub_graphs, cluster_information)



        logging.info('Identifying related clusters.')
        contig_mapper = mapper.WithinClusterMapper(self.nthreads)
        mapped_contigs, combined_contig_file = contig_mapper.map_clusters_to_each_other(cluster_folder, self.outdir)
        cluster_information, cluster_contigs_info, leftover_contigs_info, ani_info \
            = contig_mapper.collapse_related_clusters(mapped_contigs, cluster_folder,
                                                      DefaultValues.FIRST_PASS_AVA_LR_CUTOFF,
                                                      DefaultValues.FIRST_PASS_AVA_AR_CUTOFF,
                                                      DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF,
                                                      self.outdir,
                                                      cluster_information, cluster_contigs_info, leftover_contigs_info,
                                                      self.keep_related, filtered_out)


        if self.collapse_against_assembly:
            logging.info('Removing any integrated clusters.')
            assembly_mapper = mapper.AssemblyMapper(self.nthreads)
            assembly_mapped_contigs = assembly_mapper.map_clusters_against_assembly(cluster_folder, master_assembly, self.outdir, self.collapse_dir, self.minlen)

            cluster_information, cluster_contigs_info, \
                = assembly_mapper.remove_mapped_clusters(self.nthreads, assembly_mapped_contigs, cluster_folder,
                                                          DefaultValues.FIRST_PASS_AVA_LR_CUTOFF,
                                                          DefaultValues.FIRST_PASS_AVA_AR_CUTOFF,
                                                          DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF,
                                                          self.outdir,
                                                          cluster_information, cluster_contigs_info)


        logging.info('Writing summary information to file: {}/cluster_information.tsv'.format(self.results_dir))

        ''' DEBUGGING PURPOSES START'''

        cluster_information['EDGE_OCCURENCES'] = cluster_information['ID'].apply(lambda x: x1[x])
        cluster_information['TIGHTNESS'] = (cluster_information['EDGE_OCCURENCES'] + 1) / cluster_information['Total_Contigs']
        cluster_information['SIZE_MEAN'] = cluster_information['ID'].apply(lambda x: x3[x])
        cluster_information['SIZE_STD'] = cluster_information['ID'].apply(lambda x: x4[x])
        cluster_information['SIZE_CLUSTER_DEVIATION'] = cluster_information['SIZE_STD'] / cluster_information['SIZE_MEAN']

        ''' DEBUGGING PURPOSES END'''

        cluster_information.to_csv('{}/cluster_information.tsv'.format(self.results_dir), sep='\t', index=False)
        leftover_contigs_info.to_csv('{}/leftover_contigs_information.tsv'.format(self.results_dir), sep='\t', index=False)
        cluster_contigs_info.to_csv('{}/cluster_contigs_information.tsv'.format(self.results_dir), sep='\t', index=False)

        if ani_info is not None:
            logging.info('Writing information about related clusters to file: '
                         '{}/related_clusters_ANI.tsv'.format(self.results_dir, self.results_dir))
            ani_info.to_csv('{}/related_clusters_ANI.tsv'.format(self.results_dir), sep='\t', index=False)
            #ani_matrix.to_csv('{}/related_clusters_ANI_matrix.tsv'.format(self.results_dir), sep='\t')

        logging.info('')
        logging.info(
            '********************************************* RESULTS: ************************************************************')
        logging.info('Identified a total of {} clusters. Of these, {} are circular, {} are linear and {} are imperfect.'
                     .format(len(cluster_information),
                             len(cluster_information[cluster_information['Structure'] == 'Circular']),
                             len(cluster_information[cluster_information['Structure'] == 'Linear']),
                             len(cluster_information[cluster_information['Structure'] == 'Imperfect'])))
        logging.info(
            '*******************************************************************************************************************')
        logging.info('')


        logging.info('Cleaning up intermediate files.')
        # fileManager.remove_intermediates([hashfile,
        #                                   combined_assemblies,
        #                                   filteredfile,
        #                                   paf_first_out,
        #                                   '{}/{}'.format(self.outdir, DefaultValues.CIRCULAR_ALIGNMENTS_NAME),
        #                                   hashes_extracted_file])












