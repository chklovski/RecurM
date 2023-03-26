import os
import pandas as pd
import sys
import logging
import shutil
import networkx as nx
import numpy as np
import multiprocessing as mp


from recurm import fileManager
from recurm.defaultValues import DefaultValues


def pct_change(first, second):
    diff = float(second) - float(first)
    change = 0
    try:
        if diff > 0:
            change = (diff / first)
        elif diff < 0:
            diff = first - second
            change = -(diff / first)
    except ZeroDivisionError:
        return float('inf')
    return change

def determine_size_cutoffs(nthreads, assemblies, minlen, NUMBER_BINS):


    logging.info('Using {} threads to bin contig sizes into equal size bins.'.format(nthreads))

    workerQueue = mp.Queue()
    writerQueue = mp.Queue()

    for idx, params in enumerate(assemblies):
        workerQueue.put(idx)

    for _ in range(nthreads):
        workerQueue.put(None)

    processed_contigsize_dfs = mp.Manager().list()

    try:
        calcProc = []
        for _ in range(nthreads):
            calcProc.append(
                mp.Process(target=analyse_assembly_contig_size_thread,
                           args=(workerQueue, writerQueue, processed_contigsize_dfs, assemblies, minlen)))
        writeProc = mp.Process(target=report_progress_assembly_threads,
                               args=(len(assemblies), writerQueue))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put((None, None))
        writeProc.join()

    except:
        # make sure all processes are terminated
        for p in calcProc:
            p.terminate()

        writeProc.terminate()

    overall = pd.concat(processed_contigsize_dfs)
    overall = overall.sort_values(by='Size', ascending=True)
    overall.reset_index(drop=True, inplace=True)

    # split into ~equal parts:
    index_split = np.array_split(overall.index, NUMBER_BINS)

    # get top size per chunk:
    top_sizes = []
    for part in range(0, NUMBER_BINS):
        # get the last index of a size part
        ix = index_split[part][[-1]]
        # append size value to list
        top_sizes.append(overall.iloc[ix, 0].values[0])


    # now construct a list of contig names to be grouped together
    # we want reduncancy, so we want size overlap between them to be ~15%
    # which is higher than 0.9 required for cluster construction


    final_size_bins = []

    prevtop = 0

    for idx, size in enumerate(top_sizes):
        # first chunk
        if idx == 0:
            final_size_bins.append((minlen, DefaultValues.DEFAULT_MIN_BIN_SIZE))
            prevtop = DefaultValues.DEFAULT_MIN_BIN_SIZE

        else:
            # get previous top size
            prevsize = top_sizes[idx - 1]

            if size < 100000:
                multiplier = 0.35
            elif size < 250000:
                multiplier = 0.25
            elif size < 500000:
                multiplier = 0.15
            elif size < 1000000:
                multiplier = 0.1
            else:
                multiplier = 0.1

            #we want current size to be at least 10% higher than previous top size
            #we also want it to be at least 30% bigger than previous top size
            if (size > prevtop + prevtop * multiplier):
                # we want to go 11% below that (extra 1% just in case on top of 90% length ratio prerequisite for chaining contigs in clusters)
                prevsize = int(prevtop - prevtop * 0.10)
                final_size_bins.append((prevsize, size))
                prevtop = size


    logging.info('Dividing contigs into appropriate size bins with seed {}. Size distributions: {}'.format(NUMBER_BINS, final_size_bins))

    return final_size_bins

def analyse_assembly_contig_size_thread(queue_in, queue_out, processed_contigsize_dfs, assemblies, minlen):
    while True:
        idx = queue_in.get(block=True, timeout=None)
        if idx == None:
            break

        assembly = assemblies[idx]

        size, cname = [], []

        for name, seq, _ in readfq(open(assembly)):
            seq_len = len(seq)
            del seq
            if seq_len >= DefaultValues.DEFAULT_MIN_BIN_SIZE:
                size.append(seq_len)

            cname.append(name)

        if len(cname) != len(set(cname)):
            logging.error('Found duplicate header! Ensure your FASTA headers are unique within assemblies.')
            sys.exit(1)

        size_distro_df = pd.DataFrame({'Size': size})
        size_distro_df['Size'] = size_distro_df['Size'].astype(int)

        processed_contigsize_dfs.append(size_distro_df)
        queue_out.put(assembly)


def report_progress_assembly_threads(total_assemblies, queueIn):
    """Report number of processed assemblies."""

    processed = 0

    while True:
        thread = queueIn.get(block=True, timeout=None)
        if thread[0] == None:
            if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                sys.stdout.write('\n')
                sys.stdout.flush()
            break

        processed += 1

        if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
            statusStr = '    Finished processing %d of %d (%.2f%%) assembly contig sizes.' % (
                processed, total_assemblies, float(processed) * 100 / total_assemblies)
            sys.stdout.write('\r{}'.format(statusStr))
            sys.stdout.flush()

def setup_and_format_assemblies(min_len, inputfilelist, outputfolder, size_bins):

    fileManager.make_sure_path_exists(outputfolder)

    supp_assembly_files = []

    if len(size_bins) > 1:
        for floor, ceiling in size_bins:
            size_bin_folder = os.path.join(outputfolder, '{}_FLOOR_CEILING_{}_{}'.format(DefaultValues.SIZE_BIN_NAME, floor, ceiling))
            fileManager.make_sure_path_exists(size_bin_folder)
            supp_assembly_file = os.path.join(size_bin_folder, '{}_FLOOR_CEILING_{}_{}_.fasta'.format(DefaultValues.SIZE_BIN_NAME, floor, ceiling))

            supp_assembly_files.append(supp_assembly_file)


    full_bin_folder = os.path.join(outputfolder,
                                       '{}'.format(DefaultValues.MASTER_ASSEMBLY_FILE))
    fileManager.make_sure_path_exists(full_bin_folder)

    master_assembly_file = os.path.join(full_bin_folder, DefaultValues.COMBINED_ASSEMBLY_NAME)



    processed = 0

    superdict = dict()
    for entry in supp_assembly_files:
        superdict[entry] = dict()

    with open(master_assembly_file, 'w') as masterfile:

        for assembly in inputfilelist:
            for name, seq, _ in readfq(open(assembly)):
                try:
                    seq_len = len(seq)
                    if seq_len >= min_len:

                        basefile = os.path.splitext(os.path.basename(assembly))[0]
                        header = '{}{}{}{}{}'.format(basefile, DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR, seq_len, DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR, name).split(' ', 1)[0]

                        # always write to master file as we need it later regardless
                        masterfile.write('>' + header + '\n')
                        masterfile.write(seq + '\n')

                        #now go through supp_assembly_folders and check if we should write
                        for sub_assembly in supp_assembly_files:
                            floor = int(sub_assembly.split('_')[-3])
                            ceiling = int(sub_assembly.split('_')[-2])
                            if seq_len >= floor and seq_len <= ceiling:
                                superdict[sub_assembly][f'>{header}'] = seq

#                                with open (sub_assembly, 'a') as fout:
#                                    fout.write('>' + header + '\n')
#                                    fout.write(seq + '\n')
#                                fout.close()

                except Exception as e:
                    logging.error('Could not write combined assembly file: {}'.format(e))
                    sys.exit(1)

            processed += 1


            # Write the creted dicts and reset them
            for sub_assembly in supp_assembly_files:
                with open(sub_assembly, 'a') as fout:
                    for k, v in superdict[sub_assembly].items():
                    #fout.write(str(superdict[sub_assembly]))
                        fout.write(k + '\n')
                        fout.write(v + '\n')
                fout.close()

            for entry in supp_assembly_files:
                superdict[entry] = dict()


            if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                statusStr = '    Finished processing %d of %d (%.2f%%) assemblies' % (
                    processed, len(inputfilelist), float(processed) * 100 / len(inputfilelist))
                sys.stdout.write('\r{}'.format(statusStr))

                if processed == len(inputfilelist):
                    if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                        sys.stdout.write('\n')
                        sys.stdout.flush()
                else:
                    sys.stdout.flush()

        #close generated assembly files
        masterfile.close()
        #for fd in file_descriptors:
        #    fd.close()

    #re
    return supp_assembly_files, master_assembly_file



def create_final_contigs_assembly(contigs_file, main_assembly):

    final_assembly_file = os.path.join(os.path.dirname(contigs_file), DefaultValues.COMBINED_ASSEMBLY_NAME)

#    contigs = dict()
#    prev_assem = 'None'

    with open(final_assembly_file, 'w') as masterfile:

        #read in contigs names

        contig_names = pd.read_csv(contigs_file, sep='\t', names=['Contigs'])
#        contig_names['Assembly'] = contig_names['Contigs'].apply(lambda x: x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0])

        contig_dict = dict(zip(contig_names['Contigs'], contig_names.index))

        for name, seq, _ in readfq(open(os.path.abspath(main_assembly))):
            try:
                #sample = name.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0]


                # seq_len = len(seq)
                # basefile = os.path.splitext(os.path.basename(assembly))[0]
                # header = \
                # '{}{}{}{}{}'.format(basefile, DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR, seq_len,
                #                     DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR, name).split(' ',
                #                                                                               1)[0]
                #
                if name in contig_dict.keys():

                    masterfile.write('>' + name + '\n')
                    masterfile.write(seq + '\n')


            except Exception as e:
                logging.error('Could not write combined assembly file: {}'.format(e))
                sys.exit(1)
    masterfile.close()
    return final_assembly_file






def write_clusters_to_file(initial_concat_assembly, outdir, graphs, infos):
    records = {}
    master_records = {}
    cluster_out_folder = os.path.join(outdir, 'clusters')
    fileManager.make_sure_path_exists(cluster_out_folder)

    filtered_out = []

    edge_excluded = 0
    mean_excluded = 0
#    other_contigs_in_cluster = {}

    dict_edge_occur = {}
    dict_total_node = {}
    dict_mean_size = {}
    dict_std_size = {}


    for i, subgraph in enumerate(graphs):

        # identify rep contig:



        dcent = nx.degree_centrality(subgraph)
        best_degree_centrality = {k: dcent[k] for k in subgraph}
        final_selected_contig = max(best_degree_centrality, key=best_degree_centrality.get)


        # for imperfect contigs

        if infos[infos['ID'] == i]['Structure'].values[0] == 'Imperfect':
            # how many contigs does rep contig align to?
            # filter loose clusters

            edge_occurences = sum(x.count(final_selected_contig) for x in list(subgraph.edges))
            total_nodes = len(list(subgraph.nodes))


            if total_nodes < 11:  # apply n-1 calculation - rep contig needs to align to n-1 other contigs
                if edge_occurences < total_nodes - 1:
                    filtered_out.append(i)
                    #logging.debug('ID_{}: Total_edges: {}; Edge_occurences: {}'.format(i, total_edges, edge_occurences))
                    edge_excluded += 1


            elif (edge_occurences / total_nodes) < DefaultValues.TIGHTNESS_CLUSTER_CUTOFF:
                filtered_out.append(i)

            # next, assesd std size deviation from mean:
            #rep_size = final_selected_contig.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[1]

            size_array = []
            for s in list(subgraph.nodes):
                size_array.append(int(s.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[1]))

            #if abs(pct_change(float(rep_size),
            #                   float(np.mean(size_array))))
            if np.std(np.array(size_array)) / np.mean(np.array(size_array)) > DefaultValues.ALLOWED_SIZE_DEVIATION_WITHIN_CLUSTER\
                    or np.std(np.array(size_array)) > np.mean(np.array(size_array)):
                filtered_out.append(i)
                mean_excluded += 1
                #logging.debug('ID_{}: Size: {} Mean: {}'.format(i, rep_size, np.mean(size_array)))


        ''' DEBUGGING PURPOSES START'''
        edge_occurences = sum(x.count(final_selected_contig) for x in list(subgraph.edges))
        total_nodes = len(list(subgraph.edges))
        size_array = []
        for s in list(subgraph.nodes):
            size_array.append(int(s.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[1]))

        dict_edge_occur[i] = edge_occurences
        dict_total_node[i] = total_nodes
        dict_mean_size[i] = np.mean(np.array(size_array))
        dict_std_size[i] = np.std(np.array(size_array))
        ''' DEBUGGING PURPOSES END'''

        #if i not in filtered_out:
        records[str(final_selected_contig)] = i

        for name in dcent.keys():
            master_records[name] = i

    # remove ID's filtered out
    #infos = infos[~infos['ID'].isin(filtered_out)]

    infos_size = dict(zip(infos['ID'], infos['Average_Bp_Size']))
    structure = dict(zip(infos['ID'], infos['Structure']))


    for name, seq, _ in readfq(open(initial_concat_assembly)):
        if str(name) in records.keys():
            ID = records[name]

            destination = '{}/clusters/{}_cluster_ID_{}_Size_{}.fna'.format(outdir, structure[ID], ID, infos_size[ID])
            name = name.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[-1]
            write_fasta({name: seq}, destination)
        else:
            if name not in master_records.keys():
                master_records[name] = 'Not_in_cluster'

        # restore original contig
    df = pd.DataFrame(data=master_records.items(), columns=['ContigName', 'Cluster_ID'])

    df['Sample'] = df['ContigName'].apply(lambda x: x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0])
    df['Size'] = df['ContigName'].apply(lambda x: x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[1])
    df['ContigName'] = df['ContigName'].apply(lambda x: x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[2])

    cluster_contigs_info = df[df['Cluster_ID'] != 'Not_in_cluster']
    leftover_contigs_info = df[df['Cluster_ID'] == 'Not_in_cluster']

    logging.debug('Excluded due to within-cluster tightness requirements: {}'.format(edge_excluded))
    logging.debug('Excluded due to within-cluster mean requirements: {}'.format(mean_excluded))

    ''' DEBUGGING PURPOSES'''
    filtered_out = [] #TODO: revert to filtering out once proper parameters are determined

    return cluster_out_folder, cluster_contigs_info, leftover_contigs_info, filtered_out, \
           dict_edge_occur, dict_total_node, dict_mean_size, dict_std_size




def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last: break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs);  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def read_nucleotide_sequences(nucleotide_file):
    nucleotide_sequences = {}
    for name, seq, _ in readfq(open(nucleotide_file)):
        nucleotide_sequences[name] = seq
    return nucleotide_sequences


def write_fasta(seq, outputFile):
    '''write sequences to FASTA file'''
    fout = open(outputFile, 'w')

    for seqId, seq in seq.items():
        fout.write('>' + seqId + '\n')
        fout.write(seq + '\n')
    fout.close()

