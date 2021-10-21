import os
import pandas as pd
import sys
import logging
import shutil
import networkx as nx

from recurm import fileManager
from recurm.defaultValues import DefaultValues


def setup_and_format_assemblies(min_len, inputfilelist, outputfolder):
    fileManager.make_sure_path_exists(outputfolder)
    master_assembly_file = '{}/{}'.format(outputfolder, DefaultValues.COMBINED_ASSEMBLY_NAME)
    nucleotide_sequences = {}
    processed_headers = {}
    processed = 0

    with open(master_assembly_file, 'w') as fout:

        for assembly in inputfilelist:
            for name, seq, _ in readfq(open(assembly)):
                try:
                    seq_len = len(seq)
                    if seq_len >= min_len:
                        basefile = os.path.splitext(os.path.basename(assembly))[0]
                        header = '{}{}{}{}{}'.format(basefile, DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR, seq_len, DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR, name).split(' ', 1)[0]
                        if header in processed_headers.keys():
                            logging.error('Found duplicate header! Ensure your FASTA headers are unique.')
                            sys.exit(1)
                        processed_headers[header] = 1
                        fout.write('>' + header + '\n')
                        fout.write(seq + '\n')

                except Exception as e:
                    logging.error('Could not write combined assembly file: {}'.format(e))
                    sys.exit(1)

            processed += 1

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
        fout.close()
    return master_assembly_file

def write_clusters_to_file(initial_concat_assembly, outdir, graphs, infos):
    records = {}
    master_records = {}
    fileManager.make_sure_path_exists(os.path.join(outdir, 'clusters'))
    other_contigs_in_cluster = {}
    for i, subgraph in enumerate(graphs):
        dcent = nx.degree_centrality(subgraph)
        best_degree_centrality = {k: dcent[k] for k in subgraph}
        final_selected_contig = max(best_degree_centrality, key=best_degree_centrality.get)
        records[str(final_selected_contig)] = i

        for name in dcent.keys():
            master_records[name] = i


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

    df[df['Cluster_ID'] == 'Not_in_cluster'].to_csv('{}/leftover_contigs_information.tsv'.format(outdir), sep='\t', index=False)
    df[df['Cluster_ID'] != 'Not_in_cluster'].to_csv('{}/cluster_contigs_information.tsv'.format(outdir), sep='\t', index=False)





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


# def split_by_contig(inputfile):
#     nucleotide_sequences = {}
#     for name, seq, _ in readfq(open(inputfile)):
#         if name in lookupdict:
#             ID = lookupdict[name]
#
#             destination = '/scratch/microbiome/n10851381/recurm/HMP2/extracted_clusters/ID_{}/all_contigs/'.format(ID)
#             make_sure_path_exists(destination)
#             outname = '{}/{}.fna'.format(destination, name)
#             write_fasta({name: seq}, outname)
