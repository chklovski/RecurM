import sys
import logging
import subprocess
import os
from os import listdir

#set numpy to use a single thread
os.environ["NUMEXPR_MAX_THREADS"] = '1'
import warnings
import pandas as pd
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import hashlib
import multiprocessing as mp
import scipy
import networkx as nx

Hash = hashlib.sha512
MAX_HASH_PLUS_ONE = 2**(Hash().digest_size * 8)


from recurm.defaultValues import DefaultValues
from recurm import fileManager
from recurm import FASTA_manager

class Mapper():
    def __init__(self, threads):
        self.nthreads = threads
        self.__check_for_minimap2()
        self.__check_for_fastANI()

    def __check_for_minimap2(self):
        # Assume that a successful minimap2 -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['minimap2', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            logging.error("Make sure minimap2 is on your system path.")
            sys.exit(1)

    def __check_for_fastANI(self):
        # Assume that a successful fastANI -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['fastANI', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            logging.error("Make sure fastANI is on your system path.")
            sys.exit(1)

    def samplify(self, a, b):
        if a.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0] == \
                b.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0]:
            return 'SameSample'
        else:
            return 'DifferentSample'

    def hashsum(self, a, b):
        return self.str_to_probability(a) * self.str_to_probability(b)
    
    def str_to_probability(self, in_str):
        """Return a reproducible uniformly random float in the interval [0, 1) for the given string."""
        seed = in_str.encode()
        hash_digest = Hash(seed).digest()
        hash_int = int.from_bytes(hash_digest, 'big')  # Uses explicit byteorder for system-agnostic reproducibility
        return hash_int / MAX_HASH_PLUS_ONE  # Float division


    def calcLR(self, lenx, leny, alin_len, mismatch):
        if lenx > leny:
            longest = lenx
            shortest = leny
        else:
            longest = leny
            shortest = lenx

        LR = shortest / longest
        return LR

    def calcANI(self, lenx, leny, alin_len, mismatch):
        if lenx > leny:
            longest = lenx
            shortest = leny
        else:
            longest = leny
            shortest = lenx

        ANI = mismatch / alin_len
        return ANI

    def calcAR(self, lenx, leny, alin_len, mismatch):
        if lenx > leny:
            longest = lenx
            shortest = leny
        else:
            longest = leny
            shortest = lenx

        AR = alin_len / longest
        return AR

    def calcARshort(self, lenx, leny, alin_len, mismatch):
        if lenx > leny:
            longest = lenx
            shortest = leny
        else:
            longest = leny
            shortest = lenx

        ARshort = alin_len / shortest
        return ARshort


    def line_count(self, filename):
        return int(subprocess.check_output(['wc', '-l', filename]).split()[0])

    def overlap_percent(self, start1, end1, start2, end2):
        """how much does the range (start1, end1) overlap with (start2, end2)"""
        overlap = max(max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0)
        len1 = end1 - start1
        len2 = end2 - start2
        if len1 > len2:
            shortest = len2
        else:
            shortest = len1

        if overlap == 0:
            return 0
        else:
            return (overlap/shortest) * 100




class AllVsAllMapper(Mapper):


    def mapAVA(self, combined_assemblies, map_outdir):
        mapping_input = combined_assemblies
        mapping_out = os.path.join(map_outdir, DefaultValues.FIRST_PASS_MAPPING_NAME)

        try:
            cmd = "minimap2 -x ava-ont -t {} {} {} -v1 | split --additional-suffix=.paf -d --line-bytes={} - {}" \
                .format(self.nthreads, mapping_input, mapping_input, DefaultValues.PAF_CHUNK_SIZE, mapping_out)

            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished Running Minimap2')
        except Exception as e:
            logging.error('An error occured while running Minimap2: {}'.format(e))
            sys.exit(1)

        return [os.path.join(map_outdir, mapping_file) for mapping_file in listdir(map_outdir)if mapping_file.endswith('.paf')]



    def read_hashes_and_process_alignments(self, mapping_dir, hash_file, LRcutoff, ARcutoff, ANIcutoff, outfile):
        
        #read in hashes one by one;

        logging.info('Iterating through alignment data a second time to extract circular and complex alignments')

        hash_length = self.line_count(hash_file)

        hash_retrieve_params = []

        with open(hash_file) as f:

            index_buffer = []
            hash_buffer = []
            prev_sample = None

            lastline = False

            for idx, line in enumerate(f):

                index = line.split('\t')[0]
                hash = line.split('\t')[1]
                sample = line.split('\t')[2]

                if idx + 1 == hash_length:
                    lastline = True


                if prev_sample is None:
                    prev_sample = sample

                else:
                    if prev_sample != sample or lastline:  # we hit new sample or end of file

                        dup_hashes = pd.DataFrame({'OriginalIndex': index_buffer, 'Hash': hash_buffer})
                        dup_hashes = dup_hashes.astype({'OriginalIndex': 'int64', 'Hash': 'string'})

                        paf = '{}/{}.paf'.format(mapping_dir, os.path.splitext(
                            os.path.basename(prev_sample if not lastline else sample))[0])

                        map_dict = dict(zip(dup_hashes['OriginalIndex'], dup_hashes['Hash']))

                        hash_retrieve_params.append((paf, map_dict))

                        del dup_hashes

                        index_buffer = []
                        hash_buffer = []

                        prev_sample = sample

                index_buffer.append(index)
                hash_buffer.append(hash)


        # after collating indexes to retrieve, pool it and use multiprocessing for IO
        logging.info('Using {} threads to retrieve multiple alignments from alignment files.'.format(self.nthreads))


        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for idx, params in enumerate(hash_retrieve_params):
            workerQueue.put(idx)

        for _ in range(self.nthreads):
            workerQueue.put(None)


        processed_multialignment_dfs = mp.Manager().list()

        try:
            calcProc = []
            for _ in range(self.nthreads):
                calcProc.append(
                    mp.Process(target=self.__set_up_processing_sample_with_multiple_aligns,
                               args=(workerQueue, writerQueue, processed_multialignment_dfs, hash_retrieve_params, LRcutoff)))
            writeProc = mp.Process(target=self.__report_progress_paf,
                                   args=(len(hash_retrieve_params), writerQueue))

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

#        with mp.Pool(self.nthreads) as pool:
#            extracted_paf_alignments = pool.starmap(self.__process_sample_with_multiple_aligns, hash_retrieve_params)

        ''' Done processing sample '''

        extracted_paf_alignments = pd.concat(processed_multialignment_dfs)
        extracted_paf_alignments.reset_index(drop=True, inplace=True)


        extracted_aligns = extracted_paf_alignments[extracted_paf_alignments['HashCount'] == 2].copy()
        complex_aligns = extracted_paf_alignments[extracted_paf_alignments['HashCount'] > 2].copy()

        #complex aligns are fine as long as basic requirements of circular aligns are preserved

        if len(complex_aligns) > 0:
            append_to_simple_double_alignments = []

            for hashgroup in complex_aligns['Hash'].unique():
                subalign = complex_aligns[complex_aligns['Hash'] == hashgroup]

                # extract the two longest aligns and append them to be checked
                # if the two longest aren't on either side, then it'll be ignored later down the road
                subalign = subalign.sort_values(by='AR', ascending=False)
                subalign.reset_index(drop=True, inplace=True)


                append_to_simple_double_alignments.append(subalign.head(2))

            if len(append_to_simple_double_alignments) > 0:
                append_to_simple_double_alignments = pd.concat(append_to_simple_double_alignments)
                extracted_aligns = pd.concat([extracted_aligns, append_to_simple_double_alignments])



        extracted_aligns = extracted_aligns.sort_values(by='Hash')
        extracted_aligns.reset_index(drop=True, inplace=True)

        index = 0
        real_circs = []

        #Now go through the double-alignment pairs and make sure that they point in the right direction.


        while index < len(extracted_aligns):
        


            subsub = extracted_aligns.iloc[index:index + 2, :] #We've pre-filtered to 2 alignments only

            #if subsub['ANI'].values[0] < ANIcutoff or subsub['ANI'].values[1] < ANIcutoff:
            #    continue #skip
            if len(subsub[subsub['ANI'] > ANIcutoff]) == len(subsub): #all alignments pass ANI cutoff

                if subsub['AR'].values[0] + subsub['AR'].values[1] > ARcutoff:
                    # direction is conserved between alignments, ie. no palindromic alignments

                    astart = subsub[3].values[0]
                    aend = subsub[4].values[0]
                    bstart = subsub[8].values[0]
                    bend = subsub[9].values[0]
                    lena = subsub[2].values[0]
                    lenb = subsub[7].values[0]

                    astart2 = subsub[3].values[1]
                    aend2 = subsub[4].values[1]
                    bstart2 = subsub[8].values[1]
                    bend2 = subsub[9].values[1]

                    if (lena - aend) > (lena - aend2):
                        firstHalfA = True
                    else:
                        firstHalfA = False

                    if (lenb - bend) > (lenb - bend2):
                        firstHalfB = True
                    else:
                        firstHalfB = False

                    overlap1 = self.overlap_percent(astart, aend, astart2, aend2)
                    overlap2 = self.overlap_percent(bstart, bend, bstart2, bend2)



                    if len(subsub[5].unique()) == 1: # This means + to + or - to -
                    #no more than 10% overlap:
                    #TODO: check this aligns with AR to match cutoffs
                        if overlap1 < 10 and overlap2 < 10:

                            if subsub[5].values[0] == '-':
                                if firstHalfA and firstHalfB:
                                    real_circs.append(subsub)
                                elif not firstHalfA and not firstHalfB:
                                    real_circs.append(subsub)
                            #investigate other types:
                                #elif not firstHalfA and firstHalfB:
                                #    print('- dir singular: not firstHalfA and firstHalfB')
                                #    print(subsub.to_string())
                                #elif firstHalfA and not firstHalfB:
                                #    print('- dir singular: firstHalfA and not firstHalfB')
                                #    print(subsub.to_string())


                            elif subsub[5].values[0] == '+':
                                if not firstHalfA and firstHalfB:
                                    real_circs.append(subsub)
                                elif firstHalfA and not firstHalfB:
                                    real_circs.append(subsub)
                            #investigate other types:
                                #elif not firstHalfA and not firstHalfB:
                                    #print('+ dir singular: not firstHalfA and not firstHalfB')
                                    #print(subsub.to_string())
                                #elif firstHalfA and firstHalfB:
                                    #print('+ dir singular: firstHalfA and firstHalfB')
                                    #print(subsub.to_string())


                    elif len(subsub[5].unique()) == 2: # This means + to -
                        if overlap1 < 10 and overlap2 < 10:
                              if not firstHalfA and firstHalfB:
                                    real_circs.append(subsub)
                              elif firstHalfA and not firstHalfB:
                                    real_circs.append(subsub)

            index += 2


        if(len(real_circs)) > 1:
            result = pd.concat(real_circs)
    
            result.reset_index(drop=True, inplace=True)
            del result['HashCount']
            del result['Sample']
            del result['Hash']
    
            logging.info('Writing circular alignments to file.')
            result.to_csv(outfile, sep='\t', index=False, header=False)

            
            return True
        else:
            return False




    def extract_hashes_with_multiple_hits(self, infile):
        outfile = '{}_extracted'.format(infile)
        with open(outfile, 'a') as out:
            with open(infile) as f:
                prev_hash = 0
                line_buffer = []
                for line in f:
                    hash_value = line.split('\t')[1]
                    if prev_hash != hash_value:  # new hash value
                        if len(line_buffer) > 1:
                            for entry in line_buffer:
                                out.write(entry)
                        line_buffer = []
                        prev_hash = hash_value
                        line_buffer.append(line)
                    else:
                        line_buffer.append(line)
        return outfile

    def subfilter_paf_mapping_results(self, threads, paf_list, align_dir, hash_dir, chunksize, LRcutoff, ARcutoff, ANIcutoff):

        #spin off into multiple threads

        self.threads_per_bin = max(1, int(threads / len(paf_list)))
        logging.info('Processing alignment data.')
        logging.info("Data is split into {} chunks of {} each. Analysing with {} threads:"
                     .format(len(paf_list), DefaultValues.PAF_CHUNK_SIZE, threads))

        # process each bin in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for paf in paf_list:
            workerQueue.put(paf)

        for _ in range(threads):
            workerQueue.put(None)

        try:
            calcProc = []
            for _ in range(threads):
                calcProc.append(
                    mp.Process(target=self.__set_up_paf_thread, args=(workerQueue, writerQueue, align_dir, hash_dir, chunksize, LRcutoff, ARcutoff, ANIcutoff)))
            writeProc = mp.Process(target=self.__report_progress_paf, args=(len(paf_list), writerQueue))

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



        full_hash_file = os.path.join(hash_dir, DefaultValues.HASH_FILE_NAME)
        full_passed_file = os.path.join(align_dir, DefaultValues.SECOND_PASS_NAME)

        logging.info('Merging hash files and alignments kept after first pass.')
        try:
            #TODO: There has got to be a better way of doing this...
            cmd = "for i in {0}/*.hash; do awk '{{print $0 @\\t@ FILENAME;}}' $i >> {1}/{2}; rm $i; done"\
                .format(hash_dir, hash_dir, DefaultValues.HASH_FILE_NAME).replace('@', '"')
            logging.debug(cmd)
            subprocess.call(cmd, shell=True)

            cmd = "for i in {0}/*.passed; do cat $i >> {1}/{2}; rm $i; done"\
                .format(align_dir, align_dir, DefaultValues.SECOND_PASS_NAME)
            logging.debug(cmd)
            subprocess.call(cmd, shell=True)

            logging.debug('Finished concatenating hashes using awk and cat')

        except Exception as e:
            logging.error('An error occured while concatenating hashes and extracted alignments: {}'.format(e))
            sys.exit(1)

        return full_hash_file, full_passed_file



    def __set_up_paf_thread(self, queue_in, queue_out, align_dir, hash_dir, chunksize, LRcutoff, ARcutoff, ANIcutoff):

        while True:
            paf = queue_in.get(block=True, timeout=None)
            if paf == None:
                break

            basename = os.path.splitext(os.path.basename(paf))[0]
            out_align = os.path.join(align_dir, basename + '.passed')
            out_hash = os.path.join(hash_dir, basename + '.hash')


            df = pd.read_csv(paf, chunksize=chunksize,
                             names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], sep='\t')

            index = 0

            for chunk in df:

                chunk['AR'] = chunk.apply(lambda row: self.calcAR(row[2], row[7], row[11], row[10]), axis=1)
                chunk['LR'] = chunk.apply(lambda row: self.calcLR(row[2], row[7], row[11], row[10]), axis=1)
                chunk['ANI'] = chunk.apply(lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)
                chunk['Sample'] = chunk.apply(lambda row: self.samplify(row[1], row[6]), axis=1)

                chunk = chunk[chunk['Sample'] == 'DifferentSample']

                if len(chunk) > 1:
                    chunk['Hash'] = chunk.apply(lambda row: self.hashsum(row[1], row[6]), axis=1)

                    # multiple alignments need at least ~20% length coverage on a side to be a decent match
                    chunk2 = chunk[chunk['AR'] > 0.2]

                    chunk2.iloc[:, -1].to_csv(out_hash, mode='a', sep='\t', header=False)
                    del chunk['Hash']
                    del chunk['Sample']
                    del chunk2

                    chunk = chunk[chunk['AR'] > ARcutoff]
                    chunk = chunk[chunk['LR'] > LRcutoff]
                    chunk = chunk[chunk['ANI'] > ANIcutoff]

                    chunk[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'AR', 'LR', 'ANI']].to_csv(out_align, mode='a', sep='\t', header=False)

            queue_out.put(paf)


    def __set_up_processing_sample_with_multiple_aligns(self, queue_in, queue_out, align_list, param_list, LRcutoff):
        ''' Begin processing sample '''

        while True:
            idx = queue_in.get(block=True, timeout=None)
            if idx == None:
                break

            paf, map_dict = param_list[idx]
            keep_list = []

            chunksize = DefaultValues.FILTER_CHUNKSIZE

            orig_aligns = pd.read_csv(paf, chunksize=chunksize,
                                      names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                                      usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], sep='\t')

            for chunk in orig_aligns:
                sub = chunk[chunk.index.isin(map_dict.keys())]
                if len(sub) > 0:
                    keep_list.append(sub)

            extracted_aligns = pd.concat(keep_list)
            extracted_aligns.reset_index(drop=True, inplace=True)

            extracted_aligns['LR'] = extracted_aligns.apply(
                lambda row: self.calcLR(row[2], row[7], row[11], row[10]), axis=1)
            # aligned ratio may be short, but length ratio should be solid
            extracted_aligns = extracted_aligns[extracted_aligns['LR'] >= LRcutoff]
            # make sure again:
            extracted_aligns['Sample'] = extracted_aligns.apply(
                lambda row: self.samplify(row[1], row[6]), axis=1)
            extracted_aligns = extracted_aligns[extracted_aligns['Sample'] == 'DifferentSample']

            extracted_aligns = extracted_aligns.astype(
                {2: 'int64', 3: 'int64', 4: 'int64', 7: 'int64', 8: 'int64', 9: 'int64', 10: 'int64', 11: 'int64'})

            extracted_aligns['Hash'] = extracted_aligns.apply(lambda row: self.hashsum(row[1], row[6]),
                                                                              axis=1)

            hashdict = dict(extracted_aligns['Hash'].value_counts())
            extracted_aligns['HashCount'] = extracted_aligns['Hash'].apply(lambda x: hashdict[x])

            extracted_aligns['AR'] = extracted_aligns.apply(
                lambda row: self.calcAR(row[2], row[7], row[11], row[10]), axis=1)
            extracted_aligns['ANI'] = extracted_aligns.apply(
                lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)


            align_list.append(extracted_aligns)

            queue_out.put([idx])

    def __report_progress_paf(self, total_pafs, queueIn):
        """Report number of processed bins."""

        processed = 0

        while True:
            paf = queueIn.get(block=True, timeout=None)
            if paf[0] == None:
                if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                    sys.stdout.write('\n')
                    sys.stdout.flush()
                break

            processed += 1

            if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                statusStr = '    Finished processing %d of %d (%.2f%%) paf alignment chunks.' % (
                    processed, total_pafs, float(processed) * 100 / total_pafs)
                sys.stdout.write('\r{}'.format(statusStr))
                sys.stdout.flush()




class ReferenceMapper(Mapper):
    pass

class WithinClusterMapper(Mapper):

    def map_clusters_to_each_other(self, cluster_folder, outdir):
        all_contigs = fileManager.list_assembly_folder(cluster_folder, '.fna')
        combined_contig_file = os.path.join(outdir, DefaultValues.COMBINED_CONTIGS_FILE)

        with open(combined_contig_file, 'w') as fout:

            for contig in all_contigs:
                for name, seq, _ in FASTA_manager.readfq(open(contig)):
                    try:
                        basefile = os.path.splitext(os.path.basename(contig))[0]
                        #CLUSTERINFO <SEP> ORIGINALCONTIGNAME
                        header = \
                        '{}{}{}'.format(basefile, DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR, name)
                        fout.write('>' + header + '\n')
                        fout.write(seq + '\n')

                    except Exception as e:
                        logging.error('Could not write combined contigs file: {}'.format(e))
                        sys.exit(1)

        # for now we assume that we can generate a single .paf from all contigs and it won't be too big
        # this may change - to be tested

        mapping_output = '{}/{}'.format(outdir, DefaultValues.COMBINED_CONTIGS_MAPPING_FILE)
        mapping_input = combined_contig_file

        try:
            cmd = "minimap2 -x ava-ont -t {} {} {} -v1 > {}" \
                .format(self.nthreads, mapping_input, mapping_input, mapping_output)

            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished Running Minimap2')
        except Exception as e:
            logging.error('An error occured while running Minimap2: {}'.format(e))
            sys.exit(1)

        return mapping_output

    def collapse_related_clusters(self, mapping_file, cluster_folder,
                                  LRcutoff, ARcutoff, ANIcutoff, outdir,
                                  cluster_information, cluster_contigs_info, leftover_contigs_info, nocollapse):
        paf = pd.read_csv(mapping_file,
                             names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], sep='\t')

        #remove self-alignment
        paf = paf[paf[1] != paf[6]]

        paf['AR'] = paf.apply(lambda row: self.calcAR(row[2], row[7], row[11], row[10]), axis=1)
        paf['LR'] = paf.apply(lambda row: self.calcLR(row[2], row[7], row[11], row[10]), axis=1)
        paf['ANI'] = paf.apply(lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)
        paf['ARshort'] = paf.apply(lambda row: self.calcARshort(row[2], row[7], row[11], row[10]), axis=1)

        collapsee = []

        if not nocollapse:
            if len(paf[paf['ARshort'] > ARcutoff]) > 0:
    
                sel = paf[paf['ARshort'] > ARcutoff].copy()
                sel = sel[sel['ANI'] > ANIcutoff]
    
    
                #collapsee = []
                collapsed_into = []
                direction = []
    
                if len(sel) > 0:
                    logging.info('Collapsing related clusters.')
                    #indexes_for_collapse = sel.index.values
    
                    #sel = sel.sort_values(by=2, ascending=True)
                    #sel = sel.sort_values(by=7, ascending=True)
    
                    #collapse shorter cluster into longer cluster
                    sel.reset_index(drop=True, inplace=True)
    
                    for index, row in sel.iterrows():
                        if sel.at[index, 2] > sel.at[index, 7]:
                            # first element is longer
                            longer = sel.at[index, 1]
                            shorter = sel.at[index, 6]
                        else:
                            shorter = sel.at[index, 1]
                            longer = sel.at[index, 6]
    
                        shorter = shorter.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0]
                        longer = longer.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0]
    
    
                        collapsee.append(shorter)
                        collapsed_into.append(longer)
                        direction.append('-->')
    
                        #move_from = os.path.join(cluster_folder, longer + '.fna')
                        #move_to = os.path.join(cluster_folder, shorter + '.fna')
    
                        logging.debug('Collapsing {} --> {}'.format(shorter, longer))
    
                        try:
                            os.remove(os.path.join(cluster_folder, shorter + '.fna'))
                        except FileNotFoundError:
                            if shorter in collapsee:
                                continue
                            else:
                                logging.error('File not found: {}. Unable to continue.'.format(os.path.join(cluster_folder, shorter + '.fna')))
                                sys.exit(1)
    
    
                    collapsed_log = pd.DataFrame({'Collapsee': collapsee, 'Direction': direction, 'Collapsed_Into': collapsed_into})
                    collapsed_log_file = os.path.join(outdir, DefaultValues.COLLAPSED_LOG_FILE)
                    collapsed_log.to_csv(collapsed_log_file, sep='\t', index=False, header=False)
    
                    logging.info('Found and collapsed {} related clusters.'.format(len(set(collapsee))))
    
                    #update cluster information
                    remove_ids = [int(x.split('ID_')[-1].split('_')[0]) for x in collapsee]
                    new_ids = [int(x.split('ID_')[-1].split('_')[0]) for x in collapsed_into]
    
                    re_id = dict(zip(remove_ids, new_ids))
    
                    cluster_contigs_info['Cluster_ID']  = cluster_contigs_info['Cluster_ID'].\
                        apply(lambda x: 'Collapsed into {}'.format(re_id[int(x)]) if int(x) in re_id.keys() else x)
                    cluster_information = cluster_information[~cluster_information['ID'].isin(remove_ids)]
    
                    #append_to_leftovercontigs = cluster_contigs_info[cluster_contigs_info['Cluster_ID'].isin(remove_ids)].copy()
                    #cluster_contigs_info = cluster_contigs_info[~cluster_contigs_info['Cluster_ID'].isin(remove_ids)]
                    #append_to_leftovercontigs['Cluster_ID'] = new_ids
                    #cluster_contigs_info = pd.concat([cluster_contigs_info, append_to_leftovercontigs])
    
    
                # elif len(paf[paf['ARshort'] > ARcutoff]) != len(sel):
                    # Is there a difference between mapping well with low ANI and mapping decently with low ANI?
    

        # Let's set a reasonable secondary cutoff and do FASTANI on clusters passing this.

        #remove entries that have been collapsed

        paf['A'] = paf[1].apply(lambda x: x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0])
        paf['B'] = paf[6].apply(lambda x: x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0])

        rest = paf[~paf['A'].isin(collapsee)]
        rest = rest[~rest['B'].isin(collapsee)]
        rest = rest[rest['ARshort'] > DefaultValues.SECOND_PASS_CONTIGS_AR_SHORT_CUTOFF]

        #we will run FASTANI on contigs that have a (reasonably) good mapping only - must pass "SECOND_PASS_CONTIGS_AR_SHORT_CUTOFF"

        toanalyse = set(list(rest[1].values) + list(rest[6].values))
        toanalyse = ["{}/{}.fna".format(cluster_folder, x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0]) for x in toanalyse]

        FASTANI_input = os.path.join(outdir, 'FASTANI.in')
        FASTANI_output = os.path.join(outdir, 'FASTANI.out')
        pd.DataFrame({'Contigs': toanalyse}).to_csv(FASTANI_input, sep='\t', index=False, header=False)



        try:
            cmd = "fastANI -t {} --ql {} --rl {} -o {} > /dev/null 2> /dev/null" \
                .format(self.nthreads, FASTANI_input, FASTANI_input, FASTANI_output)

            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished Running fastANI')
        except Exception as e:
            logging.error('An error occured while running fastANI: {}'.format(e))
            sys.exit(1)

        os.remove(FASTANI_input)

        try:
            fastANI_results = pd.read_csv(FASTANI_output, sep='\t', names=['Query', 'Ref', 'ANI', 'BlocksMatched', 'BlocksTotal'])
        except FileNotFoundError:
            logging.warning('fastANI failed to produce output. Assuming there are no related clusters.')
            return cluster_information, cluster_contigs_info, leftover_contigs_info, None, None

        os.remove(FASTANI_output)

        #remove self-matches
        fastANI_results = fastANI_results[fastANI_results['Query'] != fastANI_results['Ref']]
        fastANI_results['Query'] = fastANI_results['Query'].apply(lambda x: x.split('/')[-1].split('.fna')[0])
        fastANI_results['Ref'] = fastANI_results['Ref'].apply(lambda x: x.split('/')[-1].split('.fna')[0])

        #remove potentially spurious results?
        #fastANI_results = fastANI_results[~(fastANI_results['BlocksTotal'] > 9) & (fastANI_results['BlocksMatched'] < 2)]

        #from duplicate entries
        fastANI_results['check_string'] = fastANI_results.apply(lambda row: ''.join(sorted([row['Query'], row['Ref']])),
                                              axis=1)
        fastANI_results = fastANI_results.drop_duplicates('check_string')


        #final_fastANI_out = os.path.join(outdir, 'clusters_relatedness.tsv')

        #logging.info('Information for related clusters written to: {}.'.format(final_fastANI_out))
        #fastANI_results[['Query', 'Ref', 'ANI']].to_csv(final_fastANI_out, sep='\t', index=False)

        # determine groups of ANI-related clusters and write to results
        G = nx.from_pandas_edgelist(fastANI_results, source='Query', target='Ref', edge_attr='ANI')

        # Build adjacency matrix
        ANI_adjacency_df = pd.DataFrame(
            nx.adjacency_matrix(G, weight='ANI').todense(),
            index=G.nodes,
            columns=G.nodes,
            dtype=float  # Compatible dtype with NaN is needed
        )
        # Overwrite the values on the diagonal
        np.fill_diagonal(ANI_adjacency_df.values, np.NaN)

        ANI_defined_groups = (list(G.subgraph(c) for c in nx.connected_components(G)))

        group_dict = {}
        for idx, group in enumerate(ANI_defined_groups):
            for entry in list(group):
                id = int(entry.split('ID_')[-1].split('_')[0])
                group_dict[id] = idx

        cluster_information['Group'] = cluster_information['ID'].apply(lambda x: 'Group_{}'.format(group_dict[x]) if x in group_dict.keys() else 'None')

        #todo: generate output (ANI and dendrogram) for EACH GROUP

        return cluster_information, cluster_contigs_info, leftover_contigs_info, fastANI_results[['Query', 'Ref', 'ANI']], ANI_adjacency_df

