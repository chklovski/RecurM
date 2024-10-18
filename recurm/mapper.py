import sys
import logging
import subprocess
import os
from os import listdir
import pickle

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


    def mapAVA(self, combined_assemblies, map_outdir, min_contig_len, fast):
        mapping_input = combined_assemblies
        mapping_out = os.path.join(map_outdir, DefaultValues.FIRST_PASS_MAPPING_NAME)
        # minimap2 options:

        #-m - minimum matching bases (effectively). Should be set to ~40% of desired minimum contig size.
        # Drastically reduces minimap2 output that needs to be parsed
        # -K and -I regulates memory usage on large datasets. -K in particular can be higher if memory > 100GB is available
        # -k (k-mer size) and -w (minimizer window) speeds up mapping at cost of more memory.

        try:


            #m - min mapping length should be ~ 0.3/0.4 of min contig length to filter out tons of tiny useless alignments
            if not fast:
                cmd = "minimap2 -x ava-ont -m {} -t {} {} {} -v1 | split --additional-suffix=.paf -d --line-bytes={} - {}" \
                    .format(min_contig_len * 0.3, self.nthreads, mapping_input, mapping_input, DefaultValues.PAF_CHUNK_SIZE, mapping_out)
                    
            else:
                cmd = "minimap2 -k16  -Xw5 -m {} -g1000 -r1000  --max-chain-skip 500 --max-chain-iter 1000 -K 5G -I 20G -2 -t {} {} {} -v1 | split --additional-suffix=.paf -d --line-bytes={} - {}" \
                    .format(min_contig_len * 0.3, self.nthreads, mapping_input, mapping_input, DefaultValues.PAF_CHUNK_SIZE, mapping_out)

            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished Running Minimap2')
        except Exception as e:
            logging.error('An error occured while running Minimap2: {}'.format(e))
            sys.exit(1)

        return [os.path.join(map_outdir, mapping_file) for mapping_file in listdir(map_outdir) if mapping_file.endswith('.paf')]



    def read_hashes_and_process_alignments(self, multialign_dir, mapping_dir, hash_file, LRcutoff, ARcutoff, ANIcutoff, outfile, short_circular, keep_inversions):


        use_short_circular = short_circular

        hashfileindex = 1

        #read in hashes one by one;

        logging.debug('Iterating through alignment data a second time to extract circular and complex alignments')

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

                        hash_save_path = os.path.join(multialign_dir, f'{hashfileindex}.pkl')
                        hashfileindex += 1

                        hash_retrieve_params.append((hash_save_path))

                        with open(hash_save_path, 'wb') as f:
                            pickle.dump([paf, map_dict], f)

                        del dup_hashes
                        del map_dict

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
                               args=(workerQueue, writerQueue, processed_multialignment_dfs, hash_retrieve_params, LRcutoff, use_short_circular)))
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

        if len(processed_multialignment_dfs) > 0:
            extracted_paf_alignments = pd.concat(processed_multialignment_dfs)
            extracted_paf_alignments.reset_index(drop=True, inplace=True)
        else:
            return False


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

        circular_labels = []

        if len(extracted_aligns) < 2:
            return False


        while index < len(extracted_aligns):



            subsub = extracted_aligns.iloc[index:index + 2, :] #We've pre-filtered to 2 alignments only

            #if subsub['ANI'].values[0] < ANIcutoff or subsub['ANI'].values[1] < ANIcutoff:
            #    continue #skip
            if len(subsub[subsub['ANI'] > ANIcutoff]) == len(subsub): #all alignments pass ANI cutoff

                passed = False
                short = False

                if (subsub['AR'].values[0] + subsub['AR'].values[1] > ARcutoff):
                    passed = True

                elif (subsub['ARshort'].values[0] + subsub['ARshort'].values[1] > ARcutoff) and use_short_circular:
                    passed = True
                    short = True

                if passed:
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

                                    if short:
                                        circular_labels.extend(('Short_Circular', 'Short_Circular'))
                                    else:
                                        circular_labels.extend(('Circular', 'Circular'))

                                elif not firstHalfA and not firstHalfB:
                                    real_circs.append(subsub)

                                    if short:
                                        circular_labels.extend(('Short_Circular', 'Short_Circular'))
                                    else:
                                        circular_labels.extend(('Circular', 'Circular'))
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

                                    if short:
                                        circular_labels.extend(('Short_Circular', 'Short_Circular'))
                                    else:
                                        circular_labels.extend(('Circular', 'Circular'))

                                elif firstHalfA and not firstHalfB:
                                    real_circs.append(subsub)

                                    if short:
                                        circular_labels.extend(('Short_Circular', 'Short_Circular'))
                                    else:
                                        circular_labels.extend(('Circular', 'Circular'))

                            #investigate other types:
                                #elif not firstHalfA and not firstHalfB:
                                    #print('+ dir singular: not firstHalfA and not firstHalfB')
                                    #print(subsub.to_string())
                                #elif firstHalfA and firstHalfB:
                                    #print('+ dir singular: firstHalfA and firstHalfB')
                                    #print(subsub.to_string())


                    elif len(subsub[5].unique()) == 2: # This means + to -
                        if keep_inversions: #this is a complex *putatively* circular alignment containing an inversion - note that it's lower confidence
                            logging.info("Including putative inversions in circular alignments.")
                            if overlap1 < 10 and overlap2 < 10:
                                  if not firstHalfA and firstHalfB:
                                      #make sure alignment starts within 100 bp of end of strand
                                      if min(astart, aend) < 101 or min(bstart, bend) < 101:
                                          if min(astart2, aend2) < 101 or min(bstart2, bend2) < 101:
                                                real_circs.append(subsub)


                                                if short:
                                                    circular_labels.extend(('Short_Circular', 'Short_Circular'))
                                                else:
                                                    circular_labels.extend(('Circular', 'Circular'))

                                  elif firstHalfA and not firstHalfB:
                                      if min(astart, aend) < 101 or min(bstart, bend) < 101:
                                          if min(astart2, aend2) < 101 or min(bstart2, bend2) < 101:
                                            real_circs.append(subsub)

                                            if short:
                                                circular_labels.extend(('Short_Circular', 'Short_Circular'))
                                            else:
                                                circular_labels.extend(('Circular', 'Circular'))


            index += 2


        if(len(real_circs)) > 1:
            result = pd.concat(real_circs)
    
            result.reset_index(drop=True, inplace=True)
            del result['HashCount']
            del result['Sample']
            del result['Hash']
            del result['ARshort']

            result['Circularity'] = circular_labels
            logging.debug('Writing circular alignments to file.')
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
        logging.debug('Processing alignment data.')
        logging.info("Alignment data is split into {} chunks of {} each. Analysing with {} threads:"
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

        if len(os.listdir(hash_dir)) == 0:
            full_hash_file = os.path.join(hash_dir, DefaultValues.HASH_FILE_NAME)
            full_passed_file = os.path.join(align_dir, DefaultValues.SECOND_PASS_NAME)
            hash_success = False
            logging.info('No hash files generated.')

            return full_hash_file, full_passed_file, hash_success


        full_hash_file = os.path.join(hash_dir, DefaultValues.HASH_FILE_NAME)
        full_passed_file = os.path.join(align_dir, DefaultValues.SECOND_PASS_NAME)

        logging.debug('Merging hash files and alignments kept after first pass.')




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

            hash_success = True

        except Exception as e:
            logging.error('An error occured while concatenating hashes and extracted alignments: {}'.format(e))
            sys.exit(1)

        return full_hash_file, full_passed_file, hash_success



    def __set_up_paf_thread(self, queue_in, queue_out, align_dir, hash_dir, chunksize, LRcutoff, ARcutoff, ANIcutoff):
    
        def filter_decider(AR, lena):
            if lena > 30000 and AR > 0.15:
                return 'YES'
            elif lena > 50000 and AR > 0.1:
                return 'YES'
            elif lena > 100000 and AR > 0.05:
                return 'YES'
            elif lena < 10000 and AR > 0.2:
                return 'YES'
            else:
                return 'NO'

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
                    # changed to account for diversity of lengths

#                    chunk['FilterAR'] = chunk.apply(lambda row: filter_decider(row['AR'], (row[2] + row[7]) / 2), axis=1)
#                    chunk2 = chunk[chunk['FilterAR'] == 'YES']                    
#                    del chunk2['FilterAR']
#                    del chunk['FilterAR']


                    chunk2 = chunk[chunk['AR'] > 0.1]
                    chunk2.iloc[:, -1].to_csv(out_hash, mode='a', sep='\t', header=False)
                    del chunk['Hash']
                    del chunk['Sample']
                    del chunk2

                    chunk = chunk[chunk['AR'] > ARcutoff]
                    chunk = chunk[chunk['LR'] > LRcutoff]
                    chunk = chunk[chunk['ANI'] > ANIcutoff]


                    chunk[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'AR', 'LR', 'ANI']].to_csv(out_align, mode='a', sep='\t', header=False)

            queue_out.put(paf)


    def __set_up_processing_sample_with_multiple_aligns(self, queue_in, queue_out, align_list, param_list, LRcutoff, use_short_circular):
        ''' Begin processing sample '''

        while True:
            idx = queue_in.get(block=True, timeout=None)
            if idx == None:
                break

            link = param_list[idx]

            with open(link, 'rb') as f:
                paf, map_dict = pickle.load(f)

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
            # make sure again:
            extracted_aligns['Sample'] = extracted_aligns.apply(
                lambda row: self.samplify(row[1], row[6]), axis=1)
            extracted_aligns = extracted_aligns[extracted_aligns['Sample'] == 'DifferentSample']

            extracted_aligns = extracted_aligns.astype(
                {2: 'int64', 3: 'int64', 4: 'int64', 7: 'int64', 8: 'int64', 9: 'int64', 10: 'int64', 11: 'int64'})


            extracted_aligns['AR'] = extracted_aligns.apply(
                lambda row: self.calcAR(row[2], row[7], row[11], row[10]), axis=1)
            extracted_aligns['ANI'] = extracted_aligns.apply(
                lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)

            extracted_aligns['ARshort'] = extracted_aligns.apply(
                lambda row: self.calcARshort(row[2], row[7], row[11], row[10]), axis=1)


            # aligned ratio may be short, but length ratio should be at least one third when using ARshort
            # or default if not using rescue_short_circular

            #if use_short_circular:
            #    extracted_aligns = extracted_aligns[extracted_aligns['LR'] >= 0.34]
            #else:
            extracted_aligns = extracted_aligns[extracted_aligns['LR'] >= LRcutoff]

            extracted_aligns = extracted_aligns[extracted_aligns['ANI'] >= DefaultValues.FIRST_PASS_AVA_ANI_CUTOFF]

            #check if we have any left
            
            if len(extracted_aligns) > 1:
            
                extracted_aligns['Hash'] = extracted_aligns.apply(lambda row: self.hashsum(row[1], row[6]),
                                                                                  axis=1)

                hashdict = dict(extracted_aligns['Hash'].value_counts())
                extracted_aligns['HashCount'] = extracted_aligns['Hash'].apply(lambda x: hashdict[x])

                extracted_aligns = extracted_aligns[extracted_aligns['HashCount'] > 1]

                if len(extracted_aligns) > 1:
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


class AssemblyMapper(Mapper):
    def map_clusters_against_assembly(self, cluster_folder, assembly_file, outdir, collapsedir, min_contig_len):


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


        #mapping_output = '{}/{}'.format(outdir, DefaultValues.ASSEMBLY_DEREP_FILE)

        try:

            #old cmd - check if dual is needed
#            cmd = "minimap2 -x ava-ont -t {} {} {} --dual=yes -v1 > {}" \
#                .format(self.nthreads, combined_contig_file, assembly_file, mapping_output)
            mapping_out = os.path.join(collapsedir, DefaultValues.COLLAPSE_MAPPING_NAME)

            cmd = "minimap2 -x ava-ont -m {} -K 5G -I 20G -2 -t {} {} {} --dual=yes -v1 | split --additional-suffix=.paf -d --line-bytes={} - {}" \
                .format(min_contig_len * 0.7, self.nthreads, combined_contig_file, assembly_file, DefaultValues.PAF_CHUNK_SIZE,
                        mapping_out)

            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished Running Minimap2')
        except Exception as e:
            logging.error('An error occured while running Minimap2: {}'.format(e))
            sys.exit(1)

       #return mapping_output
        return [os.path.join(outdir, collapsedir, mapping_file) for mapping_file in listdir(collapsedir) if
            mapping_file.endswith('.paf')]




    def remove_mapped_clusters(self, threads, paf_list, cluster_folder,
                                   LRcutoff, ARcutoff, ANIcutoff, outdir,
                                   cluster_info, cluster_contigs_info):

        #currently holding everything in memory, we should chunk and only keep good matches

        assembly_file_len = 0
        concat_df = []

        # try:
        #     with open(assembly_mapped_contigs) as f:
        #         assembly_file_len = sum(1 for line in f)
        #
        # except:
        #     logging.info('Not collapsing clusters against assembly as no valid mapping was generated.')
        #     return cluster_information, cluster_contigs_info

        if len(paf_list) < 1:
            logging.info('Not collapsing clusters against assembly as no valid mapping was generated.')
            return cluster_information, cluster_contigs_info

        else:

            # read in pafs and process them with multiple threads

            threads_per_bin = max(1, int(threads / len(paf_list)))
            logging.debug('Processing alignment data.')
            logging.info("Alignment data of contig vs assembly is split into {} chunks. Analysing with {} threads:"
                         .format(len(paf_list), threads))

            concat_df = mp.Manager().list()

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
                        mp.Process(target=self.__set_up_collapse_thread, args=(
                        workerQueue, writerQueue, concat_df, DefaultValues.FILTER_CHUNKSIZE, outdir, LRcutoff, ARcutoff, ANIcutoff)))
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


            # From here - as before

            if len(concat_df) > 0:
                paf = pd.concat(concat_df)

                cluster_info['ID'] = cluster_info['ID'].astype(int)

                paf['ID'] = paf[6].apply(lambda x: x.split('_ID_')[-1].split('_')[0])
                cluster_contig_connection = dict(zip(paf['ID'], paf[1]))


                clusters_to_remove = paf[6].unique()


                cluster_names_to_remove = [x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0] for x in clusters_to_remove]
                ids_to_remove = [int(x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0].split('_ID_')[-1].split('_')[0]) for x in cluster_names_to_remove]

                before_removal = len(cluster_info)
                removed_cluster_info = cluster_info[cluster_info['ID'].isin(ids_to_remove)]
                cluster_info = cluster_info[~cluster_info['ID'].isin(ids_to_remove)]
                after_removal = len(cluster_info)


                logging.info(f'Removed {before_removal - after_removal} clusters that were part of longer assembled contigs.')


                removed_cluster_info['Mapped_to_contig'] = removed_cluster_info['ID'].apply(lambda x: cluster_contig_connection[str(x)])

                for id in cluster_names_to_remove:
                    file = os.path.join(cluster_folder, id + '.fna')
                    os.remove(file)

                removed_cluster_info.to_csv(os.path.join(outdir, 'removed_integrated_clusters.tsv'), sep='\t', index=False)

                return cluster_info, cluster_contigs_info

            else:
                logging.info('No clusters were found to be integrated into a contig.')
                return cluster_info, cluster_contigs_info


    def __set_up_collapse_thread(self, queue_in, queue_out, align_list, chunksize, outdir, LRcutoff, ARcutoff, ANIcutoff):

        while True:
            paf = queue_in.get(block=True, timeout=None)
            if paf == None:
                break

            #paf_path = os.path.join(outdir, DefaultValues.COLLAPSE_DIR, paf)

            concat_df = []

            for chunk in pd.read_csv(paf,
                                     names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
                                     sep='\t', chunksize=chunksize):
                # remove self-alignment
                chunk = chunk[chunk[1] != chunk[6]]

                chunk['ANI'] = chunk.apply(lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)
                chunk['ARshortContig'] = (chunk[9] - chunk[8]) / chunk[7]

                chunk['AssembledContigLonger'] = chunk.apply(lambda row: 'YES' if \
                    float(row[2]) > float(
                        row[7] + row[7] * DefaultValues.ASSEMBLY_DEREP_CONTIG_PERCENT_LONGER) else 'NO', axis=1)

                chunk = chunk[chunk['AssembledContigLonger'] == 'YES']
                chunk = chunk[chunk['ANI'] > ANIcutoff]
                chunk = chunk[chunk['ARshortContig'] > ARcutoff]

                if len(chunk) > 0:
                    concat_df.append(chunk)

            if len(concat_df) > 0:
                concat_df = pd.concat(concat_df)
                align_list.append(concat_df)

            queue_out.put(paf)


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

        return mapping_output, combined_contig_file

    def collapse_related_clusters(self, mapping_file, cluster_folder,
                                  LRcutoff, ARcutoff, ANIcutoff, outdir,
                                  cluster_information, cluster_contigs_info, leftover_contigs_info,
                                  keep_related, filteredout):
        try:
            paf = pd.read_csv(mapping_file,
                                 names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], sep='\t')

            #remove self-alignment
            paf = paf[paf[1] != paf[6]]




            paf['AR'] = paf.apply(lambda row: self.calcAR(row[2], row[7], row[11], row[10]), axis=1)
            paf['LR'] = paf.apply(lambda row: self.calcLR(row[2], row[7], row[11], row[10]), axis=1)
            paf['ANI'] = paf.apply(lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)
            paf['ARshort'] = paf.apply(lambda row: self.calcARshort(row[2], row[7], row[11], row[10]), axis=1)
        except Exception as e:
            print(e)
            logging.info('Not collapsing clusters as no valid mapping was generated.')
            keep_related = True
            return cluster_information, cluster_contigs_info, leftover_contigs_info, None

        collapsee = []

        if len(paf) == 0:
            return cluster_information, cluster_contigs_info, leftover_contigs_info, None, None

        if not keep_related:
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

        # After collapsing related clusters, remove those that were previously filtered out

        #to_remove = cluster_information[cluster_information['ID'].isin(filteredout)]
        all_cluster_files = [f for f in listdir(cluster_folder) if f.endswith('.fna')]
        all_cluster_files = [os.path.join(cluster_folder, x) for x in all_cluster_files]

        tmp = {}
        for clust in all_cluster_files:
            id = int(clust.split('_ID_')[-1].split('_')[0])
            tmp[id] = clust

        to_remove = {k: tmp[k] for k in filteredout if k in tmp.keys()}

        for clust in to_remove.values():
            try:
                os.remove(clust)
            except FileNotFoundError:
                #may have been removed
                continue
        logging.info('Removed an additional {} clusters based on tightness and mean cutoffs.'.format(len(to_remove)))
        
        

        #update info
        cluster_information = cluster_information[~cluster_information['ID'].isin(filteredout)]
        new_remove_contigs = cluster_contigs_info[cluster_contigs_info['Cluster_ID'].isin(filteredout)].copy()
        cluster_contigs_info = cluster_contigs_info[~cluster_contigs_info['Cluster_ID'].isin(filteredout)]
        leftover_contigs_info = pd.concat([leftover_contigs_info, new_remove_contigs])

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

        #all_cluster_files = [f for f in listdir(cluster_folder) if f.endswith('.fna')]
        #all_cluster_files = [os.path.join(cluster_folder, x) for x in all_cluster_files]

        pd.DataFrame({'Contigs': toanalyse}).to_csv(FASTANI_input, sep='\t', index=False, header=False)


        #TODO: replace fragLen parameter with -m from main process

        try:
            cmd = "fastANI -t {} --ql {} --rl {} -o {} --fragLen 2500 --minFraction {} > /dev/null 2> /dev/null" \
                .format(self.nthreads, FASTANI_input, FASTANI_input, FASTANI_output, DefaultValues.SECOND_PASS_CONTIGS_AR_SHORT_CUTOFF)

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
            return cluster_information, cluster_contigs_info, leftover_contigs_info, None

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

        #Assume ANI is meaningless if size difference is more than 50% between two sequences

        fastANI_results['Size1'] = fastANI_results['Query'].apply(lambda x: x.split('_')[-1]).astype(int)
        fastANI_results['Size2'] = fastANI_results['Ref'].apply(lambda x: x.split('_')[-1]).astype(int)
        fastANI_results['Keep'] = fastANI_results.apply(lambda row: 'Yes' if max(row['Size1'], row['Size2']) <= min(row['Size1'], row['Size2']) * 1.5 else 'No', axis=1)

        fastANI_results = fastANI_results[fastANI_results['Keep'] == 'Yes']

        del fastANI_results['Size1']
        del fastANI_results['Size2']
        del fastANI_results['Keep']


        ''' Potential module to dereplicate circular clusters of similar size. Unclear if fastANI is reliable enough - not implemented for now 
        
        START UNIMPLEMENTED MODULE

        if not keep_related:
            # Update collapsed clusters from fastANI - this will catch circular alignments not parsed through yet 
            circular_aligns = fastANI_results[fastANI_results['ANI'] > DefaultValues.CIRCULAR_ANI_INCLUDING_CUTOFF].copy()
            circular_aligns['Circ1'] = circular_aligns['Query'].apply(lambda x: x.split('_')[0])
            circular_aligns['Circ2'] = circular_aligns['Ref'].apply(lambda x: x.split('_')[0])

            circular_aligns = circular_aligns[(circular_aligns['Circ1'] == 'Circular') | (circular_aligns['Circ2'] == 'Circular')]

            if len(circular_aligns) > 0:

                collapsee = []
                collapsed_into = []
                direction = []

                # we've generated some valid, circular alignments above threshold
                # we can iterate through and update everything
                def identify_shorter_longer(x, y, structure_x, structure_y):
                    lenx = int(x.split('_Size_')[-1].split('.')[0])
                    leny = int(y.split('_Size_')[-1].split('.')[0])

                    if lenx > leny:
                        # if we are trying to collapse a slightly shorter (up to 15%) circular into a slightly longer imperfect,
                        # do the reverse

                        if structure_x == 'Circular' and structure_y == 'Imperfect' and leny/lenx > 0.85:
                            return x, y
                        else:
                            return y, x
                    else:
                        if structure_y == 'Circular' and structure_x == 'Imperfect' and lenx/leny > 0.85:
                            return y, x
                        else:
                            return x, y

                circular_aligns[['Shorter', 'Longer']] = circular_aligns.apply(lambda row: identify_shorter_longer(row['Query'], row['Ref'], row['Circ2'], row['Circ1']),
                                                                               axis=1, result_type="expand")
                circular_aligns.reset_index(drop=True, inplace=True)
                for index, row in circular_aligns.iterrows():
                    shorter = circular_aligns.at[index, 'Shorter']
                    longer = circular_aligns.at[index, 'Longer']


                    logging.debug('Collapsing Circular clusters based on ANI: {} --> {}'.format(shorter, longer))

                    try:
                        os.remove(os.path.join(cluster_folder, shorter + '.fna'))
                    except FileNotFoundError:
                        continue #we may have already removed it before this stage.
                        #logging.error('File not found during circular ANI collapse: {}. Unable to continue.'.format(
                        #    os.path.join(cluster_folder, shorter + '.fna')))
                        #sys.exit(1)
                    collapsee.append(shorter)
                    collapsed_into.append(longer)
                    direction.append('-->')

                collapsed_log = pd.DataFrame({'Collapsee': collapsee, 'Direction': direction, 'Collapsed_Into': collapsed_into})
                collapsed_log_file = os.path.join(outdir, DefaultValues.COLLAPSED_LOG_FILE)
                # append to existing collapsed log file
                collapsed_log.to_csv(collapsed_log_file, sep='\t', index=False, header=False, mode='a' )

                logging.info('Found and collapsed {} related circular clusters based on ANI.'.format(len(set(collapsee))))

                # update cluster information
                remove_ids = [int(x.split('ID_')[-1].split('_')[0]) for x in collapsee]
                new_ids = [int(x.split('ID_')[-1].split('_')[0]) for x in collapsed_into]

                re_id = dict(zip(remove_ids, new_ids))

                cluster_contigs_info['Cluster_ID'] = cluster_contigs_info['Cluster_ID']. \
                    apply(lambda x: 'Collapsed into {}'.format(re_id[int(x)]) if (isinstance(x, int) and int(x) in re_id.keys()) else x)
                cluster_information = cluster_information[~cluster_information['ID'].isin(remove_ids)]
            del circular_aligns

        END UNIMPLEMENTED MODULE
        
        '''

        #final_fastANI_out = os.path.join(outdir, 'clusters_relatedness.tsv')

        #logging.info('Information for related clusters written to: {}.'.format(final_fastANI_out))
        #fastANI_results[['Query', 'Ref', 'ANI']].to_csv(final_fastANI_out, sep='\t', index=False)

        # determine groups of ANI-related clusters and write to results
        G = nx.from_pandas_edgelist(fastANI_results, source='Query', target='Ref', edge_attr='ANI')

        # Build adjacency matrix
        # ANI_adjacency_df = pd.DataFrame(
        #     nx.adjacency_matrix(G, weight='ANI').todense(),
        #     index=G.nodes,
        #     columns=G.nodes,
        #     dtype=float  # Compatible dtype with NaN is needed
        # )
        # # Overwrite the values on the diagonal
        #np.fill_diagonal(ANI_adjacency_df.values, np.NaN)

        ANI_defined_groups = (list(G.subgraph(c) for c in nx.connected_components(G)))

        group_dict = {}
        average_ANI_dict = {}

        for idx, group in enumerate(ANI_defined_groups):
            total_ANI = 0
            num_edges = 0
            for entry in list(group):
                id = int(entry.split('ID_')[-1].split('_')[0])
                group_dict[id] = idx
            for (node1, node2, data) in group.edges(data=True):
                total_ANI += data['ANI']
                num_edges += 1
            average_ANI_dict[idx] = total_ANI / num_edges if num_edges > 0 else 0

        cluster_information['Group'] = cluster_information['ID'].apply(lambda x: 'Group_{}'.format(group_dict[x]) if x in group_dict.keys() else 'None')
        cluster_information['Average_ANI'] = cluster_information['Group'].apply(
            lambda x: average_ANI_dict[int(x.split('_')[-1])] if x != 'None' else None)

        #todo: generate output (ANI and dendrogram) for each group


        return cluster_information, cluster_contigs_info, leftover_contigs_info, fastANI_results[['Query', 'Ref', 'ANI']]

class MergeMapping(Mapper):




    def mapMerge(self, file1, file2, map_outdir, threads):

#        mapping_input = combined_assemblies
#        mapping_out = os.path.join(map_outdir, DefaultValues.FIRST_PASS_MAPPING_NAME)
        # minimap2 options:
        try:

            #old cmd - check if dual is needed
#            cmd = "minimap2 -x ava-ont -t {} {} {} --dual=yes -v1 > {}" \
#                .format(self.nthreads, combined_contig_file, assembly_file, mapping_output)

            cmd = "minimap2 -x ava-ont {} {} --dual=yes -t {} -v1 > {}" \
                .format(file1, file2, threads, map_outdir)

            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished Running Minimap2')
        except Exception as e:
            logging.error('An error occured while running Minimap2: {}'.format(e))
            sys.exit(1)

        return map_outdir

    def process_map(self, mapfile, outdir):



        paf2 = pd.read_csv(mapfile, sep='\t',
                           names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])

        paf2['AR'] = paf2.apply(lambda row: self.calcAR(row[2], row[7], row[11], row[10]), axis=1)
        paf2['LR'] = paf2.apply(lambda row: self.calcLR(row[2], row[7], row[11], row[10]), axis=1)
        paf2['ANI'] = paf2.apply(lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)
        paf2['ARshort'] = paf2.apply(lambda row: self.calcARshort(row[2], row[7], row[11], row[10]), axis=1)
        paf2['ClustAL'] = (paf2[4] - paf2[3]) / paf2[2]  # paf.apply(lambda row: row[2] / (row[4] - row[3]), axis=1)
        paf2['ClustAL2'] = (paf2[9] - paf2[8]) / paf2[7]  # paf.apply(lambda row: row[2] / (row[4] - row[3]), axis=1)
        paf2 = paf2[paf2[1] != paf2[6]]
        paf2 = paf2[paf2['ANI'] > 0.899]


        def namer(c1, c2):
            return ''.join(sorted([c1, c2]))

        def overlap_percent(start1, end1, start2, end2):
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
                return (overlap / shortest) * 100

        paf2['namer'] = paf2.apply(lambda row: namer(row[1], row[6]), axis=1)
        paf2 = paf2[paf2['LR'] > 0.899]

        linear_subset = paf2[paf2['AR'] > 0.899].copy()

        paf2 = paf2[paf2['ClustAL2'] > 0.15]
        paf2 = paf2[paf2['ANI'] > 0.899]
        paf2 = paf2[paf2['LR'] > 0.899]
        countdict = dict(paf2['namer'].value_counts())
        paf2['namercount'] = paf2['namer'].apply(lambda x: countdict[x])
        multiples = paf2[paf2['namercount'] > 1]
        multiples = multiples.sort_values(by='namer')
        real_circs = []

        AR_THRESHOLD = 0.899

        for group in multiples['namer'].unique():
            subsub = multiples[multiples['namer'] == group]
            if subsub['namercount'].values[0] == 2:

                if (subsub['AR'].values[0] + subsub['AR'].values[1]) > AR_THRESHOLD:
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

                    overlap1 = overlap_percent(astart, aend, astart2, aend2)
                    overlap2 = overlap_percent(bstart, bend, bstart2, bend2)

                    if len(subsub[5].unique()) == 1:  # This means + to + or - to -
                        # no more than 10% overlap:
                        # TODO: check this aligns with AR to match cutoffs
                        if overlap1 < 10 and overlap2 < 10:

                            if subsub[5].values[0] == '-':
                                if firstHalfA and firstHalfB:
                                    real_circs.append(subsub)

                                elif not firstHalfA and not firstHalfB:
                                    real_circs.append(subsub)



                            elif subsub[5].values[0] == '+':
                                if not firstHalfA and firstHalfB:
                                    real_circs.append(subsub)

                                elif firstHalfA and not firstHalfB:
                                    real_circs.append(subsub)


                    elif len(subsub[5].unique()) == 2:  # This means + to -
                        if overlap1 < 10 and overlap2 < 10:
                            if not firstHalfA and firstHalfB:
                                real_circs.append(subsub)

                            elif firstHalfA and not firstHalfB:
                                real_circs.append(subsub)

        real_circs = pd.concat(real_circs)
        real_circs = real_circs.drop_duplicates(subset='namer')

        final = pd.concat([real_circs, linear_subset])

        del final['namer']
        del final['namercount']

        final = final.sort_values(by='AR', ascending=False)
        final = final.drop_duplicates(subset=6, keep='first')

        return final
