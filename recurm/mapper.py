import sys
import logging
import subprocess
import os
import pandas as pd
import numpy as np

PYTHONHASHSEED=0

from recurm.defaultValues import DefaultValues

class Mapper():
    def __init__(self, threads):
        self.nthreads = threads
        self.__check_for_minimap2()

    def __check_for_minimap2(self):
        # Assume that a successful minimap2 -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['minimap2', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            logging.error("Make sure minimap2 is on your system path.")
            sys.exit(1)

    def samplify(self, a, b):
        if a.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0] == \
                b.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0]:
            return 'SameSample'
        else:
            return 'DifferentSample'

    def hashsum(self, a, b):
        return hash(a) + hash(b)

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


    def mapAVA(self, combined_assemblies, outdir):
        mapping_output = '{}/{}'.format(outdir, DefaultValues.FIRST_PASS_MAPPING_NAME)
        mapping_input = combined_assemblies
        try:
            cmd = "minimap2 -x ava-ont -t {} {} {} --dual=yes -v1 > {}"\
                .format(self.nthreads, mapping_input, mapping_input, mapping_output)

            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished Running Minimap2')
        except Exception as e:
            logging.error('An error occured while running Minimap2: {}'.format(e))
            sys.exit(1)

        return mapping_output

    def read_hashes_and_process_alignments(self, hash_file, alignment_file, LRcutoff, ARcutoff, ANIcutoff, outfile):
        #read in hashes - this assumes hash files can be read in at once
        #roughly 10 GB RAM per TB of alignment data
        #optimize later

        chunksize = DefaultValues.FILTER_CHUNKSIZE

        dup_hashes = pd.read_csv(hash_file, sep='\t', names=['OriginalIndex', 'Hash'], dtype=str)
        dup_hashes = dup_hashes.astype({'OriginalIndex': 'int64', 'Hash': 'string'})
        #sort hashes
        indexes = dup_hashes['OriginalIndex']
        length = len(indexes)
        indexes = np.array(sorted(indexes))
        map_dict = dict(zip(dup_hashes['OriginalIndex'], dup_hashes['Hash']))

        temp_index = []
        keep_list = []

        logging.info('Iterating through alignment a second time to extract circular and complex alignments')
        orig_aligns = pd.read_csv(alignment_file, chunksize=chunksize, names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], usecols=[0,1,2,3,4,5,6,7,8,9,10], sep='\t')

        cindex = 0
        #done_chunk = False

        #Having gotten the useful indexes, run through alignment file again and pop off rows of interest
        for chunk in orig_aligns:
            cindex += len(chunk)
            # peek at first index
            if indexes[0] > cindex:
                continue  # we won't care about this chunk
            else:
                # pop off enough indexes to populate current chunksize
                i = 0
                o = 0
                while o <= cindex:
                    if i < length:
                        o = indexes[i]
                        temp_index.append(o)
                        i += 1
                    else:
                        break
                indexes = indexes[i - 1:]
                length -= len(temp_index)

            sub = chunk[chunk.index.isin(temp_index)]
            if len(sub) > 1:
                # sub['Hash'] = sub.index.to_series().map(map_dict)
                keep_list.append(sub)

            temp_index = []

        extracted_aligns = pd.concat(keep_list)
        extracted_aligns.reset_index(drop=True, inplace=True)
        del dup_hashes

        extracted_aligns = extracted_aligns.astype(
            {2: 'int64', 3: 'int64', 4: 'int64', 7: 'int64', 8: 'int64', 9: 'int64', 10: 'int64', 11: 'int64'})

        extracted_aligns['Hash'] = extracted_aligns.apply(lambda row: self.hashsum(row[1], row[6]), axis=1)

        hashdict = dict(extracted_aligns['Hash'].value_counts())
        extracted_aligns['HashCount'] = extracted_aligns['Hash'].apply(lambda x: hashdict[x])

        #TODO: do not ignore complex alignments > 2. For now ignore them.
        extracted_aligns = extracted_aligns[extracted_aligns['HashCount'] == 2]

        extracted_aligns['AR'] = extracted_aligns.apply(lambda row: self.calcAR(row[2], row[7], row[11], row[10]), axis=1)
        extracted_aligns['LR'] = extracted_aligns.apply(lambda row: self.calcLR(row[2], row[7], row[11], row[10]), axis=1)
        extracted_aligns['ANI'] = extracted_aligns.apply(lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)
        extracted_aligns['Sample'] = extracted_aligns.apply(lambda row: self.samplify(row[1], row[6]), axis=1)

        extracted_aligns = extracted_aligns[extracted_aligns['Sample'] == 'DifferentSample']

        #aligned ratio may be short, but length ratio should be solid
        extracted_aligns = extracted_aligns[extracted_aligns['LR'] >= LRcutoff]
        extracted_aligns = extracted_aligns.sort_values(by='Hash')

        index = 0
        real_circs = []

        #Now go through the double-alignment pairs and make sure that they point in the right direction.

        while index < len(extracted_aligns):

            subsub = extracted_aligns.iloc[index:index + 2, :]

            if subsub['AR'].values[0] + subsub['AR'].values[1] > ARcutoff:
                # direction is conserved between alignments, ie. no palindromic alignments
                if len(subsub[5].unique()) == 1:

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

                    #no more than 10% overlap:
                    #TODO: check this aligns with AR to match cutoffs
                    if overlap1 < 10 and overlap2 < 10:

                        if subsub[5].values[0] == '-':
                            if firstHalfA and firstHalfB:
                                real_circs.append(subsub)

                        elif subsub[5].values[0] == '+':
                            if not firstHalfA and firstHalfB:
                                real_circs.append(subsub)

            index += 2
        result = pd.concat(real_circs)
        result = result[result['ANI'] > ANIcutoff]

        result.reset_index(drop=True, inplace=True)
        del result['HashCount']
        del result['Sample']
        del result['Hash']

        logging.info('Writing circular alignments to file.')
        result.to_csv(outfile, sep='\t', index=False, header=False)




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

    def subfilter_paf_mapping_results(self, paf_file, outdir, outfile, chunksize, LRcutoff, ARcutoff, ANIcutoff):
        df = pd.read_csv(paf_file, chunksize=chunksize, names=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], sep='\t')

        file_length = self.line_count(paf_file)



        index = 0

        for chunk in df:
            chunk['AR'] = chunk.apply(lambda row: self.calcAR(row[2], row[7], row[11], row[10]), axis=1)
            chunk['LR'] = chunk.apply(lambda row: self.calcLR(row[2], row[7], row[11], row[10]), axis=1)
            chunk['ANI'] = chunk.apply(lambda row: self.calcANI(row[2], row[7], row[11], row[10]), axis=1)
            chunk['Sample'] = chunk.apply(lambda row: self.samplify(row[1], row[6]), axis=1)

            chunk = chunk[chunk['Sample'] == 'DifferentSample']

            chunk['Hash'] = chunk.apply(lambda row: self.hashsum(row[1], row[6]), axis=1)

            #multiple alignments need at least ~20% length coverage on a side to be a decent match
            chunk2 = chunk[chunk['AR'] > 0.2]


            chunk2.iloc[:, -1].to_csv('{}/{}.hashes'.format(outdir, outfile), mode='a', sep='\t', header=False)
            del chunk['Hash']
            del chunk['Sample']
            del chunk2

            chunk = chunk[chunk['AR'] > ARcutoff]
            chunk = chunk[chunk['LR'] > LRcutoff]
            chunk = chunk[chunk['ANI'] > ANIcutoff]

            chunk[[1,2,3,4,5,6,7,8,9,10,11,'AR','LR','ANI']].to_csv('{}/{}.passed'.format(outdir, outfile), mode='a', sep='\t', header=False)

            index += chunksize

            if index > file_length:
                index = file_length

            if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                statusStr = '    Finished processing {}% of initial alignment file'.\
                    format(round((float(index) / file_length) * 100), 2)
                sys.stdout.write('\r{}'.format(statusStr))

                if index == file_length:
                    if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                        sys.stdout.write('\n')
                        sys.stdout.flush()
                else:
                    sys.stdout.flush()

        return '{}/{}.hashes'.format(outdir, outfile), '{}/{}.passed'.format(outdir, outfile)

class ReferenceMapper(Mapper):
    pass

class WithinClusterMapper(Mapper):
    pass