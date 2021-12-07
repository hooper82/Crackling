import argparse
import datetime
import os
import logging
import re
import json
import dbm


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S',
)


COMPLIMENTS = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
def reverse_complement(sequence):
    """
        Returns the reverse-complement of a given sequence
    """
    return sequence.translate(COMPLIMENTS)[::-1]


def process_sequence(sequence, sequence_header):
    # Patterns for guide matching
    pattern_forward = r'(?=([ATCG]{21}GG))'
    pattern_reverse = r'(?=(CC[ACGT]{21}))'

    # New sequence deteced, process sequence
    # once for forward, once for reverse
    for pattern, strand, seqModifier in [
        [pattern_forward, '+', lambda x : x], 
        [pattern_reverse, '-', lambda x : reverse_complement(x)]
    ]:
        p = re.compile(pattern)
        for m in p.finditer(sequence):
            target23 = seqModifier(sequence[m.start() : m.start() + 23])
            yield [target23, sequence_header,  m.start(),  m.start() + 23, strand]


if __name__ == '__main__':
    # load in config
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--output_file', required=True)
    args = parser.parse_args()

    # Confirm the file exits
    assert os.path.isfile(args.input_file)
    target_file = args.input_file
    logging.info(f'{target_file} exists.')

    # Sets to keep track of Guides and sequences seen before
    candidateGuides = set()
    duplicateGuides = set()
    recordedSequences = set()


    start_time = datetime.datetime.now()

    logging.info(f'Identifying possible target sites in: {target_file}')

    # Extract guides from sequence file
    with open(args.input_file, 'r') as inFile:
        seqHeader = ''
        seq = []
        # Flag to determine whether a new seqeunce is being processed 
        newSeq = True
        for line in inFile:
            # Remove garbage from line
            line = line.strip()
            # Some lines (e.g., first line in file) can be just a line break, move to next line
            if line=='':
                continue
            # Header line, start of a new sequence
            elif line[0] == '>':
                # Process previous sequence
                if (seqHeader not in recordedSequences) and (len(seq)!=0): 
                    # Record header (If there is one)
                    if seqHeader != '':
                        recordedSequences.add(seqHeader)
                    # Concat seqence fragments
                    fullSeq = ''.join(seq)
                    # Process the full sequence
                    for guide in process_sequence(fullSeq, seqHeader):
                        # Check if guide has been seen before 
                        if guide[0] not in candidateGuides:
                            # Record guide
                            candidateGuides.add(guide[0])
                        else:
                            # Record duplicate guide
                            duplicateGuides.add(guide[0])
                # Check if seqHeader has been seen before
                if line[1:] in recordedSequences:
                    # Old Header, update newSeq to reflect we are passing over an OLD sequence
                    newSeq = False
                    seqHeader = ''
                else:
                    # New Header, update newSeq to reflect we are passing over an NEW sequence
                    newSeq = True
                    seqHeader = line[1:]
                seq = []
            # Sequence line, section of sequence
            else:
                # Check to see if a new sequence is being recorded
                if newSeq:
                    seq.append(line)

        # Process the last sequence
        # Concat seqence fragments
        fullSeq = ''.join(seq)
        for guide in process_sequence(fullSeq, seqHeader):
            # Check if guide has been seen before
            if guide[0] not in candidateGuides:
                # Record guide
                candidateGuides.add(guide[0])
            else:
                # Record duplicate guide
                duplicateGuides.add(guide[0])


    end_time = datetime.datetime.now()
    run_seconds = (end_time - start_time).total_seconds()

    candidateGuide_count = len(candidateGuides)
    duplicateGuides_count = len(duplicateGuides)
    duplicateGuides_pct = (duplicateGuides_count / candidateGuide_count) * 100
    recordedSequences_count = len(recordedSequences)

    logging.info(f'Total Time Taken : {run_seconds:.2f} seconds')
    logging.info(f'Results:')
    logging.info(f'Candidate Guides   : {candidateGuide_count:,}')
    logging.info(f'Duplicate Guides   : {duplicateGuides_count:,} ({duplicateGuides_pct:.2f} %)')
    logging.info(f'Recorded Sequences : {recordedSequences_count:,}')

    results = {
        'candidate_guides'   : sorted(list(candidateGuides)),
        'duplicate_guides'   : sorted(list(duplicateGuides)),
        'recorded_sequences' : sorted(list(recordedSequences)),
    }


    with open(args.output_file, 'w') as file:
        json.dump(results, file)

    output_file_size = os.path.getsize(args.output_file)
    
    logging.info(f'\nOutput written to {args.output_file} ({output_file_size:,} bytes)')