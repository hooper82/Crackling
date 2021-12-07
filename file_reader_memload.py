import argparse
import datetime    # Conflicts with imports in Helpers.py
import os
import logging
import re
import json

from collections import namedtuple


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S',
)


SequenceGuide = namedtuple("SequenceGuide", ['guide', 'sequence_header', 'start', 'end', 'strand'])


COMPLIMENTS = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
def reverse_complement(sequence):
    """
        Returns the reverse-complement of a given sequence
    """
    return sequence.translate(COMPLIMENTS)[::-1]


def process_sequence(sequence_header, sequence):
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
            yield SequenceGuide(target23, sequence_header,  m.start(),  m.start() + 23, strand)


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
    sequences = dict()

    start_time = datetime.datetime.now()

    logging.info(f'Identifying possible target sites in: {target_file}')

    with open(target_file, 'r') as file:
        for line in file:
            line = line.strip()

            if line[0] == '>':    # This is the header line for a new sequence. Create a new key in the sequences dict.
                header = line[1:]
                sequences[header] = []
            else:                 # Add to existing sequence.
                sequences[header].append(line)

    # Join the list of strings into a big string. This is faster than concating per line.
    for sequence_header in sequences.keys():
        sequences[sequence_header] = ''.join(sequences[sequence_header])

    # Find candidates!
    for sequence_header, sequence in sequences.items():
        for guide in process_sequence(sequence_header, sequence):

            # Check if guide has been seen before
            if guide.guide not in candidateGuides:
                # Record guide
                candidateGuides.add(guide.guide)
            else:
                # Record duplicate guide
                duplicateGuides.add(guide.guide)


    end_time = datetime.datetime.now()
    run_seconds = (end_time - start_time).total_seconds()

    candidateGuide_count = len(candidateGuides)
    duplicateGuides_count = len(duplicateGuides)
    duplicateGuides_pct = (duplicateGuides_count / candidateGuide_count) * 100

    # In the original implimentation in Crackling, there is an off-by-one error in `recordedSequences` set.
    # I belive it happened in the logic to try and read Exon Sequence files with corrupted Sequence Headers.
    # It doesn't matter. `recordedSequences` isn't used for anything in Crackling (infact, its deleted
    # in Line 251 in Crackling.py, shortly after processing the file).
    # 
    # To mirror the original implimentaiton, we're re-introducing this off-by-one error.
    recordedSequences = list(sequences.keys())
    recordedSequences.pop()            # Remove last sequence header
    recordedSequences.insert(0, "")    # Insert empty string into start of list.
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