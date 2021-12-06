import argparse
# import datetime    # Conflicts with imports in Helpers.py
import os
import logging
import tempfile
import re
from Helpers import * 


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S',
)


if __name__ == '__main__':
    # load in config
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', required=True)
    args = parser.parse_args()

    # Confirm the file exits
    assert os.path.isfile(args.file)
    target_file = args.file

    # Sets to keep track of Guides and sequences seen before
    candidateGuides = set()
    duplicateGuides = set()
    recordedSequences = set()

    def processSequence(sequence):
        # Patterns for guide matching
        pattern_forward = r'(?=([ATCG]{21}GG))'
        pattern_reverse = r'(?=(CC[ACGT]{21}))'

        # New sequence deteced, process sequence
        # once for forward, once for reverse
        for pattern, strand, seqModifier in [
            [pattern_forward, '+', lambda x : x], 
            [pattern_reverse, '-', lambda x : rc(x)]
        ]:
            p = re.compile(pattern)
            for m in p.finditer(sequence):
                target23 = seqModifier(seq[m.start() : m.start() + 23])
                yield [target23, seqHeader,  m.start(),  m.start() + 23, strand]


    logging.info(f'{target_file} exists.')

    start_time = datetime.now()

    logging.info(f'Identifying possible target sites in: {target_file}')

    # We first remove all the line breaks within a given sequence (FASTA format)
    with open(target_file, 'r') as inFile, tempfile.NamedTemporaryFile(mode='w',delete=False) as parsedFile:
        for line in inFile:
            line = line.strip()
            if line[0] == '>':
                # this is the header line for a new sequence, so we break the previous line and write the header as a new line
                parsedFile.write('\n'+line+'\n')
            else:
                # this is (part of) the sequence; we write it without line break
                parsedFile.write(line.strip())

    with open(parsedFile.name, 'r') as inFile:
        seqHeader = ''
        seq = ''
        for line in inFile:
            # Remove garbage from line
            line = line.strip()
            # Some lines (e.g., first line in file) can be just a line break, move to next line
            if line=='':
                continue
            # Header line, start of a new sequence
            elif line[0]=='>':
                # If we haven't seen the sequence OR we have found a sequence without header
                if (seqHeader not in recordedSequences) or (seqHeader=='' and seq!=''): 
                    # Record header
                    recordedSequences.add(seqHeader)
                    # Process the sequence
                    for guide in processSequence(seq):
                        # Check if guide has been seen before
                        if guide[0] not in candidateGuides:
                            # Record guide
                            candidateGuides.add(guide[0])
                            # # Record candidate guide to temp file
                            # guideBatchinator.recordEntry(guide)
                        else:
                            # Record duplicate guide
                            duplicateGuides.add(guide[0])
                # Update sequence and sequence header 
                seqHeader = line[1:]
                seq = ''
            # Sequence line, section of existing sequence
            else:
                # Append section to total sequence
                seq += line.strip()

        # Process the last sequence
        for guide in processSequence(seq):
            # Check if guide has been seen before
            if guide[0] not in candidateGuides:
                # Record guide
                candidateGuides.add(guide[0])
                # # Record candidate guide to temp file
                # guideBatchinator.recordEntry(guide)
            else:
                # Record duplicate guide
                duplicateGuides.add(guide[0])

    logging.info(f'Identified {len(candidateGuides)} possible target sites.')
    logging.info(f'\t{len(duplicateGuides)} of {len(candidateGuides)} were seen more than once.')



    end_time = datetime.now()
    run_seconds = (end_time - start_time).total_seconds()
    logging.info(f'Total Time Taken : {run_seconds:.2f} seconds')

