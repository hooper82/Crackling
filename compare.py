import argparse
import os
import json
import logging


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S',
)


def read_file(filename):
    with open(filename, 'r') as file:
        data = json.load(file)
    return data


if __name__ == '__main__':
    # load in config
    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', required=True)
    parser.add_argument('-f2', required=True)
    args = parser.parse_args()

    assert os.path.isfile(args.f1)
    assert os.path.isfile(args.f2)

    f1_data = read_file(args.f1)
    logging.info(f'Loaded {args.f1} ({os.path.getsize(args.f1):,} Bytes)')

    f2_data = read_file(args.f2)
    logging.info(f'Loaded {args.f2} ({os.path.getsize(args.f2):,} Bytes)')

    for field in ['candidate_guides', 'duplicate_guides', 'recorded_sequences']:
        if f1_data[field] == f2_data[field]:
            logging.info(f'{field} are the same in both datasets.')
        else:
            logging.warning(f'{field} are different!!')
            logging.info(f"{args.f1} : {len(f1_data[field]):,} entries")
            logging.info(f"{args.f2} : {len(f2_data[field]):,} entries")

            logging.warning(f'The following items are in {args.f1} but not {args.f2}:')
            for item in list(set(f1_data[field]) - set(f2_data[field])):
                logging.info(f"    {item}")

            logging.warning(f'The following items are in {args.f2} but not {args.f1}:')
            for item in list(set(f2_data[field]) - set(f1_data[field])):
                logging.info(f"    {item}")