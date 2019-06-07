import argparse
import os
import Levenshtein as leve
from tqdm import tqdm
from collections import Counter
import json
import sys

sys.setrecursionlimit(90000)

DEFAULT_PARAM_FILE = 'param.json'
DEFAULT_LONG_FILE = 'long.fasta'
DEFAULT_SHORT_1_FILE = 'short_1.fasta'
DEFAULT_SHORT_2_FILE = 'short_2.fasta'
DEFAULT_FIXED_LONG_FILE = 'fixed_long.fasta'
DEFAULT_MATCHES_FILE = 'matches.json'

ARGS = None
PROGRESS = {
    "now_part": 1,
    "now_done": 0,
    "total_done": 0
}
PARAM = {}

MAXDIS = 30


def parse_args():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument('DATA_DIR', type=str, help='the dataset\'s directory')
    parser.add_argument('ID', type=int, help='ID')
    ARGS = parser.parse_args()
    ARGS.PARAM_FILE = os.path.join(ARGS.DATA_DIR, DEFAULT_PARAM_FILE)
    ARGS.LONG_FILE = os.path.join(ARGS.DATA_DIR, DEFAULT_LONG_FILE)
    ARGS.SHORT_1_FILE = os.path.join(ARGS.DATA_DIR, DEFAULT_SHORT_1_FILE)
    ARGS.SHORT_2_FILE = os.path.join(ARGS.DATA_DIR, DEFAULT_SHORT_2_FILE)
    ARGS.FIXED_LONG_FILE = os.path.join(ARGS.DATA_DIR, DEFAULT_FIXED_LONG_FILE)
    ARGS.MATCHES_FILE = os.path.join(ARGS.DATA_DIR, DEFAULT_MATCHES_FILE)


def prepare_fasta_data(filename):
    content = []
    print('Load data', filename)
    with open(filename, 'r') as f:
        lines = f.readlines()
        name = 'Unknown'
        for i, line in enumerate(tqdm(lines)):
            if (i & 1) == 0:
                name = line.strip('\n')
            else:
                content.append({'name': name, 's': line.strip('\n')})
    return content


def get_comp_rev_data(content):
    tran = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return [{
        'name': data['name']+'(comp_rev)',
        's': ''.join([tran[c] for c in data['s']][::-1])
    } for data in content]


def show_match(filename, fixed_long_data, match_info):
    new_long = fixed_long_data['s']
    if match_info[0]['pos'] >= 0:
        pre = 0
    else:
        pre = -match_info[0]['pos']
    with open(filename, 'w') as f:
        f.write(' '*pre+new_long+'\n')
        for d in match_info:
            f.write((' '*(pre+d['pos']))+d['s']+'\n')
        f.write(' '*pre+new_long+'\n')


if __name__ == "__main__":
    # load ARGS
    parse_args()
    # load data
    fixed_dataset = prepare_fasta_data(ARGS.FIXED_LONG_FILE)
    match_infoset = {}
    with open(ARGS.MATCHES_FILE, 'r') as f:
        match_infoset = json.loads('['+f.read().strip('\n')[:-1]+']')

    show_match('show_matches.fasta', fixed_dataset[ARGS.ID], match_infoset[ARGS.ID])
