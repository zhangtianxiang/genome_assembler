'''
思想：
已有修复后的长传fixed_long.fasta和匹配信息matches.json
每个长串匹配了众多小串
长串做为A类点，小串作为B类点
则A B点之间建边为距离

1. 找到最小的AB距离，记录A为C
2. 不断将C上距离最小的B找出，然后找到B所在的长串与C拼接
3. 最终输出C
'''
import argparse
import os
import Levenshtein as leve
from tqdm import tqdm
from collections import Counter
import json
import sys
import math
import heapq

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

MAXDIS = 0
PART_ID = 0


def parse_args():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument('DATA_DIR', type=str, help='the dataset\'s directory')
    parser.add_argument('-fl', '--FIXED_LONG_FILE', type=str,
                        default=DEFAULT_FIXED_LONG_FILE, help='fixed long data file\'s name')
    parser.add_argument('-m', '--MATCHES_FILE', type=str,
                        default=DEFAULT_MATCHES_FILE, help='match info file \'s name')
    parser.add_argument('-p', '--PARAM_FILE', type=str,
                        default=DEFAULT_PARAM_FILE, help='the param file\'s name')
    parser.add_argument('-l', '--LONG_FILE', type=str,
                        default=DEFAULT_LONG_FILE, help='the long data file\'s name')
    parser.add_argument('-s1', '--SHORT_1_FILE', type=str,
                        default=DEFAULT_SHORT_1_FILE, help='the short_1 file\'s name')
    parser.add_argument('-s2', '--SHORT_2_FILE', type=str,
                        default=DEFAULT_SHORT_2_FILE, help='the short_2 file\'s name')
    ARGS = parser.parse_args()
    ARGS.PARAM_FILE = os.path.join(ARGS.DATA_DIR, ARGS.PARAM_FILE)
    ARGS.LONG_FILE = os.path.join(ARGS.DATA_DIR, ARGS.LONG_FILE)
    ARGS.SHORT_1_FILE = os.path.join(ARGS.DATA_DIR, ARGS.SHORT_1_FILE)
    ARGS.SHORT_2_FILE = os.path.join(ARGS.DATA_DIR, ARGS.SHORT_2_FILE)
    ARGS.FIXED_LONG_FILE = os.path.join(ARGS.DATA_DIR, ARGS.FIXED_LONG_FILE)
    ARGS.MATCHES_FILE = os.path.join(ARGS.DATA_DIR, ARGS.MATCHES_FILE)


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


if __name__ == "__main__":
    '''
    long_dna = xxx
    long_data = {
        'name': XXX
        's': long_dna
    }
    long_dataset = [long_data]
    long_datasets = [long_dataset]

    match = {
        'name': XXX
        'pos': x
        'dis': x
                's': short_dna
    }
    match_info = [match]
    match_infoset = [match_info]
    match_infosets = [match_infoset]
    '''
    # load ARGS
    parse_args()
    # load data
    with open(ARGS.PARAM_FILE, 'r') as f:
        PARAM = json.loads(f.read())
    fs1 = prepare_fasta_data(ARGS.SHORT_1_FILE)
    fs2 = prepare_fasta_data(ARGS.SHORT_2_FILE)
    fl = prepare_fasta_data(ARGS.LONG_FILE)
    fs1_rev = get_comp_rev_data(fs1)
    fs2_rev = get_comp_rev_data(fs2)
    fl_rev = get_comp_rev_data(fl)

    fixed_dataset = prepare_fasta_data(ARGS.FIXED_LONG_FILE)
    match_infoset = {}
    with open(ARGS.MATCHES_FILE, 'r') as f:
        match_infoset = json.loads('['+f.read().strip('\n')[:-1]+']')

    # update threhold
    MAXDIS = math.floor(PARAM['short_read_length'] *
                        (PARAM['short_read_error_rate']+PARAM['long_read_error_rate'])+5)

    # map_point
    name_to_Apoint = {}
    Apoint_to_name = []
    Apoint_to_data = fixed_dataset

    name_to_Bpoint = {}
    Bpoint_to_name = []
    Bpoint_to_data = []

    record_checked_A = []
    record_checked_B = []
    A_out_edges = []
    B_out_edges = []

    def get_id(short_name, short):
        sid = name_to_Bpoint.get(short_name, -1)
        if sid == -1:
            sid = len(Bpoint_to_name)
            name_to_Bpoint[short_name] = sid
            Bpoint_to_name.append(short_name)
            Bpoint_to_data.append({'name': short_name, 's': short})
            B_out_edges.append([])
            record_checked_B.append(False)
        return sid

    def encode_edge(A, B, dis, pos):
        return (dis, pos, A, B)  # 以元组存储便于比较，dis为第一关键字，pos为第二关键字

    def decode_edge(edge):
        return edge[2], edge[3], edge[0], edge[1]

    def add_edge(A, B, dis, pos):
        edge = encode_edge(A, B, dis, pos)
        A_out_edges[A].append(edge)
        B_out_edges[B].append(edge)

    # generate A point
    print('Generate A point')
    for long_data in tqdm(fixed_dataset):
        name_to_Apoint[long_data['name']] = len(Apoint_to_name)
        Apoint_to_name.append(long_data['name'])
        A_out_edges.append([])
        record_checked_A.append(False)

    # generate B point and add edges
    print('Generate B point and build graph')
    for i, match_info in enumerate(tqdm(match_infoset)):
        short_name = match_info['name']
        short = match_info['s']
        sid = get_id(short_name, short)
        add_edge(i, sid, match_info['dis'], match_info['pos'])

    # generate final DNA!
    def generate(start_Apoint):
        g = Apoint_to_data[start_Apoint]['s']
        record_checked_A[start_Apoint] = True
        h = []
        for edge in A_out_edges:
            A, B, dis, pos = decode_edge(edge)
            if record_checked_B[dis] == False:
                h.append(edge)
        heapq.heapify(h)
        offset = 0  # position 's offset
        while True:
            if len(h) == 0:
                break
            # A => mindis B
            A, B, dis, pos = decode_edge(heapq.heappop(h))
            if record_checked_B[B]:
                continue
            # merge all A (which B belongs to) to g, and update a
            for edge in B_out_edges[B]:
                A2, _, _, pos2 = decode_edge(edge)
                # try merge A2

    ans = []
    print('Final assemble')
    for i, long_data in enumerate(tqdm(fixed_dataset)):
        if record_checked_A[i] == False:
            ans.append(generate(i))
