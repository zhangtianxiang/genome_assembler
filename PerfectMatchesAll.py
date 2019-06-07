import argparse
import os
import Levenshtein as leve
from tqdm import tqdm
from collections import Counter
import json
import sys
import math

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

MAXDIS = 21
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


def match_short_data(long_data, short_data):
    long = long_data['s']
    short = short_data['s']
    mindis = -1
    minpos = 0
    for i in range(len(long)-len(short)+1):
        dis = leve.hamming(short, long[i:i+len(short)])
        if mindis == -1 or dis < mindis:
            mindis = dis
            minpos = i
    # # tail
    # for i in range(len(long)-len(short)+1, len(long)-(len(short)//2)+1):
    #     common = len(long)-i
    #     # ceil(dis*LEN_SHORT/LEN_COMMON)
    #     dis = (leve.hamming(short[:common],
    #                         long[i:i+common])*len(short)+common-1)//common
    #     if mindis == -1 or dis < mindis:
    #         mindis = dis
    #         minpos = i
    # # head
    # for i in range(-(len(short)//2), 0):
    #     common = len(short)+i
    #     # ceil(dis*LEN_SHORT/LEN_COMMON)
    #     dis = (leve.hamming(short[-i:], long[:common])
    #            * len(short)+common-1)//common
    #     if mindis == -1 or dis < mindis:
    #         mindis = dis
    #         minpos = i
    match = {
        'name': short_data['name'],
        'pos': minpos,
        'dis': mindis,
        's': short_data['s']
    }
    return match


def fix_long_data_by_matches(long_data, match_info):
    match_info = sorted(match_info, key=lambda x: x['pos'])
    long = long_data['s']
    new_match_info = []
    now_info = 0
    now_que = []
    new_long = ''
    for now_pos in range(len(long)):
        # add new match_info to que
        while now_info < len(match_info) and match_info[now_info]['pos'] <= now_pos:
            now_que.append(match_info[now_info])
            now_info += 1
        # save survived match_info to new_match_info
        while len(now_que) > 0 and now_que[0]['pos']+len(now_que[0]['s']) <= now_pos:
            new_match_info.append(now_que[0])
            now_que.pop(0)
        # count now char
        vote = {
            'A': 0,
            'T': 0,
            'G': 0,
            'C': 0
        }
        vote[long[now_pos]] = 1
        vote_max = 1
        vote_char = long[now_pos]
        for match in now_que:
            vote[match['s'][now_pos-match['pos']]] += 1
        for k, v in vote.items():
            if v > vote_max:
                vote_max = v
                vote_char = k
        new_long += vote_char
        # # remove not correct match_info
        # for i in range(len(now_que)):
        #     while i < len(now_que) and now_que[i]['s'][now_pos-now_que[i]['pos']] != vote_char:
        #         now_que.pop(i)

    # save survival match_info to new_match_info
    while len(now_que) > 0:
        new_match_info.append(now_que[0])
        now_que.pop(0)

    fixed_long_data = {
        'name': long_data['name']+'(fixed)',
        's': new_long
    }
    return fixed_long_data, new_match_info


def match_short_dataset(long_data, short_dataset):
    match_info = []
    for short_data in short_dataset:
        match = match_short_data(long_data, short_data)
        if match['dis'] < MAXDIS:
            match_info.append(match)
    fixed_data, match_info = fix_long_data_by_matches(long_data, match_info)
    return fixed_data, match_info


def get_match_info_and_save(long_dataset, short_datasets):
    global PART_ID
    PART_ID += 1
    if PROGRESS['now_part'] < PART_ID-1:  # Unknow error
        raise 'Part Error'
    elif PROGRESS['now_part'] == PART_ID-1:  # New part
        PROGRESS['now_part'] = PART_ID
        PROGRESS['now_done'] = 0
    elif PROGRESS['now_part'] == PART_ID:  # Continue
        print('Continue: have done %d/%d' %
              (PROGRESS['now_done'], len(long_dataset)))
    else:  # Done
        print('Done')
        return
    offset = PROGRESS['now_done']
    fixed_dataset = []
    match_infoset = []
    mixed_short_dataset = []
    # mix all short data
    for short_dataset in short_datasets:
        mixed_short_dataset.extend(short_dataset)
    for i, long_data in enumerate(tqdm(long_dataset)):
        if i < offset:
            continue
        fixed_data, match_info = match_short_dataset(
            long_data, mixed_short_dataset)
        save_all_and_update_progress(fixed_data, match_info)
        fixed_dataset.append(fixed_data)
        match_infoset.append(match_info)
    print('Done')


def save_all_and_update_progress(fixed_data, match_info):
    with open(ARGS.FIXED_LONG_FILE, 'a') as f:
        f.write(fixed_data['name']+'\n')
        f.write(fixed_data['s']+'\n')
    with open(ARGS.MATCHES_FILE, 'a') as f:
        f.write(json.dumps(match_info)+',\n')
    PROGRESS['now_done'] += 1
    PROGRESS['total_done'] += 1
    with open('PROGRESS.json', 'w') as f:
        f.write(json.dumps(PROGRESS))


def fix_data():
    if PROGRESS['total_done'] > 0:
        lines = []
        with open(ARGS.FIXED_LONG_FILE, 'r') as f:
            for i, line in enumerate(f.readlines()):
                if i < PROGRESS['total_done']*2:
                    lines.append(line.strip('\n'))
        assert len(lines) == PROGRESS['total_done']*2
        with open(ARGS.FIXED_LONG_FILE, 'w') as f:
            for line in lines:
                f.write(line+'\n')
        match_infoset = []
        with open(ARGS.MATCHES_FILE, 'r') as f:
            s = f.read().strip('\n')[:-1]  # 刪除逗號
            s = '['+s+']'
            match_infoset = json.loads(s)[:PROGRESS['total_done']]
        assert len(match_infoset) == PROGRESS['total_done']
        with open(ARGS.MATCHES_FILE, 'w') as f:
            for match_info in match_infoset:
                f.write(json.dumps(match_info)+',\n')
    else:
        assert PROGRESS['total_done'] == 0
        with open(ARGS.FIXED_LONG_FILE, 'w') as f:
            pass
        with open(ARGS.MATCHES_FILE, 'w') as f:
            pass


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
    # update threhold
    MAXDIS = math.floor(PARAM['short_read_length'] *
                        (PARAM['short_read_error_rate']+PARAM['long_read_error_rate'])+5)
    # load PROGRESS
    with open('PROGRESS.json', 'r') as f:
        pre = json.loads(f.read())
        if isinstance(pre, dict):
            for k, v in PROGRESS.items():
                PROGRESS[k] = pre.get(k, v)
    # fix files
    fix_data()
    # begin fix
    print('Fixing long_read')
    get_match_info_and_save(
        fl, [fs1, fs2, fs1_rev, fs2_rev])
    print('Fixing long_read(rev_comp)')
    get_match_info_and_save(
        fl_rev, [fs1, fs2, fs1_rev, fs2_rev])
