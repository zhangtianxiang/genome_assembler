'''
修复matches

多个短串匹配在这条长链上，对每个位置进行投票，得到1000个voter
所有voter，如果存在多种结果，且第二高的峰值多余最高峰值乘 (PARAM['short_read_error_rate'] * 2)
data1: 0.0 * 2 = 0
data2: 0.01 * 2 = 0.02 即允许此位置一共匹配100个，正常应该只出错2个
如果最高值的总数没有达到98%，那多出来的短串应该去匹配另外一个串，而非这个串
所以不断删去 次高值 代表的串，然后重新来一遍，直至最高值的比值能够达到98%

之后再重新求一遍voter，循环执行

这样修复每一条长链复杂度大概是匹配的短链个数*100*1000*长链个数，因为存在删除，期望应该会很低。

成功决绝了第二个点的问题达到了 9361.8
'''
import argparse
import os
import Levenshtein as leve
from tqdm import tqdm
from collections import Counter
import json
import sys
import math
import numpy as np

sys.setrecursionlimit(90000)

DEFAULT_PARAM_FILE = 'param.json'
DEFAULT_LONG_FILE = 'long.fasta'
DEFAULT_SHORT_1_FILE = 'short_1.fasta'
DEFAULT_SHORT_2_FILE = 'short_2.fasta'
DEFAULT_FIXED_LONG_FILE = 'fixed_long.fasta'
DEFAULT_MATCHES_FILE = 'matches.json'
DEFAULT_LONG_REPAIR_FILE = 'long_repair.fasta'
DEFAULT_MATCH_REPAIR_FILE = 'matches_repair.json'

ARGS = None
PARAM = {}

MINRATE = 0


def parse_args():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument('DATA_DIR', type=str, help='the dataset\'s directory')
    parser.add_argument('-fl', '--FIXED_LONG_FILE', type=str,
                        default=DEFAULT_FIXED_LONG_FILE, help='fixed long data file\'s name')
    parser.add_argument('-m', '--MATCHES_FILE', type=str,
                        default=DEFAULT_MATCHES_FILE, help='match info file \'s name')
    parser.add_argument('-mr', '--MATCH_REPAIR_FILE', type=str,
                        default=DEFAULT_MATCH_REPAIR_FILE, help='match repair info file \'s name')
    parser.add_argument('-lr', '--LONG_REPAIR_FILE', type=str,
                        default=DEFAULT_LONG_REPAIR_FILE, help='long data repair file \'s name')
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
    ARGS.MATCH_REPAIR_FILE = os.path.join(
        ARGS.DATA_DIR, ARGS.MATCH_REPAIR_FILE)
    ARGS.LONG_REPAIR_FILE = os.path.join(ARGS.DATA_DIR, ARGS.LONG_REPAIR_FILE)


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


def prepare_match_info(filename):
    match_infoset = {}
    with open(filename, 'r') as f:
        match_infoset = json.loads('['+f.read().strip('\n')[:-1]+']')
    return match_infoset


def get_comp_rev_data(content):
    tran = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return [{
        'name': data['name']+'(comp_rev)',
        's': ''.join([tran[c] for c in data['s']][::-1])
    } for data in content]


TOTAL_DELETE = 0
TOTAL_EDGE = 0

def repair(long_data, match_info):  # match_info已经按照pos排序
    global TOTAL_DELETE, TOTAL_EDGE, MINRATE
    deleted = np.zeros(len(match_info), dtype=bool)
    long = long_data['s']
    voters = [Counter() for i in range(len(long))]
    records = [{'A': [], 'T':[], 'G':[], 'C':[]} for i in range(len(long))]
    new_long = ''
    new_matches = []
    while True:
        need_fix_pos = -1
        need_fix_rate = 1.0  # 峰值的rate越低越要修复
        now_info = 0
        now_que = []
        for now_pos, (voter, record) in enumerate(zip(voters, records)):
            voter.clear()
            record['A'] = []
            record['T'] = []
            record['G'] = []
            record['C'] = []
            # add new match_info to que
            while now_info < len(match_info) and match_info[now_info]['pos'] <= now_pos:
                if deleted[now_info] == False:
                    now_que.append((match_info[now_info], now_info))
                now_info += 1
            # remove no use info
            while len(now_que) > 0 and now_que[0][0]['pos']+len(now_que[0][0]['s']) <= now_pos:
                now_que.pop(0)
            # count now char
            for match, info_pos in now_que:
                c = match['s'][now_pos-match['pos']]
                voter.update(c)
                record[c].append(info_pos)

            top = voter.most_common(1)
            if len(top) == 1:
                rate = top[0][1] / sum(voter.values())
                if need_fix_pos == -1 or need_fix_rate > rate:
                    need_fix_pos = now_pos
                    need_fix_rate = rate

        if need_fix_pos == -1:
            break
        if need_fix_rate >= MINRATE:
            break
        voter = voters[need_fix_pos]
        record = records[need_fix_pos]
        while True:
            tops = voter.most_common(2)
            if len(tops) == 1:
                break
            if tops[0][1] / sum(voter.values()) >= MINRATE:
                break
            for pos in record[tops[1][0]]:
                deleted[pos] = True
            voter[tops[1][0]] = 0

    for now_pos, voter in enumerate(voters):
        voter.update(long[now_pos])
        new_long += voter.most_common(1)[0][0]

    for i, match in enumerate(match_info):
        if deleted[i] == False:
            new_matches.append(match)
        else:
            TOTAL_DELETE += 1
        TOTAL_EDGE += 1

    repair_data = {
        'name': long_data['name']+'(repaired)',
        's': new_long
    }
    return repair_data, new_matches


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
    long_dataset = prepare_fasta_data(ARGS.LONG_FILE)
    long_comp_rev_dataset = get_comp_rev_data(long_dataset)
    long_dataset.extend(long_comp_rev_dataset)  # original long file
    match_infoset = prepare_match_info(ARGS.MATCHES_FILE)  # match infoset
    # update threhold
    MINRATE = 1.0-(PARAM['short_read_error_rate']*2)
    print('Repair long data')
    final_long_dataset = []
    final_match_infoset = []
    for long_data, match_info in zip(tqdm(long_dataset), match_infoset):
        final_long_data, final_match_info = repair(long_data, match_info)
        final_long_dataset.append(final_long_data)
        final_match_infoset.append(final_match_info)
    print('Total delete', TOTAL_DELETE,'of',TOTAL_EDGE)
    print('Write result')
    with open(ARGS.LONG_REPAIR_FILE, 'w') as f:
        for data in final_long_dataset:
            f.write(data['name']+'\n')
            f.write(data['s']+'\n')
    with open(ARGS.MATCH_REPAIR_FILE, 'w') as f:
        for info in final_match_infoset:
            f.write(json.dumps(info)+',\n')
