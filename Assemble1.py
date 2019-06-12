'''
思想：
已有修复后的长传fixed_long.fasta和匹配信息matches.json
每个长串匹配了众多小串
长串做为A类点，小串作为B类点

B类点两两之间通过A类点计算偏置的众数和总连接数
最终每个B类点两两之间存在边，有总连接数，偏置众数和数量
以偏置众数为第一关键字，总连接数为第二关键字，求最小生成树？或生成链
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
import numpy as np

sys.setrecursionlimit(90000)

DEFAULT_PARAM_FILE = 'param.json'
DEFAULT_LONG_FILE = 'long.fasta'
DEFAULT_SHORT_1_FILE = 'short_1.fasta'
DEFAULT_SHORT_2_FILE = 'short_2.fasta'
DEFAULT_FIXED_LONG_FILE = 'fixed_long.fasta'
DEFAULT_MATCHES_FILE = 'matches.json'
DEFAULT_ANS_FILE = 'contig.fasta'

ARGS = None
PARAM = {}

MAXDIS = 0

MINCNT = 0


def parse_args():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument('DATA_DIR', type=str, help='the dataset\'s directory')
    parser.add_argument('-a', '--ANS_FILE', type=str,
                        default=DEFAULT_ANS_FILE, help='ans data file\'s name')
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
    ARGS.ANS_FILE = os.path.join(ARGS.DATA_DIR, ARGS.ANS_FILE)


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


def try_merge(A, B, off):
    res = ''
    if off < 0:
        A, B = B, A
        off = -off
    common = min(len(A)-off, len(B))
    if common <= 10:  # Threshold
        return res
    error_rate = leve.hamming(A[off:off+common], B[:common])/common
    if error_rate > 0.5:  # Threshold
        return res
    res = A[:off]+B
    return res


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
        return sid

    def encode_edge(A, B, dis, pos):
        return (dis, pos, A, B)  # 以元组存储便于比较，dis为第一关键字，pos为第二关键字

    def decode_edge(edge):  # return dis, pos, A, B
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

    print('total A', len(name_to_Apoint))

    # generate B point and add edges
    print('Generate B point and build graph')
    for i, match_info in enumerate(tqdm(match_infoset)):
        for match in match_info:
            short_name = match['name']
            short = match['s']
            sid = get_id(short_name, short)
            add_edge(i, sid, match['dis'], match['pos'])

    # build Apoint graph
    offset = np.zeros((len(Apoint_to_name), len(Apoint_to_name)), dtype=int)
    count = np.zeros((len(Apoint_to_name), len(Apoint_to_name)), dtype=int)
    final_edges = []  # count,offset,u,v
    for u in range(len(Apoint_to_name)):
        A_out_edges[u] = sorted(
            A_out_edges[u], key=lambda x: x[3])  # sort by (B, )
    print('Build final graph')
    for u in tqdm(range(len(Apoint_to_name))):
        u_out_edges = A_out_edges[u]
        for v in range(u+1, len(Apoint_to_name)):  # 多个 B ？需要枚举B？
            v_out_edges = A_out_edges[v]
            counter = Counter()
            j = 0
            for i, edge in enumerate(u_out_edges):
                while j < len(v_out_edges) and v_out_edges[j][3] < edge[3]:
                    j += 1
                # same B
                k = j
                while k < len(v_out_edges) and v_out_edges[k][3] == edge[3]:
                    # pos u - pos v
                    counter.update(Counter([edge[1]-v_out_edges[k][1]]))
                    k += 1
            item = counter.most_common(1)
            if len(item) == 0:  # not match totally
                pass
            else:
                item = item[0]
                offset[u][v] = item[0]  # offset
                count[u][v] = item[1]  # count
                if count[u][v] > MINCNT:
                    final_edges.append((item[1], item[0], u, v))
    # generate final DNA!
    print('Generate final DNA')
    final_edges = sorted(final_edges, key=lambda x: x[0], reverse=True)
    print('Total edge', len(final_edges))
    fa = [-1 for i in range(len(Apoint_to_name))]
    off_to_fa_val = [0 for i in range(len(Apoint_to_name))]
    val = [long_data['s'] for long_data in fixed_dataset]

    def get_fa_and_update_off(u):
        if fa[u] == -1:
            return u
        else:  # update
            fafa = get_fa_and_update_off(fa[u])
            if fa[u] == fafa:
                return fa[u]
            off_to_fa_val[u] += off_to_fa_val[fa[u]]
            fa[u] = fafa
            return fafa

    # with open('tmp.fasta', 'w') as f:
    #     edge = final_edges[1000]
    #     cnt, off, u, v = edge
    #     print(edge)
    #     if off < 0:
    #         u, v = v, u
    #         off = -off
    #     f.write(val[u]+'\n')
    #     f.write(' '*off+val[v])

    for edge in tqdm(final_edges):
        cnt, off, u, v = edge
        fu = get_fa_and_update_off(u)
        fv = get_fa_and_update_off(v)
        if fu == fv:
            continue
        nowoff = off_to_fa_val[u]+off-off_to_fa_val[v]
        res = try_merge(val[fu], val[fv], nowoff)
        if res == '':  # merge failed
            continue
        if nowoff > 0:  # fu on left, fu is new father
            fa[fv] = fu
            off_to_fa_val[fv] = nowoff
            val[fu] = res
            val[fv] = None
        else:  # fv on left, fv is new father
            fa[fu] = fv
            off_to_fa_val[fu] = -nowoff
            val[fv] = res
            val[fu] = None

    ans = []
    for u in tqdm(range(len(Apoint_to_name))):
        if get_fa_and_update_off(u) == u:
            ans.append(val[u])
    ans = sorted(ans, key=lambda x: len(x), reverse=True)[:15]
    with open(ARGS.ANS_FILE, 'w') as f:
        for i, s in enumerate(ans):
            print('{:3d} length: {:d}'.format(i, len(s)))
            f.write('>ans_{0}/1\n{1}\n'.format(i, s))
