# -*- coding:utf-8 -*-
'''
使用DBG对结果串的两侧补全
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
import heapq

sys.setrecursionlimit(90000)

DEFAULT_PARAM_FILE = 'param.json'
DEFAULT_LONG_FILE = 'long.fasta'
DEFAULT_SHORT_1_FILE = 'short_1.fasta'
DEFAULT_SHORT_2_FILE = 'short_2.fasta'
DEFAULT_FIXED_LONG_FILE = 'fixed_long.fasta'
DEFAULT_EXTEND_FILE = 'extend_ans.fasta'

ARGS = None
DNALEN = 50


def parse_args():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument('DATA_DIR', type=str, help='the dataset\'s directory')
    parser.add_argument('ANS_FILE', type=str, help='ans data file\'s name')
    parser.add_argument('-fl', '--FIXED_LONG_FILE', type=str,
                        default=DEFAULT_FIXED_LONG_FILE, help='fixed long data file\'s name')
    parser.add_argument('-p', '--PARAM_FILE', type=str,
                        default=DEFAULT_PARAM_FILE, help='the param file\'s name')
    parser.add_argument('-l', '--LONG_FILE', type=str,
                        default=DEFAULT_LONG_FILE, help='the long data file\'s name')
    parser.add_argument('-s1', '--SHORT_1_FILE', type=str,
                        default=DEFAULT_SHORT_1_FILE, help='the short_1 file\'s name')
    parser.add_argument('-s2', '--SHORT_2_FILE', type=str,
                        default=DEFAULT_SHORT_2_FILE, help='the short_2 file\'s name')
    parser.add_argument('-e', '--EXTEND_FILE', type=str,
                        default=DEFAULT_EXTEND_FILE, help='the extended ans file\'s name')
    parser.add_argument('-done', '--DONE', type=int,
                        default=0, help='the done number of ansdna')
    ARGS = parser.parse_args()
    ARGS.PARAM_FILE = os.path.join(ARGS.DATA_DIR, ARGS.PARAM_FILE)
    ARGS.LONG_FILE = os.path.join(ARGS.DATA_DIR, ARGS.LONG_FILE)
    ARGS.SHORT_1_FILE = os.path.join(ARGS.DATA_DIR, ARGS.SHORT_1_FILE)
    ARGS.SHORT_2_FILE = os.path.join(ARGS.DATA_DIR, ARGS.SHORT_2_FILE)
    ARGS.FIXED_LONG_FILE = os.path.join(ARGS.DATA_DIR, ARGS.FIXED_LONG_FILE)
    ARGS.ANS_FILE = os.path.join(ARGS.DATA_DIR, ARGS.ANS_FILE)
    ARGS.EXTEND_FILE = os.path.join(ARGS.DATA_DIR, ARGS.EXTEND_FILE)


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


def get_comp_rev_dna(dna):
    tran = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([tran[c] for c in dna][::-1])


def get_comp_rev_dataset(dataset):
    tran = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return [{
        'name': data['name']+'(comp_rev)',
        's': get_comp_rev_dna(data['s'])
    } for data in dataset]


def generate_graph(datasets):
    point_to_dna = []
    dna_to_point = {}
    out_edges = []
    g = []

    def get_id(dna):
        point = dna_to_point.get(dna)
        if point != None:
            return point
        now = len(point_to_dna)
        dna_to_point[dna] = now
        point_to_dna.append(dna)
        out_edges.append(Counter())
        g.append([])
        return now

    def add_edge(u, v):
        pre = out_edges[u][v]
        out_edges[u].update([v])
        if pre == 0:  # Threshold
            g[u].append(v)

    def generate_from_dna(dna):
        assert len(dna) >= DNALEN
        for i in range(len(dna)-DNALEN):
            u = get_id(dna[i:i+DNALEN])
            v = get_id(dna[i+1:i+1+DNALEN])
            add_edge(u, v)

    for dataset in datasets:
        for data in tqdm(dataset):
            generate_from_dna(data['s'])

    return point_to_dna, dna_to_point, out_edges, g


def dynamic_save(i, s):
    with open(ARGS.EXTEND_FILE, 'a') as f:
        print('{:3d} length: {:d}'.format(i, len(s)))
        f.write('>ans_{0}/1\n{1}\n'.format(i, s))


if __name__ == "__main__":
    # load ARGS
    parse_args()
    # load data
    long_dataset = prepare_fasta_data(ARGS.LONG_FILE)
    rc_long_dataset = get_comp_rev_dataset(long_dataset)
    short1_dataset = prepare_fasta_data(ARGS.SHORT_1_FILE)
    rc_short1_dataset = get_comp_rev_dataset(short1_dataset)
    short2_dataset = prepare_fasta_data(ARGS.SHORT_2_FILE)
    rc_short2_dataset = get_comp_rev_dataset(short2_dataset)
    fixed_dataset = prepare_fasta_data(ARGS.FIXED_LONG_FILE)
    # rc_fixed_dataset = get_comp_rev_dataset(fixed_dataset)
    ans_dataset = prepare_fasta_data(ARGS.ANS_FILE)

    try:
        with open(ARGS.EXTEND_FILE, 'r') as f:
            pass
    except:
        with open(ARGS.EXTEND_FILE, 'w') as f:
            pass
    # get dbg points
    print('Building DB graph')
    point_to_dna, dna_to_point, out_edges, g = generate_graph(
        [short1_dataset, short2_dataset, fixed_dataset, rc_short1_dataset, rc_short2_dataset])
    total_point = len(point_to_dna)
    print('Total points', total_point)

    # delete unuse
    del long_dataset, rc_long_dataset, short1_dataset, rc_short1_dataset, short2_dataset, rc_short2_dataset, fixed_dataset #, rc_fixed_dataset
    del out_edges
    # generate sequence from a start point
    # def gen_chain(u, maxd=100000):
    #     res = []
    #     while True:
    #         top_edges = out_edges[u].most_common(1)
    #         if len(top_edges) == 0:
    #             break
    #         u = top_edges[0][0]
    #         res.append(u)
    #         if len(res) == maxd:
    #             break
    #     return res

    def gen_chain_v2(u, maxd=100000):  # 树形dp(伪)
        res = []
        # 是否访问过，访问过则认为反向边不存在，理论上最长路最多减少一半的长度
        visited = {}
        maxto = []
        orgid = []
        def get_id(u):
            ret = visited.get(u)
            if ret == None:
                newid = len(maxto)
                visited[u] = newid
                maxto.append(-1)
                orgid.append(u)
                return newid, False # new id
            else:
                return ret, True # old id
        
        def dfs(u, dep=0):
            uid, _ = get_id(u)
            maxdeepu = 1
            if dep == maxd:
                # print('Reach maxd!')
                return maxdeepu
            for v in g[u]:
                vid, exist = get_id(v)
                if exist:
                    continue
                maxdeepv = dfs(v, dep+1)
                if maxdeepv+1 > maxdeepu:
                    maxdeepu = maxdeepv+1
                    maxto[uid] = vid
            return maxdeepu

        dfs(u)
        u = maxto[visited[u]]
        while u != -1:
            res.append(orgid[u])
            u = maxto[u]
        return res

    # generate str from a start point
    def gen_str(u, maxd=100000):
        res = point_to_dna[u]
        points = gen_chain_v2(u, maxd=maxd)
        for u in points:
            res += point_to_dna[u][-1]
        return res

    # extend all ans
    def extend_tail(dna):
        # 向后extend
        TRY_STEP = 200
        assert len(dna) >= DNALEN
        step = min(TRY_STEP, len(dna)-DNALEN+1)
        max_ext = 0
        best_res = dna
        for i in range(step):
            seg = dna[-(DNALEN+i):][:DNALEN]
            st = dna_to_point.get(seg)
            if st != None:
                ext = gen_str(st, maxd = 40000)
                extlen = len(ext)-DNALEN-i
                if extlen <= max_ext:
                    # print('extend at', i, 'failed', extlen, 'continue')
                    continue
                else:
                    print('extend at', i, 'length', extlen)
                    max_ext = extlen
                    best_res = dna[:-(DNALEN+i)]+ext
                    continue
        return best_res

    # ans = []
    for i, data in enumerate(ans_dataset):
        print('Now new dna', i)
        if i < ARGS.DONE:
            print('Skip')
            continue
        pre = len(data['s'])
        print('extend tail')
        res = extend_tail(data['s'])
        res = get_comp_rev_dna(res)
        print('extend head')
        res = extend_tail(res)
        end = len(res)
        print('Pre', pre, 'Extended', end-pre)
        res = get_comp_rev_dna(res)
        dynamic_save(i, res)
        # ans.append(res)

    # ans = sorted(ans, key=lambda x: len(x), reverse=True)
    # with open(ARGS.EXTEND_FILE, 'w') as f:
    #     for i, s in enumerate(ans):
    #         print('{:3d} length: {:d}'.format(i, len(s)))
    #         f.write('>ans_{0}/1\n{1}\n'.format(i, s))
