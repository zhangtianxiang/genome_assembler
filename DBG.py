# -*- coding:utf-8 -*-
'''
最开始yy的DB图测试：效果极差
后来对算法改进并用于修复Assemble的结果，见DBGComplate
'''
import numpy as np
import sys

from collections import Counter, deque
from tqdm import tqdm

sys.setrecursionlimit(90000)
SEGLEN = 30
EPOCH = 20
ANS_FILE = 'data1/contig_dbg.fasta'
FL_FILE = 'data1/fixed_long.fasta'
S1_FILE = 'data1/short_1.fasta'
S2_FILE = 'data1/short_2.fasta'


class MyVariable:
    def __init__(self, val=0):
        self.val = val

    def value(self):
        return self.val

    def inc(self, incv=1):
        self.val += incv
        return self.val

    def dec(self, decv=1):
        self.val -= decv
        return self.val

    def reset(self, val=0):
        self.val = val


total_point = 0
total_edge = 0
total_scc = 0


def generate_graph(param, data_parts):
    global total_point
    seg_to_point = {}
    point_to_seg = []
    out_edges = []
    g = []
    rg = []  # 反图
    ind = []
    outd = []

    def get_id(seg):
        point = seg_to_point.get(seg)
        if point != None:
            return point
        now = len(point_to_seg)
        seg_to_point[seg] = now
        point_to_seg.append(seg)
        out_edges.append(Counter())
        g.append([])
        rg.append([])
        ind.append(0)
        outd.append(0)
        return now

    def add_edge(u, v):
        global total_edge
        assert u != v  # 在选取SEGLEN时不应该产生自环
        pre = out_edges[u][v]
        out_edges[u].update([v])
        if pre == 2:  # 阈值
            g[u].append(v)
            ind[v] += 1
            outd[u] += 1
            rg[v].append(u)
            total_edge = total_edge+1

    def generate_from_str(s):
        assert len(s) >= SEGLEN
        for i in range(len(s)-SEGLEN):
            u = get_id(s[i:i+SEGLEN])
            v = get_id(s[i+1:i+1+SEGLEN])
            add_edge(u, v)
            # add_edge(v, u)

    for data_part in data_parts:
        for i, line in enumerate(tqdm(data_part)):
            if (i & 1) == 0:
                continue
            if line == '':
                break
            # pre = len(line)
            line = line.strip('\n')
            # aft = len(line)
            # if pre != aft:
            #     print('pre', pre, 'aft', aft)
            generate_from_str(line)
    total_point = len(point_to_seg)
    return seg_to_point, point_to_seg, out_edges, g, rg, ind, outd


def get_SCC(g):
    global total_point, total_scc
    # 求强联通分量，tarjan算法，O(m)
    dfn = np.zeros((total_point,), dtype=int)  # dfn(i) = -1
    low = np.zeros((total_point,), dtype=int)  # low(i) = -1
    ins = np.zeros((total_point,), dtype=bool)  # instack
    sta = []  # stack
    belong = -np.ones((total_point), dtype=int)  # belong(i) = -1
    sz = [0]
    idx = MyVariable()

    def dfs(u):
        global total_scc
        nowid = idx.inc()
        dfn[u] = nowid
        low[u] = nowid
        ins[u] = True
        sta.append(u)
        for v in g[u]:
            if dfn[v] == 0:
                dfs(v)
                low[u] = min(low[u], low[v])
            elif ins[v]:
                low[u] = min(low[u], dfn[v])
        if dfn[u] == low[u]:
            total_scc += 1
            v = -1
            nowsz = 0
            while u != v:
                v = sta[-1]
                sta.pop()
                ins[v] = False
                belong[v] = total_scc
                nowsz += 1
            sz.append(nowsz)

    for u in tqdm(range(total_point)):
        if dfn[u] == 0:
            dfs(u)

    return belong, sz


def get_max_chain(g, ind):
    global total_point
    # get max_chain and remove the edges on it
    ind_rec = ind.copy()
    dp = np.ones(total_point, dtype=int)
    fr = -np.ones(total_point, dtype=int)
    max_len_p = MyVariable(0)
    q = deque()
    for u in range(total_point):
        if ind[u] == 0:
            q.append(u)
    while len(q) > 0:
        u = q[0]
        q.popleft()
        for v in g[u]:
            ind[v] -= 1
            if dp[u]+1 > dp[v]:
                dp[v] = dp[u]+1
                fr[v] = u
                if dp[v] > dp[max_len_p.value()]:
                    max_len_p.reset(v)
            if ind[v] == 0:
                q.append(v)
    
    print('Max length:', dp[max_len_p.value()])
    ind = ind_rec
    v = max_len_p.value()
    ans = [v]
    while fr[v] != -1:
        u = fr[v]
        g[u].remove(v)
        ind[v] -= 1
        ans.append(u)
        v = u
    return ans[::-1], g, ind

def seq_to_dna(seq, point_to_seg):
    ret = point_to_seg[seq[0]][:-1]
    for point in seq:
        ret += point_to_seg[point][-1]
    return ret

def run(param, long_data, short_1_data, short_2_data):
    print('processing', param.get('name'))
    assert len(short_1_data) == len(short_2_data)
    seg_to_point, point_to_seg, out_edges, g, rg, ind, outd = generate_graph(
        param, [short_1_data, short_2_data, long_data])
    print('Total points:', total_point)
    print('Total edges:', total_edge)
    belong, sz = get_SCC(g)  # 得到所有的强联通分量
    print('Total SCC:', total_scc)
    '''
    因为不存在点数大于1的SCC，所以可以直接求最长链。
    DAG求最长链，可以从入度为0的点，DP找到最长链，将其取出，将相应的边删除？？？如果存在交叉变异？
    '''
    tmpg = g.copy()
    tmpind = ind.copy()
    epoch = EPOCH
    ans = []
    for i in range(epoch):
        seq, tmpg, tmpind = get_max_chain(tmpg, tmpind)
        ans.append(seq)
    
    with open(ANS_FILE,'w') as f:
        for i in range(epoch):
            f.write('>ans_'+str(i)+'/1\n')
            f.write(seq_to_dna(ans[i], point_to_seg))
            f.write('\n')


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

if __name__ == "__main__":
    fl = prepare_fasta_data(FL_FILE)
    s1 = prepare_fasta_data(S1_FILE)
    s2 = prepare_fasta_data(S2_FILE)