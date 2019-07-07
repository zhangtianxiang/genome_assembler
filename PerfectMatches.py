'''
将A集合（短串集合）匹配到B串（长串）上，使用汉明距离，每个串找最小距离位置，并为距离设定阈值
'''
import argparse
import os
import Levenshtein as leve
from tqdm import tqdm
from collections import Counter
import math

# SHORT_LEN = 100
# LONG_LEN = 1000


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('FILE1', type=str)
    parser.add_argument('FILE2', type=str)
    parser.add_argument('OUT', type=str)
    args = parser.parse_args()
    return args


def match(short, long):
    mindis = -1
    minpos = 0
    for i in range(len(long)-len(short)+1):
        dis = leve.hamming(short, long[i:i+len(short)])
        if mindis == -1 or dis < mindis:
            mindis = dis
            minpos = i
    # tail
    for i in range(len(long)-len(short)+1, len(long)-(len(short)//2)+1):
        common = len(long)-i
        # ceil(dis*LEN_SHORT/LEN_COMMON)
        dis = (leve.hamming(short[:common],
                            long[i:i+common])*len(short)+common-1)//common
        if mindis == -1 or dis < mindis:
            mindis = dis
            minpos = i
    # head
    for i in range(-(len(short)//2), 0):
        common = len(short)+i
        # ceil(dis*LEN_SHORT/LEN_COMMON)
        dis = (leve.hamming(short[-i:], long[:common])
               * len(short)+common-1)//common
        if mindis == -1 or dis < mindis:
            mindis = dis
            minpos = i
    return mindis, minpos


def fix(long, shorts):
    print('Fixing long')
    new_shorts = []
    now_short = 0
    now_que = []
    new_long = ''
    for now_pos in tqdm(range(len(long))):
        # add new shorts to que
        while now_short < len(shorts) and shorts[now_short]['pos'] <= now_pos:
            now_que.append(shorts[now_short])
            now_short += 1
        # save survival shorts to new_shorts
        while len(now_que) > 0 and now_que[0]['pos']+len(now_que[0]['s']) <= now_pos:
            new_shorts.append(now_que[0])
            now_que.pop(0)
        # count now char
        cnt = {
            'A': 0,
            'T': 0,
            'G': 0,
            'C': 0
        }
        cnt[long[now_pos]] = 1
        cnt_max = 1
        new_char = long[now_pos]
        for short in now_que:
            cnt[short['s'][now_pos-short['pos']]] += 1
        for k, v in cnt.items():
            if v > cnt_max:
                cnt_max = v
                new_char = k
        new_long += new_char
        # remove not correct shorts
        for i in range(len(now_que)):
            while i < len(now_que) and now_que[i]['s'][now_pos-now_que[i]['pos']] != new_char:
                now_que.pop(i)

    # save survival shorts to new_shorts
    while len(now_que) > 0:
        new_shorts.append(now_que[0])
        now_que.pop(0)
    return new_long, new_shorts


def save(filename, long, new_long, shorts):
    if shorts[0]['pos'] >= 0:
        pre = 0
    else:
        pre = -shorts[0]['pos']
    with open(filename, 'w') as f:
        f.write(' '*pre+long+'\n')
        f.write(' '*pre+new_long+'\n')
        for d in shorts:
            f.write((' '*(pre+d['pos']))+d['s']+'\n')
        f.write(' '*pre+new_long+'\n')
        f.write(' '*pre+long+'\n')


if __name__ == "__main__":
    args = parse_args()
    # read data
    A = []
    B = ''
    with open(args.FILE1, 'r') as f:
        A = f.readlines()
    if len(A) & 1:
        A.pop()
    assert (len(A) & 1) == 0
    with open(args.FILE2, 'r') as f:
        B = f.read()
        B.strip('\n')
    maxdis = 14
    ans = []
    print('Matching shorts to long')
    for j, s in enumerate(tqdm(A)):
        if (j & 1) == 0:
            continue
        s = s.strip('\n')
        mindis, minpos = match(s, B)
        if mindis < maxdis:
            ans.append({'s': s, 'pos': minpos, 'dis': mindis})
    ans = sorted(ans, key=lambda x: x['pos'])
    print('Total pre', len(ans))
    new_long, new_shorts = fix(B, ans)  # use ans fix B
    print('Total final', len(new_shorts))
    # save
    save(args.OUT, B, new_long, new_shorts)
'''
O(n)汉明距离
头部尾部匹配？
ssssssssssssss
      llllllll


llllllllllllll
         ssssssssss
公共段长为x,x至少为50
汉明距离应当小于x*0.2时认为匹配
'''
