'''
对某个fasta文件求反向互补串，行行对应，可以用于对答案去重，但是最终并没有使用
'''
import argparse
import os
import Levenshtein as leve
from tqdm import tqdm
from collections import Counter
import json
import sys


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


def save_fasta(filename, content):
    with open(filename, 'w') as f:
        for data in content:
            f.write(data['name']+'\n')
            f.write(data['s']+'\n')


fa = prepare_fasta_data('data1/contigfixed.fasta')
fa_c = get_comp_rev_data(fa)
save_fasta('data1/contigfixed2.fasta', fa_c)
