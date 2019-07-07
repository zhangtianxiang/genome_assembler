'''
如果对某个数据集求匹配信息时，只求了长串的匹配信息（如为了节省时间），但是没有求长串反向互补串的匹配信息，可以通过本程序对匹配结果进行补全
原理为：正向长串匹配的所有短串进行求反向互补，匹配位置也求反向，则可以匹配到长串反向互补串上

需要在代码中指定已有多少长串的数据

输出到fixed_long_extend和matches_extend。
'''
import argparse
import os
import Levenshtein as leve
from tqdm import tqdm
from collections import Counter
import json
import sys

PRE_FIX = 'data4_2/'
LONG_READ_LEN = 1000
SHORT_READ_LEN = 100

ALL = 5000


def prepare_fasta_dataset(filename):
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


def prepare_match_infoset(filename):
    match_infoset = {}
    with open(filename, 'r') as f:
        match_infoset = json.loads('['+f.read().strip('\n')[:-1]+']')
    return match_infoset


def get_comp_rev_dataset(content):
    tran = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return [{
        'name': data['name'][:-7]+'(comp_rev)(fixed)',
        's': ''.join([tran[c] for c in data['s']][::-1])
    } for data in content]


# match_infoset=[match_info=[match]]
def get_comp_rev_match_infoset(match_infoset):
    def get_comp_rev_match_info(match_info):
        def get_comp_rev_match(match):  # match name rev，s rev_comp, pos rev
            tran = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            if match['name'][-1] == ')':
                new_name = match['name'][:-10]
            else:
                new_name = match['name']+'(comp_rev)'
            new_pos = LONG_READ_LEN-match['pos']-SHORT_READ_LEN
            new_match = {
                'name': new_name,
                'pos': new_pos,
                'dis': match['dis'],
                's': ''.join([tran[c] for c in match['s']][::-1])
            }
            return new_match
        return [get_comp_rev_match(match) for match in match_info][::-1]
    return [get_comp_rev_match_info(match_info) for match_info in match_infoset]


def save_fasta(filename, content):
    with open(filename, 'w') as f:
        for data in content:
            f.write(data['name']+'\n')
            f.write(data['s']+'\n')


def save_match(filename, match_infoset):
    with open(filename, 'w') as f:
        for match_info in match_infoset:
            f.write(json.dumps(match_info)+',\n')


fixed_dataset = prepare_fasta_dataset(PRE_FIX+'fixed_long.fasta')
match_infoset = prepare_match_infoset(PRE_FIX+'matches.json')


fixed_dataset = fixed_dataset[:ALL]
match_infoset = match_infoset[:ALL]

fixed_dataset_comp_rev = get_comp_rev_dataset(fixed_dataset)
fixed_dataset.extend(fixed_dataset_comp_rev)

match_infoset_comp_rev = get_comp_rev_match_infoset(match_infoset)
match_infoset.extend(match_infoset_comp_rev)


save_fasta(PRE_FIX+'fixed_long_extend.fasta', fixed_dataset)
save_match(PRE_FIX+'matches_extend.json', match_infoset)
