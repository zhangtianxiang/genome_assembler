import argparse
import os
import Levenshtein as leve


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('FILE1', type=str)
    parser.add_argument('FILE2', type=str)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    # read data
    A = ''
    B = ''
    with open(args.FILE1, 'r') as f:
        A = f.read()
        A.strip('\n')
    with open(args.FILE2, 'r') as f:
        B = f.read()
        B.strip('\n')
    # A in B
    assert len(A) < len(B)
    mindis = -1
    minpos = 0
    for i in range(len(B)-len(A)+1):
        dis = leve.distance(A,B[i:i+len(A)])
        if mindis == -1 or dis < mindis:
            mindis = dis
            minpos = i
    print('min distance', mindis)
    print((' '*minpos)+A)
    print(B)
'''
O(n)汉明距离
O(nlogn)jaccard
'''