import json
import argparse
import os
import importlib

DEFAULT_PARAM_FILE = 'param.json'
DEFAULT_LONG_FILE = 'long.fasta'
DEFAULT_SHORT_1_FILE = 'short_1.fasta'
DEFAULT_SHORT_2_FILE = 'short_2.fasta'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('DATA_DIR', type=str, help='the dataset\'s directory')
    parser.add_argument('ALGORITHM', type=str, help='the algorithm to run')
    parser.add_argument('-p', '--PARAM_FILE', type=str,
                        default=DEFAULT_PARAM_FILE, help='the param file\'s name')
    parser.add_argument('-l', '--LONG_FILE', type=str,
                        default=DEFAULT_LONG_FILE, help='the long data file\'s name')
    parser.add_argument('-s1', '--SHORT_1_FILE', type=str,
                        default=DEFAULT_SHORT_1_FILE, help='the short_1 file\'s name')
    parser.add_argument('-s2', '--SHORT_2_FILE', type=str,
                        default=DEFAULT_SHORT_2_FILE, help='the short_2 file\'s name')
    args = parser.parse_args()
    args.PARAM_FILE = os.path.join(args.DATA_DIR, args.PARAM_FILE)
    args.LONG_FILE = os.path.join(args.DATA_DIR, args.LONG_FILE)
    args.SHORT_1_FILE = os.path.join(args.DATA_DIR, args.SHORT_1_FILE)
    args.SHORT_2_FILE = os.path.join(args.DATA_DIR, args.SHORT_2_FILE)
    return args


if __name__ == "__main__":
    args = parse_args()
    # read data
    param = {}
    long_data = []
    short_1_data = []
    short_2_data = []
    with open(args.PARAM_FILE) as f:
        param = json.loads(f.read())
    with open(args.LONG_FILE) as f:
        long_data = f.readlines()
    with open(args.SHORT_1_FILE) as f:
        short_1_data = f.readlines()
    with open(args.SHORT_2_FILE) as f:
        short_2_data = f.readlines()
    # import algorithm
    algo = __import__(args.ALGORITHM)
    # run
    algo.run(param=param, long_data=long_data,
             short_1_data=short_1_data, short_2_data=short_2_data)
