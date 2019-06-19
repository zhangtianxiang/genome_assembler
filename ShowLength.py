import argparse


parser = argparse.ArgumentParser()
parser.add_argument('FILE', type=str)
ARGS = parser.parse_args()


def prepare_fasta_data(filename):
    content = []
    print('Load data', filename)
    with open(filename, 'r') as f:
        lines = f.readlines()
        name = 'Unknown'
        for i, line in enumerate(lines):
            if (i & 1) == 0:
                name = line.strip('\n')
            else:
                content.append({'name': name, 's': line.strip('\n')})
    return content


dataset = prepare_fasta_data(ARGS.FILE)
for data in dataset:
    print(data['name'], ':', len(data['s']))
