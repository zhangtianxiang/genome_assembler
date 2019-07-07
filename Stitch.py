'''
拼接PerfectMatchesAllMulti的结果
'''

PREFIX = 'data4/'
STITCH_LIST = [
    'fixed_long_0_1000.fasta',
    'fixed_long_1000_2000.fasta',
    'fixed_long_2000_3000.fasta',
    'fixed_long_3000_4000.fasta',
    'fixed_long_4000_5000.fasta'
]

STITCH_RESULT = 'fixed_long.json'

# STITCH_LIST = [
#     'matches_0_1000.json',
#     'matches_1000_2000.json',
#     'matches_2000_3000.json',
#     'matches_3000_4000.json',
#     'matches_4000_5000.json'
# ]

# STITCH_RESULT = 'matches.json'

lines = []

for filename in STITCH_LIST:
    with open(PREFIX+filename, 'r') as f:
        for line in f.readlines():
            line = line.strip('\n')
            if line == '':
                continue
            lines.append(line)

with open(PREFIX+STITCH_RESULT, 'w') as f:
    for line in lines:
        f.write(line+'\n')
