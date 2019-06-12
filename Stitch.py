PRE = 'data4_2/'
stitches = [
	'fixed_long.fasta',
	'fixed_long_50_1000.fasta',
	'fixed_long_1000_2000.fasta',
	'fixed_long_2000_3000.fasta',
	'fixed_long_3000_4000.fasta',
	'fixed_long_4000_5000.fasta'
	]


stitches_2 = [
	'matches.json',
	'matches_50_1000.json',
	'matches_1000_2000.json',
	'matches_2000_3000.json',
	'matches_3000_4000.json',
	'matches_4000_5000.json'
	]

lines = []

for filename in stitches_2:
	with open(PRE+filename, 'r') as f:
		for line in f.readlines():
			line = line.strip('\n');
			if line == '':
				continue
			lines.append(line)

with open(PRE+'matches_all.json', 'w') as f:
	for line in lines:
		f.write(line+'\n')