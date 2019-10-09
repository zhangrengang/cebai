import sys
from collections import OrderedDict
from xopen import xopen as open

def main(inDepth=sys.argv[1], outStat=sys.stdout):
	d_len = OrderedDict()
	d_depth0 = OrderedDict()
	for line in open(inDepth):
		temp = line.strip().split('\t')
		cid, pos, depth = temp[:3]
		depth = int(depth)
		try: d_len[cid] += 1
		except KeyError: d_len[cid] = 1
		if depth == 0:
			try: d_depth0[cid] += 1
			except KeyError: d_depth0[cid] = 1
	for cid, length in d_len.items():
		depth0 = d_depth0.get(cid, 0)
		line = [cid, length, depth0, 100.0*depth0/length]
		line = map(str, line)
		print >> outStat, '\t'.join(line)

if __name__ == '__main__':
	main()
