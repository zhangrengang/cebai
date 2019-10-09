
import sys
from math import ceil
from collections import OrderedDict
from xopen import xopen as open

def _init_bins(inFai, window_step):
	d_bins = OrderedDict()
	for line in open(inFai):
		CHR, LENGTH = line.rstrip().split('\t')[:2]
		LENGTH = int(LENGTH)
		for i in range(0, LENGTH, window_step):
			BIN = i/window_step
			try: d_bins[CHR][BIN] = []
			except KeyError: d_bins[CHR] = OrderedDict([(BIN, [])])
	return d_bins

def main(inDepth=sys.argv[1], inFai=sys.argv[2], outTab=sys.stdout, window_size=5000, window_step=2500):
	d_bins = _init_bins(inFai, window_step)
	for line in open(inDepth):
		temp = line.rstrip().split('\t')
		CHR, POS = temp[:2]
		DEPTH = temp[2:]
		POS = int(POS)
		DEPTH = map(float, DEPTH)
		first = int(ceil(1.0*(POS - window_size)/window_step))
		if first < 0:
			first = 0
		last = int(ceil(1.0*POS/window_step))
		for idx in range(first, last):
			if not d_bins[CHR][idx]:
				d_bins[CHR][idx] = DEPTH
			else:
				d_bins[CHR][idx] = [v1+v2 for v1, v2 in zip(d_bins[CHR][idx],DEPTH)]
	
	for CHR, d_bin in d_bins.items():
		for idx, sum_depth in d_bin.items():
			DEPTH = [1.0*v/window_size for v in sum_depth]
			start = idx*window_step
			end = start+window_size
			line = [CHR, start, end] + DEPTH
			line = map(str, line)
			print >> outTab, '\t'.join(line)

if __name__ == '__main__':
	main()
