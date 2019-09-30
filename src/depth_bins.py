
import sys
from math import ceil
from collections import OrderedDict
from xopen import xopen as open
def main(inDepth=sys.argv[1], outBinDeapth='mapping_bins.depth', outAvgDepth='mapping_avg.depth', window_size=5000):

	d_bins = OrderedDict()
	d_avg = OrderedDict()
	d_max_pos = OrderedDict()
	for line in open(inDepth):
		temp = line.rstrip().split()
		CHR = temp[0]
		POS = int(temp[1])
		DEPTH = map(int, temp[2:])
		TOTAL = sum(DEPTH)
		try: d_avg[CHR] += TOTAL
		except KeyError: d_avg[CHR] = TOTAL
		BIN = int(ceil(1.0 * POS / window_size))
		try: d_bins[CHR][BIN] += TOTAL
		except KeyError: 
			try: d_bins[CHR][BIN] = TOTAL
			except KeyError: d_bins[CHR] = OrderedDict([(BIN, TOTAL)])
		d_max_pos[CHR] = POS

	for k,v in d_avg.items():
		d_avg[k] = 1.0 * v / d_max_pos[k]
	for k,v in d_bins.items():
		for b, t in v.items():
			d_bins[k][b] = 1.0 * t / window_size
	f = open(outAvgDepth , 'w')
	for k,v in d_avg.items():
		print >> f, '%s\t%s\t%s' % (k,v, d_max_pos[k])
	f.close()
	f = open(outBinDeapth, 'w')
	for k,v in d_bins.items():
		for b, t in v.items():
			start = (b-1)*window_size + 1
			print >>f, '%s\t%s\t%s' % (k, start, t)
	f.close()

if __name__ == '__main__':
	main()

