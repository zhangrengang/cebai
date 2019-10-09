import sys
from collections import OrderedDict
from xopen import xopen as open

def main(inDepth=sys.argv[1], outTab=sys.stdout, maxdepth=5):
	d_bins = OrderedDict()
	for line in open(inDepth):
		temp = line.rstrip().split()
		CHR = temp[0]
		POS = int(temp[1])
		DEPTH = map(int, temp[2:])
		TOTAL = sum(DEPTH)
		if TOTAL >= maxdepth:
			TOTAL = maxdepth
		try: d_bins[CHR][TOTAL] += 1
		except KeyError:
			try: d_bins[CHR][TOTAL] = 1
			except KeyError:
				d_bins[CHR] = {TOTAL: 1}
#	f = open(outTab, 'w')
	f = outTab
	line = ['CHR'] + range(maxdepth+1)
	line = map(str, line)
	print >>f, '\t'.join(line)
	for CHR, d_count in d_bins.items():
		chr_length = sum(d_count.values())
		counts = [d_count[i] if i in d_count else 0 for i in range(maxdepth+1)]
		pecents = [100.0*v/chr_length for v in counts]
		line = [CHR] + pecents
		line = map(str, line)
		print >>f, '\t'.join(line)
#	f.close()

if __name__ == '__main__':
	main()
