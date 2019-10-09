
import sys
import numpy as np
from math import ceil
from collections import OrderedDict
from xopen import xopen as open

def main(inDepth=sys.argv[1], window_size=20000, window_step=10000, minlength=0.01):
	outBinDeapth = inDepth + '.bins'
	outplot = outBinDeapth + '.pdf'
	d_bins = OrderedDict()
	d_max_pos = OrderedDict()
	for line in open(inDepth):
		temp = line.rstrip().split()
		CHR = temp[0]
		POS = int(temp[1])
		DEPTH = map(int, temp[2:])
		TOTAL = sum(DEPTH)
		BIN_from = int(ceil(1.0 * POS / window_size))
		BIN_to = int(ceil(1.0 * (POS + window_step) / window_size))
		for BIN in range(BIN_from, BIN_to+1):
			try: d_bins[CHR][BIN] += TOTAL
			except KeyError: 
				try: d_bins[CHR][BIN] = TOTAL
				except KeyError: d_bins[CHR] = OrderedDict([(BIN, TOTAL)])
		d_max_pos[CHR] = max(d_max_pos.get(CHR,0), POS)

	for k,v in d_bins.items():
		for b, t in v.items():
			d_bins[k][b] = 1.0 * t / window_size
	f = open(outBinDeapth, 'w')
	for k,v in d_bins.items():
		for b, t in v.items():
			start = (b-1)*window_step + 1
			print >>f, '%s\t%s\t%s' % (k, start, t)
	f.close()
	last_start = 0
	Xs,Ys = [], []
	labels, label_x, vlines = [], [], []
	for (CHR, BINs) in d_bins.items():
		length = d_max_pos[CHR]
		x, y = [], []
		for BIN, depth in BINs.items():
			start = (BIN - 1)*window_size
			start += last_start
			x += [start]
			y += [depth]
		Xs += [x]
		Ys += [y]
		last_start += length
		labels += [CHR]
		label_x += [last_start - length/2]
		vlines += [last_start]
	tot_len = sum(d_max_pos.values())
	vis_labels = set([k for k, v in d_max_pos.items() if 1.0*v/tot_len >= minlength])
	height_width_ratio = sum(d_max_pos.values()) / max(d_max_pos.values())
	bin_plot(Xs,Ys, labels, label_x, vlines, height_width_ratio, outplot, vis_labels)
def bin_plot(Xs,Ys, labels, label_x, vlines, height_width_ratio, outplot, vis_labels):
	import matplotlib.pyplot as plt
	plt.figure(figsize=(5*height_width_ratio,5))
#	ny, sumy = 0, 0
	Y1s = []
	wsize = 100
	for x, y in zip(Xs,Ys):
		plt.plot(x, y)
#		ny += len(y)
#		sumy += sum(y)
		Y1s += y
		if len(x) > wsize:
			X, Y = [], []
			for i in range(len(x)):
				X += [x[i]]
				s = max(0, i-wsize/2)
				e = min(i+wsize/2, len(x))
				Y += [np.mean(y[s:e])] #[sum(y[s:e]) / (e-s)]
			plt.plot(X, Y, ls='--', lw=2, color="black")
#	ylim = sumy / ny * 2.5
	ylim = np.median(Y1s) * 2.5
	ymax = ylim
	for v in vlines:
		plt.vlines(v, 0, ymax, color="grey")
	for x, label in zip(label_x, labels):
	#	if not label.startswith('chr'):
	#		continue
		if not label in vis_labels:
			continue
		y = -ymax/30
		plt.text(x, y, label, horizontalalignment='center',verticalalignment='top',fontsize=15) #, rotation=30)
	plt.ylabel('depth')
	plt.ylim(0, ylim)
	xlim = max(vlines)
	plt.xlim(0, xlim)
	plt.savefig(outplot, bbox_inches='tight')
if __name__ == '__main__':
	main()

