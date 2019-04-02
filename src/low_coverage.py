import sys

def main(inDepth=sys.argv[1], outRegion=sys.stdout, minDepth=3, flank_length=130, min_length=5):
	LOW_POS = []
	for line in open(inDepth):
		temp = line.strip().split()
		CHR, POS, DEPTH = temp[:3]
		POS = int(POS)
		DEPTH = int(DEPTH)
		if DEPTH <= minDepth:
			LOW_POS.append((CHR, POS))
	REGIONs = pos2region(LOW_POS, flank_length=flank_length, min_length=min_length)
	for (CHR, START, END) in REGIONs:
		START = START-1
		line = [CHR, START, END]
		line = map(str, line)
		print >> outRegion, '\t'.join(line)
	
def pos2region(POSs,flank_length=0,min_length=1):
	REGIONs = []
	d_max_depth = {}
	# REGION = [[pos1, pos2, pos3],[pos4]]
	for CHR, POS in POSs:
		d_max_depth[CHR] = POS
		try:
			last_chr, last_pos = REGIONs[-1][-1]
			if CHR == last_chr and POS-last_pos <= 1 + flank_length*2:
				REGIONs[-1] += [(CHR, POS)]
			else:
				REGIONs += [[(CHR, POS)]]
		except IndexError:
			REGIONs += [[(CHR, POS)]]
#	print >>sys.stderr, REGIONs
	REGIONs2 = []
	for REGION in REGIONs:
		CHR, START = REGION[0]
		CHR, END = REGION[-1]
		if not (START-END+1) <= min_length:
			continue
		start = max(1, START-flank_length)
		end = min(d_max_depth[CHR], END+flank_length)
		REGIONs2 += [(CHR, start, end)]
#	print >>sys.stderr, REGIONs2
	REGIONs2 = complement_region(REGIONs2, d_max_depth)
	return REGIONs2

def complement_region(REGIONs, d_max_depth):
	from collections import OrderedDict
	d_region = OrderedDict()
	for (CHR, START, END) in REGIONs:
		try: d_region[CHR].append((START, END))
		except KeyError: d_region[CHR] = [(START, END)]
	REGIONs2 = []
	for CHR, REGIONs in d_region.items():
		i = 0
		START, END = REGIONs[i]
		if START > 1:
			REGIONs2 += [(CHR, 1, START-1)]
		last_start, last_end = START, END
		for START, END in REGIONs[1:]:
			if START-1 >= last_end+1:
				REGIONs2 += [(CHR, last_end+1, START-1)]
			else:
				print >>sys.stderr, [CHR, (last_start, last_end), (START, END)], 'ignored'
			last_start, last_end = START, END
		START = d_max_depth[CHR]
		if START-1 >= last_end+1:
			REGIONs2 += [(CHR, last_end+1, START-1)]
	return REGIONs2
				
if __name__ == '__main__':
	main()
