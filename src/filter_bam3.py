#coding: utf-8
'''use AS (alignment score) instead of id (identity)'''
import sys, os
import pp
import pysam
import numpy as np
import heapq

# 弃用，不并行太慢
def main(bam=sys.argv[1], outBam=sys.argv[2], min_as=120, min_cov=0.9, wsize=100, wstep=50, bins = np.linspace(0,1,101), min_depth=3):
	bamfile = pysam.AlignmentFile(bam, 'rb')
	xi, xj, xm =0,0,0
#	d_ref = {r:l for r,l in zip(bamfile.references, bamfile.lengths)}
	d_hq = {}
	xm, xn = 0,0
	for ref, length in zip(bamfile.references, bamfile.lengths):
		for start in range(0, length, wstep):
			end = start + wsize
			d_as = {rc.query_name: rc.get_tag('AS') for rc in bamfile.fetch(ref, start, end) if rc.is_proper_pair and is_covered(rc, start, end) and get_cov(rc) >= min_cov}
			if len(d_as) < min_depth:
				continue
			peak, border = find_peak_border(d_as.values(), bins=bins)
			if border < min_as:
				border = min_as
			hq_ids = set(qid for qid, id in d_id.items() if id >= border)
			print >>sys.stderr, ref, start, end, len(d_id), len(hq_ids), peak, border
			xm += len(hq_ids)
			xn += len(d_id)
			try: d_hq[ref] = d_hq[ref] | hq_ids
			except KeyError: d_hq[ref] = hq_ids
	print 'total {}, get {}'.format(xn ,xm)
	bamout = pysam.AlignmentFile(outBam, 'wb', template=bamfile)
	xi, xj, xm, xn =0,0,0,0
	for rc in bamfile:
		xi += 1
		if not rc.is_proper_pair:
			continue
		xj += 1
		cov = get_cov(rc)
		if not cov >= min_cov:
			continue
		xm += 1
		if rc.query_name in d_hq[rc.reference_name]:
			bamout.write(rc)
			xn += 1
	bamfile.close()
	bamout.close()
	print 'total {} reads, after is_proper_pair filter {}, after cov filter {}, after id filter {}'.format(xi,xj,xm, xn)
def main_pp(bam=sys.argv[1], outBam=sys.argv[2], min_as=60, min_cov=0.9, wsize=100, wstep=50, bins = np.linspace(0,151,101), min_depth=3): # 因调用pp，参数未传递
	bamfile = pysam.AlignmentFile(bam, 'rb')
	d_hq = pp_run(bam, bamfile.references, bamfile.lengths)
	bamout = pysam.AlignmentFile(outBam, 'wb', template=bamfile)
	xi, xj, xm, xn =0,0,0,0
	for rc in bamfile:
		xi += 1
		if not rc.is_proper_pair:
			continue
		xj += 1
		cov = get_cov(rc)
		if not cov >= min_cov:
			continue
		xm += 1
		if rc.query_name in d_hq[rc.reference_name]:
			bamout.write(rc)
			xn += 1
	bamfile.close()
	bamout.close()
	print 'total {} reads, after is_proper_pair filter {}, after cov filter {}, after id filter {}'.format(xi,xj,xm, xn)

def pp_run(bamfile, references, lengths, processors='autodetect'):
	ppservers = ()
	job_server = pp.Server(processors, ppservers=ppservers)
	ref_len = sorted(zip(references, lengths), key=lambda x:x[1], reverse=1)
	jobs = [(ref, job_server.submit(
							filter_id, 
							(bamfile, ref, length), 
							(find_peak_border, is_covered, get_cov, get_id), 
							("import numpy as np", "pysam", 'sys', 'heapq') ) ) \
			 for ref, length in ref_len ]
	d_hq = {}
	xm = 0
	for ref, job in jobs:
		d_hq[ref] = job()
		xm += len(d_hq[ref])
		print 'contig {} candicate {} reads'.format(ref, len(d_hq[ref]))
	print 'total candicate {} reads '.format(xm)
	return d_hq
def filter_id(bamfile, ref, length, min_cov=0.9, min_as=60, wsize=100, wstep=25, bins = np.linspace(51,151,51), min_depth=3, max_depth=100):
	bamfile = pysam.AlignmentFile(bamfile, 'rb')
	xm, xn = 0,0
	hq_idsx = set([])
	f = open('/tmp/filterbam.{}.log'.format(ref), 'w')
	for start in range(0, length, wstep):
		end = start + wsize
		d_id = {rc.query_name: rc.get_tag('AS') for rc in bamfile.fetch(ref, start, end) if rc.is_proper_pair and is_covered(rc, start, end) and get_cov(rc) >= min_cov and rc.has_tag('AS')}
		if len(d_id) < min_depth:
			continue
		peak, border = find_peak_border(d_id.values(), bins=bins, max_depth=max_depth, log=f)
		if border < min_as:
			border = min_as
		hq_items = [(qid,id) for qid, id in d_id.items() if id >= border]
		if len(hq_items) > max_depth:
			hq_items = heapq.nlargest(max_depth, hq_items, key=lambda x:x[1])
		hq_ids = set(qid for qid, id in hq_items)
		print >> f, ref, start, end, len(d_id), len(hq_items), len(hq_ids), peak, border
		xm += len(hq_ids)
		xn += len(d_id)
		hq_idsx = hq_idsx | hq_ids
#	print 'total {}, get {}'.format(xn ,xm)
	bamfile.close()
	f.close()
	return hq_idsx
	
def find_peak_border(lst, bins = np.linspace(0,1,101), max_depth=100, log=sys.stderr):
	array = np.array(lst)
	hists,bins = np.histogram(array, bins = bins)
	XY = zip(hists,bins)
	XY.reverse()
	last_hist = XY[0][0]
	print >>log, list(hists), "\t", # XY #hists, bins, XY[0][0] , XY[1][0]
	depth = last_hist
	if depth > max_depth:
		pidx = 0
		peak = 1
		bidx = 0
	elif XY[0][0] > XY[1][0] and XY[0][0] > 1: # peak at 1st pos
		pidx = 0
		peak = 1
		for i, (hist, bin) in enumerate(XY[1:], 1):
			depth += hist
			if depth > max_depth and hist <= last_hist:
				bidx = i
				break
			if hist > last_hist or hist == 0: # 
				bidx = i - 1
				break
			last_hist = hist
		else:
			bidx = i-1
	else:
		peak = 0
		for i, (hist, bin) in enumerate(XY[1:-1], 1):
			depth += hist
#			if peak == 0 and hist == 1:
#				last_hist = hist
#				continue
			if depth > max_depth:
				if peak == 0:
					pidx = i
					peak = 1
					bidx = i
					break
				elif peak == 1 and hist <= last_hist:
					bidx = i
					break
			if peak == 0 and hist > 0 and hist >= last_hist and hist > XY[i+1][0]: # peak
				pidx = i
				peak = 1
			elif peak == 1 and (hist > last_hist or hist == 0):
				bidx = i -1
				break
			last_hist = hist
		else:
			bidx = i-1
	if peak == 0:
		pidx, bidx = -1, -1
	return XY[pidx][1], XY[bidx][1]
		
def get_id(rc):
	match = sum([j for i,j in rc.cigar if i == 0])
	mismatch = rc.get_tag('NM')
	rlen = rc.query_length
	id = 1.0*(match - mismatch) / rlen
	return id
def get_cov(rc):
	rlen = rc.query_length
	hit = rlen - sum([j for i,j in rc.cigar if i in set([4,5]) ])
	cov = 1.0* hit / rlen
	return cov
def is_covered(rc, start, end):
	if rc.reference_start <= start and rc.reference_end >= end:
		return True
	else:
		return False

if __name__ == '__main__':
#	main()
	main_pp()
