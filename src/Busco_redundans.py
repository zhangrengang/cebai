import sys
import glob
#import numpypy
import numpy as np
from xopen import xopen as open

class BUSCO():
	def __init__(self, buscoDir):
		self.fullTable = glob.glob('{}/full_table_*.tsv'.format(buscoDir))[0]
	def full_table(self):
		print >>sys.stderr, 'reading', self.fullTable
		for line in open(self.fullTable):
			if line.startswith('#'):
				continue
			yield FullTableRecord(line)
	def fullTableRegions(self):
		return [rc.region for rc in self.full_table() ]
class FullTableRecord():
	def __init__(self, line):
		temp = line.rstrip().split('\t')
		self.line = temp
		try:
			Busco_id,Status,Contig,Start,End,Score,Length = temp
			Start,End,Length = map(int, (Start,End,Length))
			Score = float(Score)
		except ValueError: # Misiing
			Busco_id,Status = temp
			Contig,Start,End,Score,Length = None, None, None, None, None
		self.busco_id = Busco_id
		self.Status = Status
		self.contig = Contig
		self.start = Start
		self.end = End
		self.score = Score
		self.length = Length
		self.region = (self.contig, self.start, self.end)
	def write(self, f):
		print >>f, '\t'.join(self.line)

class Redundans():
	def __init__(self, tsv):
		self.tsv = tsv
	def __iter__(self):
		return self.parse_tsv()
	def parse_tsv(self):
		print >>sys.stderr, 'reading', self.tsv
		for line in open(self.tsv):
			if line.startswith('#'):
				continue
			yield RedundansRecord(line)
	@property
	def dict(self):
		d = {}
		for rc in self.parse_tsv():
			try: d[rc.contig] += [rc]
			except KeyError: d[rc.contig] = [rc]
		return d
class RedundansRecord():
	def __init__(self, line):
		temp = line.rstrip().split('\t')
		self.line = temp
		contig, size, target, itentity, overlap = temp
		self.contig = contig
		self.size = size
		self.target = target
		self.itentity = itentity
		self.overlap = overlap
	def write(self, f):
		print >>f, '\t'.join(self.line)

class Depth():
	def __init__(self, depthfile):
		self.depthfile = depthfile
	def __iter__(self):
		return self.parse()
	def parse(self):
		print >>sys.stderr, 'reading', self.depthfile
		for line in open(self.depthfile):
			yield DepthRecord(line)
class DepthRecord():
	def __init__(self, line):
		temp = line.rstrip().split('\t')
		self.line = temp
		temp[1:] = map(int, temp[1:])
		contig, position = temp[:2]
		self.contig = contig
		self.position = position
		self.pos = (contig, position)
		self.depths = temp[2:]
		self.depth = np.sum(self.depths)
		
	def write(self, f):
		print >>f, '\t'.join(self.line)
class DepthsRecord():
	def __init__(self, depths):
		self.depths = depths
	@property
	def mean(self):
		return round(np.mean(self.depths), 1)
	@property
	def median(self):
		return np.median(self.depths)
	@property
	def std(self):
		return round(np.std(self.depths), 1)
	@property
	def max(self):
		return np.max(self.depths)
	@property
	def min(self):
		return np.min(self.depths)
	@property
	def maxfreq(self):
		from collections import Counter
		d_count = Counter(self.depths)
		return max(d_count.items(), key=lambda x:x[1])[0]
class RegionDepth():
	def __init__(self, regions, depthfile):
		self.regions = regions
		self.depthfile = depthfile
	def __iter__(self):
		return self.parse()
	def region2pos(self, region):
		contig, start, end = region
		return [(contig, pos) for pos in range(start, end+1)]
	def regions2pos(self):
		d_pos = {}
		for region in self.regions:
			if region[0] is None:
				continue
			for pos in self.region2pos(region):
				d_pos[pos] = region
		return d_pos
	def parse(self):
		d_pos = self.regions2pos()
		d_depth = {}
		for rc in Depth(self.depthfile):
			if not rc.pos in d_pos:
				continue
			region = d_pos[rc.pos]
			try: d_depth[region] += [rc.depth]
			except KeyError: d_depth[region] = [rc.depth]
		for region in self.regions:
			try: yield region, DepthsRecord(d_depth[region])
			except KeyError: continue # Missing

def combineInfo(buscoDir, depthfile, outTsv, hetero_tsv=None):
	if hetero_tsv is None:
		d_alns = {}
	else:
		d_alns = Redundans(hetero_tsv).dict
	regions = BUSCO(buscoDir).fullTableRegions()
#	print regions
	d_depth = dict((region, rc) for region, rc in RegionDepth(regions,depthfile))
#	print d_depth
	for rc in BUSCO(buscoDir).full_table():
		region = rc.region
		try: depth = d_depth[region]
		except KeyError: continue # Missing
		depth_stat = [depth.mean, depth.std, depth.median, depth.maxfreq, depth.min, depth.max]
		try: alns = ['\t'.join(aln.line) for aln in d_alns[rc.contig]]
		except KeyError: alns = []
		line = rc.line + depth_stat + alns
		line = map(str, line)
		print >> outTsv, '\t'.join(line)

def main():
	buscoDir, depthfile = sys.argv[1:3]
	try: hetero_tsv = sys.argv[3]
	except IndexError: hetero_tsv = None
	outTsv = sys.stdout
	combineInfo(buscoDir, depthfile, outTsv, hetero_tsv)
if __name__ == '__main__':
	main()
