
import os
import sys
import numpy as np
from xopen import xopen as open

def _parse_busco_tsv(inTsv):
	d_pos = {}
	for line in open(inTsv):
		if line.startswith('#'):
			continue
		temp = line.rstrip().split('\t')
		if not temp[1] in set(['Complete', 'Duplicated']):
			continue
		Busco_id,Status,Contig,Start,End,Score,Length = temp
		Start,End = int(Start), int(End)
		for pos in range(Start, End+1):
			pos = (Contig,pos)
			try: d_pos[pos].add(Status)
			except KeyError: d_pos[pos] = set([Status])
	for pos, stats in d_pos.items():
		if len(stats) > 1:
			del d_pos[pos]
		else:
			d_pos[pos] = list(stats)[0]
	return d_pos
	
def _maxdepth(d_depth, maxdepth):
    depths = []
    for sample, d_count in d_depth.items():
        for (status, depth), freq in d_count.items():
            depths += [depth] * freq
#    mean = np.mean(depths)
 #   std = np.std(depths)
    median = np.median(depths)
    return max(int(median)*2.8, maxdepth)
#    return max(int(mean*2.5), maxdepth)
def _maxdepth2(d_depth, maxdepth):
    depths = []
    for sample, d_count in d_depth.items():
        dfps = d_count.items()
        for i in range(3, len(dfps)-3):
            if dfps[i][1] > max([dfps[i-1][1], dfps[i-2][1], dfps[i+1][1], dfps[i+2][1]]):
                depths += [dfps[i]]
    return max(max(depths, key=lambda x:x[1])[0][1] * 2.8, maxdepth)
	
def main(inDepth=sys.argv[1], prefix=sys.argv[2], inTsv=sys.argv[3], plotfmt='png', maxdepth=200):
	d_target = _parse_busco_tsv(inTsv)
	d_depth = {}
#	i = 0
	data_file = 'depth.%s.busco.data' % (prefix, )
	for line in open(inDepth):
#		i += 1
		temp = line.rstrip().split('\t')
		chr, pos = temp[:2]
		pos = int(pos)
		if not (chr, pos) in d_target:
			continue
#		if i == 1:
		status = d_target[(chr, pos)]
		nsamples = len(temp)-2
		for j in range(nsamples):
			sample = 'X%s' % (j+1,)
			depth = int(temp[j+2])
#			if depth > maxdepth:
#				depth = maxdepth
			try: d_depth[sample][(status,depth)] += 1
			except KeyError: 
				try: d_depth[sample][(status,depth)] = 1
				except KeyError: d_depth[sample] = {(status,depth): 1}
	maxdepth = _maxdepth2(d_depth, maxdepth)
	print >>sys.stderr, 'Max. Depth: %s' % (maxdepth,)
	for sample, d_count in d_depth.items():
		for (status, depth), freq in d_count.items():
			if depth > maxdepth:
				try: d_depth[sample][(status, maxdepth)] += freq
				except KeyError: d_depth[sample][(status, maxdepth)] = freq
				del d_depth[sample][(status, depth)]	
	d_status = {'Complete':'Single-copy core genes', 'Duplicated': 'Duplicated core genes'}
	f = open(data_file, 'w')
	line = ['sample', 'depth', 'freq', 'status']
	print >>f, '\t'.join(line)
	f2 = open('%s.stats' % (prefix, ), 'w')
	print >>f2, '\t'.join(['sample', 'mindepth','nbases', 'precent','mapped', 'avgdepth', 'status'])
	for sample, d_count in d_depth.items():
		for (status,depth), freq in sorted(d_count.items(), key=lambda x: x[0]):
			line = [sample, depth, freq, d_status[status] ]
			line = map(str, line)
			print >>f, '\t'.join(line)
		for xstatus in ['Complete', 'Duplicated']:
			d_count2 = dict([(depth, freq) for (status,depth), freq in d_count.items() if status == xstatus])
			total = sum(d_count2.values())
#			mapped_bases = sum(map(int, d_count.keys()))
#			avgdepth = 1.0*mapped_bases/total
			for mindepth in [1, 5, 10, 20]:
				bases = [freq for depth, freq in d_count2.items() if depth >= mindepth]
				nbases = sum(bases)
				mapped_bases = sum([depth*freq for depth, freq in d_count2.items() if depth >= mindepth])
				avgdepth = 1.0*mapped_bases/total
				percent = 100.0*nbases/ total
				line = [sample,mindepth,nbases,percent, mapped_bases,avgdepth, xstatus]
				line = map(str, line)
				print >>f2, '\t'.join(line)
	f2.close()
	f.close()
	outPlot = '%s.%s' % (prefix, plotfmt)
	r_src = '''library(ggplot2)
argv = commandArgs(TRUE)
inFile = argv[1]
outFile = argv[2]
data <- read.table(inFile, header=T, sep="\t")
extmax <- function(x){
	lower <- quantile(x, 0.7)
    maxs <- vector()
    for (i in 3:(length(x)-2)){
        if (x[i] > x[i-1] & x[i] > x[i+1] & x[i] > x[i-2] & x[i] > x[i+2] & x[i] > lower){
            maxs <- c(maxs, i)
        }
    }
    return(maxs)
}
xmax <- extmax(data$freq)
ymax <- data$freq[xmax]
#print(ymax)
ylim <- max(ymax)*1.2
xmax <- data$depth[xmax]
data2 <- data.frame(xmax=xmax,ymax=ymax)
xlim <- %s
p <- ggplot(data, aes(x=depth, y=freq, color=status)) + geom_line() + xlim(0,xlim) + geom_vline(data=data2,aes(xintercept=xmax), lty=2,col='grey') + geom_text(data=data2, aes(x=xmax, y=ymax, label=xmax), col='grey', hjust=-0.25,vjust=-0.25)
ggsave(outFile, p)
''' % (maxdepth,)
	r_src_file = 'depth.%s.busco.R' % (prefix, )
	f = open(r_src_file, 'w')
	f.write(r_src)
	f.close()
	run_cmd = 'Rscript %s %s %s' % (r_src_file, data_file, outPlot)
	os.system(run_cmd)
if __name__ == '__main__':
	main()
