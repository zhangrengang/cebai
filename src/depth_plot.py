
import os
import sys
import numpy as np
from xopen import xopen as open

def main(inDepth=sys.argv[1], prefix=sys.argv[2], plotfmt='png', maxdepth=200):
	d_depth = {}
#	i = 0
	data_file = '%s.data' % (prefix, )
	for line in open(inDepth):
#		i += 1
		temp = line.rstrip().split('\t')
#		if i == 1:
		nsamples = len(temp)-2
		for j in range(nsamples):
			sample = 'X%s' % (j+1,)
			depth = int(temp[j+2])
#			if depth > maxdepth:
#				depth = maxdepth
			try: d_depth[sample][depth] += 1
			except KeyError: 
				try: d_depth[sample][depth] = 1
				except KeyError: d_depth[sample] = {depth: 1}
	def _maxdepth(d_depth, maxdepth):
		depths = []
		for sample, d_count in d_depth.items():
			for depth, freq in d_count.items():
				depths += [depth] * freq
#		mean = np.mean(depths)
#		std = np.std(depths)
		median = np.median(depths)
		return max(int(median)*2.8, maxdepth)
#		return max(int(mean+2*std), maxdepth)
	def _maxdepth2(d_depth, maxdepth):
		depths = []
		for sample, d_count in d_depth.items():
			dfps = d_count.items()
			for i in range(3, len(dfps)-3):
				if dfps[i][1] > max([dfps[i-1][1], dfps[i-2][1], dfps[i+1][1], dfps[i+2][1]]):
					depths += [dfps[i]]
		return max(max(depths, key=lambda x:x[1])[0] * 2.8, maxdepth)
	maxdepth = _maxdepth2(d_depth, maxdepth)
	print >>sys.stderr, 'Max. Depth: %s' % (maxdepth,)
	for sample, d_count in d_depth.items():
		for depth, freq in d_count.items():
			if depth > maxdepth:
				try: d_depth[sample][maxdepth] += freq
				except KeyError: d_depth[sample][maxdepth] = freq
				del d_depth[sample][depth]
	f = open(data_file, 'w')
	line = ['sample', 'depth', 'freq']
	print >>f, '\t'.join(line)
	f2 = open('%s.stats' % (prefix, ), 'w')
	print >>f2, '\t'.join(['sample', 'mindepth','nbases', 'precent','mapped', 'avgdepth'])
	for sample, d_count in d_depth.items():
		for depth, freq in sorted(d_count.items(), key=lambda x: x[0]):
			line = [sample, depth, freq]
			line = map(str, line)
			print >>f, '\t'.join(line)
		total = sum(d_count.values())
#		mapped_bases = sum(map(int, d_count.keys()))
#		avgdepth = 1.0*mapped_bases/total
		for mindepth in [1, 2, 5, 10, 20]:
			bases = [freq for depth, freq in d_count.items() if depth >= mindepth]
			nbases = sum(bases)
			mapped_bases = sum([depth*freq for depth, freq in d_count.items() if depth >= mindepth])
			avgdepth = 1.0*mapped_bases/total
			percent = 100.0*nbases/ total
			line = [sample,mindepth,nbases,percent, mapped_bases,avgdepth]
			line = map(str, line)
			print >>f2, '\t'.join(line)
	f2.close()
	f.close()
	outPlot = '%s.%s' % (prefix, plotfmt)
	r_src = '''library(ggplot2)
argv = commandArgs(TRUE)
inFile = argv[1]
outFile = argv[2]
data <- read.table(inFile, header=T)
extmax <- function(x){
	lower <- quantile(x, 0.8)
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
data2 <- data.frame(xmax=xmax,ymax=ymax)
cov <- xmax[which(ymax==max(ymax))]
peak <- ymax[which(ymax==max(ymax))]
#print(c(cov,peak))
print(cov)
xmax <- data$depth[xmax]
data2 <- data.frame(xmax=xmax,ymax=ymax)
pois <- function(x, cov){
    dpois(x, cov, log = FALSE)
}
pois2 <- function(x, cov){
    pois(x, cov)*data$freq[cov]/pois(x, cov)[cov]
}
x <- c(1:%s)
#print(pois(x, 20))
pois.prob <- pois(x, as.integer(cov/2))+pois(x, cov)+pois(x, cov*2)+pois(x, cov*3)+pois(x, cov*4)
#poisson <- pois.prob*(peak/(pois.prob[cov]))
poisson <- pois2(x, as.integer(cov/2))+pois2(x, cov) + pois2(x, cov*2)
data3 <- data.frame(x=x, poisson=poisson)

xlim <- %s
if (length(unique(data$sample)) == 1) {
	p <- ggplot(data, aes(x=depth, y=freq)) + geom_line(color='blue')
}else{
	p <- ggplot(data, aes(x=depth, y=freq, color=sample)) + geom_line()
}
p <- p + xlim(0,xlim) + geom_vline(data=data2,aes(xintercept=xmax), lty=2,col='grey') + geom_text(data=data2, aes(x=xmax, y=ymax, label=xmax), col='grey', hjust=-0.25,vjust=-0.25) + geom_line(data=data3, aes(x=x, y=poisson), col = 'red', lty=2)
ggsave(outFile, p)
''' % (maxdepth, maxdepth)
	r_src_file = '%s.R' % (prefix, )
	f = open(r_src_file, 'w')
	f.write(r_src)
	f.close()
	run_cmd = 'Rscript %s %s %s' % (r_src_file, data_file, outPlot,)
	os.system(run_cmd)
if __name__ == '__main__':
	main()
