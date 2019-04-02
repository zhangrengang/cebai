import sys
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
def main(inSeq=sys.argv[1], inBed=sys.argv[2], outSeq=sys.stdout):
	d_bed = bed2dict(inBed)
	all_n = 0
	for rc in SeqIO.parse(inSeq, 'fasta'):
		if rc.id not in d_bed:
			SeqIO.write(rc, outSeq, 'fasta')
			continue
		regions = d_bed[rc.id]
		segments = []
		left = 0
		for start, end in sorted(regions):
			segments += [str(rc.seq[left:start])]
			segments += ['n'*(end-start)]
			left = end
		segments += [str(rc.seq[left:])]
		seq = ''.join(segments)
		d_counter = Counter(list(seq))
		num_N = d_counter['n'] + d_counter['N']
		seq = Seq(seq)
		assert len(seq) == len(rc.seq)
		rc.seq = seq
		print >> sys.stderr, rc.id, len(seq), num_N, 100.0*(len(seq)-num_N)/len(seq)
		if len(seq) == num_N:
			all_n += 1
		SeqIO.write(rc, outSeq, 'fasta')
	print >> sys.stderr, all_n, 'all N'
		
def bed2dict(inBed):
	d_bed = {}
	for line in open(inBed):
#		if line.startswith('#'):
			
		temp = line.strip().split()
		CHR, START, END = temp[:3]
		START = int(START)
		END = int(END)
		try: d_bed[CHR] += [(START, END)]
		except KeyError: d_bed[CHR] = [(START, END)]
	return d_bed

if __name__ == '__main__':
	main()
