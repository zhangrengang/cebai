import sys
from Bio import SeqIO
def seq2dict(inSeq):
	return dict([(rc.id, rc) for rc in SeqIO.parse(inSeq, 'fasta')])
def main(inSeq=sys.stdin, inFai=sys.argv[1], outSeq=sys.stdout, fill=True):
	d_seqs = seq2dict(inSeq)
	for line in open(inFai):
		temp = line.strip().split('\t')
		sid = temp[0]
		if not sid in d_seqs:
			if not fill:
				print >> sys.stderr, 'Warn: {} not found'.format(sid)
			else:
				length = int(temp[1])
				from Bio.Seq import Seq
				from Bio.SeqRecord import SeqRecord
				seq = ''.join(['N'*length])
				rc = SeqRecord(id = sid, seq = Seq(seq))
				SeqIO.write(rc, outSeq, 'fasta')
				print >> sys.stderr, 'Warn: {} not found, using N'.format(sid)
			continue
		rc = d_seqs[sid]
		SeqIO.write(rc, outSeq, 'fasta')

if __name__ == "__main__":
	main()
