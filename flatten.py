#! env python
from Bio import SeqIO
import sys
if __name__ == '__main__':
    infile=open(sys.argv[1])
    outdir=sys.argv[2]
    infas=list(SeqIO.parse(infile,'fasta'))
    def transform(s):
        return s.strip().split('|')[-1]
    for f in infas:
        SeqIO.write(f, outdir + '/'+ transform(f.id) +'.fasta', 'fasta')
    print 'done!'
