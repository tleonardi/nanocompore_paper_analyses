import pysam
import sys
from collections import Counter
bamfile = sys.argv[1]
thr = int(sys.argv[2])
outfile = sys.argv[3]
c=Counter()
with pysam.AlignmentFile(bamfile, check_sq=False, add_sam_header=True) as bam:
    for line in bam:
        c[line.reference_name]+=1

#keep = [n for n,c in c.items() if c>=thr]

with pysam.AlignmentFile(bamfile, check_sq=False, add_sam_header=True) as bam:
    with pysam.AlignmentFile(outfile, 'wb', template=bam) as outbam:
        for line in bam:
            if(not line.is_unmapped and c[line.reference_name]>thr):
                outbam.write(line)
 



