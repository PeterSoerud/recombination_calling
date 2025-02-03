import pysam
import sys

haploid_bam_in = sys.argv[1]
contig_name = sys.argv[2]

coverage = {}
bamfile = pysam.AlignmentFile(haploid_bam_in, "rb")

for pileupcolumn in bamfile.pileup(contig_name, min_base_quality = 1, min_mapping_quality = 1, truncate=True):
    if pileupcolumn.n in coverage:
        coverage[pileupcolumn.n] += 1
    else:
        coverage[pileupcolumn.n] = 1

for key,value in coverage.items():   
    print(contig_name, key, value, sep = '\t')