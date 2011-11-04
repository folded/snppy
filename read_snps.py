from snppy.snp import *
from snppy.gtf import *
import sys
import os
import glob
import collections

files = glob.glob('/u/leucegene/hiseq/*_Samples/*/transcriptome/Parsed_*/chr*.fa/snps.txt')

sample_snps = collections.defaultdict(lambda: collections.defaultdict(list))

annotation = GTFFile.read('/share/apps/CASAVA-genomes/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2011-01-27-18-25-49/Genes/genes.gtf')

for f in files:
  sample = f.split('/')[5]
  chr_snps = sample_snps[sample]
  print 'reading', f
  for snp in readSNPs_CASAVA(open(f, 'rbU')):
    chr_snps[snp.chromosome].append(snp)
