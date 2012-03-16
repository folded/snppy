import snppy.gene
import snppy.gtf
import snppy.rtree
import snppy.range
import collections

GTF_FILE = 'Homo_sapiens.GRCh37_chr1_slim.64.gtf'
# GTF_FILE = 'empty.gtf'

genes_by_chromosome = collections.defaultdict(list)
for gene in snppy.gene.genesFromGTF(snppy.gtf.GTFFile.read(open(GTF_FILE, 'rbU'))):
  genes_by_chromosome[gene.chromosome].append(gene)

chromosome_rtrees = {}
for k, v in genes_by_chromosome.iteritems():
  chromosome_rtrees[k] = snppy.rtree.RTreeNode.construct(v, 20, 20)

def readOverlappingGenes(chromosome, extents):
  rt = chromosome_rtrees.get(chromosome)
  if rt is None:
    return frozenset()
  hits = set()
  for e in extents:
    hits.update(rt.search(e))
  hits = frozenset(hits)
  return hits

def cigarToRefRange(pos, cigar):
  out = []
  for op,l in cigar:
    if op == 0:
      out.append(snppy.range.Range((pos, pos+l)))

    if op == 0 or op == 2 or op == 3:
      pos += l
      
  return out

import pysam

def transcriptAlignmentOverlap(transcript, alignment):
  return sum([ r.size[0] for r in snppy.range.intersection(sorted(transcript.exons), alignment) ])

def processIntergenicHits(hits):
  for (lo, hi), extents, read in hits:
    print 'XXX', 0.0, read[0].qname, read[1].qname

gene_counts = collections.defaultdict(int)
gene_intron_counts = collections.defaultdict(int)
transcript_counts = collections.defaultdict(int)

def processGroupedHits(hit_genes, hits):
  if len(hit_genes) > 1:
    print '* reads hit multiple genes:', [g.name for g in hit_genes], len(hits)
  elif len(hit_genes) == 1:
    print '* reads hit single gene:', [g.name for g in hit_genes], len(hits)

  hits.sort()

  for (lo, hi), extents, read in hits:
    aln_bases = sum([ e.size[0] for e in extents ])
    best = []
    for g in hit_genes:
      for t in g.transcripts:
        overlap = transcriptAlignmentOverlap(t, extents)
        if len(best) == 0 or best[0][0] < overlap:
          best = []
        best.append((overlap, read, g, t))

    frac = 0.0
    if len(best):
      if best[0][0] == 0:
        best = []
      else:
        frac = best[0][0] / float(aln_bases)

    if len(best) == 1:
      if frac == 1.0:
        gene_counts[g.name] += 1
        transcript_counts[t.name] += 1
      else:
        # hits only one transcript of one gene, but although it
        # overlaps one or more exons, it has some sequence that aligns
        # to an intron, or outside the gene.
        print 'TRN', frac, read[0].qname, read[1].qname, '->', g.name, t.name
        print '  ', extents, [ x.size[0] for x in extents ]
        print '  ', read[0]
        print '  ', read[1]


    elif len(best) > 1:
      g = list(set([ g for o, r, g, t in best ]))
      if len(g) == 1:
        if frac == 1.0:
          gene_counts[g[0].name] += 1
        else:
          # hits many transcripts of one gene equally well, but
          # although it overlaps one or more exons, it has some
          # sequence that aligns to an intron, or outside the gene.
          print 'GEN', frac, read[0].qname, read[1].qname, '->', g[0].name, '!' if len(hit_genes) > 1 else '', [ t.name for o, r, g, t in best ]
          print '  ', extents, [ x.size[0] for x in extents ]
          print '  ', read[0]
          print '  ', read[1]
      else:
        # hits multiple genes.
        # assign to one based upon the unambiguous transcription
        # evidence for the genes/transcripts to which the read pair
        # matched.
        print 'AMB', frac, read[0].qname, read[1].qname, '->', [ x.name for x in hit_genes ]
        print '  ', extents, [ x.size[0] for x in extents ]
        print '  ', read[0]
        print '  ', read[1]
    elif len(best) == 0:
      if len (hit_genes) == 1:
        # hits one gene in an intron.
        gene_intron_counts[hit_genes[0].name] += 1
      else:
        # hits more than one gene in an intron.
        print 'INT', frac, read[0].qname, read[1].qname, '->', [ x.name for x in hit_genes ]
        print '  ', extents, [ x.size[0] for x in extents ]
        print '  ', read[0]
        print '  ', read[1]

samfile = pysam.Samfile("sorted.bam", "rb")

N = 100000
C = 0

gap_counts = collections.defaultdict(int)

while 1:
  grp = collections.defaultdict(list)

  read_name = {}

  def processSingleton(r):
    # we can only do this once we've fully processed a bam file, and found the unpaired reads.
    r_hits = readOverlappingGenes(r)

  def processPair(r1, r2):
    r1_ch, r1_ex = samfile.getrname(r1.tid), cigarToRefRange(r1.pos, r1.cigar)
    r2_ch, r2_ex = samfile.getrname(r2.tid), cigarToRefRange(r2.pos, r2.cigar)
    if r1_ch != r2_ch:
      print '*', r1.qname, 'and', r2.qname, 'map to different chromosomes'
    else:
      r1_hits = readOverlappingGenes(r1_ch, r1_ex)
      r2_hits = readOverlappingGenes(r2_ch, r2_ex)

      if len(r1_hits) and len(r2_hits):
        # if there's an intersection between the genes that r1 and r2 hits, take that as the hit gene set.
        rboth_hits = tuple(sorted(r1_hits & r2_hits))
        if len(rboth_hits):
          extents = snppy.range.union(r1_ex, r2_ex)
          grp[rboth_hits].append(((extents[0].extents[0][0], extents[-1].extents[0][1]), extents, (r1, r2)))
        else:
          print '*', r1.qname, 'and', r2.qname, 'map to distinct gene sets', [ g.name for g in r1_hits ], 'and', [ g.name for g in r2_hits ]
      else:
        r_hits = tuple(sorted(r1_hits | r2_hits))
        extents = snppy.range.union(r1_ex, r2_ex)
        grp[r_hits].append(((extents[0].extents[0][0], extents[-1].extents[0][1]), extents, (r1, r2)))

  for n in xrange(N):
    read = samfile.next()
    C = C + 1

    if read.qname in read_name:
      ra = read_name[read.qname]
      rb = read
      processPair(ra, rb)

      del read_name[read.qname]
    else:
      read_name[read.qname] = read

    ex = cigarToRefRange(read.pos, read.cigar)
    for j in range(1, len(ex)):
      junction = samfile.getrname(read.tid), ex[j-1].extents[0][1], ex[j].extents[0][0]
      gap_counts[junction] += 1

  print '* reads processed:', C
  print '* reads remaining unpaired:', len(read_name)

  if frozenset() in grp:
    print '* intergenic read pairs', len(grp[frozenset()])
    processIntergenicHits(grp.pop(frozenset()))

  n_unique = 0
  n_multiple = 0
  for k, v in grp.iteritems():
    if len(k) > 1:
      n_multiple += len(v)
    else:
      n_unique += len(v)
    processGroupedHits(k, v)

  print '* read pairs mapping to a single gene', n_unique
  print '* read pairs mapping to multiple genes', n_multiple

for k, v in sorted(gap_counts.items(), key = lambda k: -k[1]):
  print '%d\t%s' % (v, k)

samfile.close()
