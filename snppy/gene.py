import range
import iter
import itertools
import sys



class Transcript(object):
  # __slots__ = (
  #   'cds',
  #   'exons',
  #   'name',
  #   'gene',
  #   'strand',
  #   'chromosome',
  #   'range'
  # )

  def __init__(self):
    pass

  def __cmp__(self, other):
    v = cmp(self.exons, other.exons)
    if v: return v
    return cmp(self.name, other.name)

  @classmethod
  def fromGTF(cls, gtf_rows):
    ob = cls()
    ob.cds = [ row.range for row in itertools.ifilter(lambda row: row.feature == 'CDS', gtf_rows) ]
    ob.exons = [ row.range for row in itertools.ifilter(lambda row: row.feature == 'exon', gtf_rows) ]
    ob.name = gtf_rows[0].transcript_id
    ob.gene = gtf_rows[0].gene_id
    ob.strand = gtf_rows[0].strand
    ob.chromosome = gtf_rows[0].seq_id
    ob.range = range.bounds(ob.exons)
    return ob



class Gene(object):
  # __slots__ = (
  #   'transcripts',
  #   'name',
  #   'strand',
  #   'chromosome',
  #   'range'
  # )

  def __init__(self):
    pass

  def __cmp__(self, other):
    v = cmp(self.range, other.range)
    if v: return v
    return cmp(self.name, other.name)

  @classmethod
  def fromGTF(cls, gtf_rows):
    ob = cls()
    ob.transcripts = sorted([ Transcript.fromGTF(rows) for rows in iter.groupBy(gtf_rows, 'transcript_id') ])
    ob.range = range.bounds([ transcript.range for transcript in ob.transcripts ])
    ob.chromosome = ob.transcripts[0].chromosome
    ob.strand = ob.transcripts[0].strand
    ob.name = ob.transcripts[0].gene
    return ob

  def decomposeTranscripts(self):
    a = []
    for i,t in enumerate(self.transcripts):
      for e in t.exons:
        a.append((e.extents[0][0], 0, i))
        a.append((e.extents[0][1], 1, i))
    a.sort()
    c = set()
    lp = a[0][0]
    tlen = 0
    ulen = 0
    n = 0
    for p, op, i in a:
      if p > lp:
        print lp, p, p - lp, len(c), c
        if len(c):
          n += 1
          tlen = tlen + (p - lp)
          if len(c) == 1:
            ulen = ulen + (p - lp)
        lp = p
      if op == 0:
        c.add(i)
      else:
        c.remove(i)

    print 'XXX %d\t%d\t%f\t%d\t%d\t' % (len(self.transcripts), n, float(ulen) / tlen, ulen, tlen)



def genesFromGTF(gtf_src):
  return [
    Gene.fromGTF(rows)
    for rows in iter.groupByConsecutive(itertools.ifilter(lambda r: r.feature in ('CDS', 'exon'), gtf_src), 'gene_id')
  ]
