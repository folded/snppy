import range
import iter
import itertools



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



def genesFromGTF(gtf_src):
  return [
    Gene.fromGTF(rows)
    for rows in iter.groupBy(itertools.ifilter(lambda r: r.feature in ('CDS', 'exon'), gtf_src), 'gene_id')
  ]