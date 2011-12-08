from snppy.range import *
import re
import os



class GFFRecord(object):
  __slots__ = (
    'seq_id', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'group'
  )

  @property
  def range(self):
    return Range((self.start, self.end))

  def length(self):
    return self.end - self.start

  def __init__(self, line):
    super(GFFRecord, self).__init__()
    line = line.rstrip('\n').split('\t', 8)
    self.seq_id = line[0]
    self.source = line[1]
    self.feature = line[2]
    self.start = int(line[3]) - 1
    self.end = int(line[4])
    self.score = None if line[5] == '.' else float(line[5])
    self.strand = line[6]
    self.frame = None if line[7] == '.' else int(line[7])
    self.group = line[8]



class GTFRecord(GFFRecord):
  __slots__ = (
    'attrs', 'gene_id', 'transcript_id'
  )

  def __init__(self, line):
    super(GTFRecord, self).__init__(line)
    self.attrs = [ attr.split(None, 1) for attr in self.group.split('; ') if len(attr.strip()) ]
    self.gene_id = self.attrs[0][1].strip('"')
    self.transcript_id = self.attrs[1][1].strip('"')



def readrows(inf, rectype):
  if type(inf) in (str, unicode):
    ext = os.path.splitext(inf)[1].lower()
    if ext == '.gz':
      import gzip
      inf = gzip.open(inf)
    elif ext == '.bz2':
      import bz2
      inf = bz2.BZ2File(inf)
    else:
      inf = open(inf)

  while 1:
    line = inf.readline()
    if line == '':
      break
    yield rectype(line)



class GFFFile(object):
  @classmethod
  def read(cls, inf):
    return readrows(inf, GFFRecord)



class GTFFile(object):
  @classmethod
  def read(cls, inf):
    return readrows(inf, GTFRecord)



__all__ = [
  'GFFRecord',
  'GTFRecord',
  'GFFFile',
  'GTFFile'
]
