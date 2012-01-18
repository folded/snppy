from snppy import range
from snppy import util
import re


class GFFRecord(object):
  __slots__ = (
    'seq_id', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'group', 'attrs'
  )

  group_re = re.compile(
    r'\s*([^;\s]+)\s+((?:[^";\s]+)|(?:"(?:[^"\\]|\\.)*)")\s*(?:;|$)'
  )

  @classmethod
  def _unquote(cls, v):
    if v.startswith('"'):
      return re.sub(r'\\(.)', lambda m: m.group(1), v[1:-1])
    else:
      return v

  @classmethod
  def _quote(cls, v):
    v = str(v)
    if re.search(r'[\s";]', v) is not None:
      v = '"' + v.replace('\\', '\\\\').replace('"', '\\"') + '"'
    return v

  @classmethod
  def fmtattrs(cls, **kw):
    return '; '.join([k + ' ' + cls._quote(v) for k, v in kw.iteritems()])

  @property
  def range(self):
    return range.Range((self.start, self.end))

  def length(self):
    return self.end - self.start

  def asStr(self, one_based = True, end_included = False):
    return '\t'.join([
      self.seq_id,
      self.source,
      self.feature,
      str(self.start + (1 if one_based else 0)),
      str(self.end + (1 if end_included else 0)),
      '.' if self.score is None else str(self.score),
      self.strand,
      '.' if self.frame is None else str(self.frame),
      self.fmtattrs(**self.attrs)])

  @classmethod
  def parse(cls, line, one_based = True, end_included = False):
    ob = cls()
    line = line.rstrip('\n').split('\t', 8)
    ob.seq_id = line[0]
    ob.source = line[1]
    ob.feature = line[2]
    ob.start = int(line[3])
    ob.end = int(line[4])
    if one_based: ob.start -= 1
    if end_included: ob.end -= 1
    ob.score = None if line[5] == '.' else float(line[5])
    ob.strand = line[6]
    ob.frame = None if line[7] == '.' else int(line[7])
    ob.group = line[8]
    ob.attrs.update([(k, cls._unquote(v)) for k, v in cls.group_re.findall(ob.group)])
    return ob

  def __init__(self):
    super(GFFRecord, self).__init__()
    self.attrs = {}
    self.seq_id = None
    self.source = None
    self.feature = None
    self.start = None
    self.end = None
    self.score = None
    self.strand = None
    self.frame = None
    self.attrs = {}
    self.group = None



class GTFRecord(GFFRecord):
  __slots__ = (
    'gene_id', 'transcript_id'
  )

  @classmethod
  def parse(cls, line, **kw):
    ob = super(GTFRecord, cls).parse(line, **kw)
    ob.gene_id = ob.attrs['gene_id']
    ob.transcript_id = ob.attrs['transcript_id']
    return ob

  def __init__(self):
    super(GTFRecord, self).__init__()



def readrows(inf, rectype, **kw):
  if type(inf) in (str, unicode):
    inf = util.open(inf, 'rbU')

  while 1:
    line = inf.readline()
    if line == '':
      break
    if line.startswith('#'):
      continue
    yield rectype.parse(line, **kw)



class GFFFile(object):
  @classmethod
  def read(cls, inf, **kw):
    return readrows(inf, GFFRecord, **kw)



class GTFFile(object):
  @classmethod
  def read(cls, inf, **kw):
    return readrows(inf, GTFRecord, **kw)


__all__ = [
  'GFFRecord',
  'GTFRecord',
  'GFFFile',
  'GTFFile'
]
