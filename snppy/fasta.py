def _blockReadFasta(inf, block_size = 64 * 1024):
  curr = ('', 0, [])
  n_read = 0

  for line in inf:
    line = line.rstrip()
    if line.startswith('>'):
      if len(curr[2]):
        blk = ''.join(curr[2])
        yield curr[0], curr[1], blk
      curr = (line[1:], 0, [])
      n_read = 0
    else:
      curr[2].append(line)
      n_read += len(line)
      while n_read >= block_size:
        blk = ''.join(curr[2])
        assert len(blk) == n_read
        yield curr[0], curr[1], blk[:block_size]
        curr = (curr[0], curr[1] + block_size, [blk[block_size:]])
        n_read -= block_size
  if len(curr[2]):
    blk = ''.join(curr[2])
    yield curr[0], curr[1], blk



class BlockQueue(object):
  def __init__(self, flank5, flank3):
    self.q5 = []
    self.c = None
    self.q3 = []

    self.flank5 = flank5
    self.flank3 = flank3
    self.n5 = 0
    self.n3 = 0
    self.seq_id = None

  def currID(self):
    return self.seq_id

  def empty(self):
    return self.seq_id is None

  def push(self, block):
    if block is not None:
      self.seq_id = block[0]
      self.q3.append(block)
      self.n3 += len(block[2])

    while (block is None and len(self.q3)) or (len(self.q3) and self.n3 - len(self.q3[0][2]) >= self.flank3):
      if self.c is not None:
        self.q5.append(self.c)
        self.n5 += len(self.c[2])
      self.c = self.q3[0]
      del self.q3[0]
      self.n3 -= len(self.c[2])
      while len(self.q5) and self.n5 - len(self.q5[0][2]) >= self.flank5:
        self.n5 -= len(self.q5[0][2])
        del self.q5[0]

      f5 = ''.join([ x[2] for x in self.q5 ])[-self.flank5:]
      f3 = ''.join([ x[2] for x in self.q3 ])[:self.flank3]
      yield self.c[0], self.c[1], f5, self.c[2], f3



def blockReadFasta(inf, block_size = 64 * 1024, flank5 = 0, flank3 = 0):
  if flank5 == 0 and flank3 == 0:
    for block in _blockReadFasta(inf, block_size):
      yield block[0], block[1], '', block[2], ''
    return

  bq = BlockQueue(flank5, flank3)

  for block in _blockReadFasta(inf, block_size):
    if not bq.empty() and bq.currID() != block[0]:
      for x in bq.push(None): yield x
      bq = BlockQueue(flank5, flank3)
    for x in bq.push(block): yield x
  for x in bq.push(None): yield x
 


__all__ = [
  'blockReadFasta'
]



import unittest

class FastaTest(unittest.TestCase):
  def testBlockReading(self):
    import cStringIO
    test = cStringIO.StringIO('>foo\nACGTGTACACATGTACACACA\n>foo2\nACGTGTACACATGTACACACA')

    for block in blockReadFasta(test, block_size = 5, flank5 = 2, flank3 = 3):
      print >>sys.stderr, block
