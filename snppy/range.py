class Range(object):
  __slots__ = ( 'extents', )

  def __str__(self):
    a, b = zip(*self.extents)
    if len(a) == 1: a, b = a[0], b[0]
    return repr(a) + '-' + repr(b)

  def __repr__(self):
    a, b = zip(*self.extents)
    if len(a) == 1: a, b = a[0], b[0]
    return '{range:' + repr(a) + '-' + repr(b) + '}'

  def __cmp__(self, other):
    for a, b in zip(self.extents, other.extents):
      v = cmp(a, b)
      if v: return v
    return 0

  def __init__(self, *extents):
    self.extents = extents

  def overlaps(self, other):
    for a, b in zip(self.extents, other.extents):
      if a[1] < b[0] or b[1] < a[0]:
        return False
    return True

  @property
  def midpoint(self):
    return tuple([ (a[0] + a[1]) / 2.0 for a in self.extents ])

  @property
  def size(self):
    return tuple([ a[1] - a[0] for a in self.extents ])


  def distance(self, other):
    def distance1(a, b):
      if a[0] >= b[1]: return a[0] - b[1]
      if b[0] >= a[1]: return b[0] - a[1]
      return 0
    return [ distance1(a, b) for a, b in zip(self.extents, other.extents) ]
    
  def overlaps(self, other):
    for a, b in zip(self.extents, other.extents):
      if b[1] < a[0] or b[0] > a[1]:
        return False
    return True

  def contains(self, other):
    for a, b in zip(self.extents, other.extents):
      if b[0] > a[0] or b[1] < a[1]:
        return False
    return True

  def isEmpty(self):
    for a in self.extents:
      if a[1] > a[0]: return False
    return True

  def intersection(self, other):
    def _intersection(a, b):
      lo = max(a[0], b[0])
      hi = max(lo, min(a[1], b[1]))
      return lo, hi

    return Range(*[ _intersection(a, b) for a, b in zip(self.extents, other.extents) ])

  @classmethod
  def bounds(cls, ranges):
    def _bounds(ranges):
      return (min([ r[0] for r in ranges ]),
              max([ r[1] for r in ranges ]))
    return Range(*[ _bounds(ranges) for ranges in zip(*[ range.extents for range in ranges ]) ])



def bounds(ranges):
  return ranges[0].bounds(ranges)



def union(r1, r2):
  a = []
  for r in r1: a.extend([ (r.extents[0][0], -1), (r.extents[0][1], +1) ])
  for r in r2: a.extend([ (r.extents[0][0], -1), (r.extents[0][1], +1) ])
  a.sort()
  out = []
  s = 0
  for p, d in a:
    s2 = s - d
    if (s == 0 and s2 == 1) or (s == 1 and s2 == 0):
      out.append(p)
    s = s2

  result = [ Range((x,y)) for x, y in zip(out[0::2], out[1::2]) ]

  return result



def intersection_1(r1, r2):
  i = j = 0
  out = []
  s = 0
  while i < len(r1)*2 or j < len(r2)*2:
    if (j == len(r2)*2 or
        (i != len(r1)*2 and
         (
          (r1[i>>1].extents[0][i&1] < r2[j>>1].extents[0][j&1]) or
          (r1[i>>1].extents[0][i&1] == r2[j>>1].extents[0][j&1] and not i&1)
         ))):
      s2 = s + (-1 if (i&1) else +1)
      p = r1[i>>1].extents[0][i&1]
      i += 1
    else:
      s2 = s + (-1 if (j&1) else +1)
      p = r2[j>>1].extents[0][j&1]
      j += 1
    if (s == 1 and s2 == 2) or (s == 2 and s2 == 1):
      out.append(p)
    s = s2

  result = [ Range((x,y)) for x, y in zip(out[0::2], out[1::2]) ]

  return result



def intersection_2(r1, r2):
  a = []
  for r in r1: a.extend([ (r.extents[0][0], -1), (r.extents[0][1], +1) ])
  for r in r2: a.extend([ (r.extents[0][0], -1), (r.extents[0][1], +1) ])
  a.sort()
  out = []
  s = 0
  for p, d in a:
    s2 = s - d
    if (s == 1 and s2 == 2) or (s == 2 and s2 == 1):
      out.append(p)
    s = s2

  result = [ Range((x,y)) for x, y in zip(out[0::2], out[1::2]) ]

  return result



def intersection_3(r1, r2):
  a1 = []
  a2 = []
  for r in r1: a1.extend([ (r.extents[0][0], -1), (r.extents[0][1], +1) ])
  for r in r2: a2.extend([ (r.extents[0][0], -1), (r.extents[0][1], +1) ])

  out = []
  s = 0
  i = j = 0
  while i < len(a1) or j < len(a2):
    if j == len(a2) or (i < len(a1) and a1[i] < a2[j]):
      p, d = a1[i]
      i += 1
    else:
      p, d = a2[j]
      j += 1
    s2 = s - d
    if (s == 1 and s2 == 2) or (s == 2 and s2 == 1):
      out.append(p)
    s = s2

  result = [ Range((x,y)) for x, y in zip(out[0::2], out[1::2]) ]

  return result




import unittest


class TestRange(unittest.TestCase):
  def test_union(self):
    r1 = [ Range((0,2)) ]
    r2 = [ Range((1,3)) ]
    print union(r1, r2)

    r1 = [ Range((0,2)) ]
    r2 = [ Range((2,4)) ]
    print union(r1, r2)

    r1 = [ Range((0,2)) ]
    r2 = [ Range((3,5)) ]
    print union(r1, r2)

    r1 = [ Range((0,2)), Range((3,5)) ]
    r2 = [ Range((1,4)) ]
    print union(r1, r2)

    r1 = [ Range((0,2)), Range((3,5)) ]
    r2 = [ Range((-2,-1)), Range((4,6)) ]
    print union(r1, r2)

    r1 = [ Range((-2,-1)), Range((0,2)), Range((3,6)) ]
    r2 = [ Range((-1,5)) ]
    print union(r1, r2)

intersection = intersection_3

__all__ = [
  'Range',
  'bounds',
  'union',
  'intersection'
]


if __name__ == '__main__':
  import random

  def test(r1, r2, rout):
    sa = set()
    sb = set()
    sout = set()
    for r in r1:
      sa.update(range(r.extents[0][0], r.extents[0][1]))
    for r in r2:
      sb.update(range(r.extents[0][0], r.extents[0][1]))
    for r in rout:
      sout.update(range(r.extents[0][0], r.extents[0][1]))
    assert sa & sb == sout

  for i in range(10000):
    r1 = sorted(random.sample(xrange(100000), 1000))
    r2 = sorted(random.sample(xrange(100000), 1000))
    r1 = [ Range((x,y)) for x, y in zip(r1[0::2], r1[1::2]) ]
    r2 = [ Range((x,y)) for x, y in zip(r2[0::2], r2[1::2]) ]
    rout = intersection_3(r1, r2)
    # test(r1, r2, rout)
