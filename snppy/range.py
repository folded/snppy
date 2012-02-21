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

  def size(self):
    return [ a[1] - a[0] for a in self.extents ]

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



__all__ = (
  'Range',
  'bounds'
)
