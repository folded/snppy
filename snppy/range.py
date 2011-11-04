class Range1D(object):
  self.__slots__ = ( 'lo', 'hi' )

  def __init__(self, lo, hi):
    self.lo = lo
    self.hi = hi

  def overlaps(self, other):
    return not (self.hi < other.lo or other.hi < self.lo)

  def intersection(self, other):
    lo = max(self.lo, other.lo)
    hi = max(lo, min(self.hi, other.hi))
    return Range1D(lo, hi)

  @classmethod
  def bounds(cls, ranges):
    return Range1D(
      min([ r.lo for r in ranges ]),
      max([ r.hi for r in ranges ])
    )



class Range2D(object):
  self.__slots__ = ( 'lo', 'hi' )

  def __init__(self, lo, hi):
    self.lo = lo
    self.hi = hi

  def overlaps(self, other):
    return not (self.hi[0] < other.lo[0] or
		self.hi[1] < other.lo[1] or
		other.hi[0] < self.lo[0] or
		other.hi[1] < self.lo[1])

  def intersection(self, other):
    lo0 = max(self.lo[0], other.lo[0])
    lo1 = max(self.lo[1], other.lo[1])
    hi0 = max(lo0, min(self.hi[0], other.hi[0]))
    hi1 = max(lo1, min(self.hi[1], other.hi[1]))
    return Range2D((lo0, lo1), (hi0, hi1))

  @classmethod
  def bounds(cls, ranges):
    return Range2D(
      (min([ r.lo[0] for r in ranges ]), min([ r.lo[1] for r in ranges ])),
      (max([ r.hi[0] for r in ranges ]), min([ r.hi[1] for r in ranges ]))
    )

def bounds(ranges):
  return ranges[0].bounds(ranges)

__all__ = (
  'Range1D',
  'Range2D',
  'bounds'
)
