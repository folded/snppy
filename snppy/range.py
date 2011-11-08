class Range(object):
  __slots__ = ( 'extents', )

  def __str__(self):
    return str(tuple([ _[0] for _ in self.extents])) + \
        '-' + \
        str(tuple([ _[1] for _ in self.extents]))

  def __init__(self, extents):
    self.extents = extents

  def overlaps(self, other):
    for a, b in zip(self.extents, other.extents):
      if a[1] < b[0] or b[1] < a[0]:
        return False
    return True

  def intersection(self, other):
    def _intersection(a, b):
      lo = max(a[0], b[0])
      hi = max(lo, min(a[1], b[1]))
      return (lo, hi)

    return Range([ _intersection(a, b) for a, b in zip(self.extents, other.extents) ])

  @classmethod
  def bounds(cls, ranges):
    def _bounds(ranges):
      return (min([ r[0] for r in ranges ]),
              max([ r[1] for r in ranges ]))
    return Range([ _bounds(ranges) for ranges in zip(*[ range.extents for range in ranges ]) ])



def bounds(ranges):
  return ranges[0].bounds(ranges)



__all__ = (
  'Range',
  'bounds'
)
