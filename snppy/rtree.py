from snppy.range import *
import itertools

class RTreeNode(object):
  __slots__ = (
    'range', 'children', 'data'
  )

  def search(self, range, out):
    if not self.range.overlaps(range): return
    out.extend(data)
    for c in self.children:
      c.search(range, out)

  def __init__(self, children = [], data = []):
    self.children = children
    self.data = data

    range_type = self.chilren[0].range.__class__ if len(self.children) else self.data[0].range.__class__

    self.range = range_type.bounds([ obj.range for obj in itertools.chain(children, data) ])

  @classmethod
  def construct(self, data):
    pass

