from snppy.range import *
import itertools
import math
import sys

class RTreeNode(object):
  __slots__ = (
    'range', 'children', 'data'
  )

  def _search(self, range, out):
    if not self.range.overlaps(range):
      return

    for d in self.data:
      if d.range.overlaps(range):
        out.append(d)

    for c in self.children:
      c._search(range, out)

  def search(self, range):
    out = []
    self._search(range, out)
    return out

  def __init__(self, children = [], data = []):
    self.children = children
    self.data = data

    range_type = self.children[0].range.__class__ if len(self.children) else self.data[0].range.__class__

    self.range = range_type.bounds([ obj.range for obj in itertools.chain(children, data) ])

  @classmethod
  def makeLeafNode(cls, obs):
    return cls(data = obs)

  @classmethod
  def makeInternalNode(cls, obs):
    return cls(children = obs)

  @classmethod
  def construct(cls, data, leaf_size, internal_size):
    ndim = len(data[0].range.extents)
    out = []
    cls.makeNodes(data, ndim, 0, 0, leaf_size, cls.makeLeafNode, out)
    while len(out) > 1:
      next = []
      cls.makeNodes(out, ndim, 0, 0, internal_size, cls.makeInternalNode, next)
      out = next
    assert len(out) == 1
    return out[0]

  @classmethod
  def makeNodes(cls, data, ndim, dim_num, dim_mask, child_size, ctor, out):
    N = len(data)

    r_best = float(N + 1)
    dim = ndim

    for i in range(ndim):
      if dim_mask & 1 << i:
        continue

      dmin = data[0].range.extents[i][0]
      dmax = data[0].range.extents[i][1]
      dsum = 0.0

      for ob in data:
        dmin = min(dmin, ob.range.extents[i][0])
        dmax = max(dmax, ob.range.extents[i][1])
        dsum = dsum + ob.range.extents[i][1] - ob.range.extents[i][0]

      if dmax > dmin:
        r = dsum / float(dmax - dmin)
      else:
        r = 0.0

      if r < r_best:
        dim = i
        r_best = r

    assert dim < ndim

    P = int(math.ceil(N / float(child_size)))
    n_parts = int(P ** (1.0 / float(ndim - dim_num)))

    data.sort(key = lambda x: (x.range.extents[dim][0] + x.range.extents[dim][1]) / 2.0)

    if n_parts == 1 or dim_num == ndim - 1:
      s = e = 0
      for i in range(P):
        e = N * (i+1) // P;
        out.append(ctor(data[s:e]))
        s = e
    else:
      s = e = 0
      for i in range(n_parts):
        e = N * (i+1) // n_parts;
        cls.makeNodes(data[s:e], ndim, dim_num + 1, dim_mask | (1 << dim), child_size, ctor, out)
        s = e


import unittest
class TestRTree(unittest.TestCase):
  def test_RTreeConstruction(self):
    class thing(object):
      def __init__(self, range):
        self.range = range
    things = [
      thing(Range(((0,1),(0,1)))),
      thing(Range(((2,3),(0,1)))),
      thing(Range(((0,1),(2,3)))),
      thing(Range(((2,3),(2,3))))
    ]
    tree = RTreeNode.construct(things, 2, 2)

    out = tree.search(Range(((0,3),(0,3))))
    assert len(out) == 4

    out = tree.search(Range(((-1,0),(-1,0))))
    assert len(out) == 1

    out = tree.search(Range(((-2,-1),(-2,-1))))
    assert len(out) == 0
