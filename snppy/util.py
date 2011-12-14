import os
import gzip
import bz2

_open = open

def open(inf, mode = 'r'):
  ext = os.path.splitext(inf)[1].lower()
  if ext == '.gz':
    inf = gzip.open(inf, mode)
  elif ext == '.bz2':
    inf = bz2.BZ2File(inf, mode)
  else:
    inf = _open(inf, mode)

  return inf

__all__ = [
  'open'
]
