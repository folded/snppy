import collections



def groupByConsecutive(iter, grouper):
  if type(grouper) is str:
    attr = grouper
    grouper = lambda obj: getattr(obj, attr)

  try:
    obj = iter.next()
    grp = grouper(obj)
    block = [ obj ]
  except StopIteration:
    return

  for obj in iter:
    newgrp = grouper(obj)
    if newgrp != grp:
      yield block
      block = [ obj ]
      grp = newgrp
    else:
      block.append(obj)

  if len(block):
    yield block

def groupBy(iter, grouper):
  if type(grouper) is str:
    attr = grouper
    grouper = lambda obj: getattr(obj, attr)
  groups = collections.defaultdict(list)
  for obj in iter:
    groups[grouper(obj)].append(obj)
  return groups.values()
