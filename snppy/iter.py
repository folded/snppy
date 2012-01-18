import collections



def groupBy(iter, grouper):
  if type(grouper) is str:
    attr = grouper
    grouper = lambda obj: getattr(obj, attr)
  groups = collections.defaultdict(list)
  for obj in iter:
    groups[grouper(obj)].append(obj)
  return groups.values()
