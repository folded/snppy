class SNP(object):
  __slots__ = (
    'chromosome', 'pos', 'depth', 'quality', 'refbase', 'genotype', 'n_acgt'
  )

  def __init__(self,
               chromosome = 'chrUn',
               pos = 0,
               depth = 0,
               quality = 0,
               refbase = 'N',
               genotype = 'NN',
               n_acgt = (0,0,0,0)):
    self.chromosome = chromosome
    self.pos = pos
    self.depth = depth
    self.quality = quality
    self.refbase = refbase
    self.genotype = genotype
    self.n_acgt = n_acgt




def readSNPs_CASAVA(inf):
  for row in inf:
    if row[0] == '#':
      continue
    row = row.split()
    s = SNP(chromosome = row[0].replace('.fa', ''),
            pos = int(row[1]),
            depth = int(row[2]),
            refbase = row[4],
            quality = int(row[5]),
            genotype = row[6],
            n_acgt = (map(int, row[-4:])))
    yield s



__all__ = (
  'SNP',
  'readSNPs_CASAVA'
)
