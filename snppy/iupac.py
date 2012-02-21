import string

__iupac_revcomp = string.maketrans('XACMGRSVTWYHKDBNxacmgrsvtwyhkdbn', 'XTGKCYSBAWRDMHVNxtgkcysbawrdmhvn')

def toIUPAC(nts):
  i = 0
  if 'A' in nts or 'a' in nts: i |= 1
  if 'C' in nts or 'c' in nts: i |= 2
  if 'G' in nts or 'g' in nts: i |= 4
  if 'T' in nts or 't' in nts: i |= 8
  return 'XACMGRSVTWYHKDBN'[i]

__iupac_to_nt = dict(
    X='',
    A='A',
    C='C',
    M='AC',
    G='G',
    R='AG',
    S='CG',
    V='ACG',
    T='T',
    W='AT',
    Y='CT',
    H='ACT',
    K='GT',
    D='AGT',
    B='CGT',
    N='ACGT'
)

def fromIUPAC(iupac):
  return __iupac_to_nt.get(iupac.upper(), 'X')

def revcompIUPAC(iupac):
  return ''.join(reversed(string.translate(iupac, __iupac_revcomp)))
