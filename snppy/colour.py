import re

# reference white: D65

Adobe_RGB_XYZ = (
  (  0.5767309,  0.2973769,  0.0270343 ),
  (  0.1855540,  0.6273491,  0.0706872 ),
  (  0.1881852,  0.0752741,  0.9911085 ),
)
        
XYZ_Adobe_RGB = (
  (  2.0413690, -0.9692660,  0.0134474 ),
  ( -0.5649464,  1.8760108, -0.1183897 ),
  ( -0.3446944,  0.0415560,  1.0154096 ),
)



Apple_RGB_XYZ = (
  (  0.4497288,  0.2446525,  0.0251848 ),
  (  0.3162486,  0.6720283,  0.1411824 ),
  (  0.1844926,  0.0833192,  0.9224628 ),

)

XYZ_Apple_RGB = (
  (  2.9515373, -1.0851093,  0.0854934 ),
  ( -1.2894116,  1.9908566, -0.2694964 ),
  ( -0.4738445,  0.0372026,  1.0912975 ),
)



sRGB_XYZ = (
  (  0.4124564,  0.2126729,  0.0193339 ),
  (  0.3575761,  0.7151522,  0.1191920 ),
  (  0.1804375,  0.0721750,  0.9503041 ),
)

XYZ_sRGB = (
  (  3.2404542, -0.9692660,  0.0556434 ),
  ( -1.5371385,  1.8760108, -0.2040259 ),
  ( -0.4985314,  0.0415560,  1.0572252 ),
)



def clip(x, x_min = 0.0, x_max = 1.0):
  return min(max(x, x_min), x_max)



def xform(X, Y, Z, M):
  R = X * M[0][0] + Y * M[1][0] + Z * M[2][0]
  G = X * M[0][1] + Y * M[1][1] + Z * M[2][1]
  B = X * M[0][2] + Y * M[1][2] + Z * M[2][2]
  return R, G, B



def XYZtosRGB(X, Y, Z):
  def T(x):
    if x > 0.0031308:
      return 1.055 * (x ** (1/2.4)) - 0.055
    else:
      return 12.92 * x

  return map(T, xform(X, Y, Z, M = XYZ_sRGB))



def sRGBtoXYZ(R, G, B):
  def Trev(x):
    if x > 0.04045:
      return ((x + 0.055)/ 1.055) ** 2.4
    else:
      return x / 12.92

  return xform(*map(Trev, (R, G, B)), M = sRGB_XYZ)



def XYZtoRGB(X, Y, Z):
  return xform(X, Y, Z, M = XYZ_Adobe_RGB)



def RGBtoXYZ(R, G, B):
  return xform(R, G, B, M = Adobe_RGB_XYZ)



def RGBtoHSV(R, G, B):
  cmin = min(R, G, B)
  cmax = max(R, G, B)
  delta = cmax - cmin

  V = cmax

  if delta == 0:
    return 0.0, 0.0, V

  S = delta / cmax

  dR = (((cmax - R) / 6) + (delta / 2)) / delta
  dG = (((cmax - G) / 6) + (delta / 2)) / delta
  dB = (((cmax - B) / 6) + (delta / 2)) / delta

  if R == cmax:
    H = dB - dG
  elif G == cmax:
    H = 1.0/3.0 + dR - dB
  else:
    H = 2.0/3.0 + dG - dR

  if H<0: H = H + 1
  if H>1: H = H - 1

  return H, S, V



def HSVtoRGB(H, S, V):
  if S == 0:
    return V, V, V

  H = H * 6.0
  H = H % 6

  I = int(H)
  p1 = V * (1 - S)
  p2 = V * (1 - S * (H - I))
  p3 = V * (1 - S * (1 - (H - I)))

  if I == 0:  return  V, p3, p1
  if I == 1:  return p2,  V, p1
  if I == 2:  return p1,  V, p3
  if I == 3:  return p1, p2,  V
  if I == 4:  return p3, p1,  V
  if I == 5:  return  V, p1, p2



def RGBtoString(R, G, B):
  R = int(clip(R) * 255 + .5)
  G = int(clip(G) * 255 + .5)
  B = int(clip(B) * 255 + .5)
  return '#%02x%02x%02x' % (R, G, B)



def stringToRGB(s):
  def H1(x):
    return int(x, 16) / 255.0
  def H2(x):
    return int(x, 16) / 15.0
  re1 = re.compile('^#([0-9A-F][0-9A-F])([0-9A-F][0-9A-F])([0-9A-F][0-9A-F])$', re.I)
  re2 = re.compile('^#([0-9A-F])([0-9A-F])([0-9A-F])$', re.I)
  m = re1.match(s)
  if m is not None:
    return map(H1, m.groups())
  m = re2.match(s)
  if m is not None:
    return map(H2, m.groups())
  return 0.0, 0.0, 0.0



import bisect
class ColourRamp(object):
  def __init__(self, control_points):
    self.c = sorted(control_points)
    self.rmin = self.c[0][0]
    self.rmax = self.c[-1][0]

  def XYZ(self, x):
    def interp(a,b,c):
      return (b-a) * c + a
    def interp3(a,b,c):
      return interp(a[0],b[0],c), interp(a[1],b[1],c), interp(a[2],b[2],c)

    if x <= self.rmin:
      x,y,z = self.c[0][1]
    elif x >= self.rmax:
      x,y,z = self.c[-1][1]
    else:
      p = bisect.bisect_left(self.c, (x, None)) - 1
      (r1,c1),(r2,c2) = self.c[p], self.c[p+1]
      x = (x - r1) / (r2 - r1)
      x,y,z = interp3(c1, c2, x)
    return x,y,z

  def RGB(self, x):
    return XYZtoRGB(*self.XYZ(x))

  @classmethod
  def equispaced(cls, range_min, range_max, *colours):
    pal = [ RGBtoXYZ(*stringToRGB(s)) for s in colours ]
    pos = [ range_min + x * float(range_max - range_min) / (len(pal)-1) for x in xrange(len(pal)) ]
    return cls(zip(pos, pal))



YlGnBu = ColourRamp.equispaced(0.0, 1.0, "#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58")
RdWtBu = ColourRamp.equispaced(0.0, 1.0, "#aa0000", "#FFFFFF", "#0000aa")
Spectral = ColourRamp.equispaced(0.0, 1.0, "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")

