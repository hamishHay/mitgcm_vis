def sh_filter(data, lons, lats, lmax = 4):
  from pyshtools.expand import SHExpandLSQ, MakeGrid2D, MakeGridDH
  import pyshtools as pysh
  import numpy as np
  cilm, chi2 = SHExpandLSQ(data, lats, lons, lmax)
  #cilm[:, 0,:] = 0
  #cilm[:, 1,:] = 0

  print(cilm)

  lons2 = np.arange(0, 361, 2.0)
  lats2 = np.arange(90, -91, -2.0)

  x, y = np.meshgrid(lons2, lats2)

  clm = pysh.SHCoeffs.from_array(cilm)
  grid = clm.expand(lat=y, lon=x)
 
  #grid.nlat = 80
  #print(grid)
  return grid #MakeGridDH(cilm, sampling=2)




if __name__=='__main__':
  print("hi")
