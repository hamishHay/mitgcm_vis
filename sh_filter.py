def sh_filter(data, lons, lats, lmax = 4):
  from pyshtools.expand import SHExpandLSQ, MakeGrid2D, MakeGridDH

  cilm, chi2 = SHExpandLSQ(data, lats, lons, lmax)
  #cilm[:, 0,:] = 0
  #cilm[:, 1,:] = 0

  return MakeGridDH(cilm,sampling=2)




if __name__=='__main__':
  print("hi")
