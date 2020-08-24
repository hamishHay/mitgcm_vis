import numpy as np
from MITgcmutils import mds
from xmitgcm import open_mdsdataset
import pyresample
import sys 
sys.path.append("/home/hay/Research/python/plotting/mitgcm_vis")

def calc_mean_EK(start=100, end=1000, step=100):

  times = np.arange(start, end+1, step, dtype=np.int)

  var_list = ['U', 'V', 'W']
  ds2 = open_mdsdataset('.', prefix=['CN', 'SN', 'XC', 'YC'], iters=100, geometry='curvilinear', ignore_unknown_vars=True)
  XC = np.array(ds2.XC)
  YC = np.array(ds2.YC)

  #EK = np.zeros(XC.shape)

  for i in range(len(times)):
    time = times[i]
    print("TIME: ", time/times[-1])
    ds = open_mdsdataset('.', prefix=var_list, iters=time, geometry='curvilinear', ignore_unknown_vars=True, default_dtype=np.float)

    U = np.array(ds.U)[0,:,:]
    V = np.array(ds.V)[0,:,:]
    W = np.array(ds.W)[0,:,:]

    U2 = np.sqrt(U*U + V*V + W*W)

    if i==0:
      EK = U2.copy()
    else:
      EK += U2


  return EK
