import numpy as np
import dask.array as da 
import time 
from MITgcmutils import mds
import h5py
import glob, os
import sys

sys.path.append("/home6/hhay/mitgcm_vis")
from project_velocity import project_velocity, interpolate_data

def get_iters(name):
  #os.chdir("/mydir")
  iters = [i.split('.')[2].lstrip("0") for i in glob.glob(name+".*.data")]
  
  try:
    iters.remove('')
  except ValueError:
    pass 

  iters = sorted([int(iter) for iter in iters])
    
  return iters


iters = get_iters("./U")
iters_tAvg = get_iters("./uVeltave")

XC = mds.rdmds('XC')
YC = mds.rdmds('YC')
RC = mds.rdmds('RC')

nx, ny = XC.shape

CS = mds.rdmds('AngleCS')
SN = mds.rdmds('AngleSN')

U = mds.rdmds('U', np.inf)



hres = 0.5 
lon = np.arange(-180, 180.01, hres) #+ hres/2
lat = np.arange(-90, 90.01, hres) #+ hres/2
lon_out, lat_out = np.meshgrid(lon, lat)

nrt, nxt, nyt = (len(RC), len(lon), len(lat))

Ui = np.zeros((len(RC), len(lat), len(lon)))
Vi = np.zeros((len(RC), len(lat), len(lon)))
Wi = np.zeros((len(RC), len(lat), len(lon)))

indx, indy, tshape = interpolate_data(U[0], XC, YC, lon_out, lat_out)[1:]




def get_data(name, iter=None):
   if iter is None:
     return mds.rdmds(name)
   else:
     return mds.rdmds(name, iter)

f = h5py.File("DATA.h5", 'w')
grp = f.create_group('DATA')
grp.create_dataset("ITER", data=np.array(iters))
grp.create_dataset("RC", data=np.array(RC))

read_list = ["T", "W"]
#print(len(iters))
#iters=iters[:1]


for dname in read_list:
  dset = grp.create_dataset(dname,(len(iters), nrt, nyt, nxt), chunks=True, dtype=np.float32)
  cnt = 0
  for iter in iters:
    start = time.time()
    data = get_data(dname, iter)#.reshape((nrt, XC.size))
    end = time.time()
    print("Load time: ", end-start)

    start = time.time()
    data_near = data[:, indx, indy].reshape((nrt,nyt,nxt))
    end = time.time()


    start = time.time()
    dset[cnt, :, :, :] = data_near
    end = time.time()
    print("h5py", end-start)
    cnt += 1

    

dset_u = grp.create_dataset("U",(len(iters), nrt, nyt, nxt), chunks=True, dtype=np.float32)
dset_v = grp.create_dataset("V",(len(iters), nrt, nyt, nxt), chunks=True, dtype=np.float32)
cnt = 0
for iter in iters:
  U = get_data("U", iter)
  V = get_data("V", iter)
    
  U, V = project_velocity(U, V, CS, SN)
 
  U_near = U[:, indx, indy].reshape((nrt,nyt,nxt))
  V_near = V[:, indx, indy].reshape((nrt,nyt,nxt))  

  dset_u[cnt, :, :, :] = U_near
  dset_v[cnt, :, :, :] = V_near
  
  cnt += 1


read_list = ['wVeltave', 'Ttave']
for dname in read_list:
  dset = grp.create_dataset(dname,(len(iters_tAvg), nrt, nyt, nxt), chunks=True, dtype=np.float32)
  cnt = 0
  for iter in iters_tAvg:
    start = time.time()
    data = get_data(dname, iter)#.reshape((nrt, XC.size))
    end = time.time()
    print("Load time: ", end-start)

    start = time.time()
    data_near = data[:, indx, indy].reshape((nrt,nyt,nxt))
    end = time.time()


    start = time.time()
    dset[cnt, :, :, :] = data_near
    end = time.time()
    print("h5py", end-start)
    cnt += 1



dset_u = grp.create_dataset("uVeltave",(len(iters_tAvg), nrt, nyt, nxt), chunks=True, dtype=np.float32)
dset_v = grp.create_dataset("vVeltave",(len(iters_tAvg), nrt, nyt, nxt), chunks=True, dtype=np.float32)
cnt = 0
for iter in iters_tAvg:
  U = get_data("uVeltave", iter)
  V = get_data("vVeltave", iter)

  U, V = project_velocity(U, V, CS, SN)

  U_near = U[:, indx, indy].reshape((nrt,nyt,nxt))
  V_near = V[:, indx, indy].reshape((nrt,nyt,nxt))

  dset_u[cnt, :, :, :] = U_near
  dset_v[cnt, :, :, :] = V_near

  cnt += 1


f.close()

