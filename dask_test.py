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
  iters.remove('')
  iters = sorted([int(iter) for iter in iters])
    
  return iters


#nx, ny = (510, 3060)
#indx = np.random.randint(0, high=nx, size=65341)
#indy = np.random.randint(0, high=ny, size=65341)

iters = get_iters("./U")

XC = mds.rdmds('XC')
YC = mds.rdmds('YC')
RC = mds.rdmds('RC')

nx, ny = XC.shape

CS = mds.rdmds('AngleCS')
SN = mds.rdmds('AngleSN')

U = mds.rdmds('U', np.inf)
V = mds.rdmds('V', np.inf)
U, V = project_velocity(U, V, CS, SN)

hres = 1.0
lon = np.arange(-180, 180.01, hres) #+ hres/2
lat = np.arange(-90, 90.01, hres) #+ hres/2
lon_out, lat_out = np.meshgrid(lon, lat)

nrt, nxt, nyt = (len(RC), len(lon), len(lat))

Ui = np.zeros((len(RC), len(lat), len(lon)))
Vi = np.zeros((len(RC), len(lat), len(lon)))
Wi = np.zeros((len(RC), len(lat), len(lon)))

indx, indy, tshape = interpolate_data(U[0], XC, YC, lon_out, lat_out)[1:]
#inds = np.unravel_index(inds, U[0].shape)
#print()
#indx, indy = np.unravel_index(inds, XC.shape )

def get_data(name, iter=None):
   if iter is None:
     return mds.rdmds(name)
   else:
     return mds.rdmds(name, iter)

f = h5py.File("DATA.h5", 'a')
grp = f.create_group('DATA')

read_list = ["T", "W"]

iters=[10]

print(indx, indy)
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
 
  #U = U.reshape((nrt, XC.size))
  #V = V.reshape((nrt, XC.size))

  U_near = U[:, indx, indy].reshape((nrt,nyt,nxt))
  V_near = V[:, indx, indy].reshape((nrt,nyt,nxt))  

  dset_u[cnt, :, :, :] = U_near
  dset_v[cnt, :, :, :] = V_near
  
  cnt += 1

f.close()

