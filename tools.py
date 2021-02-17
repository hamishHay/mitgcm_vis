# File to interpolate mitgcm data to a regular lat-lon grid
# inputs --> file_in: file name of the .data file 
#            field_name: field to perform the interpolation on
#            res: resolution in degrees of the interpolation grid 

# In addition to file_in, the following other files are required in 
# the same directory as file_in:
# XC, YC, RC, AngleCS, AngleSN

import numpy as np 
from MITgcmutils import mds
from xmitgcm import open_mdsdataset
from project_velocity import project_velocity
from interpolation import interpolate_data
import h5py

def convert_diagnostics(file_in, field_name, vectors=[], res=0.1, iter=20329):
    def load_data(name, rec=None, iter=20329):
        # if iter==np.nan, function will read all available iterations of file_in
        ds = mds.rdmds(name, itrs=iter, rec=rec, astype=np.float32)
        return ds

    field_keys = list(field_name.keys())

    field_num = len(field_keys)

    RAC = mds.rdmds('RAC')

    data_field = []
    cnt = 0
    for field in field_keys:
        # read field 
        print("Reading field " + field)
        data_field.append( load_data(file_in, rec=[ field_name[field] ], iter=iter) )
        data_field[-1][:] /= RAC
        cnt += 1

    # If there are any lat-lon vectors, rotate them to be EW-NS
    if len(vectors) > 0:
        CS = mds.rdmds('AngleCS')
        SN = mds.rdmds('AngleSN')

        for i in range(len(vectors)):
            c1 = vectors[i][0] # Vector component one
            c2 = vectors[i][1] # Vector component two
            data_field[c1], data_field[c2] = project_velocity(data_field[c1], data_field[c2], CS, SN)

    # Define target grid
    lon = np.arange(-180, 180.01, res) #+ hres/2
    lat = np.arange(-90, 90.01, res) #+ hres/2
    lon_out, lat_out = np.meshgrid(lon, lat)

    # Load input grid
    XC = mds.rdmds('XC')
    YC = mds.rdmds('YC')
    RC = mds.rdmds('RC')

    

    # Define sizes of target grid dimensions
    nrt, nxt, nyt = (len(RC), len(lon), len(lat))

    # Interpolate test field with nearest neighbour interpolation
    # Outputs the x and y indexes that map the input data onto the target grid
    # The interpolation only occurs in the x-y plane
    indx, indy, tshape = interpolate_data(data_field[0][0,:,:], XC, YC, lon_out, lat_out)[1:]


    out_file = h5py.File(file_in + ".h5", 'w')

    cnt = 0
    for field in field_keys:
        # Perform the interpolation   
        print("Interpolating field " + field)

        data_interp = data_field[cnt][:, indx, indy].reshape((nrt,nyt,nxt))

        cnt += 1
        # Save to file
        print("Saving field " + field)
        dset = out_file.create_dataset(field, data=data_interp, dtype=np.float32, chunks=True, compression="gzip")
        
    dset = out_file.create_dataset("lon", data=lon, dtype=np.float32, compression="gzip")
    dset = out_file.create_dataset("lat", data=lat, dtype=np.float32, compression="gzip")
    dset = out_file.create_dataset("RC", data=RC, dtype=np.float32, compression="gzip")
    
    out_file.close()

def convert_output(field_names, file_out, vectors=[], res=0.1, iters=20329):
    def load_data(name, rec=None, iters=20329):
        # if iter==np.nan, function will read all available iterations of file_in
        ds = mds.rdmds(name, itrs=iters, astype=np.float32)
        return ds

    field_num = len(field_names)

    data_field = []
    cnt = 0
    for field in field_names:
        # read field 
        print("Reading field " + field)
        data_field.append( load_data(field, iters=iters) )
        cnt += 1

    # If there are any lat-lon vectors, rotate them to be EW-NS
    if len(vectors) > 0:
        CS = mds.rdmds('AngleCS')
        SN = mds.rdmds('AngleSN')

        for i in range(len(vectors)):
            c1 = vectors[i][0] # Vector component one
            c2 = vectors[i][1] # Vector component two
            data_field[c1], data_field[c2] = project_velocity(data_field[c1], data_field[c2], CS, SN)

    # Define target grid
    lon = np.arange(-180, 180.01, res) #+ hres/2
    lat = np.arange(-90, 90.01, res) #+ hres/2
    lon_out, lat_out = np.meshgrid(lon, lat)

    # Load input grid
    XC = mds.rdmds('XC')
    YC = mds.rdmds('YC')
    RC = mds.rdmds('RC')

    # Define sizes of target grid dimensions
    nrt, nxt, nyt = (len(RC), len(lon), len(lat))

    # Interpolate test field with nearest neighbour interpolation
    # Outputs the x and y indexes that map the input data onto the target grid
    # The interpolation only occurs in the x-y plane
    indx, indy, tshape = interpolate_data(data_field[0][0,:,:], XC, YC, lon_out, lat_out)[1:]


    out_file = h5py.File(file_out, 'w')

    cnt = 0
    for field in field_names:
        # Perform the interpolation   
        print("Interpolating field " + field)

        data_interp = data_field[cnt][:, indx, indy].reshape((nrt,nyt,nxt))

        cnt += 1
        # Save to file
        print("Saving field " + field)
        dset = out_file.create_dataset(field, data=data_interp, dtype=np.float32, chunks=True, compression="gzip")
        
    dset = out_file.create_dataset("lon", data=lon, dtype=np.float32, compression="gzip")
    dset = out_file.create_dataset("lat", data=lat, dtype=np.float32, compression="gzip")
    dset = out_file.create_dataset("RC", data=RC, dtype=np.float32, compression="gzip")
    
    out_file.close()


if __name__=="__main__":
    # ---------------------------------------------
    # Example: 1 iteration of the diagnostics file:
    # file_in = "diagnostics_T_flux"

    # # Define the fields and their indexes
    # #                Name    :indx
    # field_names = {'ADVx_TH':0,
    #                'ADVy_TH':1,
    #                'ADVr_TH':2,
    #                'DFxE_TH':3,
    #                'DFyE_TH':4,
    #                'DFrE_TH':5,
    #                'DFrI_TH':6}

    # # Define the list of vectors, if any, that need 
    # # to be rotated to the correct coord system:
    # # Here, [0, 1] is ADVx_TH and ADVy_TH, which
    # # will be used together to perform the rotation
    # vectors = [[0, 1], [3, 4]]

    # convert_diagnostics(file_in, field_names, vectors, 0.1, iter=20329)

    # ile_in = "diagnostics_T_flux"

    # # Define the fields and their indexes
    # #                Name    :indx
    # field_names = {'UVEL':0,
    #                'VVEL':1,
    #                'WVEL':2}

    # # Define the list of vectors, if any, that need 
    # # to be rotated to the correct coord system:
    # # Here, [0, 1] is ADVx_TH and ADVy_TH, which
    # # will be used together to perform the rotation
    # vectors = [[0, 1]]

    # convert_diagnostics(file_in, field_names, vectors, 0.1, iter=20329)

    # Example 2: open raw data file
    iters = 50000
    field_names = ["U", "V", 'W']
    vectors = [[0,1]]
    file_out = "UVW_" +str(iters)+ ".h5"

    convert_output(field_names, file_out, vectors=vectors, iters=iters, res=0.5)


