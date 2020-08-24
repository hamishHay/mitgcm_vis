import numpy as np

def project_velocity(U, V, CS, SN):
    if U.ndim == 2 and V.ndim == 2:
      Umap = U*CS - V*SN
      Vmap = U*SN + V*CS 
   
    elif V.ndim == 3 and V.ndim == 3:
      Umap = U*CS[None, : :] - V*SN[None, :, :]
      Vmap = U*SN[None, :,:] + V*CS[None, :, :]

    return Umap, Vmap


def interpolate_data(data, lon1, lat1, lon2, lat2):
    import pyresample

    grid_src = pyresample.geometry.SwathDefinition(lons=lon1, lats=lat1)

    grid_target = pyresample.geometry.GridDefinition(lons=lon2,lats=lat2)


    if data.ndim==2:
        data_interp = pyresample.kd_tree.resample_nearest(grid_src, data ,grid_target,radius_of_influence=100000,fill_value=0)
        #data_interp = pyresample.kd_tree.resample_gauss(grid_src, data,
        #                                                grid_target, radius_of_influence=200000, sigmas=25000)
#       data_interp = pyresample.bilinear.resample_bilinear(data, grid_src, grid_target,
#                                                     radius=50e3, neighbours=32,
#                                                     nprocs=1, fill_value=0,
#                                                     reduce_data=True, segments=None,
#                                                     epsilon=0)
    elif data.ndim==3:
       data_interp = np.zeros( (np.shape(data)[0], np.shape(lon2)[0], np.shape(lon2)[1]) )
       for i in range(np.shape(data_interp)[0]):
          data_interp[i] =  pyresample.kd_tree.resample_nearest(grid_src, data[i] ,grid_target,radius_of_influence=100000,fill_value=0)                    

    return data_interp
