import numpy as np

def project_velocity(U, V, CS, SN):
    if U.ndim == 2 and V.ndim == 2:
      Umap = U*CS - V*SN
      Vmap = U*SN + V*CS 
   
    elif V.ndim == 3 and V.ndim == 3:
      Umap = U*CS[None, : :] - V*SN[None, :, :]
      Vmap = U*SN[None, :,:] + V*CS[None, :, :]

    return Umap, Vmap


def interpolate_data(data, lon1, lat1, lon2, lat2, inds=None):
    from scipy.spatial import cKDTree
    #import pyresample

    #grid_src = pyresample.geometry.SwathDefinition(lons=lon1, lats=lat1)

    #grid_target = pyresample.geometry.GridDefinition(lons=lon2,lats=lat2)


    #if data.ndim==2:
    #    data_interp = pyresample.kd_tree.resample_nearest(grid_src, data ,grid_target,radius_of_influence=100000,fill_value=0)
        #data_interp = pyresample.kd_tree.resample_gauss(grid_src, data,
        #                                                grid_target, radius_of_influence=200000, sigmas=25000)
#       data_interp = pyresample.bilinear.resample_bilinear(data, grid_src, grid_target,
#                                                     radius=50e3, neighbours=32,
#                                                     nprocs=1, fill_value=0,
#                                                     reduce_data=True, segments=None,
#                                                     epsilon=0)
    #elif data.ndim==3:
    
    def sph2xyz(lons, lats):
      lons = np.radians(lons)
      lats = np.radians(lats)
      
  
      x = np.cos(lons) * np.cos(lats)
      y = np.sin(lons) * np.cos(lats)
      z = np.sin(lats)
      print("making coord array")
      xyz = np.column_stack((x, y, z))
      #print(x)
      #print(xyz)
      print("coord array created")
      return xyz 

    # if no indices for nearest-neighbour interp are supplied,
    # calculate them using the KDTree algorithm
    if inds is None: 
      xyzs = sph2xyz(lon1.flatten(), lat1.flatten())
      xyzt = sph2xyz(lon2.flatten(), lat2.flatten())
      print(xyzs)

      print("Creating ckdTree")
      tree = cKDTree(xyzs)
      print("Querying ckdTree")
      d, inds = tree.query(xyzt, k = 1)
      print("Query complete")
    #rint(np.size(data))
    #data = data.compute()
    #print(type(data), type(inds), type(lon2))
    print("Getting nearest neighbour values")
    data_nearest = data.flatten()[inds].reshape(lon2.shape)
    print("Values retrieved")   
    #print(data_nearest.shape, lon2.shape, lat2.shape)
    return data_nearest, inds
