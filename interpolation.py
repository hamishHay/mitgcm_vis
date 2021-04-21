from scipy.spatial import cKDTree
import numpy as np

def interpolate_data(data, lon1, lat1, lon2, lat2, inds=None):
    
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
    # import matplotlib.pyplot as plt 

    ix, iy = np.unravel_index(inds, data.shape)
    # data_near2 = data[ix,iy].reshape(lon2.shape)
    
    return data_nearest, ix, iy, lon2.shape