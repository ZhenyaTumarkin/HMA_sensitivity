'''
script to pick points on a glacier

Input data:
RGI6 shp files
ASTER (or other) DEM, 1x1 degree tiles
Rounce2021 debris thckness  (debris shp derrived from here, alternatively could load it from Scherler 2018)
                                                                        NOTE: Scherler 2021 will not match

Outputs:
csv containing points, their weights, elevations and debris thickness
'''


###TO DO/ COULD DO
#Output data, average across the cluster, instead of taking the point value
#
#Fix pointpicking, where vectorisation fails and misses a hole 
#    (currently if point picked actually in dem hole,point is discarded, only an issue for small clusters)



import matplotlib.pyplot as plt
from sklearn import cluster

import rasterio
from rasterio import mask,merge,warp
import rasterio.features
from rasterio.io import MemoryFile

import numpy as np
import time
from shapely.geometry import Polygon,Point,shape
from shapely import maximum_inscribed_circle

import pandas as pd
import argparse
import geopandas as gpd
from geopandas import GeoDataFrame
import matplotlib.pyplot as plt
import utm
import os


t1 = time.time()




#parse rgiid input
parser = argparse.ArgumentParser()
parser.add_argument('--rgiid',type=str,required = True,help = 'input RGIId')
parser.add_argument('--outlocation',type=str,required = False,help = 'out_location')
parser.add_argument('--N_points',type=int,required = True,help = 'input num points')
args = parser.parse_args()
Id = args.rgiid
N_points = args.N_points
outlocation = args.outlocation
print(f'Glacier Id : {Id}')

inputs_path='/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'
if outlocation:
    root=f'{outlocation}'
else:
    root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'
# N_points=30   #make this an input 


edge_buffer = 0 #m   #remove this width from the edge of the glacier (avoid edge effects)
show_result = False     # plot results of point picking?


region_num = Id.split('-')[1][0:2]   # extract RGI region number from Id
Id_num = Id.split('-')[1]
region_names = {
    '13':'CentralAsia',
    '14':'SouthAsiaWest',
    '15':'SouthAsiaEast'
}
region_name= region_names[str(region_num)]


geo_df = gpd.read_file(f'{inputs_path}/{region_num}_rgi60_{region_name}/{region_num}_rgi60_{region_name}.shp')  #read region shp
debris_thick_path =f'{inputs_path}/Rounce2021_deb_thick/{region_num}/{Id_num}_hdts_m.tif'  #path for glacier debris thickness
debris_extrap_thick_path = f'{inputs_path}/Rounce2021_deb_thick/{region_num}/{Id_num}_hdts_m_extrap.tif' #path for glacier extrapolated debris thickness
path_dems = f'{inputs_path}/ASTER_DEM/'   #folder with dems



if os.path.isfile(debris_thick_path)==True:  # check which debris path to use
    debris_path = debris_thick_path
elif os.path.isfile(debris_extrap_thick_path)==True:
    debris_path = debris_extrap_thick_path
else:
    print('no debris')
    debris_path= False

######extract shp_outline of debris (from raster), and load raster 
if debris_path:
    with rasterio.open(debris_path) as src:
        debris_thickness = src.read(1, masked=True)
        #debris_thickness = np.where(debris_thickness<0.001,np.nan,debris_thickness)
        debris_transform = src.transform
        debris_th_meta = src.meta
        shape_gen = ((shape(i), j) for i, j in rasterio.features.shapes(debris_thickness, transform=src.transform))
        debris_glacier = gpd.GeoDataFrame(dict(zip(["geometry", "class"], zip(*shape_gen))), crs=src.crs)
    debris_glacier_bool = True
else:
    debris_glacier_bool = False




glacier =  geo_df[geo_df['RGIId']==Id]    # extract required glacier

#######################
# process dem, if glacier on edge join two tiles
######################
minx=glacier.bounds.minx.values.astype(int).squeeze()   # truncate at int (ie round down)
maxx = np.ceil(glacier.bounds.maxx.values).squeeze()   # round up
miny=glacier.bounds.miny.values.astype(int).squeeze()
maxy = np.ceil(glacier.bounds.maxy.values).squeeze()



######choose utm zone
rem = minx%6
utm_border = (minx-rem)/6
utm_zone_num = int(utm_border+31)
dst_crs = f'EPSG:326{utm_zone_num}'  # CRS for the required utm_zone

def mergetiles(files):   # merges two or more  tiles
    mosaic, out_transform = merge.merge(files)

    out_meta = files[0].meta.copy()       # Update metadata
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_transform
    })
    # with rasterio.open(out_name, "w", **out_meta) as dest:   # Save merged DEM
    #     dest.write(mosaic)
    return mosaic , out_transform,out_meta


if (maxx-minx ==1) & (maxy-miny==1):     # no merging required
    raster_path = f'{path_dems}ASTGTMV003_N{miny}E{minx:03d}/ASTGTMV003_N{miny}E{minx:03d}_dem.tif'
    with rasterio.open(raster_path) as src:
        region_raster = src.read(1)
        region_transform = src.transform
        region_meta = src.meta
else:
    print('need merge of nearby tiles')
    if (maxx-minx ==2) & (maxy-miny==1):   #simple edge case
        dem1 = rasterio.open(f'{path_dems}ASTGTMV003_N{miny}E{minx:03d}/ASTGTMV003_N{miny}E{minx:03d}_dem.tif')
        dem2 = rasterio.open(f'{path_dems}ASTGTMV003_N{miny}E{(minx+1):03d}/ASTGTMV003_N{miny}E{(minx+1):03d}_dem.tif')
        region_raster , out_transform,region_meta=mergetiles([dem1,dem2])
        region_raster = region_raster.squeeze()

    
    elif (maxx-minx ==1) & (maxy-miny==2):   #simple edge case
        dem1 = rasterio.open(f'{path_dems}ASTGTMV003_N{miny}E{minx:03d}/ASTGTMV003_N{miny}E{minx:03d}_dem.tif')
        dem2 = rasterio.open(f'{path_dems}ASTGTMV003_N{miny+1}E{minx:03d}/ASTGTMV003_N{miny+1}E{minx:03d}_dem.tif')
        region_raster , out_transform,region_meta=mergetiles([dem1,dem2])
        region_raster = region_raster.squeeze()
    
    elif (maxx-minx ==2) & (maxy-miny==2):  #all 4 tiles included
        dem1 = rasterio.open(f'{path_dems}ASTGTMV003_N{miny}E{minx:03d}/ASTGTMV003_N{miny}E{minx:03d}_dem.tif')
        dem2 = rasterio.open(f'{path_dems}ASTGTMV003_N{miny+1}E{minx:03d}/ASTGTMV003_N{miny+1}E{minx:03d}_dem.tif')
        dem3 = rasterio.open(f'{path_dems}ASTGTMV003_N{miny+1}E{(minx+1):03d}/ASTGTMV003_N{miny+1}E{(minx+1):03d}_dem.tif')
        dem4 = rasterio.open(f'{path_dems}ASTGTMV003_N{miny}E{(minx+1):03d}/ASTGTMV003_N{miny}E{(minx+1):03d}_dem.tif')
        region_raster , out_transform,region_meta=mergetiles([dem1, dem2,dem3,dem4])
        region_raster = region_raster.squeeze()
    else:
        print(minx,miny,maxx,maxy)
    



#############################
##now crop dem
#############################
### crop to shp file of glacier
### transform to utm coordinates (so distances become in m)
###                   Necessary later, when inscribing circle into geometry
### crop out a buffer region, so points not too close to edge
### crop out clean and debris rasters by using debris shp file

                 

glacier_buffered = glacier.to_crs(dst_crs).copy()  # create a buffer on the edge of the glacier
glacier_buffered['geometry'] = glacier_buffered.geometry.buffer(-edge_buffer)


with MemoryFile().open(**region_meta) as dataset:     #this step crops to the shape of the glacier
    dataset.write(region_raster[np.newaxis, ...])
    glacier_raster, glacier_transform = mask.mask(dataset, glacier.geometry.values, crop=True)
    glacier_meta = region_meta.copy()
    glacier_meta.update({"driver": "GTiff","height": glacier_raster.shape[1],
                        "width": glacier_raster.shape[2],"transform": glacier_transform})
    
with MemoryFile().open(**glacier_meta) as dataset:     ##### now convert to utm
    dataset.write(glacier_raster)
    utm_transform, utm_width, utm_height = rasterio.warp.calculate_default_transform(dataset.crs, dst_crs, dataset.width, dataset.height, *dataset.bounds)
    kwargs = dataset.meta.copy()
    kwargs.update({'crs': dst_crs,'transform': utm_transform,
                   'width': utm_width,'height': utm_height})
    
    utm_glacier_meta = glacier_meta.copy()
    utm_glacier_meta.update({"driver": "GTiff", "transform": utm_transform,"crs":dst_crs,'height':utm_height,'width':utm_width})
    
    with MemoryFile().open(**utm_glacier_meta) as utm_dataset:
        warp.reproject(source=rasterio.band(dataset, 1),          # actual reprojection occurs here
                destination=rasterio.band(utm_dataset, 1),
                src_transform=dataset.transform,
                src_crs=dataset.crs,
                dst_transform=utm_transform,
                dst_crs=dst_crs,
                resampling=warp.Resampling.nearest)
        whole_glacier_utm = utm_dataset.read(1)

        #print(utm_dataset.crs)


        ##now crop out buffered zone
        buffered_glacier_raster, buffered_glacier_transform = mask.mask(utm_dataset, glacier_buffered.to_crs(dst_crs).geometry.values, crop=False)
        buffered_glacier_meta = utm_dataset.meta.copy()
        buffered_glacier_meta.update({
            "driver": "GTiff",
            "height":buffered_glacier_raster.shape[1],
            "width": buffered_glacier_raster.shape[2],
            "transform": buffered_glacier_transform})
        

       






if debris_glacier_bool:

    ###read debris thickness data
    with rasterio.open(debris_path) as src:
        debris_thickness = src.read(1, masked=True)
        # deb_nan_mask = debris_thickness.mask
        # deb_nan_mask_updated = deb_nan_mask.copy()
        # deb_nan_mask_updated[debris_thickness.data <= 0.1] = True
        # debris_thickness = np.ma.array([debris_thickness.data],mask = deb_nan_mask_updated)
        #debris_thickness = np.where(debris_thickness<0.001,0,debris_thickness)
     
        debris_transform = src.transform
        debris_th_meta = src.meta

        ###reproject and resample onto constant grid the debris
        with MemoryFile().open(**debris_th_meta) as dataset:     ##### now convert to utm
            dataset.write(debris_thickness[np.newaxis, ...])
            
            utm_deb_transform, utm_deb_width, utm_deb_height = warp.calculate_default_transform(dataset.crs, dst_crs, dataset.width, dataset.height, *dataset.bounds)
            deb_th_utm = np.zeros_like(buffered_glacier_raster[0], dtype=debris_thickness.dtype)
            warp.reproject(
            source=debris_thickness,
            destination=deb_th_utm,
            src_transform=debris_transform,
            src_crs=debris_th_meta['crs'],
            dst_transform=buffered_glacier_transform,
            dst_crs=dst_crs,
            resampling=warp.Resampling.nearest
        )


    #create deb mask (where debris thickness>0)
    deb_present = deb_th_utm>0 ### change limit directly here for filtering too thin?
    deb_raster = buffered_glacier_raster.copy()
    deb_raster[0][~deb_present] = 0
    clean_raster = buffered_glacier_raster.copy()
    clean_raster[0][deb_present] = 0

    
    


    


    area_debris = len(deb_raster[deb_raster>0])
    area_clean = len(clean_raster[clean_raster>0])
    clean_weight = area_clean/(area_clean+area_debris)
    deb_weight = area_debris/(area_clean+area_debris)
    points_clean = round(clean_weight*N_points)
    points_deb = round(deb_weight*N_points)
      
    if points_deb < 1:
        points_deb = 1

    clean_transform = buffered_glacier_transform
      
        
else:
    clean_raster,clean_transform,clean_meta = buffered_glacier_raster.copy(),buffered_glacier_transform,buffered_glacier_meta.copy()
    points_clean = N_points



#################################################
#now create numpy array with inputs to the Kmeans clustering alg
#################################################
def cluster_alg(num_points,weighted_input,mask_nan,raster,land_class_weight): 
    '''
    inputs: 
    num_points = number of clusters to make
    weighted input = M*N array. M = number of points on glacier (nan_mask applied outside of function)
                                N = number of variables considered during clustering
    mask_nan = collapse 2D dem into 1D vector. value = True where points are on glacier
    raster = dem of land_cover (requred to double check point chosen is actually valid!)
    '''    


    #print(num_points)
    kmeans = cluster.KMeans(n_clusters=num_points, random_state=1).fit(weighted_input)   #Clustering occurs heree
    labels_valid = kmeans.labels_
    labels_whole =np.full((1,shape[0]*shape[1]), np.nan)

    labels_whole[mask_nan] = labels_valid
    labels_grid = labels_whole.reshape((shape[0], shape[1]))    #2D map of the different clusters
    

    ################ for each cluster pick 'center' of largest region
    mask = np.isnan(labels_grid) ==False
    shapes = rasterio.features.shapes(labels_grid,mask = mask,transform = clean_transform) # vectorise the labels
    dict_polygons = {'class_list' : [],
                    'area_list' : [],}
    polygons_list = []
    for polygon_all,val in shapes:   #iterate adding the polygons to a dictionary
        polygon = Polygon(polygon_all['coordinates'][0])
        polygons_list.append(polygon)
        dict_polygons['area_list'].append(polygon.area)
        dict_polygons['class_list'].append(val)
    df_polygons=pd.DataFrame.from_dict(dict_polygons)
    mask_max_area = df_polygons.groupby('class_list')['area_list'].idxmax()   #in each cluster type, find the largest polygon by area
    #print(mask_max_area)
    polygons_max = np.array(polygons_list)[mask_max_area]    ##select the largest polygon in each cluster
    # cluster_total_area = df_polygons.groupby('class_list')['area_list'].sum() 
    # weight = np.array(cluster_total_area/total_area)   #find weight
    weights_new = np.array(np.bincount(labels_valid, minlength=num_points))

    x_ind = np.full(num_points,0)
    y_ind = np.full(num_points,0)
    for i,polygon in enumerate(polygons_max):
        inscribed_circ = maximum_inscribed_circle(polygon)   #inscribe circle to find widest region of cluster
        x_tmp= (inscribed_circ.coords[0][0] - clean_transform.c)/clean_transform.a
        y_tmp= (clean_transform.f - inscribed_circ.coords[0][1] )/clean_transform.a
        # print(x_tmp,y_tmp)
        # print(inscribed_circ.coords[0])
        #print((inscribed_circ.coords[0][1]),(inscribed_circ.coords[0][0]))
        if raster[int(y_tmp),int(x_tmp)]!=0:  #check x,y xhosen actually in glacier boundary
            
            x_ind[i] = x_tmp
            y_ind[i] = y_tmp
        else:
            print('point_discarded')
            
    
    mask_zero = (x_ind!=0)&(y_ind!=0)
    weight_masked = (weights_new)[mask_zero]
    normalised_weight = weight_masked/np.sum(weight_masked)*land_class_weight
    return x_ind[mask_zero],y_ind[mask_zero],labels_grid,normalised_weight

    
   




x_s = np.arange(1,buffered_glacier_raster.shape[2]+1,1)    # np.arange max val is val-1 (hence need +1)
y_s = np.arange(1,buffered_glacier_raster.shape[1]+1,1) 

xx,yy = np.meshgrid(x_s,y_s)  # create 2D arrays of coordinates
shape = np.shape(xx)      # example shape of glacier grid  (from UTM transform but before cropping the buffered region)
X = xx.reshape(1, shape[0]*shape[1])  #collapse to 1D arrays
Y = yy.reshape(1, shape[0]*shape[1])

total_area = np.shape(whole_glacier_utm[whole_glacier_utm>0])[0]   #find area of glacier



#create weighted input for clean ice
elevation = clean_raster[0].reshape(1, shape[0]*shape[1])
mask_nan = elevation>0   #1D mask for on glacier points
elevation_mask,X_mask,Y_mask = elevation[mask_nan],X[mask_nan],Y[mask_nan]
non_scaled_matrix = np.array([elevation_mask,X_mask,Y_mask])
scaled_matrix = np.zeros_like(non_scaled_matrix).astype(float)
for i in range (len(non_scaled_matrix)):    #rescale all of the inputs so they are between 0 and 1 (and have range 1)
    scaled_matrix[i]= (non_scaled_matrix[i] - np.nanmin(non_scaled_matrix[i]))/(np.nanmax(non_scaled_matrix[i])-np.nanmin(non_scaled_matrix[i]))

weights = [5,1,1]  # elevation, x,y
weighted_input_clean = np.multiply(scaled_matrix.T,np.array([weights]))


x_ind_clean,y_ind_clean,labels_clean,weight_clean = cluster_alg(points_clean,weighted_input_clean,mask_nan,clean_raster[0],clean_weight)

x_clean = utm_transform.c + x_ind_clean*clean_transform.a   #deal with python's/rasterio's way of handling copordinates
y_clean = utm_transform.f - y_ind_clean*clean_transform.a
df_clean = pd.DataFrame({'x':x_clean,'y':y_clean,})    #create dataframe for output
df_clean['deb_bool'] = np.zeros_like(df_clean['x']).astype(int)    
df_clean['weight'] = weight_clean
df_clean['deb_thickness_m'] = np.zeros_like(df_clean['x']).astype(int)
df_clean['elev_m'] = clean_raster[0,y_ind_clean,x_ind_clean]
df_clean['lat'],df_clean['lon'] = utm.to_latlon(x_clean,y_clean, utm_zone_num, 'N') 

#####now debris covered, includes debris thickness!
if debris_glacier_bool: 
    debris_thickness_1d = deb_th_utm.reshape(1, shape[0]*shape[1])
    elevation = deb_raster[0].reshape(1, shape[0]*shape[1])
    mask_nan = (debris_thickness_1d>0)&(elevation>0)
    elevation_mask,X_mask,Y_mask,deb_thick_mask = elevation[mask_nan],X[mask_nan],Y[mask_nan],debris_thickness_1d[mask_nan]
    non_scaled_matrix = np.array([elevation_mask,X_mask,Y_mask,deb_thick_mask])
    scaled_matrix = np.zeros_like(non_scaled_matrix).astype(float)
    for i in range (len(non_scaled_matrix)):  #rescale to between 0,1
        range = np.nanmax(non_scaled_matrix[i])-np.nanmin(non_scaled_matrix[i])
        if range !=0:
            scaled_matrix[i]= (non_scaled_matrix[i] - np.nanmin(non_scaled_matrix[i]))/(np.nanmax(non_scaled_matrix[i])-np.nanmin(non_scaled_matrix[i]))
        else:
            scaled_matrix[i] = np.zeros_like(non_scaled_matrix[i])
 
    weights = [5,1,1,2]  # elevation, x,y, deb_thickness


    weighted_input_deb = np.multiply(scaled_matrix.T,np.array([weights]))
    x_ind_deb,y_ind_deb,labels_deb,weight_deb = cluster_alg(points_deb,weighted_input_deb,mask_nan,deb_raster[0],deb_weight)

    x_deb= utm_transform.c + x_ind_deb*clean_transform.a
    y_deb = utm_transform.f - y_ind_deb*clean_transform.a
    df_deb = pd.DataFrame({'x':x_deb,'y':y_deb,})
    df_deb['deb_bool'] = np.ones_like(df_deb['x']).astype(int)
    df_deb['deb_thickness_m'] = deb_th_utm[y_ind_deb,x_ind_deb]
    df_deb['elev_m'] = deb_raster[0,y_ind_deb,x_ind_deb]
    df_deb['lat'],df_deb['lon'] = utm.to_latlon(x_deb,y_deb, utm_zone_num, 'N') 
   
    df_deb['weight'] = weight_deb

    print(df_deb)
    df_out = pd.concat([df_clean,df_deb])

else:
    df_out = df_clean


##############check points lie in RGI6 outline (probably unnecessary due to check in cluster alg)
points = [Point(xi, yi) for xi, yi in zip(df_out['x'], df_out['y'])]
mask = np.array([glacier.to_crs(dst_crs).geometry.contains(pt) for pt in points])
df_out = df_out[mask.squeeze()]
if any(mask) == False:
    print('points ouside of region')




############normalise weights (in case  of clusters being removed)
tot_weight = df_out.weight.sum()
df_out['weight'] = df_out['weight']/tot_weight
# df_out = df_out[df_out['weight']>1/(N_points*10)]  #remove points that contributr less than 0.5% to final MB
# tot_weight = df_out.weight.sum()
# df_out['weight'] = df_out['weight']/tot_weight   #re-normalise again (to add to 1 after points removed)

if show_result:  # plot results
    fig,ax = plt.subplots(1,3,figsize=(12,16))
    deb_raster = np.where(deb_raster<10,np.nan,deb_raster)
    buffered_glacier_raster = np.where(buffered_glacier_raster <10,np.nan,buffered_glacier_raster )
    ax[0].imshow(whole_glacier_utm,cmap='Blues')
    ax[0].imshow(deb_raster[0],cmap = 'Greens')
    levels = np.arange(0, 8875, 100)
    contours = ax[0].contour(buffered_glacier_raster[0], levels=levels, colors='black', linewidths=0.5)
    ax[1].imshow(labels_clean,cmap = 'Blues')
    ax[1].imshow(labels_deb,cmap = 'Greens',vmin=-2)
    ax[2].imshow(deb_th_utm,cmap = 'Greens')

    #now show the points
    for a in ax:
        a.scatter(x_ind_clean,y_ind_clean,c = 'red')
        a.scatter(x_ind_deb,y_ind_deb,c = 'orange')
    ax[0].set_title('elevation')
    ax[1].set_title('clusters')
    ax[2].set_title('debris thickness')
    
    plt.show()

#print(scaled_matrix)



    # cmap = plt.cm.tab20c
    # color = cmap(i)
    #print(color)


#     plt.scatter(inscribed_circ[0],inscribed_circ[1],c='red')
# plt.imshow(labels_grid, cmap='Blues')
# plt.show()



#####write results to file
try:
    
    file_name=f'{root}/Forcing_data/{Id}/coords_out_{Id}.csv'
    file =  open(file_name,'w')
    file.write(f'UTM used,{dst_crs}\n')
    file.close()
    df_out.to_csv(file_name, mode='a',header = True,  index=False)
except FileNotFoundError:
    print('out_file_not_in_list (extract_points_kmeans.py)')

print(f'Glacier points picked in {time.time()-t1} seconds')






















