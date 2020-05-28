# -*- coding: utf-8 -*-
"""
@author: Nadia Bloemendaal, nadia.bloemendaal@vu.nl

For more information, please see 
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. 
Generation of a global synthetic tropical cyclone hazard dataset using STORM. 
Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

This is the STORM module for simulation of the TC genesis location

Copyright (C) 2020 Nadia Bloemendaal. All versions released under the GNU General Public License v3.0
"""
#==============================================================================
# Variables used in this code:
    #grid_copy: copy of grid, the 1x1 deg matrix with genesis counts
    #no_storms: number of storms per year
    #lon_genesis_list,lat_genesis_list: lists of lon and lat points of genesis
    #value: value in copy_grid
    #weighted_list_index: list of indices (weighted) corresponding to the value in grid_copy. Example: if the value is 5, the correspondng index is listed 5 times.
    #idx: randomly chosen index of weighted_list_index
    #row,col: row and column of grid_copy corresponding to the chosen index
    #lat_pert,lon_pert: random value between 0 and 1, indicating the genesis location. This is NOT the actual longitude latitude, but indexed based on the coordinates of 'grid'
#==============================================================================

import numpy as np
import random
import os
import sys
dir_path=os.path.dirname(os.path.realpath(sys.argv[0]))
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

from SELECT_BASIN import Basins_WMO
def Check_EP_formation(lat,lon): 
    """
    Check if formation is in Eastern Pacific (this should be inhibited if basin==NA)
    Parameters
    ----------
    lat : latitude coordinate of genesis.
    lon : longitude coordinate of genesis

    Returns
    -------
    l : 1=yes (formation in EP),0=no (no formation in EP).

    """    
    if lat<=60. and lon<260.:
        l=1
    elif lat<=17.5 and lon<270.:
        l=1
    elif lat<=15. and lon<275.:
        l=1
    elif lat<=10. and lon<276.:
        l=1
    elif lat<=9. and lon<290.:
        l=1
    else:
        l=0
    return l
    
    
def Check_NA_formation(lat,lon): 
    """
    Check if formation is in North Atlantic (this should be inhibited if basin==EP)
    Parameters
    ----------
    lat : latitude coordinate of genesis
    lon : longitude coordinate of genesis

    Returns
    -------
    l : 1=yes (formation in NA) 0=no (no formation in NA).

    """
    if lat<=60. and lat>17.5 and lon>260.:
        l=1
    elif lat<=17.5 and lat>15. and lon>270.:
        l=1
    elif lat<=15. and lat>10 and lon>275.:
        l=1
    elif lat<=10. and lon>276.:
        l=1
    else:
        l=0
    return l

def Check_if_landfall(lat,lon,basin,land_mask):
    """
    Parameters
    ----------
    lat : latitude position of TC
    lon : longitude position of TC
    lat1 : upper left corner latitude coordinate of basin
    lon0 : upper left corner longitude coordinate of basin
    land_mask : land-sea mask

    Returns
    -------
    l : 0=no landfall, 1=landfall

    """
    s,monthdummy,lat0_WMO,lat1_WMO,lon0_WMO,lon1_WMO=Basins_WMO(basin)
    
    x=int(10*(lon-lon0_WMO))
    y=int(10*(lat1_WMO-lat))
    l=land_mask[y,x]   
    return l

def Startingpoint(no_storms,monthlist,basin):
    """
    This function samples the genesis locations of every TC in a given year

    Parameters
    ----------
    no_storms : number of TCs in given year
    monthlist : months in which TCs were formed.
    basin : basin.

    Returns
    -------
    lon_coordinates : list of longitude coordinates of genesis locations
    lat_coordinates : list of latitude coordinates of genesis locations

    """
    basins=['EP','NA','NI','SI','SP','WP']
    basin_name = dict(zip(basins,[0,1,2,3,4,5]))
    idx=basin_name[basin]
    
    lon_coordinates=[]
    lat_coordinates=[]
    weighted_list_index=[]
    weighted_list=[]
    
        
    s,monthdummy,lat0,lat1,lon0,lon1=Basins_WMO(basin)
  
    land_mask=np.loadtxt(os.path.join(dir_path,'Land_ocean_mask_'+str(basin)+'.txt'))

    for month in monthlist:
        
        grid_copy=np.loadtxt(os.path.join(dir_path,'GRID_GENESIS_MATRIX_'+str(idx)+'_'+str(month)+'.txt'))
    
    
        grid_copy=np.array(grid_copy)
        grid_copy=np.round(grid_copy,1)   
        #==============================================================================
        # Make a list with weighted averages. The corresponding grid-index is calculated as len(col)*row_index+col_index
        #============================================================================== 
    
        for i in range(0,len(grid_copy[:,0])):
            for j in range(0,len(grid_copy[0,:])):
                if grid_copy[i,j]>-1:
                    value=int(10*grid_copy[i,j])
                else:
                    value=0
                if value>0.:
                    for k in range(0,value):
                        weighted_list.append(value)
                        weighted_list_index.append(i*(len(grid_copy[0,:])-1)+j)    
                                
        #==============================================================================
        # The starting longitude and latitude coordinates 
        #==============================================================================
        var=0
        while var==0:
            idx0=random.choice(weighted_list_index)

            row=int(np.floor(idx0/(len(grid_copy[0,:])-1)))
            col=int(idx0%(len(grid_copy[0,:])-1))
            lat_pert=random.uniform(0,0.94) #take 0.94 to make sure the randomly selected point is still inside the grid box after rounding off. 
            lon_pert=random.uniform(0,0.94)
            lon=lon0+round(col+lon_pert,1)
            lat=lat1-round(row+lat_pert,1)
            
            
            if lon<lon1 and lat<lat1:
                check=Check_if_landfall(lat,lon,basin,land_mask) #make sure the coordinate isn't on land
                
                if basin=='EP':
                    check=Check_NA_formation(lat,lon)
                if basin=='NA':
                    check=Check_EP_formation(lat,lon)
                    
                if check==0:            
                    lon_coordinates.append(lon)
                    lat_coordinates.append(lat)
                    var=1  
                else:
                    var=0
    
        if len(lon_coordinates)==no_storms:       
            return lon_coordinates,lat_coordinates
