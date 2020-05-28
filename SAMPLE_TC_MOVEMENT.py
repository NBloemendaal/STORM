# -*- coding: utf-8 -*-
"""
@author: Nadia Bloemendaal, nadia.bloemendaal@vu.nl

For more information, please see 
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. 
Generation of a global synthetic tropical cyclone hazard dataset using STORM. 
Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

This is the STORM module for simulation of the TC track

Copyright (C) 2020 Nadia Bloemendaal. All versions released under the GNU General Public License v3.0
"""
import numpy as np
from SELECT_BASIN import Basins_WMO
import os
import sys
dir_path=os.path.dirname(os.path.realpath(sys.argv[0]))
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

def find_lat_index_bins(basin,lat):
    """
    find index of latitude bin in coefficients list

    Parameters
    ----------
    basin : basin.
    lat : latitude position.

    Returns
    -------
    latindex : index of bin.

    """
    s,monthdummy,lat0,lat1,lon0,lon1=Basins_WMO(basin)
    base=5
    latindex=np.floor(float(lat-lat0)/base)
    return latindex

def LAT_JAMES_MASON(dlat,lat,a,b,c):
    """
    Parameters
    ----------
    dlat : backward change in latitude (dlat0, lat[i]-lat[i-1]).
    lat : latitude.
    a,b,c : coefficients

    Returns
    -------
    y : forward change in latitude (dlat1, lat[i+1]-lat[i]).

    """
    y=a+b*dlat+c/lat
    return y

def LON_JAMES_MASON(dlon,a,b):
    """
    Parameters
    ----------
    dlon : backward change in longitude (dlon0, lon[i]-lon[i-1])
    a,b : coefficients

    Returns
    -------
    y : forward change in longitude (dlon1, lon[i+1]-lon[i])

    """
    y=a+b*dlon
    return y

def Check_if_landfall(lat,lon,lat1,lon0,land_mask):
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
    x_coord=int(10*(lon-lon0))
    y_coord=int(10*(lat1-lat))
    l=land_mask[y_coord,x_coord]  

    return l
  
def TC_movement(lon_genesis_list,lat_genesis_list,basin): 
    """
    Parameters
    ----------
    lon_genesis_list : list of longitudinal positions of genesis in a year
    lat_genesis_list : list of latitudinal positions of genesis in a year
    basin : basin
    
    Returns
    ---------
    latall : all latitude positions of the eye of TC for every TC in a year
    lonall : all longitude positions of the eye of TC for every TC in a year
    landfallall : landfall (0=no 1=yes) along the track for every TC in a year
    """

    basins=['EP','NA','NI','SI','SP','WP']
    basin_name = dict(zip(basins,[0,1,2,3,4,5]))
    idx=basin_name[basin]
    
    constants=np.loadtxt(os.path.join(__location__,'JM_LONLATBINS_'+str(idx)+'.txt'))
    
    land_mask=np.loadtxt(os.path.join(__location__,'Land_ocean_mask_'+str(basin)+'.txt'))
      
    
    s,monthdummy,lat0,lat1,lon0,lon1=Basins_WMO(basin)
    latall=[]
    lonall=[]
    landfallall=[]
   
    for lat_genesis,lon_genesis in zip(lat_genesis_list,lon_genesis_list):                 
        #load data for longitude/latitude
        latlijst=[]
        lonlijst=[] 
        landfalllijst=[]
        lat=lat_genesis
        lon=lon_genesis
        latlijst.append(lat)
        lonlijst.append(lon)
        landfall=Check_if_landfall(lat,lon,lat1,lon0,land_mask) #1=landfall 0=no landfall  
        landfalllijst.append(landfall)
        
        var=0 #var is the 'stop-parameter'. The track generation ends when the storm moves over land or out of the basin
        while var==0:
            ind=int(find_lat_index_bins(basin,lat))
            #constants values for latitude/longitude   
            [a0,a1,b0,b1,b2,Elatmu,Elatstd,Elonmu,Elonstd,Dlat0mu,Dlat0std,Dlon0mu,Dlon0std]=constants[ind]    
            
            if len(latlijst)==1: #if this is the first time step after genesis, we need to sample the first change in lon/lat/pressure
                dlat0=np.random.normal(Dlat0mu,Dlat0std,1)
                dlon0=np.random.normal(Dlon0mu,Dlon0std,1)              
            
            dlat1=LAT_JAMES_MASON(dlat0,lat,b0,b1,b2)
            if basin=='SP' or basin=='SI':
                if lat>-10.:                    
                    dlat1=float(dlat1-np.abs(np.random.normal(Elatmu,Elatstd)))
                else:
                    dlat1=float(dlat1+np.random.normal(Elatmu,Elatstd))
                    
            else:
                if lat<10.:
                    dlat1=float(dlat1+np.abs(np.random.normal(Elatmu,Elatstd)))
                else:
                    dlat1=float(dlat1+np.random.normal(Elatmu,Elatstd))
                    

            dlon1=LON_JAMES_MASON(dlon0,a0,a1)
            epsilon=np.random.normal(Elonmu,Elonstd)
            dlon1=float(dlon1+epsilon)
            if np.abs(lat)>=45:
              if dlon1<0.:
                dlon1=0
                            
            lat=round(dlat1+lat,1)
            lon=round(dlon1+lon,1)
            
            
            dlat0=dlat1 
            dlon0=dlon1
            
                    
                                      
            if lat<=lat1-0.1 and lat>lat0 and lon<=lon1-0.1 and lon>lon0: #if storm is still inside the domain. The 0.1's are added to make sure we don't end up at the upper/rightmost edge
                latlijst.append(lat)
                lonlijst.append(lon) 
                landfall=Check_if_landfall(lat,lon,lat1,lon0,land_mask) #1=landfall 0=no landfall  
                landfalllijst.append(landfall)
                    
            else: #the storm has moved out of the domain. Stop the algorithm
                var=1             
                latall.append(latlijst) #all locations of all storms (seperate entries) are contained in these two lists
                lonall.append(lonlijst)
                landfallall.append(landfalllijst)
            

    return(latall,lonall,landfallall)           

