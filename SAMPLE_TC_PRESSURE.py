# -*- coding: utf-8 -*-
"""
@author: Nadia Bloemendaal, nadia.bloemendaal@vu.nl

For more information, please see 
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. 
Generation of a global synthetic tropical cyclone hazard dataset using STORM. 
Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

This is the STORM module for simulation of the TC pressure

Copyright (C) 2020 Nadia Bloemendaal. All versions released under the [] license
"""

import numpy as np
from SELECT_BASIN import Basins_WMO
from math import radians, cos, sin, asin, sqrt
from SAMPLE_RMAX import Add_Rmax
import math
import sys
import os
dir_path=os.path.dirname(os.path.realpath(sys.argv[0]))
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
import random

def Calculate_Vmax(Penv,Pc,coef):
    """
    Function to convert pressure to vmax. The equation is based on the empirical wind-pressure relationship 
    (Harper 2002, Atkinson and Holliday 1977)
        
    Input:
        Penv: Environmental pressure (in hPa)
        Pc: Central/minimum pressure in the TC
        a,b: coefficients of empirical wind/pressure relationship
        
    Output:
        Vmax10: 10-minute mean maximum sustained wind speed of the TC (m/s)
    """
    [a,b]=coef
    Vmax10=a*(Penv-Pc)**b

    return Vmax10

def Calculate_Pressure(Vmax10,Penv,coef):
    """
    Convert Vmax to Pressure following the empirical wind-pressure relationship (Harper 2002, Atkinson and Holliday 1977)
    
    Input: 
        Vmax: 10-min mean maximum wind speed in m/s
        Penv: environmental pressure (hPa)
        a,b: coefficients. See Atkinson_Holliday_wind_pressure_relationship.py
    
    Returns:
        Pc: central pressure in the eye
    
    """
    [a,b]=coef
    Pc=Penv-(Vmax10/a)**(1./b)
    
    return Pc

                  
def TC_Category(V):
    """
    Find the category on the Saffir-Simpson Hurricane Wind Scale

    Parameters
    ----------
    V : max wind speed (m/s).

    Returns
    -------
    cat : category (0-5, 0=Tropical Storm).

    """
    if V>=18. and V<33.:
        cat=0
    elif V>=33. and V<43.:
        cat=1
    elif V>=43. and V<50.:
        cat=2
    elif V>=50 and V<58.:
        cat=3
    elif V>=58. and V<70.:
        cat=4
    elif V>=70.:
        cat=5
    
    else:
        cat=-1
    
    return cat

def find_index_pressure(basin,lat,lon,lat0,lon0,lon1):
    """
    Find the index for the coefficient list corresponding to the lon/lat position of the TC

    Parameters
    ----------
    basin : basin.
    lat : latitude position of TC.
    lon : longitude position of TC.
    lat0 : upper left corner latitude of basin.
    lon0 : upper left corner longitude of basin.
    lon1 : upper right corner longitude of basin.

    Returns
    -------
    ind : index.

    """
    base=5
    latindex=np.floor(float(lat-lat0)/base)
    lonindex=np.floor(float(lon-lon0)/base)
    maxlon=(lon1-lon0)/5.
    ind=latindex*maxlon+lonindex
    return ind
    
def PRESSURE_JAMES_MASON(dp,pres,a,b,c,d,mpi):
    """
    Function to calculate the change in pressure

    Parameters
    ----------
    dp : backward change in pressure (dp0, pressure[i]-pressure[i-1]).
    pres : pressure (hPa).
    a,b,c,d : coefficients.
    mpi : mpi in hPa.

    Returns
    -------
    y : forward change in pressure (dp1, pressure[i+1]-pressure[i]).

    """
    if pres<mpi:
        presmpi=0
    else:
        presmpi=pres-mpi
    y=a+b*dp+c*np.exp(-d*presmpi)
    return y

def haversine(lat1,lon1,lat2,lon2):
    """
    function to calculate the distance between two coordinates

    Parameters
    ----------
    lat1 : latitude point 1.
    lon1 : longitude point 1.
    lat2 : latitude point 2.
    lon2 : longitude point 2.

    Returns
    -------
    km : distance in km.

    """
    lon1,lat1,lon2,lat2=map(radians,[lon1,lat1,lon2,lat2])
    dlon=abs(lon1-lon2)
    dlat=abs(lat2-lat1)
    A1=sin(dlat/2)**2.+cos(lat1)*cos(lat2)*sin(dlon/2)**2.
    C2=2.*asin(sqrt(A1))
    r=6371.
    km=C2*r
    
    return km

def decay_after_landfall(lat_landfall,lon_landfall,latlijst,lonlijst,p,coef,Penv):
    """
    Function to calculate the decay after landfall. From Kaplan&DeMaria 1995
    
    Input:
        prev_lat,prev_lon: previous latitude and longitude (one time step before landfall)
        lat_landfall, lon_landfall: latitude and longitude coordinate at landfall
        latlijst,lonlijst: TC track (latitude and longitude list)
        p: pressure at landfall (hPa)
        
        **This is all needed for the wind-pressure relationships**
        coefS: set of coefficients to calculate S
        Penv: Environmental pressure (in hPa)
        C: Forward speed in kt
    
    Output: 
        pressure_decay: central pressure evolution after landfall (in hPa)
        wind_decay: maximum wind speed evolution after landfall (in m/s)
    """
    
    #Coefficients from Kaplan & DeMaria
    #wind is calculated in knots
    C1=0.0109 #kth-2
    D1=-0.0503 #kth-2
    R=0.9
    t0=150
    alpha=0.095 #h-1
    vb=26.7 #kt at R=0.9
    
    v0=Calculate_Vmax(Penv,p,coef)  #wind speed in m/s 
    
    wind_decay=[]
    pressure_decay=[]
    pressure_decay.append(p)
    wind_decay.append(v0)
    
    v0=v0/0.5144444444 #wind speed at landfall, in kt
    D0=1. # km
    
    v=v0
    t=3
    j=1
    pres_landfall=p
    
    while v>35 or j<len(latlijst): #While the storm hasn't dissipated (wind speed lower
        #than 18m/s or, equivalently, 35 kt) or moved out of the basin, proceed
        
        #Distance needs to be greater than 1. So this means that we are going to look at moments AFTER landfall
        try:
            D=haversine(lat_landfall,lon_landfall,latlijst[j],lonlijst[j])
            #D is given in km 
            
            if D==0.: #storm is stationairy at the landfall location
                pressure_decay.append(pres_landfall)
                wind_decay.append(v0*0.5144444) #v in m/s
                j=j+1
                t=t+3
                
            if D>1:
                M=C1*t*(t0-t)
                
                b_KM=D1*t*(t0-t)
                
                C_KM=M*np.log(D/D0)+b_KM
                
                v=vb+(R*v0-vb)*np.exp(-alpha*t)-C_KM # v in kt
                            
                pres_landfall=Calculate_Pressure(v*0.514444,Penv,coef) #v in m/s
                                
                pres_landfall=round(pres_landfall,1)
                
                pressure_decay.append(pres_landfall)
                wind_decay.append(v*0.514444) #v in m/s
                
                if v*0.51444<18.:
                    return pressure_decay, wind_decay
                else:
                    t=t+3 #we have 3-hourly data, so time after landfall = t=t+3
                    j=j+1 #index of remainder of lat/lon list after landfall 
            else:
                v=-100.
        except IndexError: #in this case, the storm has moved out of the basin
            v=-100.
     
    return pressure_decay,wind_decay   

def distance_from_coast(lon,lat,fpath,degree_in_km=111.12):
    """
    Calculate the distance from coast

    Parameters
    ----------
    lon : longitude position of TC.
    lat : latitude position of TC.
    fpath : land/sea mask.
    degree_in_km : The default is 111.12.

    Returns
    -------
    mindist : distance to coast in km.

    """
    if lon>180:
        lon=lon-360.

    D = np.load(fpath,encoding = 'latin1').tolist()

    lons,lats = D['lons'],D['lats']

    dists = np.sqrt((lons-lon)**2+(lats-lat)**2)
    mindist=np.min(dists)*degree_in_km
    return mindist  

def add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,lijst,TC_data,idx):
    """
    Add parameters to the TC data list when TC is dissipated/moved out of basin

    Parameters
    ----------
    pressure_list : array of pressure (hPa).
    wind_list : array of wind (m/s).
    latfull : array of latitude coordinates.
    lonfull : array of longitude coordinates.
    year : year of TC occurrence.
    storm_number : TC storm number.
    month : month of TC occurrence.
    basin : basin.
    landfallfull : array of landfall (0=no 1=yes).
    lijst : dummy indicating the duration of the TC.
    TC_data : existing array of TC data to which will be appended.
    idx : basin idx.

    Returns
    -------
    TC_data : array of TC data.

    """
    rmax_list=Add_Rmax(pressure_list)
    
    x=min(len(landfallfull),len(lijst))    
             
    for l in range(0,x):
        if landfallfull[l]==1.:
            distance=0
        else:
            distance=distance_from_coast(lonfull[l],latfull[l],(os.path.join(dir_path,'coastal_basemap_data.npy')))
            
        category=TC_Category(wind_list[l])
        
        TC_data.append([year,month,storm_number,l,idx,latfull[l],lonfull[l],pressure_list[l],wind_list[l],rmax_list[l],category,landfallfull[l],distance])
    

    return TC_data

def TC_pressure(basin,latlist,lonlist,landfalllist,year,storms,monthlist,TC_data):  
    """
    Calculate TC pressure

    Parameters
    ----------
    basin : basin.
    latlist : array of TC track latitude positions.
    lonlist : array of TC track longitude positions.
    landfalllist : array of TC landfall (0=no 1=yes).
    year : year
    storms : number of storms.
    monthlist : months of TC occurrence.
    TC_data : array of TC data.

    Returns
    -------
    TC_data : array of TC data + new TCs

    """

    basin_name = dict(zip(['EP','NA','NI','SI','SP','WP'],[0,1,2,3,4,5]))
    
    idx=basin_name[basin]

    latidx_penv=np.linspace(90,-90,721)
    lonidx_penv=np.linspace(0,359.75,1440)    
    
    JM_pressure=np.load(os.path.join(__location__,'COEFFICIENTS_JM_PRESSURE.npy')).item()
            
    Genpres=np.load(os.path.join(__location__,'DP0_PRES_GENESIS.npy')).item()
            
    WPR_coefficients=np.load(os.path.join(__location__,'COEFFICIENTS_WPR_PER_MONTH.npy')).item()
    
    Genwind=np.load(os.path.join(__location__,'GENESIS_WIND.npy')).item()

    intlist=[5,3,2,5,5,5]
    
    int_thres=intlist[idx]
       
    s,monthdummy,lat0,lat1,lon0,lon1=Basins_WMO(basin)
    
    wind_threshold=18. #if vmax<18, the storm is a tropical depression and we stop tracking it.

    for storm_number,month,latfull,lonfull,landfallfull in zip(range(0,int(storms)),monthlist,latlist,lonlist,landfalllist):
        i=0
        vmax=0
        count=0
        p=np.nan
        
        #This is the full MSLP field, with lat0=90 deg, lat1=-90 deg, lon0=0 deg, lon1=359.75 deg. len(lat)=721, len(lon)=1440
        Penv_field=np.loadtxt(os.path.join(dir_path,'Monthly_mean_MSLP_'+str(month)+'.txt'))
               
        constants_pressure=JM_pressure[idx][month]
        constants_pressure=np.array(constants_pressure)
        
        coef=WPR_coefficients[idx][month]
        coef=np.array(coef)
        
        p_threshold=min(constants_pressure[:,6])-10.
        
        EP=Genpres[idx][month]
        

        while i<len(latfull): 
            lat,lon,landfall=latfull[i],lonfull[i],landfallfull[i]
            
            lat_dummy=np.abs(latidx_penv-lat).argmin()
            lon_dummy=np.abs(lonidx_penv-lon).argmin()
            
            Penv=Penv_field[lat_dummy,lon_dummy]           
            
            if lat0<=lat<=lat1 and lon0<=lon<=lon1: #make sure we're inside the basin
                
                if ((p<p_threshold) | math.isnan(p)): #something went wrong. start again
                    i=0
                    vmax=0
                    
                if i==0: 
                    vmax=random.choice(Genwind[idx][month])
                    p=Calculate_Pressure(vmax,Penv,coef)
                    
                    
                    pressure_list=[]    
                    wind_list=[]
                    #at genesis, we need to sample the genesis pressure and dp1. This is done basin-wide:
                    
                    [Pmu,Pstd,DP0mu,DP0std,dpmin,dpmax]=EP
                    dp0=np.random.normal(DP0mu,DP0std)

                    dp1=-1.*np.abs(dp0)

                    pressure_list.append(p) 
                    
                    wind_list.append(vmax)
                    
                    i=i+1

                #next: check if the storm makes landfall. In that case, we move to the dissipation-formula                
                if landfall==1: #landfall --> use Kaplan and DeMaria Formula for dissipation
                    if ((p<p_threshold) | math.isnan(p)):
                        print('Landfall',p,p_threshold)
                        i=0
                        vmax=0
                    	
                     
                    elif vmax<wind_threshold or p>Penv: #The storm makes landfall as a tropical depression
                        TC_data=add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,pressure_list,TC_data,idx)
                        i=1000000000000000 
                     
                    else:
                       #calculate the landfall pressure
                        
                      ind=int(find_index_pressure(basin,lat,lon,lat0,lon0,lon1)) #find index for pressure
                      
                      [c0,c1,c2,c3,EPmu,EPstd,mpi]=constants_pressure[ind]
                      y=PRESSURE_JAMES_MASON(dp1,p,c0,c1,c2,c3,mpi) 
                      epsilon=np.random.normal(EPmu,EPstd)
                      dp0=float(y+epsilon)
                      
                      while dp0<dpmin: #if more intensification than seen in the underlying dataset
                          if y-dpmin>EPmu-2.*EPstd: #epsilon should be resampled
                              epsilon=np.random.normal(EPmu,EPstd)
                              dp0=y+epsilon
                          else: #y is already smaller than dpmin. 
                              dp0=np.random.uniform(dpmin,0)
                      
                      while dp0>dpmax: #if more weakening than seen in the underlying dataset
                          if y-dpmax<EPmu+2.*EPstd:
                              epsilon=np.random.normal(EPmu,EPstd)
                              dp0=y+epsilon
                          else:
                              dp0=np.random.uniform(0,dpmax)      
                      
                      if p<mpi:#if pressure has dropped below mpi
                          if dp0<0: #if intensification
                              if count<2: #if intensification has been going on for less than 2 time steps
                                  count=count+1
                              else:
                                  dp0=abs(dp0)
                                      
                      else: #the storm is above mpi
                          count=0  
  
                      p=round(dp0+p,1)
                      dp1=dp0
                      
                      if vmax<wind_threshold or p>Penv: #The storm is no longer a tropical storm
                          #print('Dissipated',len(pressure_list),len(landfallfull))
                          TC_data=add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,pressure_list,TC_data,idx)
                          i=10000000000000000000000000000000
                        
                      else:
                          pressure_list.append(p)
                          
                          vmax=Calculate_Vmax(Penv,p,coef)
                          vmax=round(vmax,1)
                          wind_list.append(vmax)
  
                    if any(c<1 for c in landfallfull[i:]): #check whether the storm moves back over the ocean
                        
                        check_move_ocean=i+np.where(np.array(landfallfull[i:])==0.)[0][0]
                        #storm moves back over open ocean: apply decay function for i till check_move_ocean
                        
                        if check_move_ocean>i+3: #if this is not the case, we're crossing a very small island and no decay function should be used then 
                            
                            decay_pressure, decay_wind = decay_after_landfall(lat,lon,latfull[i:i+check_move_ocean],lonfull[i:i+check_move_ocean],p,coef,Penv)
                            for d in range(len(decay_pressure)):
                                pressure_list.append(decay_pressure[d])
                                wind_list.append(decay_wind[d])

                        #if the storm has decayed before moving back over the ocean:    
                        if wind_list[-1]<wind_threshold:
                            TC_data=add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,pressure_list,TC_data,idx)
                                                        
                            i=10000000000000000000000000.
                            
                        #if the storm has not decayed:
                        else:   
                            dp1=pressure_list[-1]-pressure_list[-2] 
                            p=pressure_list[-1]
                            i=check_move_ocean
                            
                    else: #the storm does not move back over open ocean, so use the decay function until the storm has dissipated
                        decay_pressure, decay_wind = decay_after_landfall(lat,lon,latfull[i:],lonfull[i:],p,coef,Penv)
                        for d in range(len(decay_pressure)):
                            pressure_list.append(decay_pressure[d])
                            wind_list.append(decay_wind[d])
                        
                        #print('Decayed over land',len(pressure_list),len(landfallfull))                         
                        TC_data=add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,pressure_list,TC_data,idx)
                            
                        i=1000000000000
                   
                else: #no landfall
                    if ((p<p_threshold) | math.isnan(p)):
                        print('No landfall',p,p_threshold)
                        i=0
                        vmax=0
                        
                    elif vmax<wind_threshold or p>Penv and i>3: #The storm is no longer a tropical storm
                        #print('Dissipated',len(pressure_list),len(landfallfull))
                        TC_data=add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,pressure_list,TC_data,idx)
                        
                        i=1000000000000000                       
                        
                    else: #apply James-Mason formula to find next change in pressure
                        ind=int(find_index_pressure(basin,lat,lon,lat0,lon0,lon1)) #find index for pressure
                        
                        [c0,c1,c2,c3,EPmu,EPstd,mpi]=constants_pressure[ind] 
                        
                        y=PRESSURE_JAMES_MASON(dp1,p,c0,c1,c2,c3,mpi) 
                        
                        epsilon=np.random.normal(EPmu,EPstd)

                        dp0=float(y+epsilon)  

                        
                        while dp0<dpmin: #if more intensification than seen in the underlying dataset
                            if y-dpmin>EPmu-2.*EPstd: #epsilon should be resampled
                                epsilon=np.random.normal(EPmu,EPstd)
                                dp0=y+epsilon
                            else: #y is already smaller than dpmin. 
                                dp0=np.random.uniform(dpmin,0)
                        
                        while dp0>dpmax: #if more weakening than seen in the underlying dataset
                            if y-dpmax<EPmu+2.*EPstd:
                                epsilon=np.random.normal(EPmu,EPstd)
                                dp0=y+epsilon
                            else:
                                dp0=np.random.uniform(0,dpmax)  
                                
                        if p<mpi:#if pressure has dropped below mpi
                            if dp0<0: #if intensification
                                if count<2: #if intensification has been going on for less than 2 time steps
                                    count=count+1
                                else:
                                    dp0=abs(dp0)

                        else: #the storm is above mpi
                            count=0  
                        
                        if i<int_thres:
                          dp0=-1.*np.abs(dp0)
                        p=round(dp0+p,1)
                        dp1=dp0
                        
                        if vmax<wind_threshold or p>Penv: #The storm is no longer a tropical storm
                          #print('Dissipated',len(pressure_list),len(landfallfull))
                          TC_data=add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,pressure_list,TC_data,idx)
                          i=10000000000000000000000000000000
                        
                        else:
                          pressure_list.append(p)
                          
                          vmax=Calculate_Vmax(Penv,p,coef)
                          vmax=round(vmax,1)
                          wind_list.append(vmax)
                        
                          i=i+1 
                        
            else: #we are outside the basin. Move on to the next storm
                TC_data=add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,pressure_list,TC_data,idx)                
                i=100000000000000000.
                
        if i==len(latfull):
            TC_data=add_parameters_to_TC_data(pressure_list,wind_list,latfull,lonfull,year,storm_number,month,basin,landfallfull,pressure_list,TC_data,idx)
           
    return(TC_data)




