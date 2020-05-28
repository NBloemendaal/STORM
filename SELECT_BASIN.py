# -*- coding: utf-8 -*-
"""
@author: Nadia Bloemendaal, nadia.bloemendaal@vu.nl

For more information, please see 
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. 
Generation of a global synthetic tropical cyclone hazard dataset using STORM. 
Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

This is the STORM module for simulation of genesis month, frequency, and basin boundaries

Copyright (C) 2020 Nadia Bloemendaal. All versions released under the [] license
"""
import numpy as np
import random
import os
import sys
dir_path=os.path.dirname(os.path.realpath(sys.argv[0]))
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

#Basin indices: 
# 0 = EP = Eastern Pacific
# 1 = NA = North Atlantic
# 2 = NI = North Indian
# 3 = SI = South Indian
# 4 = SP = South Pacific
# 5 = WP = Western Pacific

def Genesis_month(idx,storms):
    """
    Sample the genesis months for every TC
    Parameters
    ----------
    idx : basin index (0=EP 1=NA 2=NI 3=SI 4=SP 5=WP).
    storms : number of TCs.

    Returns
    -------
    monthall : list of all genesis months.

    """
    monthlist=np.load(os.path.join(__location__,'GENESIS_MONTHS.npy')).item()
    
    monthall=[]
    for i in range(0,storms):
        monthall.append(np.random.choice(monthlist[idx]))
    
    return monthall


    
def Storms(idx): 
    """
    Sample the number of TC formations in a given year

    Parameters
    ----------
    idx : basin index (0=EP 1=NA 2=NI 3=SI 4=SP 5=WP).

    Returns
    -------
    s : number of storms.

    """
    mu_list=np.loadtxt(os.path.join(__location__,'POISSON_GENESIS_PARAMETERS.txt'))
    #mu_list has the shape [EP,NA,NI,SI,SP,WP]
    
    mu=mu_list[idx]

    poisson=np.random.poisson(mu,10000)
    s=random.choice(poisson)
    return s

def Basins_WMO(basin):
    """
    Basin definitions

    Parameters
    ----------
    basin : basin.

    Returns
    -------
    s : number of storms.
    month : list of genesis months.
    lat0 : lower left corner latitude.
    lat1 : upper right corner latitude.
    lon0 : lower left corner longitude.
    lon1 : upper right corner longitude.

    """
    #We follow the basin definitions from the IBTrACS dataset, but with lat boundaries set at 60 N/S
    #The ENP/AO border will be defined in the algorithm later. 
    
    basins=['EP','NA','NI','SI','SP','WP']
    basin_name = dict(zip(basins,[0,1,2,3,4,5]))
    idx=basin_name[basin]
    
    s=Storms(idx)
    
    month=Genesis_month(idx,s)
  
    if idx==0: #Eastern Pacific
        lat0,lat1,lon0,lon1=5,60,180,285
    if idx==1: #North Atlantic
        lat0,lat1,lon0,lon1=5,60,255,359
    if idx==2: #North Indian
        lat0,lat1,lon0,lon1=5,60,30,100
    if idx==3: #South Indian
        lat0,lat1,lon0,lon1=-60,-5,10,135
    if idx==4: #South Pacific
        lat0,lat1,lon0,lon1=-60,-5,135,240
    if idx==5: #Western Pacific
        lat0,lat1,lon0,lon1=5,60,100,180
        
    return s,month,lat0,lat1,lon0,lon1 


