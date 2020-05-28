# -*- coding: utf-8 -*-
"""
@author: Nadia Bloemendaal, nadia.bloemendaal@vu.nl

For more information, please see 
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. 
Generation of a global synthetic tropical cyclone hazard dataset using STORM. 
Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

This is the STORM model master program

Copyright (C) 2020 Nadia Bloemendaal. All versions released under the [] license.
"""

import numpy as np

#Custom made modules
from SELECT_BASIN import Basins_WMO

from SAMPLE_STARTING_POINT import Startingpoint
from SAMPLE_TC_MOVEMENT import TC_movement
from SAMPLE_TC_PRESSURE import TC_pressure

import os
import sys
dir_path=os.path.dirname(os.path.realpath(sys.argv[0]))
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

import time
start_time=time.time()

#==============================================================================
# Step 1: Define basin and number of years to run
#==============================================================================
#please set basin (EP,NA,NI,SI,SP,WP)
basin='EP'
loop=0 #ranges between 0 and 9 to simulate slices of 1000 years

total_years=1000 #set the total number of years you'd like to simulate

TC_data=[] #This list is composed of: [year,storm number,lat,lon,pressure,wind,rmax,category,Holland B parameter,precipitation,landfall flag]
#==============================================================================
#     Step 2: load grid with weighted genesis counts
#==============================================================================
for year in range(0,total_years):
    storms_per_year,genesis_month,lat0,lat1,lon0,lon1=Basins_WMO(basin) 

    if storms_per_year>0:
            #==============================================================================
            # Step 3: Generate (list of) genesis locations
            #==============================================================================
            lon_genesis_list,lat_genesis_list=Startingpoint(storms_per_year,genesis_month,basin) 
                            
            #==============================================================================
            # Step 4: Generate initial conditions    
            #==============================================================================
            latlist,lonlist,landfalllist=TC_movement(lon_genesis_list,lat_genesis_list,basin)
                    
            TC_data=TC_pressure(basin,latlist,lonlist,landfalllist,year,storms_per_year,genesis_month,TC_data)

TC_data=np.array(TC_data)

np.savetxt(os.path.join(__location__,'STORM_DATA_IBTRACS_'+str(basin)+'_'+str(total_years)+'_YEARS_'+str(loop)+'.txt'),TC_data,fmt='%5s',delimiter=',')