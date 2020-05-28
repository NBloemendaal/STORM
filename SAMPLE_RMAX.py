# -*- coding: utf-8 -*-
"""
@author: Nadia Bloemendaal, nadia.bloemendaal@vu.nl

For more information, please see 
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. 
Generation of a global synthetic tropical cyclone hazard dataset using STORM. 
Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

This is the STORM module for the calculation of radius to maximum winds

Copyright (C) 2020 Nadia Bloemendaal. All versions released under the GNU General Public License v3.0.
"""

import numpy as np
import os
import sys
dir_path=os.path.dirname(os.path.realpath(sys.argv[0]))
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

def sample_rmax(p,rmax_pres):
    if p>940.:
        r=np.random.choice(rmax_pres[2],1)
    elif p<=940. and p>920.:
        r=np.random.choice(rmax_pres[1],1)
    else:
        r=np.random.choice(rmax_pres[0],1)
    return float(r)


def Add_Rmax(pressure):    
    rmax_pres=np.load(os.path.join(__location__,'RMAX_PRESSURE.npy'),allow_pickle=True).item()
        
    #sample rmax at genesis
    rmaxlist=[]
    rgenesis=sample_rmax(pressure[0],rmax_pres)
    rmaxlist.append(rgenesis)   
     
    #sample rmax at minimum pressure    
    p_min=np.min(pressure)
    ind=pressure.index(p_min)
    ind=int(ind)
    
    rmin=sample_rmax(p_min,rmax_pres)
        
    if rmin<rgenesis and ind>0: #if the new radius is smaller than the old one, AND the min pressure is not at genesis
        for i in range(1,ind+1):
            radius=i*(rmin-rmaxlist[0])/ind+rmaxlist[0]
            rmaxlist.append(radius)
    else: #the sampled radius is larger OR min pressure is at genesis
        for i in range(1,ind+1):
            rmaxlist.append(rgenesis)
     
    rind=rmaxlist[-1]
    #sample radius at dissipation
    rdis=sample_rmax(pressure[-1],rmax_pres)
    if rdis>rind: #if radius is larger than the last added radius to the rmax_list:
        for i in range(ind+1,len(pressure)):
            radius=(rdis-rind)/(len(pressure)-1-ind)*i+rdis-(len(pressure)-1)*(rdis-rind)/(len(pressure)-1-ind)
            rmaxlist.append(radius)
    else:
        for i in range(ind+1,len(pressure)):
            rmaxlist.append(rind)
    
    return rmaxlist
