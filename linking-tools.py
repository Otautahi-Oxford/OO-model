"""
Adds tools for providing orbits/rates to test orbit/shift'n'stack methods to 
help with searching for ISOs in surveys

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


import usefulFunctions
import orbitSampler


def get_orbs(N, rh, delta_r, cone_lon, cone_lat, cone_halfangle):
    
    Ntrial = int(1e6) # gives about 100 in cone
    Nsampled = 0
    orbs = pd.DataFrame([])
    while Nsampled<N:
        print(Nsampled)
        trialorbs = orbitSampler.sampleOrbits(Ntrial, rh+delta_r/2, 0, 0)
        # MJD range just 0-0 as we want static population at one time
        
        nu = usefulFunctions.t2nu(0, trialorbs)
        orb_r_vecs = usefulFunctions.nuOrb2r(nu, trialorbs)
        orb_lon, orb_lat, orb_rh = usefulFunctions.vec2spherical(orb_r_vecs)
            
        trialmask = ((rh-delta_r/2<=orb_rh)
                &(orb_rh<=rh+delta_r/2)
                &(usefulFunctions.angleBetween(orb_lon,orb_lat,cone_lon,cone_lat)<=cone_halfangle)
                )
        orbs = pd.concat([orbs, trialorbs])
        Nsampled += np.count_nonzero(trialmask)
    return orbs.iloc[:N]


def get_rates(orbs, rh, delta_r, cone_lon, cone_lat, cone_halfangle):    
    delta_t = 1 #days?
    
    nu1 = usefulFunctions.t2nu(0, orbs)
    nu2 = usefulFunctions.t2nu(delta_t, orbs)
    
    vec1 = usefulFunctions.nuOrb2r(nu1, orbs)
    vec2 = usefulFunctions.nuOrb2r(nu2, orbs)
    
    delta_vec = (vec2-vec1)
    print("median delta in au", np.median(np.sqrt(delta_vec[...,0]**2 + delta_vec[...,1]**2 + delta_vec[...,2]**2)))
    
    # Adding Earth's motion - assumes opposition
    Earth_lons = cone_lon + np.array([-1,1])*2*3.14/365.25 * delta_t/2
    Earth_vec1 = np.array([[np.cos(Earth_lons[0]),np.sin(Earth_lons[0]),0]])
    Earth_vec2 = np.array([[np.cos(Earth_lons[1]),np.sin(Earth_lons[1]),0]])
    
    lon1, lat1, range1 = usefulFunctions.vec2spherical(vec1 - Earth_vec1)
    lon2, lat2, range2 = usefulFunctions.vec2spherical(vec2 - Earth_vec2)
    
    # rate_loncoslat = np.cos(cone_lat)*(lon2 - lon1)/delta_t
    rate_loncoslat = np.where(lon2>lon1, 1,-1)*usefulFunctions.angleBetween(lon1, cone_lat,
                                                  lon2, cone_lat)/delta_t
    # above prevents 0/360 bug
    rate_lat = (lat2 - lat1)/delta_t
    
    rate_range = (range2 - range1)/delta_t
    
    return rate_range, rate_loncoslat, rate_lat


def plot_rates(N, rh, cone_lon,cone_lat, suptext=None):
    
    cone_halfangle = 1.75*3.14/180 # LSST field
        
    delta_r = 5 
    
    orbs = get_orbs(N, rh, delta_r, cone_lon, cone_lat, cone_halfangle)
    
    rate_range, rate_loncoslat, rate_lat = get_rates(orbs, rh, delta_r, cone_lon, cone_lat, cone_halfangle)
    
    cov = np.cov([rate_range, rate_loncoslat, rate_lat])
    print('rate_loncoslat width = ', np.sqrt(cov[1,1])*180/3.14, ' deg/day')
    print('rate_lat width = ', np.sqrt(cov[2,2])*180/3.14, ' deg/day')
    
    plotwidth = 5*(cov[1,1]*cov[2,2])**0.25
    
    Earthv = 3.14*2*1/365.25
    rhv = np.sqrt(887.13/rh)
    
    asPerHour = True
    fig,axs = plt.subplots(ncols=2, figsize=[8,4])
    ax=axs[0]
    ax.scatter(rate_loncoslat*180/3.14 *(60**2/24 if asPerHour else 1), 
               rate_lat*180/3.14*(60**2/24 if asPerHour else 1),
               s=0.1)
    ax.scatter((Earthv-rhv)/(rh-1), 0, marker='x', c='r', s=20)
    ax.set_xlim((np.mean(rate_loncoslat)
                   +plotwidth*np.array([-1,1])
                   )*180/3.14*(60**2/24 if asPerHour else 1))
    ax.set_ylim((np.mean(rate_lat)
                   +plotwidth*np.array([-1,1])
                   )*180/3.14*(60**2/24 if asPerHour else 1))
    ax.set_xlabel('longitude rate cos(latitude) / ' + ('arcsec/hour' if asPerHour else 'deg/day'))
    ax.set_ylabel('latitude rate / ' + ('arcsec/hour' if asPerHour else 'deg/day'))
    ax=axs[1]
    ax.hist(rate_range, bins=np.linspace(*(np.mean(rate_range)+
                                      3*np.sqrt(cov[0,0])*np.array([-1,1])), 31)
            )
    ax.set_xlabel('range rate / au/day')
    fig.suptitle(f'{suptext}rh={rh:.0f} au, cone_lon={cone_lon*180/3.14:.0f} deg, cone_lat={cone_lat*180/3.14:.0f} deg')



#rates example
cosmos_lon, cosmos_lat = 150*3.14/180, -10*3.14/180 #cosmos field

# plot_rates(10,40, cosmos_lon, cosmos_lat, suptext='COSMOS FIELD\n')

N = int(1e3)

for rh in [10,40,100]:
    plot_rates(N, rh, cosmos_lon, cosmos_lat, suptext='COSMOS FIELD\n')









