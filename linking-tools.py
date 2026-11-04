"""
Adds tools for providing orbits/rates to test orbit/shift'n'stack methods so
that ISOs can be searched for in data


To give to THOR (does test orbits in N-bod code, provide distribution entering sphere)

and Luis Slaazar-Manzano for LSST DDF shift'n'stack rates at r

and Matt Holman/ari heinze Heliolink/poindexter



Further thoughts: For rates, even if normal shift and stack doesn't use distance,
rates will change with distance and therefore mag limit, so do rates as functio of distance.
May need to do some maths?

For TNO distances, can probably do from persepctive of Sun rather than from Earth?
Not entirely, as parallax motion matters, but could do cone from Sun then add 
add Earth's motion effect afterwards?

More maths - 
TO get rates at given rh, take observable sphere with radius rh, as increasing
beyond this will only add objects with q>rh and these never reach expected distance
Still very inefficient

Can I increase sampling into cone+right distance?
Can't restrict to just perihelia at rh as q<<rh can still then pass through 
cone at rh

If there alternate method, eg. pick position, trace orbits passing thoguh that 
point, evaluate velocity distribution at vinfs they end up at?






TODO

Rates along cones, for LSST DDFs. To be most useful 
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
    
    lon1, lat1, rh1 = usefulFunctions.vec2spherical(vec1)
    lon2, lat2, rh2 = usefulFunctions.vec2spherical(vec2)
    
    # rate_loncoslat = np.cos(cone_lat)*(lon2 - lon1)/delta_t
    rate_loncoslat = np.where(lon2>lon1, 1,-1)*usefulFunctions.angleBetween(lon1, cone_lat,
                                                  lon2, cone_lat)/delta_t
    # above prevents 0/360 bug
    rate_lat = (lat2 - lat1)/delta_t
    
    rate_rh = (rh2 - rh1)/delta_t
    
    return rate_rh, rate_loncoslat, rate_lat


def plot_rates(N, rh, cone_lon,cone_lat, suptext=None):
    
    cone_halfangle = 1.75*3.14/180 # LSST field
        
    delta_r = 5 
    
    orbs = get_orbs(N, rh, delta_r, cone_lon, cone_lat, cone_halfangle)
    
    rate_rh, rate_loncoslat, rate_lat = get_rates(orbs, rh, delta_r, cone_lon, cone_lat, cone_halfangle)
    
    cov = np.cov([rate_rh, rate_loncoslat, rate_lat])
    print('rate_loncoslat width = ', np.sqrt(cov[1,1])*180/3.14, ' deg/day')
    print('rate_lat width = ', np.sqrt(cov[2,2])*180/3.14, ' deg/day')
    
    plotwidth = 5*(cov[1,1]*cov[2,2])**0.25
    
    asPerHour = True
    fig,axs = plt.subplots(ncols=2, figsize=[8,4])
    ax=axs[0]
    ax.scatter(rate_loncoslat*180/3.14 *(60**2/24 if asPerHour else 1), 
               rate_lat*180/3.14*(60**2/24 if asPerHour else 1),
               s=0.1)
    ax.set_xlim((np.mean(rate_loncoslat)
                   +plotwidth*np.array([-1,1])
                   )*180/3.14*(60**2/24 if asPerHour else 1))
    ax.set_ylim((np.mean(rate_lat)
                   +plotwidth*np.array([-1,1])
                   )*180/3.14*(60**2/24 if asPerHour else 1))
    ax.set_xlabel('longitude rate cos(latitude) / ' + ('arcsec/hour' if asPerHour else 'deg/day'))
    ax.set_ylabel('latitude rate / ' + ('arcsec/hour' if asPerHour else 'deg/day'))
    ax=axs[1]
    ax.hist(rate_rh, bins=np.linspace(*(np.mean(rate_rh)+
                                      3*np.sqrt(cov[0,0])*np.array([-1,1])), 31)
            )
    ax.set_xlabel('rh rate / au/day')
    fig.suptitle(f'{suptext}rh={rh:.0f} au, cone_lon={cone_lon*180/3.14:.0f} deg, cone_lat={cone_lat*180/3.14:.0f} deg')



#rates example
cosmos_lon, cosmos_lat = 150*3.14/180, -10*3.14/180 #cosmos field

N = int(1e3)

for rh in [10,40,100]:
    plot_rates(N, rh, cosmos_lon, cosmos_lat, suptext='COSMOS FIELD\n')

