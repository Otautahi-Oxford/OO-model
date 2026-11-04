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


cone_lon, cone_lat, cone_halfangle = 150*3.14/180, -10*3.14/180, 1.75*3.14/180

# Ra dec of cosmos = 

N = int(1e7)
r_obs = 20
delta_r = 5


MJD_start = 0
MJD_end = 0 # just want static population


def get_orbs(N, r_obs, delta_r, cone_lon, cone_lat, cone_halfangle):
    
    Ntrial = int(1e6) # gives about 100 in cone
    Nsampled = 0
    orbs = pd.DataFrame([])
    while Nsampled<N:
        print(Nsampled)
        trialorbs = orbitSampler.sampleOrbits(Ntrial, r_obs, MJD_start, MJD_end)
        
        nu = usefulFunctions.t2nu(0, trialorbs)
        r_vecs = usefulFunctions.nuOrb2r(nu, trialorbs)
        lon, lat, rh = usefulFunctions.vec2spherical(r_vecs)
            
        trialmask = ((r_obs-delta_r<=rh)
                &(rh<=r_obs)
                &(usefulFunctions.angleBetween(lon,lat,cone_lon,cone_lat)<=cone_halfangle)
                )
        orbs = pd.concat([orbs, trialorbs])
        Nsampled += np.count_nonzero(trialmask)
    return orbs.iloc[:N]


def get_rates(orbs, r_obs, delta_r, cone_lon, cone_lat, cone_halfangle):    
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
    
    return rate_loncoslat, rate_lat


orbs = get_orbs(1000, r_obs, delta_r, cone_lon,cone_lat, cone_halfangle)

rate_loncoslat, rate_lat = get_rates(orbs, r_obs, delta_r, cone_lon, cone_lat, cone_halfangle)

cov = np.cov(rate_loncoslat, rate_lat)
print('rate_loncoslat width = ', np.sqrt(cov[0,0])*180/3.14, ' deg/day')
print('rate_lat width = ', np.sqrt(cov[1,1])*180/3.14, ' deg/day')

asPerHour = False
fig,ax = plt.subplots()
ax.scatter(rate_loncoslat*180/3.14 *(60**2/24 if asPerHour else 1), 
           rate_lat*180/3.14*(60**2/24 if asPerHour else 1),
           s=0.1)
ax.set_xlabel('longitude rate cos(latitude) / ' + ('arsec/hour' if asPerHour else 'deg/day'))
ax.set_ylabel('latitude rate / ' + ('arsec/hour' if asPerHour else 'deg/day'))
ax.set_title(f'r_obs={r_obs:.1f} au, cone_lon={cone_lon*180/3.14:.1f} deg, cone_lat={cone_lat*180/3.14:.1f} deg')



