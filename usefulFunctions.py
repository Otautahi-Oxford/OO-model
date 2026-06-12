
"""
A selection of useful functions for orbit calculations. Assumes an orbit
is stored like row in output of orbitSampler.py, i.e.: 
orb['sma'] is the semi-major axis length, 
orb['ecc'] is the eccentricty

nu is the true anomaly, an angular position on the orbit, equal to the 
angle between position of object and its perihelion point, as measured from the Sun.
Negative before perihelion, positive after.

"""

import numpy as np

sidereal_year = 365.256363 # for mu
mu = 4*np.pi**2/(sidereal_year**2) # gravitational parameter, AU^3/JD^2  



def angleBetween(lon1,lat1,lon2,lat2):
    """
    angle between two points on sphere, uses haversine formula for stability
    https://en.wikipedia.org/wiki/Haversine_formula
    """
    dlon = lon2-lon1
    dlat = lat2-lat1
    # 0.9999* would prevent nan due to numerical error when angle=180deg, add if needed
    return 2*np.arcsin(np.sqrt(np.sin(dlat/2)**2 + np.cos(lon1)*np.cos(lat2)*np.sin(dlon/2)**2))
    

def t2nu(t, orb):
    """
    Solves Kepler eqution, based on Farnoccia 2013, fast! 
    
    t can be input as an array of MJDs, this calculates true anomaly (nu) at
    those MJDs
    """
    
    M = np.sqrt(mu/(-orb['sma']**3))*(t - orb['tau']) #mean anomaly
    F = (M/np.sqrt(M**2+0.0001))*np.minimum(np.arcsinh(np.abs(M)/(orb['ecc']-1)),
                                            1+np.log((orb['ecc']+2*np.abs(M))/(np.exp(1)*orb['ecc']-2)))
    # The slightly weird np.sqrt(M**2+0.0001) is to get abs(M) that doesn't go to zero to avoid nans
    
    def deltaKepler(Fk):
        return -(orb['ecc']*np.sinh(Fk) - Fk - M)/(orb['ecc']*np.cosh(Fk) - 1)
    
    maxiter = 100 # 100 is more than enough, at tol=1e-6 highest k I saw was 4
    tol = 1e-6 # rad 
    for k in range(maxiter):
        dF = deltaKepler(F)
        F += dF
        if not np.any(np.abs(dF)>=tol):
            # contemplate speed up of not doing check every iteration.
            # by doing multiple dFs per loop
            break
    # print(k)
    return 2*np.arctan(np.sqrt((orb['ecc']+1)/(orb['ecc']-1))*np.tanh(F/2))


def nuOrb2r(nu, orb):
    """
    Calculates the 3D helicentric position at true anomalies nu 
    (can be input as np array)
    
    I think this only works for hyperbolics
    
    
    Calculates heliocentric distance and argument angle, then performs two 
    rotations, by raan and by inc 
    
    Column stack ensure's correct shape
    """
    dist = -orb['sma']*(orb['ecc']**2-1)/(1+orb['ecc']*np.cos(nu))
    assert np.all(dist>0) # catches when nu is out of range
    arg = orb['argp']+nu
    return dist.reshape(-1,1)*np.column_stack([np.cos(orb['raan'])*np.cos(arg)
              -np.sin(orb['raan'])*np.cos(orb['inc'])*np.sin(arg),
              np.sin(orb['raan'])*np.cos(arg)
              +np.cos(orb['raan'])*np.cos(orb['inc'])*np.sin(arg),
              np.sin(orb['inc'])*np.sin(arg)])


    
    
    
#### EXAMPLE ####

# import orbitSampler

# orb = orbitSampler.sampleOrbits(1,5, 61000, 61500).iloc[0]

# MJDs = 61250 + np.array([-300, -150, 0, 150, 300])
# print("\nFor MJDs: ", MJDs)
# nu = t2nu(MJDs, orb)
# print("\nTrue anomalies are (in radians): ", nu)
# positions = nuOrb2r(nu, orb)
# print("\nwith positions in Cartesian Ecliptic coords: \n", positions)




