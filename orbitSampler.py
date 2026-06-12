import numpy as np
from scipy.optimize import newton
import astropy.coordinates as coords
import astropy.units as u
import pandas as pd


import matplotlib as mpl
cmap2d = mpl.colormaps['Blues']

"""
Implements orbit sampling method of Dorsey et al. (2025) Appendix A.
http://doi.org/10.3847/PSJ/adf8ca
https://ui.adsabs.harvard.edu/abs/2025PSJ.....6..214D/abstract
https://arxiv.org/abs/2502.16741

Velocities in velDist.csv are in km/s. For internal consistency, in this code
all units in AU, days, AU/day, and radians. These are the units in the output 
of sampleOrbits(), apart from the U,V, and W components of the pre-encounter 
velocity, which are converted back into km/s. Orbital angles are relative to 
Ecliptic plane.

Commonly used parameters:
A = mu/vinf**2 is the absolute value of an ISO's (negative) semimajoraxis
Bmax = np.sqrt(r_obs**2 + 2*A*r_obs) is the maximum impact parameter an ISO
    with velocity vinf can have and still pass within sphere of radius r_obs


sampleOrbits() is the main function.



"""


gaiadf = pd.read_csv('velDist.csv', index_col=0)


# globals

sidereal_year = 365.256363 # for mu
tropical_year = 365.2422 # for survey length

mu = 4*np.pi**2/(sidereal_year**2) # gravitational parameter, AU^3/JD^2  
vconv = 86400/(149597870.7) # velocity unit conversion, (AU per JD)/(km per s)

vinf_table = vconv*np.sqrt(gaiadf['U']**2 + gaiadf['V']**2 + gaiadf['W']**2).to_numpy()
# AU/JD
# NB separating vinf at this point allows gaiadf to be kept as-is 
# (velocity units of km/s) and everything that follows as AU/day and as np arrays

beta=1
ISO_weights = (gaiadf['oneOverESF']*10**(beta*gaiadf['mh'])).to_numpy()
ISO_weights /= ISO_weights


def calc_number(r_obs, T, n0):
    """Calculates the number of ISOs that spend any time within an observable 
    sphere of radius r_obs (au) sphere over survey of length T (days), accounting 
    for gravitational focussing. Simply proportional to the unfocussed 
    density n0. Equal to the integral of Equation (A5) in Dorsey et al. (2025)
    """
    A = mu/vinf_table**2 #AU
    Bmax = np.sqrt(r_obs**2+2*A*r_obs) #AU
    weights = ISO_weights
    weights /= weights.sum() # must be normalised
    number = 2*np.pi*n0*(weights*(
        2*Bmax**3/3 + A**3*np.log(1+(r_obs+Bmax)/A) - (r_obs+A)*A*Bmax + vinf_table*T*Bmax**2/2)).sum(axis=-1)
    return number


def drawVels(N, r_obs, T):
    """
    Draws N samples from the velocity distribution of ISOs passing through a 
    sphere of radius r_obs (AU) over time T (days), given by Equation (A5)
    in Dorsey et al. (2025).
    Returns array of indices corresponding to rows of gaiadf.
    """
    A = mu/vinf_table**2 #AU
    Bmax = np.sqrt(r_obs**2+2*A*r_obs) #AU
    Ap = A/Bmax
    
    marginalised_weights = ISO_weights*Bmax**3*(2/3 
                            - Ap*(np.sqrt(Ap**2+1)- Ap**2*np.log((1+np.sqrt(1+Ap**2))/Ap))
                            + vinf_table*T/(2*Bmax))
    # weights that give U,V,W distribution of velocities marginalised over B, phi and tau
    csw = np.append([0], np.cumsum(marginalised_weights/marginalised_weights.sum()))
    randos = np.random.uniform(size=N) # uniform draws from [0,1]
    indx = np.digitize(randos, csw)-1
    return indx


def Bcdf(x, A, r_obs, T, offset=0.0):
    """
    Velocity-dependent cumulative distribtuion function of impact parameter B, 
    equal to Eq (A4) divided by Eq (A5) in Dorsey et al. (2025).
    Inputs are x=B/Bmax, A = mu/vinf**2
    Offset is for root finding in drawBs().
    """
    vinf = np.sqrt(mu/A)
    Bmax = np.sqrt(r_obs**2+2*A*r_obs)
    Ap = A/Bmax
    cdf = (
        2*(1-(1-x**2)**(3/2))/3 - Ap*(
            (Ap**2+x**2)*np.log((np.sqrt(1-x**2)+np.sqrt(1+Ap**2))/np.sqrt(Ap**2+x**2))
            - Ap**2*np.log((1+np.sqrt(1+Ap**2))/Ap) + np.sqrt(Ap**2+1)*(1-np.sqrt(1-x**2))
            )
        + vinf*T*x**2/(2*Bmax)
          )/(2/3 - Ap*(np.sqrt(Ap**2+1)-Ap**2*np.log((1+np.sqrt(1+Ap**2))/Ap)) + vinf*T/(2*Bmax)
        )
    return cdf - offset
    

def drawBs(indx, r_obs, T):
    """
    Input array of indices corresponding to velocities in gaiadf (i.e. output of 
    drawVels). Returns array of corresponding impact parameters B drawn from 
    the velocity-dependent cumulative distribution function Bcdf(), with method
    of Dorsey et al. (2025) Appendix A.
    
    Inverts cdf by root-finding.
    
    #TODO INVESTIGATE
    Exactly one input datapoint always makes it fail, i=57894 (old table) which is also the lowest speed and
    quite highly weighted. Could just cut it, but why does it fail?
    A = 10751
    Bmax = 415
    Maybe it is in finite difference formula it goes funny? For now, just deleted it from table
    
    #TODO INVESTIGATE
    Approxiately 1/100000 get result=True but x=nan.
    This happends when rando is too high (> about 0.99999), so just rerunning them
    """
    A = mu/vinf_table[indx]**2
    Bmax = np.sqrt(r_obs**2+2*r_obs*A)
    randos = np.random.uniform(0, 0.9999, size=len(indx)) # uniform draws from [0,1]
    # upperlimit of 0.9999 stops weird error
    x, *results = newton(Bcdf, args=(A,r_obs,T,randos), x0 = randos, full_output=True)
    # remember Bcdf takes input as x=B/Bmax, not B
    
    #checking
    # fails = np.nonzero(~results[0])[0]
    assert np.all(results[0])
    i_nans = np.nonzero(np.isnan(x))[0]
    print(f'\nRerunning nans: {i_nans}')
    for i in i_nans:
        print(i)
        print(A[i])
        print(randos[i])
        print(results[0][i])
        print(results[1][i])
        newrando = np.random.uniform(0,1)
        x[i], *newresult = newton(Bcdf, args=(A[i],r_obs,T,newrando), x0 = newrando, full_output=True)
        print(x[i], newrando, newresult)
    
    assert np.all(~np.isnan(x))
    assert np.all(x>=0)
    assert np.all(x<1)
    return Bmax*x


def t_res_func(x, A, r_obs):
    """
    Residence time an ISO spends within observable sphere of radius r_obs.
    Equal to Equation (A1) of Dorsey et al. (2025)
    Inputs are x=B/Bmax, A = mu/vinf**2 
    """
    vinf = np.sqrt(mu/A)
    Bmax = np.sqrt(r_obs**2+2*A*r_obs)
    Ap = A/Bmax
    ts = (2*Bmax/vinf)*(np.sqrt(1-x**2) - Ap*np.log((np.sqrt(1-x**2)+np.sqrt(1+Ap**2))/np.sqrt(Ap**2+x**2)))
    return ts


def v2lonlat(U,V,W):
    """
    gets ECLIPTIC longitude and latitude of heading from GALACTIC velocity vector
    """
    l = np.arctan2(V,U)
    b = np.arctan(W/np.sqrt(U**2 + V**2))
    crds = coords.SkyCoord(l=l*u.rad,b=b*u.rad, frame='galactic').barycentricmeanecliptic
    return (crds.lon.to(u.rad).value, crds.lat.to(u.rad).value)


def sampleOrbits(N, r_obs, MJD_start, MJD_end, savepath=None):
    """
    Makes N draws from the distribution of ISOs passing within an observable 
    sphere of radius r_obs within given MJD range using method of Dorsey et al. 
    (2025) Appendix A.
    """
    T_survey = MJD_end - MJD_start
    indx = drawVels(N, r_obs, T_survey)
    vinf = vinf_table[indx]
    A = mu/vinf**2 #useful quantity for integrals, equal to -1*semi-major axis (therefore is positive)
    Bmax = np.sqrt(r_obs**2+2*A*r_obs)
    B = drawBs(indx, r_obs, T_survey) # impact parameter, AU
    phi = np.random.uniform(0, 2*np.pi, N) # measured from ECLIPTIC +ve z-axis, anticlockwise
    t_res = t_res_func(B/Bmax, A, r_obs) # residence times, Julian days
    
    tau = MJD_start + np.random.uniform(-t_res/2, T_survey + t_res/2)
    # MJD of perihelion
    
    #### Orbital Elements ####
    #now converts my parameters (U, V, W, B, phi) 
    # to classical orbit elements (a, e, i, RAAN, argp)
    
    sma = -A # semi-major axis, AU
    lon, lat = v2lonlat(gaiadf['U'].to_numpy()[indx],
                        gaiadf['V'].to_numpy()[indx],
                        gaiadf['W'].to_numpy()[indx])
    
    v_unit_vec = np.array([np.cos(lon)*np.cos(lat),
                           np.sin(lon)*np.cos(lat),
                                       np.sin(lat)]).transpose(1,0)
    B_unit_vec = np.array([-np.cos(lon)*np.sin(lat)*np.cos(phi) -np.sin(lon)*np.sin(phi),
                           -np.sin(lon)*np.sin(lat)*np.cos(phi) +np.cos(lon)*np.sin(phi),
                           np.cos(lat)*np.cos(phi)]).transpose(1,0)
    # transposes here make 0th axis the different samples, and the 1st axis
    #the element of the vector. This is so matrix functions work
    print(f'CHECK: {np.abs(np.sum(B_unit_vec*v_unit_vec,axis=1)).max()} should be 0')
    #checks perpendicular
    
    h_unit_vec = np.cross(B_unit_vec,v_unit_vec)
    print(f'CHECK: {np.abs(np.sum(h_unit_vec**2, axis=1)-1).max()} should be 0')
    #checks normalised
    
    N_unit_vec = np.cross(np.array([0,0,1]), h_unit_vec) # ascending node vector
    N_unit_vec /= np.sqrt(np.sum(N_unit_vec**2, axis=1)).reshape(-1,1)
    print(f'CHECK: {np.abs(np.sum(N_unit_vec**2, axis=1)-1).max()} should be 0')
    #checks normalised
    
    ecc_vec = np.cross(v_unit_vec,np.cross(B_unit_vec,v_unit_vec))*(B/A).reshape(-1,1) + v_unit_vec
    ecc = np.sqrt(np.sum(ecc_vec**2, axis=1))
    ecc_unit_vec = ecc_vec/ecc.reshape(-1,1)
    print(f'CHECK: {ecc.min()} should be >=1')
    print(f'CHECK: {np.abs(np.sqrt(ecc**2-1)*A - B).max()} should be 0')
    # checks equality of all
    print(f'CHECK: {np.abs(np.sum(ecc_unit_vec**2, axis=1)-1).max()} should be 0')
    #checks normalised
    
    raan = np.where(N_unit_vec[:,1]>=0, 
                         np.arccos(N_unit_vec[:,0]),
                         2*np.pi-np.arccos(N_unit_vec[:,0]))
    argp = np.where(ecc_unit_vec[:,2]>=0, 
                       np.arccos(np.sum(N_unit_vec*ecc_unit_vec,axis=1)),
                       2*np.pi - np.arccos(np.sum(N_unit_vec*ecc_unit_vec,axis=1)))
    inc = np.arccos(h_unit_vec[:,2]) 
    
    perihel = -A + np.sqrt(A**2+B**2)

    #### Output ####
    
    d = {'vinf': vinf, # speed of ISO at infinity - in AU/day !
         'B': B, # speed of ISO at infinity - in AU/day !
         'perihel': perihel, # perihelion distance - AU
         'ecc': ecc, # eccentricity
         'inc': inc, # inclination - radians because I'm a theorist
         'raan': raan, # RA of the ascending node - radians
         'argp': argp, # argument of perihelion - radians
         'tau': tau, # MJD of perihelion
         'sma': sma, # semi-major axis length (negative as hyperbolae) - AU
         'phi': phi, # position angle of impact parameter
         't_res': t_res, # residence time of observable sphere - days
         'U': gaiadf['U'].to_numpy()[indx],
         'V': gaiadf['V'].to_numpy()[indx],
         'W': gaiadf['W'].to_numpy()[indx], 
         # components of pre-encounter velocity relative to Galactic frame - KM/S!
         'source_id': gaiadf['source_id'].to_numpy()[indx], # ID of corresponding star in gaia
         }
    orbtabledf = pd.DataFrame(data=d)
    
    if savepath is not None:
        orbtabledf.to_csv(savepath)
    
    return orbtabledf
    
