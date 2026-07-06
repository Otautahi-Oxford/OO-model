# The Ōtautahi–Oxford Model
A data-driven model of the velocity distribution of interstellar objects (ISOs) in the solar neighbourhood, combined with a fast analytical method of sampling their trajectories through the solar system.

## The Velocity Distribution 


Interstellar objects disperse away from parent stars quickly, mixing into the Galaxy-spanning population beyond the possiblity of backtracing beyond a few Myr. However, ISOs are are ejected relatively slowly, and are subject to the same Galactic dynamics as stars, meaning the ISO population as a whole will share the kinematics of the Milky Way's stellar population. 

This means we can use nearby stars to predict the properties of the ISOs we observe passing through the solar system **_despite the fact that these stars are not the parent stars of those ISOs_**. The resulting velocity distribution is complex and non-Gaussian, with over- and under-densities introduced by orbital resonances with the Milky Way's spiral arms and bar. 
![Galactic in-plane velocity disribution with three known ISO velocities](https://github.com/Otautahi-Oxford/OO-model/blob/main/_readme_files/UV.png)

Our prediction of the ISO velocity distribution is based on a sample of ~200,000 stars within 200pc of the Sun observed by *Gaia*, debiased from survey selection effects. It is stored here in *velDist.csv* in the form of a weighted list of 3D velocities in km/s relative to the Sun in the Galactic frame (*U*, *V*, *W*). Two weightings are given: 

- **oneOverESF**, which gives the weighting which debiases the sample to give an accurate estimate of the velocity distribution of the *sine morte* stellar population, and

- **ISO_weights**, which is equal to **oneOverESF** with an additional factor proportional to the metal mass fraction of each star, 10<sup>\[Fe/H\]</sup>. 

**ISO_weights** is our suggested weighting for the ISO velocity distribution, as it accounts for our preliminary expectation that higher-metallicity stars will contribute more ISOs than lower-metallicty stars based on similar correlations in exoplanet occurrence rates. However, alternate weightings can and should be experimented with as we discover more ISOs to see which fits the observed distribution best.

The remaining columns of *velDist.csv* give some of the properties of each star in the sample, including the *Gaia* source ID for crossmatching. For more details see 
[Hopkins et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025AJ....169...78H/abstract).

![Distribution of ISO radiants with three known ISOs](https://github.com/Otautahi-Oxford/OO-model/blob/main/_readme_files/radiants.png)


## Orbit Sampling Method

The above velocity distribution gives the kinematics of ISOs in the solar neighbourhood but before entering the solar system, far from the gravitational influence of the Sun. However, in the inner solar system where ISOs are observable, the ISO population is subject to complex velocity-dependent effects such as gravitational focussing.

*orbitSampler.py* and its function *sampleOrbits* implements the analytical trajectory sampling method described in appendix A of [Dorsey et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025PSJ.....6..214D/abstract), allowing the sampling of any number of ISO trajectories passing through a sphere over the course of a finite-length survey. 
It assumes the above *Gaia*-based velocity distribution, and assumes ISO move on hyperbolae centred on the solar system barycentre.
This method is fast and efficient as it does not require inefficient rejection sampling or slow orbit integration to capture the effect of gravitational focussing.

![Sample of orbits passing through inner solar system](https://github.com/Otautahi-Oxford/OO-model/blob/main/_readme_files/oneDirectionOrbs.png)


# Cite 

[Hopkins et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025AJ....169...78H/abstract) for the velocity distribution.

[Dorsey et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025PSJ.....6..214D/abstract) (appendix A) for the orbit sampling method. 


