# The Ōtautahi–Oxford Model
A data-driven model of the velocity distribution of interstellar objects (ISOs) in the solar neighbourhood, combined with a fast analytical method of sampling their trajectories through the solar system.

## The Velocity Distribution 

![Galactic in-plane velocity disribution with three known ISO velocities](https://github.com/Otautahi-Oxford/OO-model/blob/main/_readme_files/UV.png)
![Distribution of ISO radiants with three known ISOs](https://github.com/Otautahi-Oxford/OO-model/blob/main/_readme_files/radiants.png)

Interstellar objects disperse away from parent stars quickly, mixing into the Galaxy-spanning population beyong the possiblity of backtracing beyond a few Myr. However, ISOs are are ejected relatively slowly, and are subject to the same Galactic dynamics as stars, meaning the ISO population as a whole will share the kinematics of the Milky Way's stellar population. 

This means we can use nearby stars to predict the properties of the ISOs passing through the solar system **__despite the fact that these stars are not the parent stars of those ISOs__**. The resulting velocity distribution is complex and non-Gaussian, with over- and under-densities introduced by orbital resonances with the Milky Way's spiral arms and bar. 

Our prediction of the ISO velocity distribution is based on a sample of ~200,000 stars within 200pc of the Sun observed by *Gaia*, debiased from survey selection effects. It is stored here in *velDist.csv* in the form of a weighted list of 3D velocities in km/s relative to the Sun in the Galactic frame (*U*, *V*, *W*). Two weightings are given: 

- **oneOverESF**, which gives the weighting which debiases the sample to give an accurate estimate of the velocity distribution of the *sine morte* stellar population, and

- **ISO_weights**, which is equal to *oneOverESF* with an additional factor proportional to the metal mass fraction of each star, 10<sup>\[M/H\]</sup>. 

*ISO_weights* is our suggested weighting for the ISO velocity distribution, as it accounts for our preliminary expectation that higher-metallicity stars will contribute more ISOs than lower-metallicty stars, based on similar correlations in exoplanet occurrence rates. However, alternate weightings can and should be experimented with as we discover more ISOs to see which fits observations best!

The remaining columns of *velDist.csv* give some of the properties of each star in the sample. For more details see 
[Hopkins et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025AJ....169...78H/abstract).


## Orbit Sampling Method

![Sample of orbits passing through inner Solar system](https://github.com/Otautahi-Oxford/OO-model/blob/main/_readme_files/oneDirectionOrbs.png)



# Cite 

[Hopkins et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025AJ....169...78H/abstract) for the velocity distribution.

[Dorsey et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025PSJ.....6..214D/abstract) (appendix A) for the orbit sampling method. 


