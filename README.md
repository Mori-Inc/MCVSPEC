# X-ray Spectral Models for Magnetic CV's

There are three versions of MCVSPEC, our spectral model for X-ray emissions from magnetic cataclysmic variables. 

The common input parameters are: reflectOn (switch parameter to determine whether to model reflection from WD surface; 1 if yes), M (WD mass), f (fractional accretion area), luminosity (L, in 10^33 ergs/s), Z_wd (WD surface abundance relative to solar), Z (accretion column abundance relative to solar), and cos i (inclination angle of reflecting surface). Depending on source type and spin equilibrium, there's an additional parameter used to estimate magnetospheric radius.

POLARSPEC: model for polars. Accepts B (magnetic field in 10^6 G)

IPSPEC: model for IP's, assuming spin equilibrium. Accepts P_spin (spin period in s)

DEQSPEC: model for IP's, assuming spin disequilibrium. Accepts R_m/R_wd (magnetospheric radius ratio)
