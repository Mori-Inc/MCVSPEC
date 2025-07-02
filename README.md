# X-ray Spectral Models for Magnetic CV's

This is MCVSPEC, our spectral model for modelling X-ray emissions from magnetic cataclysmic variables (Bridges et al. in prep).

MCVSPEC has three versions. The common input parameters are:
| Variable       | Units                 | Description                                                           |
|----------------|-----------------------|-----------------------------------------------------------------------|
| `reflectOn`    | 0, 1, or 2             | Switch to determine if and how reflection from WD surface is modeled (0 = no reflection; 1 = use reflect at varying heights along the column, which is more time-consuming; 2 = use reflect once at shock height)|
| `M`            | Solar masses (M☉)     | Mass of the white dwarf (WD)                                           |
| `f`            | Dimensionless (0 to 1) | Fractional accretion area of the WD                                    |
| `L`            | 10^33 ergs/s           | Luminosity of the WD                                                   |
| `Z_wd`         | Relative to solar      | WD surface abundance                                  |
| `Z`            | Relative to solar      | Accretion column abundance                            |
| `cos i`        | Dimensionless (0 to 1) | Cosine of the inclination angle of the reflecting surface              |

Depending on source type and spin equilibrium, there's an additional parameter used to estimate magnetospheric radius.

POLARSPEC: model for polars. Accepts B (magnetic field in 10^6 G)

IPSPEC: model for IP's, assuming spin equilibrium. Accepts P_spin (spin period in s)

DEQSPEC: model for IP's, assuming spin disequilibrium. Accepts R_m/R_wd (magnetospheric radius ratio)
