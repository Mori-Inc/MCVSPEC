# MCVSPEC: An X-Ray Spectral Model for Magnetic Cataclysmic Variables

## Introduction
MCVSPEC is a model for post-shock accretion flow in magnetic Cataclysmic Variables (mCVs).
It solves for the full thermal profile of the post-shock accretion column (PSAC) and can produce an X-Ray spectrum from that profile. Currently MCVSPEC supports two interfaces: python and xspec. These interfaces are not perfectly equivalent (in particular the python interface provides more flexible accsess to the thermal profiles while the xspec interface provides a more accurate spectrum) so it is recommended that both are installed where possible.

## Dependencies
* CMake
* Python
### XSPEC Interface
* [HEASOFT](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/)
### Python Interface
* [Pybind11](https://pybind11.readthedocs.io/en/stable/)
* [NumPy](https://numpy.org)
* [Astropy](https://docs.astropy.org/en/stable/index.html)
* [PyAtomDB](https://atomdb.readthedocs.io/en/master/)

## Installation

### Clone Repository (Or Download Release)
We will assume that the repository is located in `/path/to/MCVSPEC`. You will also want to create a build directory (`/path/to/MCVSPEC/build` is fine, but it can be anywhere)

### Build MCVSPEC
From your build directory execute

`cmake /path/to/MCVSPEC` (optionally specify an installation directory with `-DCMAKE_INSTALL_PREFIX`)

By default, both the python and xspec interface will be built. This requires that cmake can find `HEASoft` and `pybind11` respectively. `HEASoft` is found through `$HEADAS` environmental variable. If you would like to turn off either of the interfaces it can be done with the cmake variables BUILD_PYBINDINGS and
BUILD_XSMODEL.

After cmake configures your build, simply run `make` and `make install`.

`cmake` will execute the `HEASoft` utility `initpackage` supplying your install directory with the necessary build files.
`make` will execute `hmake` producing the mcvspec xspec library.

### Initialize Package
If installing pyMCVSPEC simply add `/path/to/install/pymcvspec` to your `$PYTHONPATH`
For the XSPEC interface one can run, from the xspec prompt, `lmod mcvspec /path/to/install/xspec` to load the mcvspec models into their xspec session. If you would like to configure XSPEC to load mcvspec on start you can always add `load /path/to/install/xspec/libmcvspec.dylib` the `global_customize.tcl` file in `HEASoft` (described in the "Customizing system-wide" section of the XSPEC manual)

## Usage.

## XSPEC Interface

MCVSPEC contains a total of eight models. All eight models have the following paramters:

| Variable          | Units                  | Description                                                           |
|-------------------|------------------------|-----------------------------------------------------------------------|
| `M`               | Solar masses (M☉)      | Mass of the white dwarf (WD)                                          |
| *`Mag Var`*       | Model dependant        | Variable to set the magnetic field strength of the WD                 |
| *`Accretion Var`* | Model dependant        | Variable to set the accretion rate of the WD                          |
| *`Area Var`*      | Model dependant        | Variable to set the accretion area of the WD                          | 
| `abund`           | Relative to solar      | Accretion column metalicity                                           |
| `cosAngle`        | Dimensionless (0 to 1) | Cosine of the inclination angle of the PSAC                           |
| `distance`        | parsecs                | distance to source (used only for flux normalization)                 |
| `reflectOn`       | 0, 1                   | Toggle reflection                                                     |

Four of the models are designed for fitting polars (`polarspec`) and the other four, for intermediate polars (`ipspec`). For a polar the *`Mag Var`* is simply `B`, the sufrace magnetic field of the the WD in MG. For intermediate polars we assume that the B-field is not known but can be deduced from the spin period of the IP. Under the assumption of spin-equilibirum, a magnetospheric radius (R<sub>m</sub>)can be deduced from the WD spin period. The Magnetic field is then inferred from this radius and the accretion rate of the IP. The *`Mag Var`* is then the spin period (`Pspin`) in seconds. A second variable, `corotRatio` exists which is the ratio of the magnetoshperic radius to the cortation radius. This permits the user to choose an arbitrary magnetospheric radius should they so desire. 

The "default" models: `polarspec` and `ipspec` use the bolometric luminosity `L`, measured in units of 10<sup>33</sup> ergs/s, to deduce the total mass accretion rate and the fractional accretion area (`f`) to deduce the accretion area. Three additional variants of these models exist which are prepended with a combination of "d" (for "m Dot") and "a" (for area) which allow the user to, alternatively, directly specify the specific accretion rate, with parameter `mdot` measured in g cm<sup>-2</sup> s<sup>-1</sup>, and area, with paramter `area` measured in 10<sup>15</sup> cm<sup>2</sup>.

Internally, all eight models produce the same underlying object, an instance of the `Cataclysmic_Variable` class. The only difference between the eight models is the mechanisms by which the parameters of this class are determined. 

## Python Interface

MCVSPECs python interface utilizes `AtomDB`s python api to generate the primary thermal X-Ray spectrum but has no implementation for X-Ray reflection at this time.

To load the python module:
```python
import mcvspec
import astropy.units as u
# create an instace of a cataclysmic variable
mass = 0.63*u.M_sun
b_field = 13.6*u.MG
accretion_rate = 1*u.g/u.cm**2/u.s
accretion_area = 1e15*u.cm**2
abund = 0.5
incl = 60*u.deg
dist = 160*u.pc
my_source = mcvspec.cataclysmic_variable(mass, b_field, accretion_rate, 
    accretion_area, metalicity=abund, orbital_inclination=incl
    distance=dist)

# get information about the thermal profile
my_source.altitude
my_source.electron_temperature
my_source.ion_density

# get spectrum
e_bins = np.arange(3., 79., 0.4)*u.keV
flux = my_source.spectrum(e_bins)
```
Note that there is no reflection implemented in python yet.

The python interface has 3 classes: `cataclysmic_variable`, `polar`, and `intermediate_polar` which mirror the construction of the XSPEC model. `polar` and `intermediate_polar` inherit from the `cataclysmic_variable` class and use observable paramters in their construction while `cataclysmic_variable` uses physical paramters of the PSAC.

## The WD Mass-Radius Relationship

MCVSPEC uses a numerical solution to the WD mass-radius relationship following the model of [Hamada and Salpeter (1961)](https://ui.adsabs.harvard.edu/abs/1961ApJ...134..683H/abstract). This solution, at discrete points, is stored in a file, `data/mass_radius.txt`, and is interpolated to estimate the mass-radius relationship. This file is generated with the python script `data/wd_mass_radius.py` and then converted into a C++ header during the build process. Hamada and Salpeter's model depends on the white dwarf composition which is set with the parameters $\mu=\frac{A}{Z}$ and $Z$. The default `data/mass_radius.txt` was generated for $\mu=2$ and $Z=6$. To generate a different mass-radius relationship simply execute the script with `python wd_mass_radius.py -Z <your_z> -mu <your_mu>` and then rebuild the software. One could also test their own mass-radius relationship by generating an identically formatted `data/mass_radius.txt` and rebuilding.
