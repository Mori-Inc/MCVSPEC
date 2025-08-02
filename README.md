# MCVSPEC: An X-Ray Spectral Model for Magnetic Cataclysmic Variables

## Introduction
MCVSPEC is a model for the post-shock accretion flow in magnetic Cataclysmic Variables (mCVs).
It solves for the full thermal profile of the post-shock accretion column (PSAC) and can produce an X-Ray spectrum from that profile. Currently MCVSPEC supports two interfaces: python and xspec. These interfaces are not totally equivalent (in particular the python interface provides more flexible accsess to the thermal profiles while the xspec interface provides a more accurate spectrum) so it is recommended that both are installed where possible.

## XSPEC Interface

MCVSPEC contains two models: `polarspec` and `ipsepc` which are used for their respective class of objects. Both models share the common inputs:

| Variable       | Units                  | Description                                                           |
|----------------|------------------------|-----------------------------------------------------------------------|
| `reflectOn`    | 0, 1                   | Toggle reflection                                                     |
| `M`            | Solar masses (M☉)      | Mass of the white dwarf (WD)                                          |
| `L`            | 10^33 ergs/s           | Luminosity of the WD                                                  |
| `Z`            | Relative to solar      | Accretion column abundance                                            |
| `cos i`        | Dimensionless (0 to 1) | Cosine of the inclination angle of the reflecting surface             |
| `areaScal`     | Dimensionless          | Column cross-sectional scaling exponent A~(1+x/R)^n                   |
| `distance`     | parsecs                | distance to source (used only for flux normalization)                 |

The magnetic field strength is determined differently for each class of mCV.

| Model Name | Magnetic Variable | Units    | Description                                 |
|------------|-------------------|----------|---------------------------------------------|
| polarspec  | `B`               | MG       | Surface magnetic field of WD                |
| ipspec     | `Pspin`           | s        | Spin period of WD                           |

`ipsepc` is our spin-equilibrium model for intermediate polars, it computes the magnetospheric radius from the spin period (and accretion rate) by setting the ram pressure equal to the magnetic pressure. It has one additional parameter: `CoRotRatio`, which is the ratio of the magnetoshperic radius to the co-rotation radius for additional flexibility.

Additionally, the accretion area can be specified in one of two ways. Either (in the case of `ipsepc` and `polarspec`) the fractional accretion area (`f`) can be specified or (in the case of `aripsepc` and `arpolarspec`)the accretion area (A) in units of 10<sup>15</sup> cm<sup>2</sup>.

## Python Interface

MCVSPECs python interface utilizes `AtomDB`s python api to generate the primary thermal X-Ray spectrum but has no implementation for X-Ray reflection at this time.

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

By default both the python and xspec interface will be built which require that cmake can find `HEASoft` and `pybind11` respectively. `HEASoft` is found through `$HEADAS` environmental variable. If you would like to turn off either of the interfaces it can be done with the cmake variables BUILD_PYBINDINGS and
BUILD_XSMODEL.

After cmake configures your build simply run `make` and `make install`.

`cmake` will execute the `HEASoft` utility `initpackage` supplying your install directory with the necessary build files.
`make` will execute `hmake` producing the mcvspec xspec library.

### Initialize Package
If installing pyMCVSPEC simply add `/path/to/install/pymcvspec` to your `$PYTHONPATH`
For the XSPEC interface one can run, from the xspec prompt, `lmod mcvspec /path/to/install/xspec` to load the mcvspec models into their xspec session. If you would like to configure XSPEC to load mcvspec on start you can always add `load /path/to/install/xspec/libmcvspec.dylib` the `global_customize.tcl` file in `HEASoft` (described in the "Customizing system-wide" section of the XSPEC manual)

## Usage

### Load the mcvspec models in xspec
Once loaded into XSPEC you can then choose between `polarspec`, `ipsepc`, `apolarspec`, or `aipspec` depending on your source.

To load the python module:
```python
import mcvspec
import astropy.units as u
# create an instace of a polar
mass = 0.63*u.M_sun
b_field = 13.6*u.MG
metalicity = 0.5
luminosity = 0.6e33*u.erg/u.s
fractional_area = 1e-3
cos_incl_angle = 0.5
src_distance = 160*u.pc
my_source = mcvspec.polar(mass, b_field, metalicity,
    luminosity, fractional_area, cos_incl_angle,
    src_distance)

# get information about the thermal profile
my_source.altitude
my_source.electron_temperature
my_source.ion_density

# get spectrum
e_bins = np.arange(3., 79., 0.4)*u.keV
flux = my_source.spectrum(e_bins)
```
Note that there is no reflection implemented in python yet
