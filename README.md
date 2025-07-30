# MCVSPEC: An X-Ray Spectral Model for Magnetic Cataclysmic Variables

## Introduction
MCVSPEC is a model for the post-shock accretion flow in magnetic Cataclysmic Variables (mCVs).
It solves for the full thermal profile of the post-shock accretion column (PSAC) and produces an X-Ray spectrum by convolving the thermal profile with the optically-thin thermal bremsstrahlung spectrum produced by APEC. Additioanlly, MCVSPEC is capable of reflecting the PSAC column off the WD surface.

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

## Dependencies
* [HEASOFT](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/)
* CMake
* Python (recommended)
* [Pybind11](https://pybind11.readthedocs.io/en/stable/) (recommended)

## Installation

### Clone Repository (Or Download Release)
We will assume that the repository is located in `/path/to/MCVSPEC`. You will also want to create a build directory (`/path/to/MCVSPEC/build` is fine, but it can be anywhere)

### Build MCVSPEC
`HEASoft` must be initialized and (optionally but highly recommended) the virtual environment you wish to install pyMCVSPEC to should be activated.

From your build directory execute

`cmake /path/to/MCVSPEC` (optionally specify an installation directory with `-DCMAKE_INSTALL_PREFIX`)

`make`

`make install`

`cmake` will execute the `HEASoft` utility `initpackage` supplying your install directory with the necessary build files.
`make` will execute `hmake` producing the mcvspec xspec library.

### Initialize Package
If installing pyMCVSPEC simply add `/path/to/install/pymcvspec` to your `$PYTHONPATH`

## Usage

### Load the mcvspec models in xspec
From the xspec prompt run `lmod mcvspec /path/to/install/xspec` to load the mcvspec models. You can then choose between `polarspec`, `ipsepc`, or `deqspec` depending on your source.

To load the python module:
```python
import mcvspec
# create an instace of a polar
my_source = mcvspec.polar(mass, b_field, metalicity,
    luminosity, fractional_area, cos_incl_angle,
    src_distance, reflection)
# solve for the accretion column profile (int argument is max number of itterations for shock height determination)
my_source.execute(10000)
# get information about the thermal profile
my_source.altitude()
my_source.electron_temperature()
my_source.ion_density()
```
