# MCVSPEC: An X-Ray Spectral Model for Magnetic Cataclysmic Variables

## Introduction
MCVSPEC is a model for the post-shock accretion flow in magnetic Cataclysmic Variables (mCVs).
It solves the hydrodynamic equations for the flow of material in the post shock region of the mCV accretion column producing a thermal profile of the column.
The final X-Ray spectrum is produced from the sum of APEC models and their reflected components.

MCVSPEC contains three models for different classes of mCVs. All 3 use the following input parameters are:
| Variable       | Units                 | Description                                                           |
|----------------|-----------------------|-----------------------------------------------------------------------|
| `reflectOn`    | 0, 1, or 2             | Switch to determine if and how reflection from WD surface is modeled (0 = no reflection; 1 = use reflect once at shock height; 2 = varying heights along the column, which is more time-consuming)|
| `M`            | Solar masses (M☉)     | Mass of the white dwarf (WD)                                           |
| `f`            | Dimensionless (0 to 1) | Fractional accretion area of the WD                                    |
| `L`            | 10^33 ergs/s           | Luminosity of the WD                                                   |
| `Z`            | Relative to solar      | Accretion column abundance                            |
| `cos i`        | Dimensionless (0 to 1) | Cosine of the inclination angle of the reflecting surface              |
| `distance`     | parsecs                | distance to source (used only for flux normalization)              |

The magnetic field strength is determined differently for each class of mCV.


| Model Name | Magnetic Variable | Units    | Description | Object Description                     |
|------------|-------------------|----------|-------------|--------------------------|
| polarspec  | `B`               | MG       | Surface magnetic field of WD | polar                                  |
| ipspec     | `Pspin`           | s        | Spin period of WD | intermediate polar in spin equilibirum |
| deqspec    | `Rm/R_wd`         | unitless | Ratio of magnetospheric radius to WD radius | intermediate polar not in spin equilibirum |


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
