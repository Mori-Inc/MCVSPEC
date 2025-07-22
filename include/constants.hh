#pragma once

#include <valarray>
using std::valarray;

// physical constants
const double pi = 3.14159265358979323846264338327950;
const double light_speed = 2.99792458e10; // cm/s
const double electron_mass = 9.109383713928e-28; // grams
const double proton_mass = 1.67262192595e-24; // grams
const double solar_mass = 1.989100e+33; // grams
const double solar_radius = 6.9599e10; // cm
const double grav_const = 6.672590e-8; // dyn cm^2 g^-2
const double boltz_const = 1.380658e-16; // erg/K
const double boltz_const_kev = 8.617333262e-8; // keV/K
const double planck_const = 6.62607015e-27; // erg s
const double hbar = planck_const/(2*pi);
const double fine_structure_constant = 7.2973525643e-3;
const valarray<int> atomic_charge = {1,2,6,7,8,10,12,13,14,16,18,20,26,28}; // charges of elements in abundances array
const valarray<double> atomic_mass = {1.007975,4.002602,12.0106,14.006855,15.9994,20.17976,24.3055,
                                26.98153843,28.085,32.0675,39.8775,40.0784,55.8452,58.69344};
// conversion factors
const double erg_to_kev = 6.241509074461e8;
const double amu_to_g =  1.6605390689252e-24; // mass of amu in grams
const double pc_to_cm = 3.0856775814913673e18;
// constants of the model
const double gaunt_factor = 1.2;
const double wd_mol_mass = 2.; // white dwarf molar mass
const double shock_ratio = 4.;
const double alpha = 2.0; // cyclotron cooling rate pressure exponent
const double beta = 3.85; // cyclotron cooling rate density exponent
