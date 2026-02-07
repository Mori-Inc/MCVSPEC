#pragma once

#include <vector>
using std::vector;

// physical constants
constexpr double pi = 3.14159265358979323846264338327950;
constexpr double c = 2.99792458e10; // cm/s
constexpr double m_e = 9.109383713928e-28; // grams
constexpr double m_p = 1.67262192595e-24; // grams
constexpr double m_sol = 1.989100e+33; // grams
constexpr double r_sol = 6.9599e10; // cm
constexpr double grav_const = 6.672590e-8; // dyn cm^2 g^-2
constexpr double k_b = 1.380658e-16; // erg/K
constexpr double k_b_keV = 8.617333262e-8; // keV/K
constexpr double planck_const = 6.62607015e-27; // erg s
constexpr double hbar = planck_const/(2*pi);
constexpr double alpha = 7.2973525643e-3;
constexpr vector<double> atomic_charge = {1,2,6,7,8,10,12,13,14,16,18,20,26,28}; // charges of elements in abundances array
constexpr vector<double> atomic_mass = {1.007975,4.002602,12.0106,14.006855,15.9994,20.17976,24.3055,
                                26.98153843,28.085,32.0675,39.8775,40.0784,55.8452,58.69344};
constexpr double bremss_coeff = sqrt(512*pi/(27.*m_e*m_e*m_e))*alpha*alpha*alpha*hbar*hbar; // cgs bremms constant for hydrogen plasma
constexpr double cyclotron_coeff = 8.07e-2/(k_b*k_b);
constexpr double exchange_coeff = 4*alpha*alpha*hbar*hbar*c*c*sqrt(2*pi*m_e);
constexpr double coulomb_log_coeff = 2*m_e/(pi*alpha*c*hbar*hbar*hbar);

// conversion factors
constexpr double erg_to_kev = 6.241509074461e8;
constexpr double amu_to_g =  1.6605390689252e-24; // mass of amu in grams
constexpr double pc_to_cm = 3.0856775814913673e18;
constexpr double ryd_to_erg = 2.1798723611035845e-11; // rydberg energy unit in ergs
// constants of the model
constexpr double kT_grid_spacing = 0.5; // keV
constexpr double altitude_grid_spacing = 0.1; // fractional
