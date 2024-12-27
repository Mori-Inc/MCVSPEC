#ifndef MCVSPEC_H
#define MCVSPEC_H

#include <cmath>
#include <iostream>
#include <valarray>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include <xsTypes.h>
#include <funcWrappers.h>

#include "mcvspec_namespaces.hh"

using std::string;
using namespace mcvspec_physical_constants;
using namespace mcvspec_model_constants;
using namespace mcvspec;

double Calculate_Magnetic_Field();
double Calculate_White_Dwarf_Radius();
double Calculate_Accretion_Rate();
double Estimate_Shock_Height();
double Calculate_B_Free_Shock_Height();
double Calculate_Shock_Temperature();
double Calculate_Electron_Density();
double Calculate_Epsilon_Zero();
double Epsilon_Diff(double, void*);
double Epsilon_Diff_Derivative(double, void*);
void Epsilon_Diff_and_Derivative(double, void*, double*, double*);
double Root_Finder(double, int, double);
int Normalized_Position_Derivative(double, const double[], double[], void*);
void Runge_Kutta(double, double);
void Shock_Height_Shooting(double, int);
void MCVspec_Spectrum(const RealArray&, const int, RealArray&, const string&);

#endif
