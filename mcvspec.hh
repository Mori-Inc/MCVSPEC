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

using std::string;
using std::function;
using std::valarray;
using std::cout;
using std::cerr;
using std::endl;

const double pi = 3.14159265358979323846264338327950; // pi
const double solar_mass = 1.989100e+33; // solar mass in grams
const double solar_radius = 6.9599e10; // solar radius in cm
const double grav_const = 6.672590e-8; // Newton's gravitational constant in cgs units
const double wd_mol_mass = 2.0; // mean molecular mass of White dwarf
const double col_mol_mass = 0.615; // mean molecular mass of acretion column (see Cropper 1999)
const double bremss_const = 6.99e16; // bremsstrahlung constant (see Saxton: this coefficient gives the flux norm close to Suleimanov's model for the same WD mass input - 18% difference)
const double boltz_const = 1.380658e-16; // Boltzmann constant in cgs
const double boltz_const_kev = 8.617333262e-8;
const double hydrg_mass = 1.672623e-24; // Mass of hydrogen in grams
const double alpha = 2.0; // Thermodynamic Constant
const double beta = 3.85; // Thermodynamic Constant
const double proton_molar_mass = 1.0072764665789; // molar mass of a proton
const double electron_molar_mass = 5.485799090441e-4; // molar mass of an electron
const double mass_limit = 5.816*solar_mass/(pow(wd_mol_mass,2.0));
const double helium_ratio = (2*col_mol_mass-electron_molar_mass-proton_molar_mass)/(2*electron_molar_mass+4*proton_molar_mass-6*col_mol_mass);
const double electron_ion_ratio = (1+2*helium_ratio)/(1+helium_ratio);
const double apec_norm = 2.62511E-34; // (1/4)*(1 km/1 kpc)^2 normalization needed for apec
const double shock_velocity = 0.25; // downstream velocity at shock, normalized to upstream free fall velocity at shock
const double initial_veloicty = 0; // velocity at WD surface, normalized to upstream velocity at shock
const double initial_height = 0;
const double absolute_err = 1e-3;
const double relative_err = 1e-2;

// variables for user input
inline double mass, b_field, p_spin, luminosity, col_abund, wd_abund, mag_ratio, cos_incl;
inline int reflection_sel;

inline double fractional_area, accretion_area, accretion_rate, specific_accretion;
inline double shock_height, velocity_at_shock, shock_temperature, shock_electron_dens, b_free_shock_height;;
inline double free_fall_velocity, wd_radius, mag_radius;
inline double epsilon_shock, epsilon_zero; // ratio of brems cooling time to cytclotron cooling time

const int max_num_grid_points = 200;
inline int num_grid_points;
inline valarray<double> density(max_num_grid_points);
inline valarray<double> pressure(max_num_grid_points);
inline valarray<double> temperature(max_num_grid_points);
inline valarray<double> altitude(max_num_grid_points);
inline valarray<double> velocity(max_num_grid_points);
inline valarray<double> electron_dens(max_num_grid_points);
inline valarray<double> offset_altitude(max_num_grid_points);

inline double Calculate_White_Dwarf_Radius(double m){
    // TODO: update to solve WD equation of state numerically
    return solar_radius*(0.0225/wd_mol_mass)*sqrt(1.0-pow(m/mass_limit,4./3.))/cbrt(m/mass_limit);
}
inline double Calculate_Accretion_Rate(double m, double luminosity, double r_wd, double r_m){
    return luminosity/(grav_const*m*(1./r_wd - 1./r_m));
}
inline double Calculate_Accretion_Rate(double m, double luminosity, double r_wd){
    return luminosity/(grav_const*m*(1./r_wd));
}
inline double Calculate_Magnetic_Field(double m, double accretion_rate, double r_wd, double r_m){
    return sqrt(32.*accretion_rate)*pow(grav_const*m, 1./4.)*pow(r_m,7./4.)*pow(r_wd,-3);
}
inline double Root_Finder(function<double(double, void*)> func, function<double(double, void*)> derivative, void* args, double estimate, int max_itter, double tolerance){
    double x, x_new;
    int i = 0;
    while(i < max_itter && x > tolerance){
        x -= func(x, args)/derivative(x, args);
        i++;
    }
    return x;
}
inline double Epsilon_Diff(double eps_s, void* e_zero){
    double eps_zero = *(double *)e_zero;
    return eps_zero*pow(1.+eps_s, 0.425) - eps_s;
}
inline double Epsilon_Diff_Derivative(double eps_s, void* e_zero){
    double eps_zero = *(double *)e_zero;
    return 0.425*eps_zero*pow(1. + eps_s, -0.575) - 1.;
}
inline double Calculate_B_Free_Shock_Height(double v_freefall, double m_dot){
    // TODO: replace url with doi
    double integral = (39.*sqrt(3.) - 20*pi)/96.; // value of integral from EQ 7a of Wu 1994: https://articles.adsabs.harvard.edu/pdf/1994ApJ...426..664W
    return pow(v_freefall,3.)*integral/(2.*bremss_const*m_dot);
}
inline double Calculate_Epsilon(double b_field, double shock_temp, double shock_e_dens, double shock_height){
    return 9.1e-3*pow(b_field/10.,2.85)*pow(shock_temp/1e8,2.)*pow(shock_e_dens/1e16,-1.85)*pow(shock_height/1e7,-0.85);
}
inline double Calculate_Shock_Temperature(double radius, double shock_height, double mag_radius){
    return (3./8.)*grav_const*mass*col_mol_mass*hydrg_mass*(1/(radius+shock_height) - 1/mag_radius)/boltz_const;
}
inline double Calculate_Electron_Density(double accretion_rate, double v_free_fall){
    return (accretion_rate/v_free_fall)*shock_velocity*electron_ion_ratio/(hydrg_mass*col_mol_mass);
}


inline int Normalized_Position_Derivative(double velocity_normed, const double position_normed[], double pos_derivative[], void* eps_s){
    double epsilon_s = *(double *)eps_s;
    double prefactor = pow(free_fall_velocity,3)/(2.*bremss_const*specific_accretion);
    pos_derivative[0] = prefactor*(pow(velocity_normed,1.5)*(5.0-8.*velocity_normed)/sqrt(1-velocity_normed)/
                        (1.+epsilon_s*pow(3.,-alpha)*pow(4.,alpha+beta)*
                        pow(1.0-velocity_normed,alpha)*pow(velocity_normed,beta)));
    return GSL_SUCCESS;
}
inline void Runge_Kutta(double initial_veloicty, double initial_height){
    gsl_odeiv2_system accretion_column = {Normalized_Position_Derivative, nullptr, 1, &epsilon_shock};
    const gsl_odeiv2_step_type *ode_type = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *ode_stepper = gsl_odeiv2_step_alloc(ode_type, 1);
    gsl_odeiv2_control *ode_control = gsl_odeiv2_control_y_new(absolute_err, relative_err);
    gsl_odeiv2_evolve *ode_evolver = gsl_odeiv2_evolve_alloc(1);
    double vel = initial_veloicty;
    double height[1] = {initial_height};
    int n = 0;
    double initial_step = absolute_err;
    velocity[0] = vel*free_fall_velocity;
    altitude[0] = height[0]*shock_height;
    while (vel < shock_velocity){
        int status = gsl_odeiv2_evolve_apply(ode_evolver, ode_control, ode_stepper, &accretion_column,
                                             &vel, shock_velocity ,&initial_step, height);
        if(status != GSL_SUCCESS){
            break;
        }
        velocity[n] = vel*free_fall_velocity;
        altitude[n] = height[0]*shock_height;
        n += 1;
    }
    gsl_odeiv2_evolve_free(ode_evolver);
    gsl_odeiv2_control_free(ode_control);
    gsl_odeiv2_step_free(ode_stepper);

    for (int i = 0; i < n; i++){
        density[i] = specific_accretion/velocity[i];
        pressure[i] = specific_accretion*(free_fall_velocity - velocity[i]);
        temperature[i] = shock_temperature*velocity[i]*(free_fall_velocity-velocity[i])*(16./(3.*pow(free_fall_velocity, 2)));
        electron_dens[i] = density[i]*electron_ion_ratio/(hydrg_mass*col_mol_mass);
    }
    num_grid_points =  n;
}
inline void Shock_Height_Shooting(double tolerance, int max_itter){
    double old_height, error;

    error = 1e100;
    int i = 0;

    Runge_Kutta(initial_veloicty, initial_height);
    shock_height = altitude[num_grid_points-1];
    old_height = altitude[num_grid_points-1];

    while(i < max_itter && error > tolerance){
        free_fall_velocity = sqrt(2.*grav_const*mass*(1/(wd_radius+shock_height) - 1/mag_radius));
        epsilon_shock = Calculate_Epsilon(b_field, shock_temperature, shock_electron_dens, shock_height);

        Runge_Kutta(initial_veloicty, initial_height);

        shock_height = altitude[num_grid_points-1];
        if (wd_radius + shock_height > mag_radius) {
            shock_height = old_height;
            cerr << "Accretion column height has exceded magnetospheric radius" << endl;
        }
        error = abs((old_height-shock_height)/old_height);
        old_height = shock_height;
        i++;
    }
}

inline void MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string){

    int n = flux.size();

    double segment_height = altitude[0];

    valarray<double> flux_from_layer(n);
    valarray<double> flux_error(n);

    valarray<double> apec_parameters(3);
    valarray<double> refl_parameters(5);

    refl_parameters[0] = 1-sqrt(1.0-1.0/pow(1+altitude[num_grid_points-1]/wd_radius,2));
    refl_parameters[1] = 0.;
    refl_parameters[2] = wd_abund;
    refl_parameters[3] = wd_abund;
    refl_parameters[4] = cos_incl;

    for(int i=0; i<num_grid_points; i++){
        for(int j=0; j<n; j++){
            flux[j] += flux_from_layer[j];
            flux_from_layer[j] = 0;
        }

        apec_parameters[0] = temperature[i]*boltz_const_kev;
        apec_parameters[1] = col_abund;
        apec_parameters[2] = 0.0;

        if (apec_parameters[0] > 64.0){
            CXX_bremss(energy, apec_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }
        else{
            CXX_apec(energy, apec_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }

        // normalizes spectrum appropriately
        if(i!=0){
            segment_height = altitude[i]-altitude[i-1];
        }
        for(int j=0; j<n; j++){
            flux_from_layer[j] *= apec_norm*segment_height*pow(electron_dens[i]*1e-7,2)/electron_ion_ratio;
        }

        // apply reflect to each slice
        if (reflection_sel == 2){
            refl_parameters[0] = 1-sqrt(1.0-1.0/pow(1+altitude[i]/wd_radius,2));
            CXX_reflect(energy, refl_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }
    }

    for(int j=0; j<n; j++){
        flux[j] += flux_from_layer[j];
        flux_from_layer[j] = 0;
    }

    if(reflection_sel==1){
        CXX_reflect(energy, refl_parameters, spectrum_num, flux, flux_error, init_string);
    }
}

#endif
