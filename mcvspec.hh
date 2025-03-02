#ifndef MCVSPEC_H
#define MCVSPEC_H

#include <cmath>
#include <iostream>
#include <valarray>
#include <algorithm>
#include <vector>
#include <xsTypes.h>
#include <funcWrappers.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

#include "tableau.hh"

using std::string;
using std::function;
using std::valarray;
using std::vector;
using std::find;
using std::abs;
using std::max;
using std::min;
using std::pow;
using std::cout;
using std::endl;

// physical constants
const double pi = 3.14159265358979323846264338327950; // pi
const double electron_mass = 9.109383713928e-28;
const double proton_mass = 1.67262192595e-24;
const double solar_mass = 1.989100e+33; // solar mass in grams
const double solar_radius = 6.9599e10; // solar radius in cm
const double grav_const = 6.672590e-8; // Newton's gravitational constant in cgs units
const double boltz_const = 1.380658e-16; // Boltzmann constant in cgs
const double boltz_const_kev = 8.617333262e-8;
const double planck_const = 6.62607015e-27;
const double fine_structure_constant = 7.2973525643e-3;
// conversion factors
const double erg_to_kev = 6.241509074461e8;
const double amu_to_g =  1.6605390689252e-24; // mass of amu in grams
const double pc_to_cm = 3.0856775814913673e18;
// constants of the model
const double shock_ratio = 4.;
const double alpha = 2.0; // Thermodynamic Constant
const double beta = 3.85; // Thermodynamic Constant
const double gaunt_factor = 1.2;
inline double bremss_const = sqrt(2*pi/3)*4*planck_const*planck_const*pow(fine_structure_constant,3)*gaunt_factor/
                             (3*pow(electron_mass,3)*pi*pi*sqrt(2)*pow(proton_mass/electron_mass + 1., 1.5)); // See Rybicki and Lightman
inline double epsilon_const = 1;
inline double kt_const = 3./16.;
inline double electron_density_const = 1;
inline double avg_atomic_charge;
inline double abundances[14] = {1.,0,0,0,0,0,0,0,0,0,0,0,0,0};
const int atomic_charge[14] = {1,2,6,7,8,10,12,13,14,16,18,20,26,28}; // charges of elements in abundances array
const double atomic_mass[14] = {1.007975,4.002602,12.0106,14.006855,15.9994,20.17976,24.3055,
                                26.98153843,28.085,32.0675,39.8775,40.0784,55.8452,58.69344};

const double wd_mol_mass = 2.;
const double mass_limit = 5.816*solar_mass/(pow(wd_mol_mass,2.0));
const double initial_veloicty = 1/shock_ratio; // velocity at WD surface, normalized to upstream velocity at shock
const double initial_height = 1.;
const double absolute_err = 1e-8;
const double relative_err = 1e-6;
const double max_itter = 10000;

// variables for user input
inline double mass, b_field, p_spin, luminosity, col_abund, wd_abund, mag_ratio, cos_incl, source_distance;
inline int reflection_sel;

inline double fractional_area, accretion_area, accretion_rate, specific_accretion;
inline double shock_height, velocity_at_shock, shock_temperature, shock_electron_dens, b_free_shock_height;;
inline double free_fall_velocity, wd_radius, mag_radius;
inline double epsilon_shock; // ratio of brems cooling time to cytclotron cooling time

inline vector<double> density;
inline vector<double> pressure;
inline vector<double> temperature;
inline vector<double> altitude;
inline vector<double> velocity;
inline vector<double> electron_dens;
inline vector<double> offset_altitude;

inline void Set_Abundances(double metalicity){
    abundances[0] = FunctionUtility::getAbundance(atomic_charge[0]);
    abundances[1] = FunctionUtility::getAbundance(atomic_charge[1]);
    double norm = abundances[0]+abundances[1];
    for(int i = 2; i < 14; i++){
        abundances[i] = metalicity*FunctionUtility::getAbundance(atomic_charge[i]);
        norm += metalicity*abundances[i];
    }
    double avg_ion_mass = 0;
    avg_atomic_charge = 0;
    double avg_charge_squared = 0;
    for(int i = 0; i < 14; i++){
        avg_ion_mass += (abundances[i]/norm)*atomic_mass[i];
        avg_atomic_charge += (abundances[i]/norm)*atomic_charge[i];
        avg_charge_squared += (abundances[i]/norm)*atomic_charge[i]*atomic_charge[i];
    }
    avg_ion_mass *= amu_to_g;
    bremss_const = sqrt(2.*pi/3.)*gaunt_factor*4*planck_const*planck_const*pow(fine_structure_constant,3);
    bremss_const /= 3*pow(electron_mass, 3)*pi*pi;
    bremss_const *= avg_charge_squared*avg_atomic_charge/sqrt(avg_atomic_charge+1);
    bremss_const /= sqrt(pow(avg_atomic_charge + avg_ion_mass/electron_mass,3));
    epsilon_const = (7.62e-2/gaunt_factor)*pow(10.,0.025)*(pow(shock_ratio-1.,2)/pow(shock_ratio,5.85))*(pow(electron_mass, 3.85)/pow(boltz_const, 2));
    epsilon_const *= (avg_atomic_charge/avg_charge_squared)*pow(avg_atomic_charge/(1.+avg_atomic_charge),2);
    epsilon_const *= pow((avg_atomic_charge + avg_ion_mass/electron_mass)/avg_atomic_charge, 3.85);
    kt_const = electron_mass*(avg_atomic_charge + avg_ion_mass/electron_mass)/(1.+avg_atomic_charge);
    electron_density_const = avg_atomic_charge/(avg_atomic_charge + avg_ion_mass/electron_mass);
}


inline double Calculate_White_Dwarf_Radius(double m){
    // TODO: update to solve WD equation of state numerically
    return solar_radius*(0.0225/wd_mol_mass)*sqrt(1.0-pow(m/mass_limit,4./3.))/cbrt(m/mass_limit);
}
inline double Calculate_Accretion_Rate(double m, double luminosity, double r_wd){
    return luminosity*r_wd/(grav_const*m);
}
inline double Calculate_Accretion_Rate(double m, double luminosity, double r_wd, double r_m){
    return luminosity/(grav_const*m*(1./r_wd - 1./r_m));
}
inline double Calculate_Magnetic_Field(double m, double accretion_rate, double r_wd, double r_m){
    return sqrt(32.*accretion_rate)*pow(grav_const*m, 1./4.)*pow(r_m,7./4.)*pow(r_wd,-3);
}
inline double Calculate_Free_Fall_Velocity(double m, double r_wd, double h_s){
    return sqrt(2*grav_const*m/(r_wd+h_s));
}
inline double Calculate_Free_Fall_Velocity(double m, double r_wd, double h_s, double r_m){
    return sqrt(2*grav_const*m*(1./(r_wd+h_s) - 1./r_m));
}
inline double Calculate_Shock_Density(double m_dot, double vel){
    return m_dot/vel;
}
inline double Calculate_Epsilon(double v_freefall){
    double density = Calculate_Shock_Density(specific_accretion, v_freefall);
    return epsilon_const*pow(accretion_area,-17./40.)*pow(b_field, 57./20.)*pow(density,-37./20.)*pow(v_freefall,4.);
}
inline double Calculate_B_Free_Shock_Height(double v_freefall, double m_dot){
    double integral = (39.*sqrt(3.) - 20*pi)/96.; // value of integral from EQ 7a of Wu 1994 DOI: 10.1086/174103
    return pow(v_freefall,3.)*integral/(2.*bremss_const*m_dot);
}

inline double Normalized_Position_Derivative(double vel, double pos, void* eps_s){
    double epsilon_s = *(double *)eps_s;
    double prefactor = pow(free_fall_velocity,3)/(2.*shock_height*bremss_const*specific_accretion);
    double cyclotron = 1/(1+epsilon_s*pow(4.,alpha+beta)*pow(3.,-alpha)*pow(1-vel,alpha)*pow(vel,beta));
    return prefactor*vel*vel*(5.0-8.0*vel)/sqrt(vel*(1-vel))*cyclotron;
}

inline void Dormand_Prince(function<double(double, double, void*)> func, void* args, vector<double>* t, vector<double>* y, double t_bound, vector<double>* t_eval, vector<double>* y_eval, double abs_err, double rel_err, int max_itter){
    double k[n_stages+1];
    k[0] = func((*t)[0], (*y)[0],args);
    double q[order-1];
    y_eval->resize(t_eval->size());

    double dir = (0. < (t_bound-(*t)[0])) - ((t_bound-(*t)[0]) < 0.); // direction of integration

    double tol = abs_err + rel_err*abs((*y)[0]);
    double f_0 = k[0];
    double h_0 = 1e-2*abs((*y)[0])/abs(f_0);
    double f_1 = func((*t)[0]+dir*h_0, (*y)[0]+dir*h_0*f_0, args);
    double delta = abs(f_1-f_0)/(tol*h_0);
    double h_1 = pow(1e-2/max(delta,abs(f_0/tol)),0.2);

    double h = min(1e2*h_0,h_1);

    double y_new, t_new, err, factor, dy, sigma;
    bool step_failed = false;

    while(dir*(t_bound-t->back()) > 0 && t->size() < max_itter){
        h = min(h,dir*(t_bound-t->back()));
        t_new = t->back()+dir*h;

        for(int i = 1; i<n_stages; i++){
            dy = 0;
            for(int j = 0; j<i; j++){
                dy += a[i][j]*k[j];
            }
            k[i] = func(t->back()+c[i]*dir*h, y->back()+dir*h*dy, args);
        }
        y_new = 0;
        err = 0;
        for(int i = 0; i<n_stages; i++){
            y_new += b[i]*k[i];
            err += e[i]*k[i];
        }
        y_new = y->back() + dir*h*y_new;
        k[n_stages] = func(t_new, y_new, args);
        tol = abs_err + rel_err*max(abs(y->back()), abs(y_new));
        err = h*abs(err+e[n_stages]*k[n_stages])/tol;

        if(err < 1.){
            if (err == 0){
                factor = 10.;
            }
            else if(step_failed){
                factor = min(0.9*pow(err,-0.2), 1.);
            }
            else{
                factor = min(0.9*pow(err,-0.2), 10.);
            }
            step_failed = false;
        }
        else if(err >= 1.){
            step_failed = true;
            h *= max(0.9*pow(err,-0.2), 0.2);
            continue;
        }
        else{
            // if err is nan
            h *= 0.2;
            step_failed = true;
            continue;
        }

        if(t_eval->size() > 0){
            for(int i = 0; i<order-1; i++){
                q[i] = 0.;
                for(int j = 0; j<n_stages+1; j++){
                    q[i] += k[j]*p[j][i];
                }
            }
            for(double t_interp:*t_eval){
                if((dir*(t_new-t_interp)<0)|(dir*(t_interp-t->back())<0)){
                    continue;
                }
                sigma = (t_interp-t->back())/(dir*h);;
                dy = 0.;
                for(int i = 0; i < order-1; i++){
                    dy += q[i]*pow(sigma,i+1);
                }
                dy *= dir*h;
                int ind = find(t_eval->begin(), t_eval->end(), t_interp) - t_eval->begin();
                (*y_eval)[ind] = y->back() + dy;
            }
        }

        h *= factor;
        y->push_back(y_new);
        t->push_back(t_new);
        k[0] = k[n_stages];
    }

    if (t->size() == max_itter){
        cout << "Warning steps exceded max_itter, completion is not guaranteed. final (t,y) = (" << t->back() << ", " << y->back() << ')' << endl;
    }
}

inline void Shock_Height_Shooting(double tolerance, int max_itter, bool is_ip){
    double error = 1e100;
    int itters = 0;
    vector<double> norm_velocity;
    vector<double> norm_height;
    vector<double> vel_eval;
    vector<double> pos_eval;
    double slope;

    while(itters < max_itter && error > tolerance){
        norm_velocity = {initial_veloicty};
        norm_height = {initial_height};
        Dormand_Prince(Normalized_Position_Derivative, &epsilon_shock, &norm_velocity, &norm_height, 1e-2*absolute_err, &vel_eval, &pos_eval,
            absolute_err, relative_err, max_itter);
        slope = Normalized_Position_Derivative(norm_velocity.back(), norm_height.back(), &epsilon_shock);
        error = abs(norm_height.back() - slope*norm_velocity.back());

        shock_height *= 1-norm_height.back();
        free_fall_velocity = Calculate_Free_Fall_Velocity(mass, wd_radius, shock_height);
        if(is_ip){
            free_fall_velocity = Calculate_Free_Fall_Velocity(mass, wd_radius, shock_height, mag_radius);
        }
        epsilon_shock = Calculate_Epsilon(free_fall_velocity);
        itters++;
    }
    norm_velocity = {initial_veloicty};
    norm_height = {initial_height};
    double kT = erg_to_kev*((shock_ratio-1)/(shock_ratio*shock_ratio))*kt_const*free_fall_velocity*free_fall_velocity;
    while(kT > 0.){
        vel_eval.push_back((1.0-sqrt(1.0-4.*kT/(erg_to_kev*kt_const*free_fall_velocity*free_fall_velocity)))/2);
        kT -= 1.;
    }
    Dormand_Prince(Normalized_Position_Derivative, &epsilon_shock, &norm_velocity, &norm_height, 1e-2*absolute_err, &vel_eval, &pos_eval,
        absolute_err, relative_err, max_itter);
    velocity.resize(vel_eval.size());
    altitude.resize(vel_eval.size());
    density.resize(vel_eval.size());
    pressure.resize(vel_eval.size());
    temperature.resize(vel_eval.size());
    electron_dens.resize(vel_eval.size());
    for (int i = 0; i < vel_eval.size(); i++){
        velocity[i] = vel_eval[i]*free_fall_velocity;
        altitude[i] = pos_eval[i]*shock_height;
        density[i] = specific_accretion/velocity[i];
        pressure[i] = specific_accretion*free_fall_velocity*(1 - vel_eval[i]);
        temperature[i] = kt_const*free_fall_velocity*free_fall_velocity*vel_eval[i]*(1-vel_eval[i])/boltz_const;
        electron_dens[i] = electron_density_const*density[i]/electron_mass;
    }
}

inline void MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string){

    int n = flux.size();

    double segment_height = 0;

    valarray<double> flux_from_layer(n);
    valarray<double> flux_error(n);

    valarray<double> apec_parameters(3);
    valarray<double> refl_parameters(5);

    refl_parameters[0] = 1-sqrt(1.0-1.0/pow(1+altitude.back()/wd_radius,2));
    refl_parameters[1] = 0.;
    refl_parameters[2] = wd_abund;
    refl_parameters[3] = wd_abund;
    refl_parameters[4] = cos_incl;

    for(int i=0; i<altitude.size(); i++){
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
        if(i==0){
            segment_height = abs(altitude[i+1]-altitude[i])/2.;
        }
        else if(i==altitude.size()-1){
            segment_height = abs(altitude[i]-altitude[i-1])/2.;
        }
        else{
            segment_height = 0.5*abs(altitude[i+1]-altitude[i-1]);
        }
        for(int j=0; j<n; j++){
            flux_from_layer[j] *= (accretion_area*segment_height*electron_dens[i]*electron_dens[i]*1e-14)/avg_atomic_charge;
            flux_from_layer[j] /= (4*pi*source_distance*source_distance);
        }

        // apply reflect to each slice
        if (reflection_sel == 1){
            refl_parameters[0] = 1-sqrt(1.0-1.0/pow(1+altitude[i]/wd_radius,2));
            CXX_reflect(energy, refl_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }
    }

    for(int j=0; j<n; j++){
        flux[j] += flux_from_layer[j];
        flux_from_layer[j] = 0;
    }

    if(reflection_sel == 2){
        CXX_reflect(energy, refl_parameters, spectrum_num, flux, flux_error, init_string);
    }
}

#endif
