#ifndef EOS_H
#define EOS_H

#include "constants.hh"
#include "integration.hh"

const double pressure_constant = pi*pow(electron_mass,4)*pow(light_speed,5)/(3*pow(planck_const,3));
const double density_constant = (8.*pi/3.)*wd_mol_mass*amu_to_g*pow(electron_mass*light_speed/planck_const,3);


valarray<double> Chandrasekhar_White_Dwarf_Equation(double, valarray<double>, void*);

inline double Radius_Shooting(double mass, int max_itter){
    int itters = 0;
    double y_0 = 2;
    double solved_mass = 2*solar_mass;
    Integrator white_dwarf(Chandrasekhar_White_Dwarf_Equation, 2);
    vector<double> eta;
    valarray<double> bounds;

    while(abs(1 - solved_mass/mass) > relative_err && itters < max_itter){
        bounds = {1e-12 + 1./y_0, 1e3};
        white_dwarf.Integrate(&y_0, 1e-32, 1e5, {1.,0.}, bounds);
        eta = white_dwarf.t;
        solved_mass = (-4*pi/(density_constant*density_constant))*pow((2*pressure_constant/(pi*grav_const)),3./2.)*(eta.back()*eta.back())*white_dwarf.y.back()[1];
        y_0 *= cbrt(mass/solved_mass);
        itters++;
    }
    bounds = {1e-10 + 1./y_0, 1e3};
    white_dwarf.Integrate(&y_0, 1e-32, 1e5, {1.,0.}, bounds);
    eta = white_dwarf.t;
    double radius = sqrt(2*pressure_constant/(pi*grav_const*density_constant*density_constant))*(eta.back()/y_0);
    return radius;
}

#endif
