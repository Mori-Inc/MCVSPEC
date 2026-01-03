#pragma once

#include "constants.hh"
#include "dipole.hh"
#include "integration.hh"

class Cataclysmic_Variable{
    protected:
        // input properties
        const double mass, radius, b_field, inverse_mag_radius, corotation_ratio, distance;
        const double accretion_rate, accretion_area, metalicity, pressure_ratio, incl_angle;
        double shock_height;
        valarray<double> abundances; // fractional abundance of elements in accretion column
        // derived column properties
        double avg_ion_mass, avg_atomic_charge, density_const, exchange_const, bremss_const, cyclotron_const;
        // boundary conditions
        double w_s, x_s, v_s, pe_s, s_s;
        // thermal profile
        valarray<double> velocity, altitude, electron_temperature, ion_temperature, electron_density, ion_density, electron_pressure, total_pressure, volume;
        // utilities
        const int refl;
        Dipole geometry;
        double upper_bound, lower_bound, upper_landing, lower_landing;
        // unit conversion
        const double length_conv, vel_conv, time_conv, volume_conv, mass_conv, energy_conv, density_conv;


    public:
        Cataclysmic_Variable(double,double,double,double,double,double,double,double,double,double,int);

        void Bracket_Shock_Height();
        void Shock_Height_Shooting();
        void Build_Column_Profile();
        void Print_Properties();
        void Update_Shock_Height(double);
        double Get_Landing_Altitude();
        static double Get_Radius(double);
        static double Get_Accretion_Rate(double, double, double, double);

    protected:
        virtual void Set_Abundances() = 0;
        void Set_Cooling_Constants();
        void Guess_Shock_Height();
        void Flow_Equation(double,const valarray<double>&, valarray<double>&) const;
        Integrator<Cataclysmic_Variable, &Cataclysmic_Variable::Flow_Equation> accretion_column{*this};
};
