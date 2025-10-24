#pragma once

#include "constants.hh"
#include "dipole.hh"
#include "integration.hh"

class Cataclysmic_Variable{
    protected:
        // input white dwarf properties
        double mass, radius, b_field, inverse_mag_radius, corotation_ratio, distance;
        // input column properties
        double free_fall_speed, scaled_mdot;
        double accretion_rate, accretion_area, metalicity, shock_height, shock_area, pressure_ratio, incl_angle;
        valarray<double> abundances; // fractional abundance of elements in accretion column
        // derived column properties
        double avg_ion_mass, avg_atomic_charge;
        double density_const, coulomb_log_const, exchange_const, bremss_const, cyclotron_const;
        double shock_mdot;
        // thermal profile
        valarray<double> velocity, altitude, electron_temperature, ion_temperature, electron_density, ion_density, electron_pressure, total_pressure, volume;
        // utilities
        Dipole geometry;
        int refl;
        double upper_bound, lower_bound;


    public:
        Cataclysmic_Variable(double,double,double,double,double,double,double,double,double,double,int);

        void Shock_Height_Shooting();
        void Build_Column_Profile();
        void Print_Properties();
        void Update_Shock_Height(double);
        double Get_Landing_Altitude();
        static double Get_Radius(double);
        static double Get_Accretion_Rate(double, double, double, double);

    protected:
        virtual void Set_Abundances(double) = 0;
        void Set_Cooling_Constants();
        void Guess_Shock_Height();
        void Flow_Equation(double,const valarray<double>&, valarray<double>&) const;
        Integrator<Cataclysmic_Variable, &Cataclysmic_Variable::Flow_Equation> accretion_column{*this};
};
