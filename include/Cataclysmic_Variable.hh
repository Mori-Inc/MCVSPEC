#pragma once

#include "constants.hh"
#include "integration.hh"
#include <valarray>
#include <xsTypes.h>
#include <funcWrappers.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

class Cataclysmic_Variable{
    protected:
        // input white dwarf properties
        double mass, b_field, radius, inverse_mag_radius, distance;
        double non_dim_radius;
        // input column properties
        double accretion_rate, accretion_area, metalicity, shock_height, shock_speed, pressure_ratio, incl_angle, area_exponent;
        valarray<double> abundances; // fractional abundance of elements in accretion column
        // derived column properties
        double avg_ion_mass, avg_atomic_charge;
        double density_const, force_const, cooling_ratio_const, coulomb_log_const, exchange_const, bremss_const;
        double cooling_ratio, shock_mdot;
        // thermal profile
        valarray<double> velocity, altitude, electron_temperature, ion_temperature, electron_density, ion_density, electron_pressure, total_pressure;
        // utilities
        int refl;
        Integrator accretion_column;

    public:
        Cataclysmic_Variable(double,double,double,double,double,double,double,int); // base constructor
        Cataclysmic_Variable(double,double,double,double,double,double,double,double,int); // polars
        Cataclysmic_Variable(double,double,double,double,double,double,double,int,double); // ips
        virtual void Set_Abundances(double);
        void Set_Accretion_Rate(double);
        void Set_Inverse_Mag_Radius(double);
        valarray<double> Flow_Equation(double vel, valarray<double> pos);
        void Shock_Height_Shooting();
        void Build_Column_Profile();
        virtual void MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string);
        void Print_Properties();

    protected:
        void Set_Cooling_Constants();
        void Set_Radius();
        void Guess_Shock_Height();
        void Update_Shock_Height(double);
        double Get_Landing_Altitude();
};
