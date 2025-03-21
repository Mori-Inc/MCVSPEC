#ifndef CV_H
#define CV_H

#include "constants.hh"
#include "integration.hh"
#include <xsTypes.h>
#include <funcWrappers.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

class Cataclysmic_Variable{
    protected:
        // input white dwarf properties
        double mass, b_field, radius, inverse_mag_radius, distance;
        // input column properties
        double accretion_rate, accretion_area, metalicity, shock_height, pre_shock_speed, pressure_ratio, incl_angle;
        valarray<double> abundances, charge; // fractional abundance of elements in accretion column
        // derived column properties
        double cooling_ratio, bremss_constant, exchange_constant, avg_ion_mass, avg_atomic_charge;
        double thermal_constant, cooling_ratio_const, kt_const, electron_density_const;
        // thermal profile
        valarray<double> altitude, electron_temperature, ion_temperature, electron_density, ion_density;
        // utilities
        int refl;
        Integrator accretion_column;

    public:
        Cataclysmic_Variable(double,double,double,double,double,double,double,int); // polars
        Cataclysmic_Variable(double,double,double,double,double,double,int,double); // ips
        void Set_Abundances(double);
        void Set_Accretion_Rate(double);
        void Set_Inverse_Mag_Radius(double);
        static valarray<double> Flow_Equation(double vel, valarray<double> pos, void* my_class_instance);
        void Shock_Height_Shooting(int);
        void MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string);
        void Print_Properties();

    private:
        void Set_Cooling_Constants();
        static valarray<double> Chandrasekhar_White_Dwarf_Equation(double, valarray<double>, void*);
        void Radius_Shooting(int);
        void Set_Pre_Shock_Speed(int);
        void Set_Cooling_Ratio();
};

#endif
