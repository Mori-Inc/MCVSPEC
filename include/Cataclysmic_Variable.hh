#pragma once

#include "constants.hh"
#include "dipole.hh"
#include "integration.hh"

class Cataclysmic_Variable{
    protected:
        static constexpr size_t n_dim = 4;
        // input properties
        const double mass, radius, b_field, inverse_mag_radius, corotation_ratio, distance;
        const double accretion_rate, accretion_area, metalicity, pressure_ratio, incl_angle;
        double shock_height;
        vector<double> abundances; // fractional abundance of elements in accretion column
        // derived column properties
        double avg_ion_mass, avg_atomic_charge, density_const, exchange_const, bremss_const, cyclotron_const;
        // boundary conditions
        double w_s, x_s, v_s, pe_s, s_s;
        // thermal profile
        vector<double> altitude, volume, velocity, density, total_pressure, electron_pressure, electron_density, electron_temperature, ion_temperature;
        // utilities
        const int refl;
        Dipole geometry;
        // unit conversion
        const double length_conv, vel_conv, time_conv, volume_conv, mass_conv, energy_conv, density_conv;
        // if solution found
        bool valid_solution = true;


    public:
        Cataclysmic_Variable(double,double,double,double,double,double,double,double,double,double,double,double,int);

        void Bracket_Shock_Position(double&,double&,double&,double&);
        void Determine_Shock_Position();
        template <typename func>
        void Build_Grid(func,const vector<double>&,vector<State<n_dim>>&);
        void Build_Column_Profile();
        void Print_Properties();
        void Update_Shock_Position(double);
        double Get_Landing_Altitude(double);
        static double Get_Radius(double);
        static double Get_Accretion_Rate(double, double, double, double);
        void Flow_Equation(double,const State<n_dim>&, State<n_dim>&) const;

    protected:
        virtual void Set_Abundances() = 0;
        void Set_Cooling_Constants();
        struct Diff_EQ {
            Cataclysmic_Variable& self;
            void operator()(double t, const State<n_dim>& y, State<n_dim>& dydt) const {
              self.Flow_Equation(t, y, dydt);
            }
        };
        Integrator<4, Diff_EQ> accretion_column;
};
