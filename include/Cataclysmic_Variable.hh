#pragma once

#include "dipole.hh"
#include "integration.hh"
#include <vector>

double Mass_to_Radius(double);
double Luminosity_to_Accretion_Rate(double, double, double, double);

class Cataclysmic_Variable{
    protected:
        static constexpr size_t n_dim = 4;
        static constexpr size_t n_grid_vars = 3;
        // input properties
        const double mass, radius, b_field, inverse_mag_radius, corotation_ratio, distance;
        const double accretion_rate, accretion_area, metallicity, pressure_ratio, incl_angle;
        std::vector<double> abundances; // fractional abundance of elements in accretion column
        // derived column properties
        double avg_ion_mass, avg_atomic_charge, mass_to_number_density, exchange_const, bremss_const, cyclotron_const;
        // thermal profile
        std::vector<double> altitude, volume, velocity, density, total_pressure, electron_pressure, electron_density, electron_temperature, ion_temperature;
        // shock boundary
        State<n_dim> shock_boundary;
        double shock_entropy;
        // utilities
        const int refl;
        Dipole geometry;
        // unit conversion
        const double length_conv, vel_conv, time_conv, volume_conv, mass_conv, energy_conv, density_conv, pressure_conv;
        // if solution found
        bool valid_solution = true;


    public:
        Cataclysmic_Variable(double,double,double,double,double,double,double,double,double,double,double,double,int);
        virtual ~Cataclysmic_Variable() = default;

        void Flow_Equation(double,const State<n_dim>&, State<n_dim>&) const;

        void Find_Shock_Position();
        void Build_Column_Profile();

        void Print_Properties() const;

    protected:
        virtual void Set_Abundances() = 0;
        void Set_Cooling_Constants();
        int Bracket_Shock_Position(double&,double&,double&,double&) const;
        template <typename func>
        void Build_Grid(func,const State<n_grid_vars>&,std::vector<State<n_dim>>&) const;
        void Compute_Shock_Bound(double, State<n_dim>&, double&) const;
        double Landing_Altitude(double) const;
        struct Diff_EQ {
            Cataclysmic_Variable& self;
            void operator()(double t, const State<n_dim>& y, State<n_dim>& dydt) const {
              self.Flow_Equation(t, y, dydt);
            }
        };
        const double abs_err = 1e-8;
        const double rel_err = 1e-6;
        mutable Integrator<4, Diff_EQ> accretion_column;
};
