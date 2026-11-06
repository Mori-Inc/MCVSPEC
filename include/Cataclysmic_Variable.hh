#pragma once

#include "dipole.hh"
#include "integration.hh"
#include <vector>

struct White_Dwarf{
    double mass;
    double radius;
    double b_field;
    double cos_inclination;
    double inverse_mag_radius;
    double corotation_radius;
    double distance;
    White_Dwarf(double m, double r, double b, double cosi, double irm, double cort, double d):
        mass(m), radius(r), b_field(b), cos_inclination(cosi), inverse_mag_radius(irm), corotation_radius(cort), distance(d)
    {};
    White_Dwarf() = default;
};

struct Accretion_Column{
    double accretion_rate;
    double accretion_area;
    double metallicity;
    double shock_pressure_ratio;
    double sin_mag_colat;
    Accretion_Column(double mdot, double a, double m, double delt, double sinb):
        accretion_rate(mdot), accretion_area(a), metallicity(m), shock_pressure_ratio(delt), sin_mag_colat(sinb)
    {};
    Accretion_Column() = default;
};

struct Tolerance{
    double absolute_error;
    double relative_error;
    double kT_grid_spacing;
    double altitude_grid_spacing;
    Tolerance(double abserr, double relerr, double dkT, double dz):
        absolute_error(abserr), relative_error(relerr), kT_grid_spacing(dkT), altitude_grid_spacing(dz)
    {};
    Tolerance() = default;
};

double Mass_to_Radius(double);
double Luminosity_to_Accretion_Rate(double, White_Dwarf);

class Cataclysmic_Variable{
    protected:
        static constexpr size_t n_dim = 4;
        static constexpr size_t n_grid_vars = 3;
        // input properties
        const White_Dwarf white_dwarf;
        const Accretion_Column accretion_column;
        const Tolerance error_control;
        const Dipole geometry;
        // unit conversion (cgs values of MCVSPEC nd unit system)
        const double length_conv, vel_conv, accretion_rate_conv;
        const double time_conv, mass_conv, volume_conv, energy_conv, density_conv, pressure_conv;

        // derived column properties
        std::vector<double> abundances; // fractional abundance of elements in accretion column
        double avg_ion_mass, avg_atomic_charge, mass_to_number_density, exchange_const, bremss_const, cyclotron_const;
        // thermal profile
        std::vector<double> position, altitude, volume_element, velocity, density, total_pressure, electron_pressure, electron_density, electron_temperature, ion_temperature;
        // shock boundary condition
        State<n_dim> shock_boundary;
        double shock_entropy;
        // if solution found
        bool valid_solution = true;


    public:
        Cataclysmic_Variable(White_Dwarf, Accretion_Column, Tolerance);
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
        mutable Integrator<4, Diff_EQ> integrator;
};
