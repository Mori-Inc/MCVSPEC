#pragma once

#include "integration.hh"
#include "Cataclysmic_Variable.hh"
#include <vector>

class Saxton_CV{ // no dipole, two-temp, gravity, cyclotron
    public:
        static constexpr size_t n_dim = 3;
        static constexpr size_t n_grid_vars = 3;
    protected:
        // input properties
        const White_Dwarf white_dwarf;
        const Accretion_Column accretion_column;
        const Tolerance error_control;
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
        double shock_speed;
        // if solution found
        bool valid_solution = true;


    public:
        Saxton_CV(White_Dwarf, Accretion_Column, Tolerance);

        void Flow_Equation(double,const State<n_dim>&, State<n_dim>&) const;

        void Find_Shock_Position();
        void Build_Column_Profile();

        void Print_Properties() const;

        int Solve_Profile();

    protected:
        void Set_Abundances();
        void Set_Cooling_Constants();
        int Bracket_Shock_Position(double&,double&,double&,double&) const;
        template <typename func>
        void Build_Grid(func,const State<n_grid_vars>&,std::vector<State<n_dim>>&, std::vector<double>& ) const;
        void Compute_Shock_Bound(double, State<n_dim>&, double&) const;
        double Landing_Altitude(double) const;
        struct Diff_EQ {
            Saxton_CV& self;
            void operator()(double t, const State<n_dim>& y, State<n_dim>& dydt) const {
              self.Flow_Equation(t, y, dydt);
            }
        };
        mutable Integrator<n_dim, Diff_EQ> integrator;

    public:
        const double Get_Mass() const {return white_dwarf.mass;}
        const double Get_B_Field() const {return white_dwarf.b_field;}
        const double Get_inv_Mag_Radius() const {return white_dwarf.inverse_mag_radius;}
        const double Get_Corotation_Radius() const {return white_dwarf.corotation_radius;}
        const double Get_Distance() const {return white_dwarf.distance;}
        const double Get_Mdot() const {return accretion_column.accretion_rate;}
        const double Get_Area() const {return accretion_column.accretion_area;}
        const double Get_Abund() const {return accretion_column.metallicity;}
        const double Get_Shock_Ratio() const {return accretion_column.shock_pressure_ratio;}
        const double Get_Cos_Inclination_Angle() const {return white_dwarf.cos_inclination;}
        const double Get_mBar() const {return avg_ion_mass;}
        const double Get_ZBar() const {return avg_atomic_charge;}
        const double Get_Dens_to_ne() const {return mass_to_number_density;}
        const double Get_Exch_Const() const {return exchange_const;}
        const double Get_Bremss_Const() const {return bremss_const;}
        const double Get_Cycl_Const() const {return cyclotron_const;}
        const double Get_Length_Conv() const {return length_conv;}
        const double Get_Mass_Conv() const {return mass_conv;}
        const double Get_Time_Conv() const {return time_conv;}
        const double Get_Vel_Conv() const {return vel_conv;}
        const double Get_Vol_Conv() const {return volume_conv;}
        const double Get_Energy_Conv() const {return energy_conv;}
        const double Get_Density_Conv() const {return density_conv;}
        const double Get_Column_Coord() const {return 0;}
        const std::vector<double>& Get_Position() const {return position;}
        const std::vector<double>& Get_Altitude() const {return altitude;}
        const std::vector<double>& Get_Volume_Element() const {return volume_element;}
        const std::vector<double>& Get_Velocity() const {return velocity;}
        const std::vector<double>& Get_Density() const {return density;}
        const std::vector<double>& Get_Pressure() const {return total_pressure;}
        const std::vector<double>& Get_Electron_Pressure() const {return electron_pressure;}
        const std::vector<double>& Get_Electron_Density() const {return electron_density;}
        const std::vector<double>& Get_Electron_Temperature() const {return electron_temperature;}
        const std::vector<double>& Get_Ion_Temperature() const {return ion_temperature;}
        const std::vector<double>& Get_Abundance() const {return abundances;}
        State<n_dim> state{}, deriv{};
};

class Cropper_CV{ // gravity, cyclotron
    public:
        static constexpr size_t n_dim = 2;
        static constexpr size_t n_grid_vars = 2;
    protected:
        // input properties
        const White_Dwarf white_dwarf;
        const Accretion_Column accretion_column;
        const Tolerance error_control;
        // unit conversion (cgs values of MCVSPEC nd unit system)
        const double length_conv, vel_conv, accretion_rate_conv;
        const double time_conv, mass_conv, volume_conv, energy_conv, density_conv, pressure_conv;

        // derived column properties
        std::vector<double> abundances; // fractional abundance of elements in accretion column
        double avg_ion_mass, avg_atomic_charge, mass_to_number_density, bremss_const, cyclotron_const;
        // thermal profile
        std::vector<double> position, altitude, volume_element, velocity, density, total_pressure, electron_pressure, electron_density, electron_temperature, ion_temperature;
        // shock boundary condition
        State<n_dim> shock_boundary;
        double shock_speed;
        // if solution found
        bool valid_solution = true;


    public:
        Cropper_CV(White_Dwarf, Accretion_Column, Tolerance);

        void Flow_Equation(double,const State<n_dim>&, State<n_dim>&) const;

        void Find_Shock_Position();
        void Build_Column_Profile();

        void Print_Properties() const;

        int Solve_Profile();

    protected:
        void Set_Abundances();
        void Set_Cooling_Constants();
        int Bracket_Shock_Position(double&,double&,double&,double&) const;
        template <typename func>
        void Build_Grid(func,const State<n_grid_vars>&,std::vector<State<n_dim>>&, std::vector<double>& ) const;
        void Compute_Shock_Bound(double, State<n_dim>&, double&) const;
        double Landing_Altitude(double) const;
        struct Diff_EQ {
            Cropper_CV& self;
            void operator()(double t, const State<n_dim>& y, State<n_dim>& dydt) const {
              self.Flow_Equation(t, y, dydt);
            }
        };
        mutable Integrator<n_dim, Diff_EQ> integrator;

    public:
        const double Get_Mass() const {return white_dwarf.mass;}
        const double Get_B_Field() const {return white_dwarf.b_field;}
        const double Get_inv_Mag_Radius() const {return white_dwarf.inverse_mag_radius;}
        const double Get_Corotation_Radius() const {return white_dwarf.corotation_radius;}
        const double Get_Distance() const {return white_dwarf.distance;}
        const double Get_Mdot() const {return accretion_column.accretion_rate;}
        const double Get_Area() const {return accretion_column.accretion_area;}
        const double Get_Abund() const {return accretion_column.metallicity;}
        const double Get_Shock_Ratio() const {return accretion_column.shock_pressure_ratio;}
        const double Get_Cos_Inclination_Angle() const {return white_dwarf.cos_inclination;}
        const double Get_mBar() const {return avg_ion_mass;}
        const double Get_ZBar() const {return avg_atomic_charge;}
        const double Get_Dens_to_ne() const {return mass_to_number_density;}
        const double Get_Exch_Const() const {return 0;}
        const double Get_Bremss_Const() const {return bremss_const;}
        const double Get_Cycl_Const() const {return cyclotron_const;}
        const double Get_Length_Conv() const {return length_conv;}
        const double Get_Mass_Conv() const {return mass_conv;}
        const double Get_Time_Conv() const {return time_conv;}
        const double Get_Vel_Conv() const {return vel_conv;}
        const double Get_Vol_Conv() const {return volume_conv;}
        const double Get_Energy_Conv() const {return energy_conv;}
        const double Get_Density_Conv() const {return density_conv;}
        const double Get_Column_Coord() const {return 0;}
        const std::vector<double>& Get_Position() const {return position;}
        const std::vector<double>& Get_Altitude() const {return altitude;}
        const std::vector<double>& Get_Volume_Element() const {return volume_element;}
        const std::vector<double>& Get_Velocity() const {return velocity;}
        const std::vector<double>& Get_Density() const {return density;}
        const std::vector<double>& Get_Pressure() const {return total_pressure;}
        const std::vector<double>& Get_Electron_Pressure() const {return electron_pressure;}
        const std::vector<double>& Get_Electron_Density() const {return electron_density;}
        const std::vector<double>& Get_Electron_Temperature() const {return electron_temperature;}
        const std::vector<double>& Get_Ion_Temperature() const {return ion_temperature;}
        const std::vector<double>& Get_Abundance() const {return abundances;}
        State<n_dim> state{}, deriv{};
};


















class Wu_CV{ // gravity, cyclotron
    public:
        static constexpr size_t n_dim = 1;
        static constexpr size_t n_grid_vars = 2;
    protected:
        // input properties
        const White_Dwarf white_dwarf;
        const Accretion_Column accretion_column;
        const Tolerance error_control;
        // unit conversion (cgs values of MCVSPEC nd unit system)
        const double length_conv, vel_conv, accretion_rate_conv;
        const double time_conv, mass_conv, volume_conv, energy_conv, density_conv, pressure_conv;

        // derived column properties
        std::vector<double> abundances; // fractional abundance of elements in accretion column
        double avg_ion_mass, avg_atomic_charge, mass_to_number_density, bremss_const, cyclotron_const;
        // thermal profile
        std::vector<double> position, altitude, volume_element, velocity, density, total_pressure, electron_pressure, electron_density, electron_temperature, ion_temperature;
        // shock boundary condition
        State<n_dim> shock_boundary;
        double shock_speed;
        // if solution found
        bool valid_solution = true;


    public:
        Wu_CV(White_Dwarf, Accretion_Column, Tolerance);

        void Flow_Equation(double,const State<n_dim>&, State<n_dim>&) const;

        void Find_Shock_Position();
        void Build_Column_Profile();

        void Print_Properties() const;

        int Solve_Profile();

        double vff;

    protected:
        void Set_Abundances();
        void Set_Cooling_Constants();
        int Bracket_Shock_Position(double&,double&,double&,double&);
        template <typename func>
        void Build_Grid(func,const State<n_grid_vars>&,std::vector<State<n_dim>>&, std::vector<double>& ) const;
        void Compute_Shock_Bound(double, State<n_dim>&, double&);
        double Landing_Altitude(double);
        struct Diff_EQ {
            Wu_CV& self;
            void operator()(double t, const State<n_dim>& y, State<n_dim>& dydt) const {
              self.Flow_Equation(t, y, dydt);
            }
        };
        mutable Integrator<n_dim, Diff_EQ> integrator;

    public:
        const double Get_Mass() const {return white_dwarf.mass;}
        const double Get_B_Field() const {return white_dwarf.b_field;}
        const double Get_inv_Mag_Radius() const {return white_dwarf.inverse_mag_radius;}
        const double Get_Corotation_Radius() const {return white_dwarf.corotation_radius;}
        const double Get_Distance() const {return white_dwarf.distance;}
        const double Get_Mdot() const {return accretion_column.accretion_rate;}
        const double Get_Area() const {return accretion_column.accretion_area;}
        const double Get_Abund() const {return accretion_column.metallicity;}
        const double Get_Shock_Ratio() const {return accretion_column.shock_pressure_ratio;}
        const double Get_Cos_Inclination_Angle() const {return white_dwarf.cos_inclination;}
        const double Get_mBar() const {return avg_ion_mass;}
        const double Get_ZBar() const {return avg_atomic_charge;}
        const double Get_Dens_to_ne() const {return mass_to_number_density;}
        const double Get_Exch_Const() const {return 0;}
        const double Get_Bremss_Const() const {return bremss_const;}
        const double Get_Cycl_Const() const {return cyclotron_const;}
        const double Get_Length_Conv() const {return length_conv;}
        const double Get_Mass_Conv() const {return mass_conv;}
        const double Get_Time_Conv() const {return time_conv;}
        const double Get_Vel_Conv() const {return vel_conv;}
        const double Get_Vol_Conv() const {return volume_conv;}
        const double Get_Energy_Conv() const {return energy_conv;}
        const double Get_Density_Conv() const {return density_conv;}
        const double Get_Column_Coord() const {return 0;}
        const std::vector<double>& Get_Position() const {return position;}
        const std::vector<double>& Get_Altitude() const {return altitude;}
        const std::vector<double>& Get_Volume_Element() const {return volume_element;}
        const std::vector<double>& Get_Velocity() const {return velocity;}
        const std::vector<double>& Get_Density() const {return density;}
        const std::vector<double>& Get_Pressure() const {return total_pressure;}
        const std::vector<double>& Get_Electron_Pressure() const {return electron_pressure;}
        const std::vector<double>& Get_Electron_Density() const {return electron_density;}
        const std::vector<double>& Get_Electron_Temperature() const {return electron_temperature;}
        const std::vector<double>& Get_Ion_Temperature() const {return ion_temperature;}
        const std::vector<double>& Get_Abundance() const {return abundances;}
        State<n_dim> state{}, deriv{};
};
