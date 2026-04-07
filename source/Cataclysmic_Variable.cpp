#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "integration.hh"
#include "mass_radius.hh"
#include "gaunt.hh"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

double Luminosity_to_Accretion_Rate(double luminosity, White_Dwarf wd){
    double accretion_rate = luminosity/(grav_const*wd.mass*((1./wd.radius) - wd.inverse_mag_radius));
    return accretion_rate;
}

double Mass_to_Radius(double mass){
    int left_ind = 0;
    int i = mass_radius_length/2;
    int right_ind = mass_radius_length-1;
    while(right_ind-left_ind > 1){
        i = left_ind + (right_ind-left_ind)/2;
        if(mass>white_dwarf_mass[i]){
            left_ind = i;
        }
        else{
            right_ind = i;
        }
    }
    double delta_r = white_dwarf_radius[right_ind]-white_dwarf_radius[left_ind];
    double delta_m = white_dwarf_mass[right_ind]-white_dwarf_mass[left_ind];
    double radius = white_dwarf_radius[left_ind] + (delta_r/delta_m)*(mass-white_dwarf_mass[left_ind]);
    return radius;
}

Cataclysmic_Variable::Cataclysmic_Variable(White_Dwarf wd, Accretion_Column col, Tolerance tol):
    white_dwarf(wd), accretion_column(col), error_control(tol), geometry(accretion_column.sin_mag_colat*accretion_column.sin_mag_colat),
    length_conv(white_dwarf.radius), vel_conv(sqrt(2*grav_const*white_dwarf.mass/white_dwarf.radius)), accretion_rate_conv(geometry.a_0*col.accretion_rate),
    time_conv(length_conv/vel_conv), mass_conv(accretion_rate_conv*length_conv*length_conv*time_conv), volume_conv(length_conv*length_conv*length_conv),
    energy_conv(mass_conv*vel_conv*vel_conv), density_conv(mass_conv/volume_conv), pressure_conv(energy_conv/volume_conv),
    integrator(Diff_EQ{*this},error_control.absolute_error,error_control.relative_error)
{}

void Cataclysmic_Variable::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
    // constants related to column composition
    avg_ion_mass = 0;
    avg_atomic_charge = 0;
    double avg_charge_squared = 0;
    double avg_charge_sqr_over_mass = 0;

    for(size_t i=0; i<abundances.size(); i++){
        avg_ion_mass += abundances[i]*atomic_mass[i]*amu_to_g;
        avg_atomic_charge += abundances[i]*atomic_charge[i];
        avg_charge_squared += abundances[i]*atomic_charge[i]*atomic_charge[i];
        avg_charge_sqr_over_mass += abundances[i]*atomic_charge[i]*atomic_charge[i]/(atomic_mass[i]*amu_to_g);
    }

    mass_to_number_density = avg_atomic_charge/(avg_ion_mass + avg_atomic_charge*m_e);
    const double j_e = mass_to_number_density*accretion_column.accretion_rate*geometry.a_0;
    const double b_sqr = white_dwarf.b_field*white_dwarf.b_field;
    bremss_const = bremss_coeff*(avg_charge_squared/avg_atomic_charge)*sqrt(j_e*j_e*j_e);
    bremss_const *= sqrt(time_conv*time_conv*time_conv*mass_conv)*length_conv/energy_conv;
    cyclotron_const = cyclotron_coeff*(avg_atomic_charge/avg_charge_squared);
    cyclotron_const *= b_sqr*pow(b_sqr/(geometry.a_0*accretion_column.accretion_area), 0.425)/(j_e*j_e);
    cyclotron_const *= pow(geometry.a_0*geometry.a_0/j_e, 1.85);
    cyclotron_const *= pressure_conv*pressure_conv*pow(vel_conv,3.85);
    exchange_const = exchange_coeff*avg_charge_sqr_over_mass*sqrt(j_e*j_e*j_e*j_e*j_e);
    exchange_const *= length_conv*length_conv*time_conv*time_conv*time_conv*sqrt(time_conv/energy_conv)/energy_conv;
}

void Cataclysmic_Variable::Compute_Shock_Bound(double w_s, State<n_dim>& bound, double& s_s) const{
    double r_s, proj_r_w, convergance, scale_factors[3];
    geometry.update_coordinates(w_s, r_s, proj_r_w, convergance, scale_factors);
    double mdot = 1./(scale_factors[0]*scale_factors[2]);
    double vff = sqrt(1./r_s - length_conv*white_dwarf.inverse_mag_radius);
    bound[0] = w_s;
    bound[1] = vff;
    bound[2] = vff*0.25;
    bound[3] = (accretion_column.shock_pressure_ratio/(1+accretion_column.shock_pressure_ratio))*mdot*0.75*vff;
    s_s = mdot*0.75*vff*pow(0.25*vff/mdot, 5./3.);
}

void Cataclysmic_Variable::Flow_Equation(double entropy,const State<n_dim>& state, State<n_dim>& derivs) const{
    const double& s = entropy;
    const double& w = state[0];
    const double& x = state[1];
    const double& v = state[2];
    const double& pe = state[3];

    double r, proj_r_w, convergance, scale_factors[3];
    geometry.update_coordinates(w, r, proj_r_w, convergance, scale_factors);

    const double& h = scale_factors[0];
    const double area = scale_factors[1]*scale_factors[2];
    const double mdot = 1./area;
    const double p = mdot*(x-v);


    const double ne_cgs = mass_to_number_density*density_conv*mdot/v;
    const double kT_cgs = pressure_conv*pe/ne_cgs;
    const double gff = gaunt::gaunt_factor(kT_cgs);
    const double coulomb_log = 0.5*log(coulomb_log_coeff*kT_cgs*kT_cgs/ne_cgs);

    const double cyc = (cyclotron_const/gff)*pe*pe*v*v*v*h*pow(v*v/h, 0.425);
    const double rad = bremss_const*gff*sqrt(pe/(h*h*h*v*v*v))*(1+cyc);
    const double exch = (exchange_const*coulomb_log/(bremss_const*gff))*(avg_atomic_charge*(p-pe)/pe - 1.0)/(pe*v*h); // ratio of exchange to radiation
    const double grav = 0.5*proj_r_w/(r*r)/rad;
    const double conv = v*convergance/(h*rad);

    double dw_ds = -1.5*p*v/(s*h*rad);
    double dx_ds = -1.5*(p/s)*(p*conv/mdot - grav);
    double dv_ds = 1.5*(p*v/s)*(2.0/mdot + 5*p*conv/mdot -3*grav)/(5*x - 8*v);
    double dp_ds = (p/s)*(1.0 - exch + 5*pe*(1.5*grav - 1.0/mdot - v*conv)/(5*x - 8*v));

    derivs[0] = dw_ds;
    derivs[1] = dx_ds;
    derivs[2] = dv_ds;
    derivs[3] = dp_ds;
}

double Cataclysmic_Variable::Landing_Altitude(double w_s) const {
    // return signed distance from WD surface in w
    double s_s;
    State<n_dim> y;
    Compute_Shock_Bound(w_s, y, s_s);
    double s = s_s;
    integrator.Initialize(s, 0, y);
    while(s/s_s > 1e-2){
         integrator.Step(s, y);
    }
    double error = 1;
    State<n_dim> slope{};
    Flow_Equation(s,  y, slope);
    double s_prev = s;
    double dw_prev = slope[0];

    while(error > error_control.absolute_error){
        integrator.Step(s, y);
        Flow_Equation(s,  y, slope);
        error = 0.5*std::abs((slope[0]-dw_prev)/(s-s_prev))*s*s; //difference between linear and quadratic extrapolation on w
        s_prev = s;
        dw_prev = slope[0];
    }
    integrator.Step(s, y);
    Flow_Equation(s,  y, slope);
    double landing = y[0] - slope[0]*s;
    return landing - geometry.w_0;
}

// Bracket the shock height by minimizing the landing coordinate
// (which is maximizing the landing altitude since w ~ -r close to the surface and at the pole)
// exit if
// 1. a point is found with w_l < w_0, in which case we have bracketed our solution
// 2. a minimum is found with w_l > w_0 in which case no solution exists
int Cataclysmic_Variable::Bracket_Shock_Position(double& upper_bound, double& lower_bound, double& upper_landing, double& lower_landing) const {
    upper_bound = geometry.w_0;
    upper_landing = Landing_Altitude(upper_bound);

    double r=1.;
    double dr = 0.01;
    auto r_to_w = [this](double r){ return sqrt(1-geometry.u*r)/(r*r); };
    double samples[3] = {Landing_Altitude(r_to_w(r-dr)),
                        upper_landing,
                        Landing_Altitude(r_to_w(r+dr))};

    double dwl_drs[2] = {(samples[1]-samples[0])/dr, (samples[2]-samples[1])/dr};
    double step = 0;
    while(dwl_drs[0]*dwl_drs[1] > 0){ // while minima not bounded
        if(samples[2]<0){
            lower_bound = r_to_w(r+dr);
            lower_landing = samples[2];
            return 0;
        }
        upper_bound = r_to_w(r+dr);
        upper_landing = samples[2];
        step = 0.5*dr*(samples[0]-samples[2])/(samples[0]-2*samples[1]+samples[2]);
        r += std::min(1.,step);
        samples[0] = Landing_Altitude(r_to_w(r-dr));
        samples[1] = Landing_Altitude(r_to_w(r));
        samples[2] = Landing_Altitude(r_to_w(r+dr));
        dwl_drs[0] = (samples[1]-samples[0])/dr;
        dwl_drs[1] = (samples[2]-samples[1])/dr;
    }

    if(samples[1] > 0){ // if minima > 0
        return 1;
    }

    lower_bound = r_to_w(r);
    lower_landing = samples[1];
    return 0;
}

void Cataclysmic_Variable::Find_Shock_Position(){
    double upper_bound, lower_bound;
    double upper_landing, lower_landing;
    int success = Bracket_Shock_Position(upper_bound, lower_bound, upper_landing, lower_landing);
    // if no minimum skip
    if(success != 0){
        valid_solution = false;
        return;
    }
    const double k1 = 0.2/(upper_bound-lower_bound);
    const double n0 = 1;
    double nmax = log2((upper_bound-lower_bound)/(2*error_control.absolute_error)) + n0;
    int i=0;
    double new_bound, new_altitude, midpoint, regula_falsi, truncation, projection, dir;
    while(upper_bound-lower_bound > error_control.absolute_error){
        midpoint = (upper_bound+lower_bound)/2;
        regula_falsi = (upper_landing*lower_bound - lower_landing*upper_bound)/(lower_landing-upper_landing);
        dir = (0. < (midpoint-regula_falsi)) - ((midpoint-regula_falsi) < 0.);
        truncation = k1*(upper_bound-lower_bound)*(upper_bound-lower_bound); // k2 = 2

        if(truncation <= std::abs(midpoint-regula_falsi)){
            new_bound = regula_falsi + dir*truncation;
        }
        else{
            new_bound = midpoint;
        }

        projection = error_control.absolute_error*(pow(2,nmax-i)) - (upper_bound-lower_bound)/2;
        if(std::abs(new_bound-midpoint) > projection){
            new_bound = midpoint - dir*projection;
        }

        new_altitude = Landing_Altitude(new_bound);
        if(new_altitude>0){
            upper_bound = new_bound;
            upper_landing = new_altitude;
        }
        else if(new_altitude<0){
            lower_bound = new_bound;
            lower_landing = new_altitude;
        }
        else{
            upper_bound = new_bound;
            lower_bound = new_bound;
        }
        i++;
    }
    State<n_dim> bound{};
    double s_s;
    Compute_Shock_Bound(0.5*(upper_bound+lower_bound), bound, s_s);
    shock_boundary = bound;
    shock_entropy = s_s;
}

template <typename func>
void Cataclysmic_Variable::Build_Grid(func grid_func, const State<n_grid_vars>& grid_spacing, std::vector<State<n_dim>>& grid) const {
    if(!valid_solution){
        return;
    }
    constexpr int n_segments=16;
    // integration variables
    double t = shock_entropy;
    State<n_dim> y = shock_boundary;
    integrator.Initialize(t, 0, y);

    // grid variables
    State<n_grid_vars> grid_vars{};
    grid_func(t,y,grid_vars);
    grid.reserve(size_t(std::abs(grid_vars[0]/grid_spacing[0]))+1);
    grid.push_back(y);

    // segment variables
    double t_left, t_right, t_cross;
    State<n_dim> y_left, y_right;
    State<n_grid_vars> grid_left, grid_right;
    y_right = y;
    grid_right = grid_vars;

    auto crossing_found = [&grid_spacing, &grid_vars, &grid_left, &grid_right, &t_cross](const double t_left, const double t_right){
        double delta_t, distance_right, distance_left;
        double min_dt = t_right-t_left;
        bool found=false;
        for(size_t i=0; i<n_grid_vars; ++i){
            distance_right = std::abs((grid_right[i]-grid_vars[i])/grid_spacing[i]);
            distance_left = std::abs((grid_left[i]-grid_vars[i])/grid_spacing[i]);
            delta_t = (t_right-t_left)*(1-distance_left)/(distance_right-distance_left);
            if(distance_right>1 && std::abs(delta_t)<std::abs(min_dt)){
                found=true;
                min_dt = delta_t;
            }
        }
        if(!found){
            return false;
        }
        t_cross = t_left + min_dt;
        return true;
    };

    while(grid_vars[0] > grid_spacing[0] && t > error_control.absolute_error){
        double t0 = t;
        integrator.Dense_Step(t,y);
        double dt = (t-t0)/double(n_segments);
        for(size_t seg=0; seg<n_segments; seg++){
            t_left = t0 + dt*seg;
            t_right = t0 + dt*(seg+1);
            y_left = y_right;
            integrator.Interpolate(t_right, y_right);
            grid_left = grid_right;
            grid_func(t_right, y_right, grid_right);
            while(crossing_found(t_left, t_right)){
                t_left = t_cross;
                integrator.Interpolate(t_left, y_left);
                grid_func(t_left, y_left, grid_left);
                grid.push_back(y_left);
                grid_vars = grid_left;
            }
        }
    }
}

void Cataclysmic_Variable::Build_Column_Profile(){
    if(!valid_solution){
        return;
    }
    const double dkTe = (error_control.kT_grid_spacing/erg_to_kev)*mass_to_number_density/(vel_conv*vel_conv);
    const double dkTi = dkTe/avg_atomic_charge;
    double r, proj, conv, scale_factors[3];
    geometry.update_coordinates(shock_boundary[0], r, proj, conv, scale_factors);
    const State<n_grid_vars> grid_spacing = {dkTi, dkTe, error_control.altitude_grid_spacing*(r-1)};

    auto kT_e = [](const State<n_dim>& y, double scale_factors[3]){
        return y[2]*y[3]*scale_factors[0]*scale_factors[2];
    };
    auto kT_i = [](const State<n_dim>& y, double scale_factors[3]){
        return y[2]*(y[1]-y[2]-y[3]*scale_factors[0]*scale_factors[2]);
    };
    auto grid_func = [this, kT_e, kT_i, &r, &proj, &conv, &scale_factors](const double& t, const State<n_dim>& y, State<n_grid_vars>& vars){
        geometry.update_coordinates(y[0], r, proj, conv, scale_factors);
        vars[0] = kT_i(y, scale_factors);
        vars[1] = kT_e(y, scale_factors);
        vars[2] = r;
    };

    std::vector<State<n_dim>> grid;
    Build_Grid(grid_func, grid_spacing, grid);

    int n_points = grid.size();
    position.resize(n_points);
    volume_element.resize(n_points);
    velocity.resize(n_points);
    altitude.resize(n_points);
    total_pressure.resize(n_points);
    electron_pressure.resize(n_points);
    electron_density.resize(n_points);
    density.resize(n_points);
    electron_temperature.resize(n_points);
    ion_temperature.resize(n_points);

    double mdot;

    for(size_t i=0; i<n_points; i++){
        geometry.update_coordinates(grid[i][0], r, proj, conv, scale_factors);
        position[i] = grid[i][0];
        volume_element[i] = (accretion_column.accretion_area*length_conv/geometry.a_0)*scale_factors[0]*scale_factors[1]*scale_factors[2];
        altitude[i] = length_conv*(r-1);
        velocity[i] = vel_conv*grid[i][2];
        mdot = 1./(scale_factors[0]*scale_factors[2]);
        total_pressure[i] = pressure_conv*mdot*(grid[i][1]-grid[i][2]);
        electron_pressure[i] = pressure_conv*grid[i][3];
        density[i] = density_conv*mdot/grid[i][2];
        electron_density[i] = mass_to_number_density*density[i];
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*(total_pressure[i]-electron_pressure[i])/(electron_density[i]/avg_atomic_charge);
    }
}

void Cataclysmic_Variable::Print_Properties() const{
    using std::cout;
    using std::endl;
    if(altitude.size() < 1){
        return;
    }
    cout << "===================================================" << endl;
    cout << "                   mCV Properties                  " << endl;
    cout << "===================================================" << endl;
    cout << " mass:               " << white_dwarf.mass/m_sol << " M_solar" << endl;
    cout << " radius:             " << white_dwarf.radius/r_sol << " R_solar" << endl;
    cout << " B_field:            " << white_dwarf.b_field/1e6 << " MG" << endl;
    if(white_dwarf.inverse_mag_radius != 0){
        cout << " R_m/R:              " << (1./white_dwarf.inverse_mag_radius)/white_dwarf.radius << endl;
    }
    cout << " accretion rate:     " << accretion_column.accretion_rate*accretion_column.accretion_area << " g/s" << endl;
    cout << " accretion rate:     " << density[0]*velocity[0] << " --> " <<  density[density.size()-1]*velocity[velocity.size()-1] << " g/cm2/s" << endl;
    cout << " shock height:       " << altitude[0]/white_dwarf.radius << " (h/R_wd)" << endl;
    cout << " shock height:       " << altitude[0] << " cm" << endl;
    cout << " shock temperature:  " << electron_temperature[0] << " keV" << endl;
    cout << " density:            " <<  electron_density[0] << " --> " << electron_density[electron_density.size()-1] <<  " e-/cm3" << endl;
}
