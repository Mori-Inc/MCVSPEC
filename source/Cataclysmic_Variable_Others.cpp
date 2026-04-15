#include <Cataclysmic_Variable.hh>
#include <Cataclysmic_Variable_Others.hh>
#include "gaunt.hh"
#include <iostream>

Saxton_CV::Saxton_CV(White_Dwarf wd, Accretion_Column col, Tolerance tol):
    white_dwarf(wd), accretion_column(col), error_control(tol),
    length_conv(white_dwarf.radius), vel_conv(sqrt(2*grav_const*white_dwarf.mass/white_dwarf.radius)), accretion_rate_conv(col.accretion_rate),
    time_conv(length_conv/vel_conv), mass_conv(accretion_rate_conv*length_conv*length_conv*time_conv), volume_conv(length_conv*length_conv*length_conv),
    energy_conv(mass_conv*vel_conv*vel_conv), density_conv(mass_conv/volume_conv), pressure_conv(energy_conv/volume_conv),
    integrator(Diff_EQ{*this},error_control.absolute_error,error_control.relative_error)
{
    Set_Abundances();
}

void Saxton_CV::Set_Abundances(){
    abundances.resize(n_elements);
    abundances = {1.00e+00, 9.77e-02, 3.63e-04, 1.12e-04, 8.51e-04, 1.23e-04,
                  3.80e-05, 2.95e-06, 3.55e-05, 1.62e-05, 3.63e-06, 2.29e-06,
                  4.68e-05, 1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
    double total= abundances[0]+abundances[1];
    for(size_t i=2; i<abundances.size(); i++){
        abundances[i] *= accretion_column.metallicity;
        total += abundances[i];
    }
    std::transform(abundances.begin(),abundances.end(),abundances.begin(),[total](double x) {return x/total;});
    Set_Cooling_Constants();
}

int Saxton_CV::Solve_Profile(){
    Find_Shock_Position();
    Build_Column_Profile();
    if(!valid_solution){
        return -1;
    }
    return 0;
}

void Saxton_CV::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
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
    const double j_e = mass_to_number_density*accretion_column.accretion_rate;
    const double b_sqr = white_dwarf.b_field*white_dwarf.b_field;
    bremss_const = bremss_coeff*(avg_charge_squared/avg_atomic_charge)*sqrt(j_e*j_e*j_e);
    bremss_const *= sqrt(time_conv*time_conv*time_conv*mass_conv)*length_conv/energy_conv;
    cyclotron_const = cyclotron_coeff*(avg_atomic_charge/avg_charge_squared);
    cyclotron_const *= b_sqr*pow(b_sqr/(accretion_column.accretion_area), 0.425)/(j_e*j_e);
    cyclotron_const *= pow(j_e, -1.85);
    cyclotron_const *= pressure_conv*pressure_conv*pow(vel_conv,3.85);
    exchange_const = exchange_coeff*avg_charge_sqr_over_mass*sqrt(j_e*j_e*j_e*j_e*j_e);
    exchange_const *= length_conv*length_conv*time_conv*time_conv*time_conv*sqrt(time_conv/energy_conv)/energy_conv;
}

void Saxton_CV::Compute_Shock_Bound(double r_s, State<n_dim>& bound, double& v_s) const{
    double vff = sqrt(1./r_s - length_conv*white_dwarf.inverse_mag_radius);
    bound[0] = r_s;
    bound[1] = vff;
    bound[2] = (accretion_column.shock_pressure_ratio/(1+accretion_column.shock_pressure_ratio))*0.75*vff;
    v_s = 0.25*vff;
}

void Saxton_CV::Flow_Equation(double velocity,const State<n_dim>& state, State<n_dim>& derivs) const{
    const double& v = velocity;
    const double& r = state[0];
    const double& x = state[1];
    const double& pe = state[2];

    const double ne_cgs = mass_to_number_density*density_conv/v;
    const double kT_cgs = pressure_conv*pe/ne_cgs;
    const double gff = gaunt::gaunt_factor(kT_cgs);
    const double coulomb_log = 0.5*log(coulomb_log_coeff*kT_cgs*kT_cgs/ne_cgs);

    const double grav = 0.5/(r*r);
    const double cyc = (cyclotron_const/gff)*pe*pe*v*v*v*pow(v, 0.85);
    const double rad = bremss_const*gff*sqrt(pe/(v*v*v))*(1+cyc);
    const double p = x-v;
    const double exch = (exchange_const*coulomb_log/(bremss_const*gff))*(avg_atomic_charge*(p-pe)/pe - 1.0)/(pe*v); // ratio of exchange to radiation

    double dr_dv = (5*x - 8*v)/(2*rad + 3*grav);
    double dx_dv = (8 - 5*x/v)/(3 + 2*rad/grav);
    double dp_dv = 2*(rad-exch)*dr_dv/(3*v) - 5*pe/(3*v);

    derivs[0] = dr_dv;
    derivs[1] = dx_dv;
    derivs[2] = dp_dv;
}

double Saxton_CV::Landing_Altitude(double r_s) const {
    // return signed distance from WD surface in w
    double v_s;
    State<n_dim> y;
    Compute_Shock_Bound(r_s, y, v_s);
    double v = v_s;
    integrator.Initialize(v, 0, y);
    while(v/v_s > 1e-2){
         integrator.Step(v, y);
    }
    double error = 1;
    State<n_dim> slope{};
    Flow_Equation(v,  y, slope);
    double v_prev = v;
    double dr_prev = slope[0];

    while(error > error_control.absolute_error){
        integrator.Step(v, y);
        Flow_Equation(v,  y, slope);
        error = 0.5*std::abs((slope[0]-dr_prev)/(v-v_prev))*v*v; //difference between linear and quadratic extrapolation on w
        v_prev = v;
        dr_prev = slope[0];
    }
    integrator.Step(v, y);
    Flow_Equation(v,  y, slope);
    double landing = y[0] - slope[0]*v;
    return landing - 1.0;
}

// Bracket the shock height by minimizing the landing coordinate
// (which is maximizing the landing altitude since w ~ -r close to the surface and at the pole)
// exit if
// 1. a point is found with w_l < w_0, in which case we have bracketed our solution
// 2. a minimum is found with w_l > w_0 in which case no solution exists
int Saxton_CV::Bracket_Shock_Position(double& upper_bound, double& lower_bound, double& upper_landing, double& lower_landing) const {
    lower_bound = 1.0;
    lower_landing = Landing_Altitude(lower_bound);

    double r=1.;
    double dr = 0.01;
    double samples[3] = {Landing_Altitude(r-dr),
                        lower_landing,
                        Landing_Altitude(r+dr)};

    double drl_drs[2] = {(samples[1]-samples[0])/dr, (samples[2]-samples[1])/dr};
    double step = 0;
    while(drl_drs[0]*drl_drs[1] > 0){ // while minima not bounded
        if(samples[2]>0){
            upper_bound = r+dr;
            upper_landing = samples[2];
            return 0;
        }
        lower_bound = r+dr;
        lower_landing = samples[2];
        step = 0.5*dr*(samples[0]-samples[2])/(samples[0]-2*samples[1]+samples[2]);
        r += std::min(1.,std::abs(step));
        samples[0] = Landing_Altitude(r-dr);
        samples[1] = Landing_Altitude(r);
        samples[2] = Landing_Altitude(r+dr);
        drl_drs[0] = (samples[1]-samples[0])/dr;
        drl_drs[1] = (samples[2]-samples[1])/dr;
    }

    if(samples[1] < 0){ // if minima > 0
        return 1;
    }

    upper_bound = r;
    upper_landing = samples[1];
    return 0;
}

void Saxton_CV::Find_Shock_Position(){
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
    double v_s;
    Compute_Shock_Bound(0.5*(upper_bound+lower_bound), bound, v_s);
    shock_boundary = bound;
    shock_speed = v_s;
}

template <typename func>
void Saxton_CV::Build_Grid(func grid_func, const State<n_grid_vars>& grid_spacing, std::vector<State<n_dim>>& grid, std::vector<double>& t_grid) const {
    if(!valid_solution){
        return;
    }
    constexpr int n_segments=16;
    // integration variables
    double t = shock_speed;
    State<n_dim> y = shock_boundary;
    integrator.Initialize(t, 0, y);

    // grid variables
    State<n_grid_vars> grid_vars{};
    grid_func(t,y,grid_vars);
    grid.reserve(size_t(std::abs(grid_vars[0]/grid_spacing[0]))+1);
    t_grid.reserve(size_t(std::abs(grid_vars[0]/grid_spacing[0]))+1);
    grid.push_back(y);
    t_grid.push_back(t);

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
                t_grid.push_back(t_left);
                grid_vars = grid_left;
            }
        }
    }
}

void Saxton_CV::Build_Column_Profile(){
    if(!valid_solution){
        return;
    }
    const double dkTe = (error_control.kT_grid_spacing/erg_to_kev)*mass_to_number_density/(vel_conv*vel_conv);
    const double dkTi = dkTe/avg_atomic_charge;

    const State<n_grid_vars> grid_spacing = {dkTi, dkTe, error_control.altitude_grid_spacing*(shock_boundary[0]-1)};

    auto kT_e = [](const double t, const State<n_dim>& y){
        return t*y[2];
    };
    auto kT_i = [](const double t, const State<n_dim>& y){
        return t*(y[1]-t-y[2]);
    };
    auto grid_func = [this, kT_e, kT_i](const double& t, const State<n_dim>& y, State<n_grid_vars>& vars){
        vars[0] = kT_i(t,y);
        vars[1] = kT_e(t,y);
        vars[2] = y[0];
    };

    std::vector<State<n_dim>> grid;
    std::vector<double> t_grid;
    Build_Grid(grid_func, grid_spacing, grid, t_grid);

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

    for(size_t i=0; i<n_points; i++){
        position[i] = grid[i][0];
        volume_element[i] = accretion_column.accretion_area*length_conv;
        altitude[i] = length_conv*(grid[i][0]-1);
        velocity[i] = vel_conv*t_grid[i];
        total_pressure[i] = pressure_conv*(grid[i][1]-t_grid[i]);
        electron_pressure[i] = pressure_conv*grid[i][2];
        density[i] = density_conv/t_grid[i];
        electron_density[i] = mass_to_number_density*density[i];
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*(total_pressure[i]-electron_pressure[i])/(electron_density[i]/avg_atomic_charge);
    }
}





















Cropper_CV::Cropper_CV(White_Dwarf wd, Accretion_Column col, Tolerance tol):
    white_dwarf(wd), accretion_column(col), error_control(tol),
    length_conv(white_dwarf.radius), vel_conv(sqrt(2*grav_const*white_dwarf.mass/white_dwarf.radius)), accretion_rate_conv(col.accretion_rate),
    time_conv(length_conv/vel_conv), mass_conv(accretion_rate_conv*length_conv*length_conv*time_conv), volume_conv(length_conv*length_conv*length_conv),
    energy_conv(mass_conv*vel_conv*vel_conv), density_conv(mass_conv/volume_conv), pressure_conv(energy_conv/volume_conv),
    integrator(Diff_EQ{*this},error_control.absolute_error,error_control.relative_error)
{
    Set_Abundances();
}

void Cropper_CV::Set_Abundances(){
    abundances.resize(n_elements);
    abundances = {1.00e+00, 9.77e-02, 3.63e-04, 1.12e-04, 8.51e-04, 1.23e-04,
                  3.80e-05, 2.95e-06, 3.55e-05, 1.62e-05, 3.63e-06, 2.29e-06,
                  4.68e-05, 1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
    double total= abundances[0]+abundances[1];
    for(size_t i=2; i<abundances.size(); i++){
        abundances[i] *= accretion_column.metallicity;
        total += abundances[i];
    }
    std::transform(abundances.begin(),abundances.end(),abundances.begin(),[total](double x) {return x/total;});
    Set_Cooling_Constants();
}

int Cropper_CV::Solve_Profile(){
    Find_Shock_Position();
    Build_Column_Profile();
    if(!valid_solution){
        return -1;
    }
    return 0;
}

void Cropper_CV::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
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
    const double j_e = mass_to_number_density*accretion_column.accretion_rate;
    const double b_sqr = white_dwarf.b_field*white_dwarf.b_field;
    bremss_const = bremss_coeff*(avg_charge_squared/avg_atomic_charge)*sqrt(j_e*j_e*j_e);
    bremss_const *= sqrt(time_conv*time_conv*time_conv*mass_conv)*length_conv/energy_conv;
    cyclotron_const = cyclotron_coeff*(avg_atomic_charge/avg_charge_squared);
    cyclotron_const *= b_sqr*pow(b_sqr/(accretion_column.accretion_area), 0.425)/(j_e*j_e);
    cyclotron_const *= pow(j_e, -1.85);
    cyclotron_const *= pressure_conv*pressure_conv*pow(vel_conv,3.85);
}

void Cropper_CV::Compute_Shock_Bound(double r_s, State<n_dim>& bound, double& v_s) const{
    double vff = sqrt(1./r_s - length_conv*white_dwarf.inverse_mag_radius);
    bound[0] = r_s;
    bound[1] = vff;
    v_s = 0.25*vff;
}

void Cropper_CV::Flow_Equation(double velocity,const State<n_dim>& state, State<n_dim>& derivs) const{
    const double& v = velocity;
    const double& r = state[0];
    const double& x = state[1];

    const double p = x-v;
    const double ne_cgs = mass_to_number_density*density_conv/v;
    const double kT_cgs = pressure_conv*(avg_atomic_charge/(avg_atomic_charge+1))*p/ne_cgs;
    const double gff = gaunt::gaunt_factor(kT_cgs);

    const double grav = 0.5/(r*r);
    const double cyc = (cyclotron_const/gff)*p*p*v*v*v*pow(v, 0.85);
    const double rad = bremss_const*gff*sqrt(p/(v*v*v))*(1+cyc);

    double dr_dv = (5*x - 8*v)/(2*rad + 3*grav);
    double dx_dv = (8 - 5*x/v)/(3 + 2*rad/grav);

    derivs[0] = dr_dv;
    derivs[1] = dx_dv;
}

double Cropper_CV::Landing_Altitude(double r_s) const {
    // return signed distance from WD surface in w
    double v_s;
    State<n_dim> y;
    Compute_Shock_Bound(r_s, y, v_s);
    double v = v_s;
    integrator.Initialize(v, 0, y);
    while(v/v_s > 1e-2){
         integrator.Step(v, y);
    }
    double error = 1;
    State<n_dim> slope{};
    Flow_Equation(v,  y, slope);
    double v_prev = v;
    double dr_prev = slope[0];

    while(error > error_control.absolute_error){
        integrator.Step(v, y);
        Flow_Equation(v,  y, slope);
        error = 0.5*std::abs((slope[0]-dr_prev)/(v-v_prev))*v*v; //difference between linear and quadratic extrapolation on w
        v_prev = v;
        dr_prev = slope[0];
    }
    integrator.Step(v, y);
    Flow_Equation(v,  y, slope);
    double landing = y[0] - slope[0]*v;
    return landing - 1.0;
}

// Bracket the shock height by minimizing the landing coordinate
// (which is maximizing the landing altitude since w ~ -r close to the surface and at the pole)
// exit if
// 1. a point is found with w_l < w_0, in which case we have bracketed our solution
// 2. a minimum is found with w_l > w_0 in which case no solution exists
int Cropper_CV::Bracket_Shock_Position(double& upper_bound, double& lower_bound, double& upper_landing, double& lower_landing) const {
    lower_bound = 1.0;
    lower_landing = Landing_Altitude(lower_bound);

    double r=1.;
    double dr = 0.01;
    double samples[3] = {Landing_Altitude(r-dr),
                        lower_landing,
                        Landing_Altitude(r+dr)};

    double drl_drs[2] = {(samples[1]-samples[0])/dr, (samples[2]-samples[1])/dr};
    double step = 0;
    while(drl_drs[0]*drl_drs[1] > 0){ // while minima not bounded
        if(samples[2]>0){
            upper_bound = r+dr;
            upper_landing = samples[2];
            return 0;
        }
        lower_bound = r+dr;
        lower_landing = samples[2];
        step = 0.5*dr*(samples[0]-samples[2])/(samples[0]-2*samples[1]+samples[2]);
        r += std::min(1.,std::abs(step));
        samples[0] = Landing_Altitude(r-dr);
        samples[1] = Landing_Altitude(r);
        samples[2] = Landing_Altitude(r+dr);
        drl_drs[0] = (samples[1]-samples[0])/dr;
        drl_drs[1] = (samples[2]-samples[1])/dr;
    }

    if(samples[1] < 0){ // if minima > 0
        return 1;
    }

    upper_bound = r;
    upper_landing = samples[1];
    return 0;
}

void Cropper_CV::Find_Shock_Position(){
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
    double v_s;
    Compute_Shock_Bound(0.5*(upper_bound+lower_bound), bound, v_s);
    shock_boundary = bound;
    shock_speed = v_s;
}

template <typename func>
void Cropper_CV::Build_Grid(func grid_func, const State<n_grid_vars>& grid_spacing, std::vector<State<n_dim>>& grid, std::vector<double>& t_grid) const {
    if(!valid_solution){
        return;
    }
    constexpr int n_segments=16;
    // integration variables
    double t = shock_speed;
    State<n_dim> y = shock_boundary;
    integrator.Initialize(t, 0, y);

    // grid variables
    State<n_grid_vars> grid_vars{};
    grid_func(t,y,grid_vars);
    grid.reserve(size_t(std::abs(grid_vars[0]/grid_spacing[0]))+1);
    t_grid.reserve(size_t(std::abs(grid_vars[0]/grid_spacing[0]))+1);
    grid.push_back(y);
    t_grid.push_back(t);

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
                t_grid.push_back(t_left);
                grid_vars = grid_left;
            }
        }
    }
}

void Cropper_CV::Build_Column_Profile(){
    if(!valid_solution){
        return;
    }
    const double dkTe = (avg_atomic_charge/(avg_atomic_charge+1))*(error_control.kT_grid_spacing/erg_to_kev)*mass_to_number_density/(vel_conv*vel_conv);

    const State<n_grid_vars> grid_spacing = {dkTe, error_control.altitude_grid_spacing*(shock_boundary[0]-1)};

    auto kT_e = [](const double t, const State<n_dim>& y){
        return t*(y[1]-t);
    };
    auto grid_func = [this, kT_e](const double& t, const State<n_dim>& y, State<n_grid_vars>& vars){
        vars[0] = kT_e(t,y);
        vars[1] = y[0];
    };

    std::vector<State<n_dim>> grid;
    std::vector<double> t_grid;
    Build_Grid(grid_func, grid_spacing, grid, t_grid);

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

    for(size_t i=0; i<n_points; i++){
        position[i] = grid[i][0];
        volume_element[i] = accretion_column.accretion_area*length_conv;
        altitude[i] = length_conv*(grid[i][0]-1);
        velocity[i] = vel_conv*t_grid[i];
        total_pressure[i] = pressure_conv*(grid[i][1]-t_grid[i]);
        electron_pressure[i] = total_pressure[i]/(1.0 + 1./avg_atomic_charge);
        density[i] = density_conv/t_grid[i];
        electron_density[i] = mass_to_number_density*density[i];
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*(total_pressure[i]-electron_pressure[i])/(electron_density[i]/avg_atomic_charge);
    }
}



















Wu_CV::Wu_CV(White_Dwarf wd, Accretion_Column col, Tolerance tol):
    white_dwarf(wd), accretion_column(col), error_control(tol),
    length_conv(white_dwarf.radius), vel_conv(sqrt(2*grav_const*white_dwarf.mass/white_dwarf.radius)), accretion_rate_conv(col.accretion_rate),
    time_conv(length_conv/vel_conv), mass_conv(accretion_rate_conv*length_conv*length_conv*time_conv), volume_conv(length_conv*length_conv*length_conv),
    energy_conv(mass_conv*vel_conv*vel_conv), density_conv(mass_conv/volume_conv), pressure_conv(energy_conv/volume_conv),
    integrator(Diff_EQ{*this},error_control.absolute_error,error_control.relative_error)
{
    Set_Abundances();
}

void Wu_CV::Set_Abundances(){
    abundances.resize(n_elements);
    abundances = {1.00e+00, 9.77e-02, 3.63e-04, 1.12e-04, 8.51e-04, 1.23e-04,
                  3.80e-05, 2.95e-06, 3.55e-05, 1.62e-05, 3.63e-06, 2.29e-06,
                  4.68e-05, 1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
    double total= abundances[0]+abundances[1];
    for(size_t i=2; i<abundances.size(); i++){
        abundances[i] *= accretion_column.metallicity;
        total += abundances[i];
    }
    std::transform(abundances.begin(),abundances.end(),abundances.begin(),[total](double x) {return x/total;});
    Set_Cooling_Constants();
}

int Wu_CV::Solve_Profile(){
    Find_Shock_Position();
    Build_Column_Profile();
    if(!valid_solution){
        return -1;
    }
    return 0;
}

void Wu_CV::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
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
    const double j_e = mass_to_number_density*accretion_column.accretion_rate;
    const double b_sqr = white_dwarf.b_field*white_dwarf.b_field;
    bremss_const = bremss_coeff*(avg_charge_squared/avg_atomic_charge)*sqrt(j_e*j_e*j_e);
    bremss_const *= sqrt(time_conv*time_conv*time_conv*mass_conv)*length_conv/energy_conv;
    cyclotron_const = cyclotron_coeff*(avg_atomic_charge/avg_charge_squared);
    cyclotron_const *= b_sqr*pow(b_sqr/(accretion_column.accretion_area), 0.425)/(j_e*j_e);
    cyclotron_const *= pow(j_e, -1.85);
    cyclotron_const *= pressure_conv*pressure_conv*pow(vel_conv,3.85);
}

void Wu_CV::Compute_Shock_Bound(double r_s, State<n_dim>& bound, double& v_s){
    vff = sqrt(1./r_s - length_conv*white_dwarf.inverse_mag_radius);
    bound[0] = r_s;
    v_s = 0.25*vff;
}

void Wu_CV::Flow_Equation(double velocity,const State<n_dim>& state, State<n_dim>& derivs) const{
    const double& v = velocity;
    const double& r = state[0];

    const double p = vff-v;
    const double ne_cgs = mass_to_number_density*density_conv/v;
    const double kT_cgs = pressure_conv*(avg_atomic_charge/(avg_atomic_charge+1))*p/ne_cgs;
    const double gff = gaunt::gaunt_factor(kT_cgs);

    const double cyc = (cyclotron_const/gff)*p*p*v*v*v*pow(v, 0.85);
    const double rad = bremss_const*gff*sqrt(p/(v*v*v))*(1+cyc);

    double dr_dv = (5*vff-8*v)/(2*rad);

    derivs[0] = dr_dv;
}

double Wu_CV::Landing_Altitude(double r_s){
    // return signed distance from WD surface in w
    double v_s;
    State<n_dim> y;
    Compute_Shock_Bound(r_s, y, v_s);
    double v = v_s;
    integrator.Initialize(v, 0, y);
    while(v/v_s > 1e-2){
         integrator.Step(v, y);
    }
    double error = 1;
    State<n_dim> slope{};
    Flow_Equation(v,  y, slope);
    double v_prev = v;
    double dr_prev = slope[0];

    while(error > error_control.absolute_error){
        integrator.Step(v, y);
        Flow_Equation(v,  y, slope);
        error = 0.5*std::abs((slope[0]-dr_prev)/(v-v_prev))*v*v; //difference between linear and quadratic extrapolation on w
        v_prev = v;
        dr_prev = slope[0];
    }
    integrator.Step(v, y);
    Flow_Equation(v,  y, slope);
    double landing = y[0] - slope[0]*v;
    return landing - 1.0;
}

// Bracket the shock height by minimizing the landing coordinate
// (which is maximizing the landing altitude since w ~ -r close to the surface and at the pole)
// exit if
// 1. a point is found with w_l < w_0, in which case we have bracketed our solution
// 2. a minimum is found with w_l > w_0 in which case no solution exists
int Wu_CV::Bracket_Shock_Position(double& upper_bound, double& lower_bound, double& upper_landing, double& lower_landing){
    lower_bound = 1.0;
    lower_landing = Landing_Altitude(lower_bound);

    double r=1.;
    double dr = 0.01;
    double samples[3] = {Landing_Altitude(r-dr),
                        lower_landing,
                        Landing_Altitude(r+dr)};

    double drl_drs[2] = {(samples[1]-samples[0])/dr, (samples[2]-samples[1])/dr};
    double step = 0;
    while(drl_drs[0]*drl_drs[1] > 0){ // while minima not bounded
        if(samples[2]>0){
            upper_bound = r+dr;
            upper_landing = samples[2];
            return 0;
        }
        lower_bound = r+dr;
        lower_landing = samples[2];
        step = 0.5*dr*(samples[0]-samples[2])/(samples[0]-2*samples[1]+samples[2]);
        r += std::min(1.,std::abs(step));
        samples[0] = Landing_Altitude(r-dr);
        samples[1] = Landing_Altitude(r);
        samples[2] = Landing_Altitude(r+dr);
        drl_drs[0] = (samples[1]-samples[0])/dr;
        drl_drs[1] = (samples[2]-samples[1])/dr;
    }

    if(samples[1] < 0){ // if minima > 0
        return 1;
    }

    upper_bound = r;
    upper_landing = samples[1];
    return 0;
}

void Wu_CV::Find_Shock_Position(){
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
    double v_s;
    Compute_Shock_Bound(0.5*(upper_bound+lower_bound), bound, v_s);
    shock_boundary = bound;
    shock_speed = v_s;
    vff = 4*shock_speed;
}

template <typename func>
void Wu_CV::Build_Grid(func grid_func, const State<n_grid_vars>& grid_spacing, std::vector<State<n_dim>>& grid, std::vector<double>& t_grid) const {
    if(!valid_solution){
        return;
    }
    constexpr int n_segments=16;
    // integration variables
    double t = shock_speed;
    State<n_dim> y = shock_boundary;
    integrator.Initialize(t, 0, y);

    // grid variables
    State<n_grid_vars> grid_vars{};
    grid_func(t,y,grid_vars);
    grid.reserve(size_t(std::abs(grid_vars[0]/grid_spacing[0]))+1);
    t_grid.reserve(size_t(std::abs(grid_vars[0]/grid_spacing[0]))+1);
    grid.push_back(y);
    t_grid.push_back(t);

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
                t_grid.push_back(t_left);
                grid_vars = grid_left;
            }
        }
    }
}

void Wu_CV::Build_Column_Profile(){
    if(!valid_solution){
        return;
    }
    const double dkTe = (avg_atomic_charge/(avg_atomic_charge+1))*(error_control.kT_grid_spacing/erg_to_kev)*mass_to_number_density/(vel_conv*vel_conv);

    const State<n_grid_vars> grid_spacing = {dkTe, error_control.altitude_grid_spacing*(shock_boundary[0]-1)};

    auto kT_e = [&](const double t, const State<n_dim>& y){
        return t*(vff-t);
    };
    auto grid_func = [this, kT_e](const double& t, const State<n_dim>& y, State<n_grid_vars>& vars){
        vars[0] = kT_e(t,y);
        vars[1] = y[0];
    };

    std::vector<State<n_dim>> grid;
    std::vector<double> t_grid;
    Build_Grid(grid_func, grid_spacing, grid, t_grid);

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

    for(size_t i=0; i<n_points; i++){
        position[i] = grid[i][0];
        volume_element[i] = accretion_column.accretion_area*length_conv;
        altitude[i] = length_conv*(grid[i][0]-1);
        velocity[i] = vel_conv*t_grid[i];
        total_pressure[i] = pressure_conv*(vff-t_grid[i]);
        electron_pressure[i] = total_pressure[i]/(1.0 + 1./avg_atomic_charge);
        density[i] = density_conv/t_grid[i];
        electron_density[i] = mass_to_number_density*density[i];
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*(total_pressure[i]-electron_pressure[i])/(electron_density[i]/avg_atomic_charge);
    }
}
