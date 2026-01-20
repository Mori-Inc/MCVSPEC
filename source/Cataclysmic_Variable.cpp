#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "integration.hh"
#include "mass_radius.hh"
#include "gaunt.hh"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::cerr;
using std::abs;
using std::vector;

static double previous_shock_height = 0;
static int f_evals = 0;
Cataclysmic_Variable::Cataclysmic_Variable(double m, double r, double b, double mdot, double area, double inv_r_m, double corot_rat, double abund, double theta, double dist, int reflection):
    mass(m), radius(r), b_field(b),  inverse_mag_radius(inv_r_m), corotation_ratio(corot_rat), distance(dist), accretion_rate(mdot), accretion_area(area), metalicity(abund),
    pressure_ratio(.75), incl_angle(theta), refl(reflection), geometry(sin(1*pi/180)*sin(1*pi/180)),
    length_conv(radius), vel_conv(sqrt(2*grav_const*mass/radius)), time_conv(length_conv/vel_conv), volume_conv(length_conv*length_conv*length_conv),
    mass_conv((geometry.a_0*accretion_rate/accretion_area)*volume_conv/vel_conv),
    energy_conv(mass_conv*vel_conv*vel_conv), density_conv(mass_conv/volume_conv)
{}

void Cataclysmic_Variable::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
    // constants related to column composition
    avg_ion_mass = 0;
    avg_atomic_charge = 0;
    double avg_charge_squared = 0;
    double avg_charge_sqr_over_mass = 0;

    for(uint i=0; i<abundances.size(); i++){
        avg_ion_mass += abundances[i]*atomic_mass[i]*amu_to_g;
        avg_atomic_charge += abundances[i]*atomic_charge[i];
        avg_charge_squared += abundances[i]*atomic_charge[i]*atomic_charge[i];
        avg_charge_sqr_over_mass += abundances[i]*atomic_charge[i]*atomic_charge[i]/(atomic_mass[i]*amu_to_g);
    }

    density_const = avg_atomic_charge/(1 + m_e*avg_atomic_charge/avg_ion_mass);
    bremss_const = bremss_coeff*(avg_charge_squared/avg_atomic_charge)*pow(density_const/avg_ion_mass,1.5);
    bremss_const /= energy_conv*length_conv*length_conv/(mass_conv*mass_conv);
    cyclotron_const = cyclotron_coeff*(avg_atomic_charge/avg_charge_squared)*pow(density_const/avg_ion_mass,-3.85);
    cyclotron_const *= pow(b_field/sqrt(4-3*geometry.u), 2.85)*pow(accretion_area/geometry.a_0,-0.425);
    cyclotron_const *= vel_conv*vel_conv*vel_conv*vel_conv/pow(density_conv, 1.85);
    exchange_const = exchange_coeff*avg_charge_sqr_over_mass*pow(density_const/avg_ion_mass, 2.5);
    exchange_const /= energy_conv*energy_conv*length_conv*length_conv/(mass_conv*mass_conv*mass_conv);
}

double Cataclysmic_Variable::Get_Accretion_Rate(double luminosity, double mass, double radius, double inverse_mag_radius){
    double accretion_rate = luminosity/(grav_const*mass*((1./radius) - inverse_mag_radius));
    return accretion_rate;
}

double Cataclysmic_Variable::Get_Radius(double mass){
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

void Cataclysmic_Variable::Update_Shock_Position(double shock_pos){
    w_s = shock_pos;
    double r_s, proj_r_w, convergance, metric[3];
    geometry.update_coordinates(w_s, r_s, proj_r_w, convergance, metric);
    shock_height = (r_s-1)*radius;
    double mdot = 1./(metric[0]*metric[2]);

    double vff = sqrt(1./r_s - radius*inverse_mag_radius);

    x_s = vff;
    v_s = vff/4;
    pe_s = (pressure_ratio/(pressure_ratio+1))*mdot*(x_s-v_s);
    s_s = mdot*(x_s-v_s)*pow(v_s/mdot, 5./3.);
}

void Cataclysmic_Variable::Flow_Equation(double entropy,const vector<double>& state, vector<double>& derivs) const{
    const double& s = entropy;
    const double& w = state[0];
    const double& x = state[1];
    const double& v = state[2];
    const double& pe = state[3];

    double r, proj_r_w, convergance, metric[3];
    geometry.update_coordinates(w, r, proj_r_w, convergance, metric);

    const double area = metric[0]*metric[2];
    const double mdot = 1./area;
    const double p = mdot*(x-v);
    const double dens = mdot/v;
    const double dens3 = dens*dens*dens;
    const double dens5 = dens3*dens*dens;
    const double chi = (1+avg_atomic_charge)/avg_atomic_charge;
    const double b_scale = sqrt(4-3*geometry.u*r)/(r*r*r);

    const double ne_cgs = (density_conv*dens)*density_const/avg_ion_mass;
    const double kT_cgs = (energy_conv/volume_conv)*pe/ne_cgs;
    const double gff = gaunt::gaunt_factor(kT_cgs);
    const double coulomb_log = 0.5*log(coulomb_log_coeff*kT_cgs*kT_cgs/ne_cgs);

    const double grav = -0.5*proj_r_w/(r*r);
    const double cyc = (cyclotron_const/gff)*pe*pe*pow(b_scale/dens, 2.85)/(dens*pow(area,0.425));
    const double rad = bremss_const*gff*sqrt(pe*dens3)*(1+cyc);
    const double exch = exchange_const*coulomb_log*sqrt(dens5/pe)*(p/pe - chi);
    const double geom = convergance*v*(x-v)/metric[1];

    const double common_factor = p/s;

    double dw_ds = -1.5*common_factor*v/(metric[1]*rad);
    double dx_ds = -1.5*common_factor*(grav + geom)/rad;
    double dv_ds = (1.5*common_factor*v/(5*x - 8*v))*(2./mdot + 3*grav/rad + 5*geom/rad);
    double dp_ds = common_factor*(1 - exch/rad - 2.5*(pe/(5*x-8*v))*(2./mdot + 3*grav/rad + 3*(v/(x-v))*geom/rad));

    derivs[0] = dw_ds;
    derivs[1] = dx_ds;
    derivs[2] = dv_ds;
    derivs[3] = dp_ds;
    f_evals++;
}

double Cataclysmic_Variable::Get_Landing_Altitude(double w_s){
    // return signed distance from WD surface in w
    Update_Shock_Position(w_s);
    double s = s_s;
    vector<double> y = {w_s, x_s, v_s, pe_s};
    accretion_column.Initialize(s, 0, y);
    while(s/s_s > 1e-2){
         accretion_column.Step(s, y);
    }
    double error = 1;
    vector<double> slope(4);
    Flow_Equation(s,  y, slope);
    double s_prev = s;
    double dw_prev = slope[0];

    while(error > 1e-8){
        accretion_column.Step(s, y);
        Flow_Equation(s,  y, slope);
        error = 0.5*abs((slope[0]-dw_prev)/(s-s_prev))*s*s; //difference between linear and quadratic extroplation on w
        s_prev = s;
        dw_prev = slope[0];
    }
    accretion_column.Step(s, y);
    Flow_Equation(s,  y, slope);
    double landing = y[0] - slope[0]*s;
    return landing - geometry.w_0;
}

// Bracket the shock height by minimizing the landing coordinate
// (which is maximizing the landing altitude since w ~ -r close to the surface and at the pole)
// exit if
// 1. a point is found with w_l < w_0, in which case we have bracketed our solution
// 2. a minimum is found with w_l > w_0 in which case no solution exists
void Cataclysmic_Variable::Bracket_Shock_Position(double& upper_bound, double& lower_bound, double& upper_landing, double& lower_landing){
    upper_bound = geometry.w_0;
    upper_landing = Get_Landing_Altitude(upper_bound);
    if(upper_landing < 0){
        cerr << "Error: column is inverted?" << endl;
    }

    double logw = -log(geometry.w_0);
    double dlw = 0.01;
    double samples[3] = {Get_Landing_Altitude(1./exp(logw-dlw)),
                        upper_landing,
                        Get_Landing_Altitude(1./exp(logw+dlw))};
    double dwl_dlw[2] = {(samples[1]-samples[0])/dlw, (samples[2]-samples[1])/dlw};
    double step = 0;
    while(dwl_dlw[0]*dwl_dlw[1] > 0){ // while minima not bounded
        if(samples[2]<0){
            lower_bound = 1./exp(logw+dlw);
            lower_landing = samples[2];
            return;
        }
        upper_bound = 1./exp(logw+dlw);
        upper_landing = samples[2];
        step = 0.5*dlw*(samples[0]-samples[2])/(samples[0]-2*samples[1]+samples[2]);
        logw += step;
        samples[0] = Get_Landing_Altitude(1./exp(logw-dlw));
        samples[1] = Get_Landing_Altitude(1./exp(logw));
        samples[2] = Get_Landing_Altitude(1./exp(logw+dlw));
        dwl_dlw[0] = (samples[1]-samples[0])/dlw;
        dwl_dlw[1] = (samples[2]-samples[1])/dlw;
    }
    if(samples[1] > 0){ // if minima > 0
        cerr << "No valid solutions" << endl;
    }
    lower_bound = 1./exp(logw);
    lower_landing = samples[1];
}

void Cataclysmic_Variable::Determine_Shock_Position(){
    double upper_bound, lower_bound;
    double upper_landing, lower_landing;
    Bracket_Shock_Position(upper_bound, lower_bound, upper_landing, lower_landing);
    double k1 = 0.2/(upper_bound-lower_bound);
    double n0 = 1;
    double nmax = log2((upper_bound-lower_bound)/(2*h_s_tolerance)) + n0;
    int i=0;
    double new_bound, new_altitude, midpoint, regula_falsi, truncation, projection, dir;
    while(upper_bound-lower_bound > h_s_tolerance){
        midpoint = (upper_bound+lower_bound)/2;
        regula_falsi = (upper_landing*lower_bound - lower_landing*upper_bound)/(lower_landing-upper_landing);
        dir = (0. < (midpoint-regula_falsi)) - ((midpoint-regula_falsi) < 0.);
        truncation = k1*(upper_bound-lower_bound)*(upper_bound-lower_bound); // k2 = 2

        if(truncation <= abs(midpoint-regula_falsi)){
            new_bound = regula_falsi + dir*truncation;
        }
        else{
            new_bound = midpoint;
        }

        projection = h_s_tolerance*(pow(2,nmax-i)) - (upper_bound-lower_bound)/2;
        if(abs(new_bound-midpoint) > projection){
            new_bound = midpoint - dir*projection;
        }

        new_altitude = Get_Landing_Altitude(new_bound);
        if(new_altitude>0){
            upper_bound = new_bound;
            upper_landing = new_altitude;
        }
        else if(new_altitude<0){
            lower_bound = new_bound;
            lower_landing = new_altitude;
        }
        else{
            upper_bound = new_altitude;
            lower_bound = new_altitude;
        }
        i++;
    }
    Update_Shock_Position((upper_bound+lower_bound)/2);
    previous_shock_height = w_s;
}

template <typename func>
void Cataclysmic_Variable::Build_Grid(func grid_func, const vector<double>& grid_spacing, vector<vector<double>>& grid){
    const int n_segments = 16; // number of subdivisions to make between each RK step to search for roots
    double t = s_s;
    vector<double> y = {w_s, x_s, v_s, pe_s};
    double dt;
    vector<double> last_grid_vars(grid_spacing.size());
    vector<vector<double>> current_grid_vars(2,vector<double>(grid_spacing.size()));
    vector<double> distance(2);
    double t_interp, t_grid;
    vector<double> y_interp(y.size());


    accretion_column.Initialize(t, 0, y);
    grid.push_back(y);
    grid_func(t,y,last_grid_vars);
    current_grid_vars[1] = last_grid_vars;

    while(last_grid_vars[0] > grid_spacing[0]){
        t_interp = t;
        y_interp = y;
        accretion_column.Dense_Step(t, y);

        dt = (t-t_interp)/n_segments;
        for(int i=0; i<n_segments; i++){
            current_grid_vars[0] = current_grid_vars[1];
            distance[0] = distance[1];
            t_interp += dt;
            accretion_column.Interpolate(t_interp, y_interp);
            grid_func(t_interp,y_interp,current_grid_vars[1]);
            for(int j=0; j<grid_spacing.size(); j++){
                distance[1] = abs((current_grid_vars[1][j]-last_grid_vars[j])/grid_spacing[j]);
                if(distance[1]>1){
                    t_grid = t_interp + (1-distance[1])*dt/(distance[1]-distance[0]);
                    accretion_column.Interpolate(t_grid, y_interp);
                    grid_func(t_grid,y_interp,last_grid_vars);
                    grid.push_back(y_interp);
                }
            }
        }
    }
}

void Cataclysmic_Variable::Build_Column_Profile(){

    const double dkTe = (kT_grid_spacing/erg_to_kev)*density_const/(avg_ion_mass*vel_conv*vel_conv);
    const double dkTi = dkTe/avg_atomic_charge;
    const vector<double> grid_spacing = {dkTe, dkTi, altitude_grid_spacing*shock_height};

    double r, proj, conv, metric[3];
    vector<double> dy_dt(4);
    auto kT_e = [](const vector<double>& y, double metric[3]){
        return y[2]*y[3]*metric[0]*metric[2];
    };
    auto kT_i = [](const vector<double>& y, double metric[3]){
        return y[2]*(y[1]-y[2]-y[3]*metric[0]*metric[2]);
    };
    auto grid_func = [this, kT_e, kT_i, &r, &proj, &conv, &metric, &dy_dt](const double& t, const vector<double>& y, vector<double>& vars){
        geometry.update_coordinates(y[0], r, proj, conv, metric);
        Flow_Equation(t,y,dy_dt);

        vars[0] = kT_e(y, metric);
        vars[1] = kT_i(y, metric);
        vars[2] = r;
    };

    vector<vector<double>> grid;
    Build_Grid(grid_func, grid_spacing, grid);

    int n_points = grid.size();
    velocity.resize(n_points);
    altitude.resize(n_points);
    total_pressure.resize(n_points);
    electron_pressure.resize(n_points);
    electron_density.resize(n_points);
    density.resize(n_points);
    electron_temperature.resize(n_points);
    ion_temperature.resize(n_points);
    volume.resize(n_points);

    double mdot, a, b;

    for(uint i=0; i<n_points; i++){
        geometry.update_coordinates(grid[i][0], r, proj, conv, metric);
        altitude[i] = length_conv*(r-1);
        velocity[i] = vel_conv*grid[i][2];
        mdot = 1./(metric[0]*metric[2]);
        total_pressure[i] = (energy_conv/volume_conv)*mdot*(grid[i][1]-grid[i][2]);
        electron_pressure[i] = (energy_conv/volume_conv)*grid[i][3];
        density[i] = density_conv*mdot/grid[i][2];
        electron_density[i] = (density_const/avg_ion_mass)*density[i];
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*(total_pressure[i]-electron_pressure[i])/(electron_density[i]/avg_atomic_charge);
        if(i==0){
            a = grid[i][0];
        }
        else{
            a = (grid[i-1][0] + grid[i][0])/2;
        }
        if(i==n_points-1){
            b = grid[i][0];
        }
        else{
            b = (grid[i][0] + grid[i+1][0])/2;
        }
        geometry.update_coordinates(a, r, proj, conv, metric);
        volume[i] = metric[0]*metric[1]*metric[2];
        geometry.update_coordinates(b, r, proj, conv, metric);
        volume[i] += metric[0]*metric[1]*metric[2];
        geometry.update_coordinates((a+b)/2, r, proj, conv, metric);
        volume[i] += 4*metric[0]*metric[1]*metric[2];
        volume[i] *= length_conv*length_conv*length_conv*(b - a)/6;
    }
}

void Cataclysmic_Variable::Print_Properties(){
    cout << "===================================================" << endl;
    cout << "                   mCV Properties                  " << endl;
    cout << "===================================================" << endl;
    cout << " mass:               " << mass/m_sol << " M_solar" << endl;
    cout << " radius:             " << radius/r_sol << " R_solar" << endl;
    cout << " B_field:            " << b_field/1e6 << " MG" << endl;
    if(inverse_mag_radius != 0){
        cout << " R_m/R:              " << (1./inverse_mag_radius)/radius << endl;
    }
    cout << " accretion rate:     " << accretion_rate << " g/s" << endl;
    cout << " accretion rate:     " << density[0]*velocity[0] << " --> " <<  density[density.size()-1]*velocity[velocity.size()-1] << " g/cm2/s" << endl;
    cout << " shock height:       " << shock_height/radius << " (h/R_wd)" << endl;
    cout << " shock height:       " << shock_height << " cm" << endl;
    cout << " shock temperature:  " << electron_temperature[0] << " keV" << endl;
    cout << " density:            " <<  electron_density[0] << " --> " << electron_density[electron_density.size()-1] <<  " e-/cm3" << endl;
    cout << " fevals: " << f_evals << endl;
}
