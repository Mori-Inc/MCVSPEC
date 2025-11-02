#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "integration.hh"
#include "mass_radius.hh"
#include "gaunt.hh"
#include <cmath>
#include <iostream>
#include <valarray>

using std::cout;
using std::endl;
using std::abs;

static double previous_shock_height = 0;
Cataclysmic_Variable::Cataclysmic_Variable(double m, double r, double b, double mdot, double inv_r_m, double corot_rat, double area, double theta, double dist, int reflection):
    mass(m), radius(r), b_field(b),  inverse_mag_radius(inv_r_m), corotation_ratio(corot_rat), distance(dist), accretion_rate(mdot), accretion_area(area), pressure_ratio(.75), incl_angle(theta), refl(reflection), geometry(sin(10*pi/180)*sin(10*pi/180))
{
    if(inverse_mag_radius>0){
        b_field = sqrt(32*accretion_rate*sqrt(grav_const*mass/pow(inverse_mag_radius,7)))/(radius*radius*radius);
    }
}

void Cataclysmic_Variable::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
    avg_ion_mass = (abundances*atomic_mass).sum()*amu_to_g;
    avg_atomic_charge = (abundances*atomic_charge).sum();
    double avg_charge_squared = (abundances*atomic_charge*atomic_charge).sum();
    double avg_charge_sqr_over_mass = (abundances*atomic_charge*atomic_charge/atomic_mass).sum()/amu_to_g;
    density_const = avg_atomic_charge/(1 + m_e*avg_atomic_charge/avg_ion_mass);
    double r, dr_dw, proj_r_w, convergance, metric[3];
    double w = sqrt(1-geometry.u);
    geometry.update_coordinates(w, r, dr_dw, proj_r_w, convergance, metric);
    area_scale = accretion_area/(metric[0]*metric[2]);
    scaled_mdot = accretion_rate/area_scale;
    free_fall_speed = sqrt(2*grav_const*mass/radius);

    coulomb_log_const = 2*m_e/(pi*alpha*c*hbar*hbar*hbar);
    exchange_const = 4*(alpha*hbar*c)*(alpha*hbar*c)*sqrt(2*pi*m_e*pow((density_const/avg_ion_mass),5))*avg_charge_sqr_over_mass;
    exchange_const *= scaled_mdot*radius/(free_fall_speed*free_fall_speed*free_fall_speed*free_fall_speed*free_fall_speed);
    bremss_const = sqrt(512*pi/(27*m_e*m_e*m_e))*alpha*alpha*alpha*hbar*hbar;
    bremss_const *= (avg_charge_squared/avg_atomic_charge)*sqrt(pow(density_const/avg_ion_mass,3));
    bremss_const *= scaled_mdot*radius/(free_fall_speed*free_fall_speed*free_fall_speed);
    cyclotron_const = 8.07e-2*(avg_atomic_charge/(avg_charge_squared*k_b*k_b))*pow(density_const/avg_ion_mass,-3.85);
    cyclotron_const *= pow(scaled_mdot,-1.85)*pow(free_fall_speed,5.85);
    cyclotron_const *= pow(area_scale,-0.425)*pow(b_field/sqrt(4-3*geometry.u), 2.85);
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

void Cataclysmic_Variable::Guess_Shock_Height(){
    if(previous_shock_height != 0){
        Update_Shock_Height(previous_shock_height);
    }
    else{
        const double integral = (39.*sqrt(3.) - 20*pi)/96.; // value of integral from EQ 7a of Wu 1994 DOI: 10.1086/174103
        const double shock_speed = sqrt(2*grav_const*mass*((1./radius) - inverse_mag_radius));
        const double c_br = bremss_const*free_fall_speed*free_fall_speed*free_fall_speed/(scaled_mdot*radius);
        double h_s = (pow(shock_speed,3.)*integral*accretion_area/(2*c_br*accretion_rate))/radius;
        Update_Shock_Height(geometry.w_0-h_s);
        h_s = (pow(v_s*free_fall_speed,3.)*integral*accretion_area/(2*c_br*accretion_rate))/radius;
        Update_Shock_Height(geometry.w_0-h_s);
    }
}

void Cataclysmic_Variable::Update_Shock_Height(double h_s){
    w_s = h_s;
    double r_s, dr_dw, proj_r_w, convergance, metric[3];
    geometry.update_coordinates(w_s, r_s, dr_dw, proj_r_w, convergance, metric);
    shock_height = (r_s-1)*radius;;

    v_s = 1./r_s;
    if(inverse_mag_radius > 0){
        const double ir_m = radius*inverse_mag_radius;
        const double r_c = 1./(corotation_ratio*ir_m);
        const double reduction_factor = 1. - 0.5/pow(1./corotation_ratio + 0.53/ir_m - 0.38, 3.7 + 0.44/ir_m);
        v_s -= ir_m + 0.5*(1./(ir_m*ir_m) - h_s*h_s)/(r_c*r_c*r_c);
        v_s *= reduction_factor;
    }
    v_s = 0.25*sqrt(v_s);

    const double mdot_s = 1./(metric[0]*metric[2]);
    p_s = 3*mdot_s*v_s;
    pe_s = p_s/(1 + 1/pressure_ratio);
    s_s = p_s*pow(v_s/mdot_s, 5./3.);
}

void Cataclysmic_Variable::Flow_Equation(double entropy,const valarray<double>& pos_vel_pres_epres, valarray<double>& derivs) const{
    const double& s = entropy;
    const double& w = pos_vel_pres_epres[0];
    const double& v = pos_vel_pres_epres[1];
    const double& p = pos_vel_pres_epres[2];
    const double& pe = pos_vel_pres_epres[3];
    double r, dr_dw, proj_r_w, convergance, metric[3];
    geometry.update_coordinates(w, r, dr_dw, proj_r_w, convergance, metric);
    const double area = metric[0]*metric[2];

    const double mdot = 1./area;
    const double dens = mdot/v;
    const double dens3 = dens*dens*dens;
    const double dens5 = dens*dens*dens*dens*dens;
    const double chi = (1+avg_atomic_charge)/avg_atomic_charge;
    const double b_scale = sqrt(4-3*geometry.u*r)/(r*r*r);

    const double ne = (scaled_mdot/free_fall_speed)*density_const*dens/avg_ion_mass; // cgs
    const double kT = scaled_mdot*free_fall_speed*pe/ne; // cgs
    const double gff = gaunt::gaunt_factor(kT);
    const double coulomb_log = 0.5*log(coulomb_log_const*kT*kT/ne);

    const double grav = mdot*proj_r_w/(2*r*r);
    const double cyc = (cyclotron_const/gff)*pe*pe*pow(b_scale/dens, 2.85)/(dens*pow(area,0.425));
    const double rad = bremss_const*gff*sqrt(pe*dens3)*(1+cyc);
    const double exch = exchange_const*coulomb_log*sqrt(dens5/pe)*(p/pe - chi);

    double dw_ds = -1.5*v*p/(metric[1]*s*rad);
    double dv_ds = (3*v*p/s)*(1. + 1.5*grav/rad + 2.5*p*v*convergance/(metric[1]*rad));
    double dp_ds = -1.5*p*grav/(s*rad) - mdot*dv_ds;
    double dpe_ds = (p/s)*(1-exch/rad) + 2.5*v*p*pe*convergance/(s*metric[1]*rad) - 5*pe*dv_ds/(3*v);

    derivs[0] = dw_ds;
    derivs[1] = dv_ds;
    derivs[2] = dp_ds;
    derivs[3] = dpe_ds;
}

double Cataclysmic_Variable::Get_Landing_Altitude(){
    // return signed distance from WD surface in w
    double s = s_s;
    valarray<double> y = {w_s, v_s, p_s, pe_s};
    accretion_column.Initialize(s, 0, y);
    while(y[1]/v_s > 1e-2){
        accretion_column.Step(s, y);
    }
    double error = 1;
    valarray<double> slope(4);
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

void Cataclysmic_Variable::Bracket_Shock_Height(){
    double w = sqrt(1-geometry.u);
    double r, dr_dw, proj_r_w, convergance, metric[3];
    geometry.update_coordinates(w, r, dr_dw, proj_r_w, convergance, metric);

    upper_bound = w + 1e-6/dr_dw; // one bound is small perturbation above surface
    Update_Shock_Height(upper_bound);
    upper_landing = Get_Landing_Altitude();
    if(upper_landing < 0){
        while(upper_landing < 0){
            lower_bound = upper_bound;
            lower_landing = upper_landing;
            upper_bound += 1e-6/dr_dw;
            Update_Shock_Height(upper_bound);
            upper_landing = Get_Landing_Altitude();
        }
        return;
    }
    double d = 4;
    lower_bound = w/d;
    Update_Shock_Height(lower_bound);
    lower_landing = Get_Landing_Altitude();
    if(lower_landing > 0){
        while(lower_landing > 0){
            upper_bound = lower_bound;
            upper_landing = lower_landing;
            d += 1;
            lower_bound = w/d;
            Update_Shock_Height(lower_bound);
            lower_landing = Get_Landing_Altitude();
        }
        return;
    }
}

void Cataclysmic_Variable::Shock_Height_Shooting(){
    Bracket_Shock_Height();
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

        Update_Shock_Height(new_bound);
        new_altitude = Get_Landing_Altitude();
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
    Update_Shock_Height((upper_bound+lower_bound)/2);
    previous_shock_height = shock_height;
}

void Cataclysmic_Variable::Build_Column_Profile(){
    // generate a velocity grid that is has a spacing of ~ kT_grid_spacing
    // interpolate between each point in the RK grid to find a set of velocities to evaluate our integral at

    // determine de-dimensionalized grid size
    const double kTe_const = erg_to_kev*free_fall_speed*free_fall_speed*avg_ion_mass/density_const;
    const double kTi_const = avg_atomic_charge*kTe_const;
    const double dkTe = kT_grid_spacing/kTe_const;
    const double dkTi = kT_grid_spacing/kTi_const;
    const valarray<double> grid_spacing = {altitude_grid_spacing, dkTe, dkTi};

    // vars for integration
    double s = entropy_boundary;
    valarray<double> y = {w_s, v_s, p_s, pe_s};
    accretion_column.Initialize(s, 0, y);
    double s_old = s;
    double ds;

    // vars for interval splitting
    valarray<double> y_mid(4), grid_vars(3), grid_vars_l(3), grid_vars_r(3), crossing(3);
    y_mid = y;
    const double& x = y_mid[0];
    const double& v = y_mid[1];
    const double& p = y_mid[2];
    const double& pe = y_mid[3];
    double r, dr_dw, proj_r_w, convergance, metric[3];
    geometry.update_coordinates(y[0], r, dr_dw, proj_r_w, convergance, metric);
    double area = metric[0]*metric[2];
    grid_vars_r = {x, pe*v*area, (p-pe)*v*area};
    grid_vars = {x, pe*v*area, (p-pe)*v*area};
    // vars for root finding
    double s_high, s_low, s_mid;
    valarray<double> root_vars(3);

    // grid
    vector<valarray<double>> grid;
    grid.push_back(y);

    int seg;
    bool root_found=false;
    while(s>1e-8 && grid_vars[1] > 0.5*dkTe){
        s_old = s;
        accretion_column.Dense_Step(s, y);
        ds = (s_old-s)/16.;
        seg=0;

        while(seg<16){
            grid_vars_l = grid_vars_r;

            s_mid = s_old - ds*(seg+1);
            accretion_column.Interpolate(s_mid, y_mid);
            geometry.update_coordinates(y_mid[0], r, dr_dw, proj_r_w, convergance, metric);
            area = metric[0]*metric[2];
            grid_vars_r = {x, pe*v*area, (p-pe)*v*area};

            crossing = (abs(grid_vars_l-grid_vars)/grid_spacing - 1)*(abs(grid_vars_r-grid_vars)/grid_spacing - 1);

            for(uint i=0; i<crossing.size(); i++){
                if(crossing[i]<0){
                    root_found = true;
                    s_high = s_old - ds*seg;
                    s_low = s_old - ds*(seg+1);
                    while(s_high-s_low > 1e-8){
                        s_mid = (s_high+s_low)/2;
                        accretion_column.Interpolate(s_mid, y_mid);
                        geometry.update_coordinates(y_mid[0], r, dr_dw, proj_r_w, convergance, metric);
                        area = metric[0]*metric[2];
                        root_vars = {x, pe*v*area, (p-pe)*v*area};
                        if(abs(root_vars[i]-grid_vars[i])/grid_spacing[i] - 1 > 0){
                            s_low = s_mid;
                        }
                        else{
                            s_high = s_mid;
                        }
                    }
                    accretion_column.Interpolate(s_high, y_mid);
                    geometry.update_coordinates(y_mid[0], r, dr_dw, proj_r_w, convergance, metric);
                    area = metric[0]*metric[2];
                    root_vars = {x, pe*v*area, (p-pe)*v*area};
                    crossing = (abs(grid_vars_l-grid_vars)/grid_spacing - 1)*(abs(root_vars-grid_vars)/grid_spacing - 1); // update crossing with the new
                }
            }// after checking crossing I will have found the earliest grid point in a given segment
            seg++;
            // update the grid
            if(root_found){
                grid_vars = root_vars;
                grid.push_back(y_mid);
                root_found=false;
                seg--; // repeate search on segment in case
                s_mid = s_high-1e-8;
                accretion_column.Interpolate(s_mid, y_mid);
                geometry.update_coordinates(y_mid[0], r, dr_dw, proj_r_w, convergance, metric);
                area = metric[0]*metric[2];
                grid_vars_r = {x, pe*v*area, (p-pe)*v*area}; // shift left bound to just after our previous root
            }
        }
    }
    int n_points = grid.size();
    velocity.resize(n_points);
    altitude.resize(n_points);
    total_pressure.resize(n_points);
    electron_pressure.resize(n_points);
    electron_density.resize(n_points);
    ion_density.resize(n_points);
    electron_temperature.resize(n_points);
    ion_temperature.resize(n_points);
    volume.resize(n_points);

    double mdot, a, b;

    for(uint i=0; i<n_points; i++){
        geometry.update_coordinates(grid[i][0], r, dr_dw, proj_r_w, convergance, metric);
        altitude[i] = radius*(r-1);
        velocity[i] = grid[i][1]*free_fall_speed;
        total_pressure[i] = scaled_mdot*free_fall_speed*grid[i][2];
        electron_pressure[i] = scaled_mdot*free_fall_speed*grid[i][3];
        mdot = scaled_mdot/(metric[0]*metric[2]);
        electron_density[i] = (mdot/velocity[i])*density_const/avg_ion_mass;
        ion_density[i] = electron_density[i]/avg_atomic_charge;
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*(total_pressure[i]-electron_pressure[i])/ion_density[i];

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
        geometry.update_coordinates(a, r, dr_dw, proj_r_w, convergance, metric);
        volume[i] = metric[0]*metric[1]*metric[2];
        geometry.update_coordinates(b, r, dr_dw, proj_r_w, convergance, metric);
        volume[i] += metric[0]*metric[1]*metric[2];
        geometry.update_coordinates((a+b)/2, r, dr_dw, proj_r_w, convergance, metric);
        volume[i] += 4*metric[0]*metric[1]*metric[2];
        volume[i] *= radius*area_scale*(b - a)/6;
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
    cout << " accretion rate:     " << shock_mdot << "-->" << accretion_rate/accretion_area << " g/cm2/s" << endl;
    cout << " shock height:       " << shock_height/radius << " (h/R_wd)" << endl;
    cout << " shock temperature:  " << electron_temperature[0] << " keV" << endl;
}
