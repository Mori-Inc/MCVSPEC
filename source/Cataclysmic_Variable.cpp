#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "integration.hh"
#include "mass_radius.hh"
#include <cmath>
#include <iostream>
#include <XSFunctions/Utilities/FunctionUtility.h>

using std::cout;
using std::endl;
using std::abs;
using std::ceil;

static double previous_shock_height = 0;

valarray<double> Flow_Equation_Wrapper(double t, valarray<double> y, void* cv_instance){
    return ((Cataclysmic_Variable*)(cv_instance))->Flow_Equation(t, y);
}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double r, double b, double mdot, double inv_r_m, double area, double theta, double n, double dist, int reflection):
    mass(m), radius(r), b_field(b),  inverse_mag_radius(inv_r_m), distance(dist), accretion_rate(mdot), accretion_area(area), pressure_ratio(.75), incl_angle(theta), area_exponent(n),  refl(reflection),
    accretion_column(Flow_Equation_Wrapper, this, 3)
{
    if(inverse_mag_radius>0){
        b_field = sqrt(32*accretion_rate*sqrt(grav_const*mass/pow(inverse_mag_radius,7)))/(radius*radius*radius);
    }
}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double r, double b, double mdot, double inv_r_m, double metals, double area, double theta, double n, double dist, int reflection):
    mass(m), radius(r), b_field(b),  inverse_mag_radius(inv_r_m), distance(dist), accretion_rate(mdot), accretion_area(area), metalicity(metals), pressure_ratio(.75), incl_angle(theta), area_exponent(n),  refl(reflection),
    accretion_column(Flow_Equation_Wrapper, this, 3)
{
    if(inverse_mag_radius>0){
        b_field = sqrt(32*accretion_rate*sqrt(grav_const*mass/pow(inverse_mag_radius,7)))/(radius*radius*radius);
    }
    Set_Abundances(metalicity);
    Guess_Shock_Height();
}

void Cataclysmic_Variable::Set_Abundances(double metalicity){
    abundances.resize(atomic_charge.size());
    abundances[0] = FunctionUtility::getAbundance(atomic_charge[0]);

    abundances[1] = FunctionUtility::getAbundance(atomic_charge[1]);
    for(int i = 2; i < 14; i++){
        abundances[i] = metalicity*FunctionUtility::getAbundance(atomic_charge[i]);
    }
    abundances = abundances/abundances.sum();
    Set_Cooling_Constants();
}

void Cataclysmic_Variable::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
    avg_ion_mass = (abundances*atomic_mass).sum()*amu_to_g;
    avg_atomic_charge = (abundances*atomic_charge).sum();
    double avg_charge_squared = (abundances*atomic_charge*atomic_charge).sum();
    double avg_charge_sqr_over_mass = (abundances*atomic_charge*atomic_charge/atomic_mass).sum()/amu_to_g;
    density_const = avg_atomic_charge/(1 + m_e*avg_atomic_charge/avg_ion_mass);

    force_const = grav_const*mass/(radius*radius);
    cooling_ratio_const = 8.07e-2*avg_atomic_charge*pow(b_field, 2.85)*pow(avg_ion_mass/density_const,3.85);
    cooling_ratio_const /= gaunt_factor*avg_charge_squared*k_b*k_b*pow(accretion_area,0.425);
    coulomb_log_const = 0.5*log(2*m_e/(pi*alpha*c)) + 1.5*log(avg_ion_mass/(hbar*density_const));
    exchange_const = 4*(alpha*hbar*c)*(alpha*hbar*c)*sqrt(2*pi*m_e*pow((density_const/avg_ion_mass),5))*avg_charge_sqr_over_mass;
    bremss_const = sqrt(512*pi/(27*m_e*m_e*m_e))*alpha*alpha*alpha*hbar*hbar*gaunt_factor;
    bremss_const *= (avg_charge_squared/avg_atomic_charge)*sqrt(pow(density_const/avg_ion_mass,3));
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
        double integral = (39.*sqrt(3.) - 20*pi)/96.; // value of integral from EQ 7a of Wu 1994 DOI: 10.1086/174103
        shock_speed = sqrt(2*grav_const*mass*((1./radius) - inverse_mag_radius));
        Update_Shock_Height(pow(shock_speed,3.)*integral*accretion_area/(2*bremss_const*accretion_rate));
    }
}

void Cataclysmic_Variable::Update_Shock_Height(double h_s){
    shock_height = h_s;
    shock_speed = sqrt(2*grav_const*mass*((1./(radius+shock_height)) - inverse_mag_radius));
    shock_mdot = accretion_rate/(accretion_area*pow(1+shock_height/radius, area_exponent));
    non_dim_radius = radius/shock_height;
    cooling_ratio = cooling_ratio_const*pow(shock_speed,5.85)/pow(shock_mdot, 1.85);
}

valarray<double> Cataclysmic_Variable::Flow_Equation(double vel, valarray<double> pos_pres_epres){
    double pos = pos_pres_epres[0];
    double press = pos_pres_epres[1];
    double e_press = pos_pres_epres[2];

    double mdot = pow((1+non_dim_radius)/(pos+non_dim_radius), area_exponent);
    double dens = mdot/vel;
    double coulomb_log = coulomb_log_const + 2.5*log(shock_speed) - 0.5*log(shock_mdot) + 0.5*log(e_press*e_press/(dens*dens*dens));

    double gravity = force_const*dens/((1+pos/non_dim_radius)*(1+pos/non_dim_radius));
    gravity *= (shock_height/(shock_speed*shock_speed));
    double exchange = exchange_const*coulomb_log*sqrt(dens*dens*dens*dens*dens/e_press)*(press/e_press - ((1+avg_atomic_charge)/avg_atomic_charge));
    exchange *= (shock_mdot*shock_height/(shock_speed*shock_speed*shock_speed));
    double radiation = bremss_const*sqrt(dens*dens*dens*e_press);
    radiation *= 1 + cooling_ratio*e_press*e_press*pow(dens,-3.85)*pow(1+pos/non_dim_radius,-8.55-0.425*area_exponent);
    radiation *= (shock_mdot*shock_height/(shock_speed*shock_speed*shock_speed));

    double dpos_dvel = (5*press - 3*mdot*vel)/(2*radiation + 3*vel*gravity - 5*press*vel*area_exponent/(non_dim_radius+pos));
    double dpress_dvel = -gravity*dpos_dvel - mdot;
    double depress_dvel = 2*(radiation - exchange/(shock_speed*shock_speed))*dpos_dvel/(3*vel);
    depress_dvel -= 5*e_press*((1/vel) + area_exponent*dpos_dvel/(non_dim_radius+pos))/3.;

    return {dpos_dvel,dpress_dvel,depress_dvel};
}

double Cataclysmic_Variable::Get_Landing_Altitude(double cutoff_alt){
    accretion_column.Integrate(this, 0.25, 1e-4, {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))}, cutoff_alt, 0);
    double slope = Flow_Equation(accretion_column.t.back(),  accretion_column.y.back())[0];
    return accretion_column.y.back()[0] - slope*accretion_column.t.back();
}

double Cataclysmic_Variable::Get_Landing_Altitude(){
    accretion_column.Integrate(this, 0.25, 1e-4, {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))});
    double slope = Flow_Equation(accretion_column.t.back(),  accretion_column.y.back())[0];
    return accretion_column.y.back()[0] - slope*accretion_column.t.back();
}

void Cataclysmic_Variable::Bracket_Shock_Height(double integration_limit){
    double height_lim = integration_limit/shock_height;
    Update_Shock_Height((upper_bound+lower_bound)/2);
    double x_final = Get_Landing_Altitude(height_lim);
    // double check that we bracket the problem
    if(x_final<0){ // double check that upper bound is an upper bound
        lower_bound = shock_height;
        Update_Shock_Height(upper_bound);
        x_final = Get_Landing_Altitude(height_lim);
        while(x_final<0){
            upper_bound = 2*upper_bound - lower_bound;
            lower_bound = (upper_bound+lower_bound)/2;
            Update_Shock_Height(upper_bound);
            x_final = Get_Landing_Altitude(height_lim);
        }
    }
    else{ // double check that the lower bound is a lower bound
        upper_bound = shock_height;
        Update_Shock_Height(lower_bound);
        x_final = Get_Landing_Altitude(height_lim);
        while(x_final>0){
            lower_bound = 2*lower_bound - upper_bound;
            upper_bound = (upper_bound+lower_bound)/2;
            Update_Shock_Height(lower_bound);
            x_final = Get_Landing_Altitude(height_lim);
        }
    }
    if(integration_limit>0){
        while(upper_bound-lower_bound > 10*integration_limit){
            Update_Shock_Height((upper_bound+lower_bound)/2);
            x_final = Get_Landing_Altitude(height_lim);
            if(x_final < 0){
                lower_bound = shock_height;
            }
            else{
                upper_bound = shock_height;
            }
        }
    }
    else{
        while(x_final>1e-8){
            Update_Shock_Height((upper_bound+lower_bound)/2);
            x_final = Get_Landing_Altitude();
            if(x_final < 0){
                lower_bound = shock_height;
            }
            else{
                upper_bound = shock_height;
            }
        }
    }
}

void Cataclysmic_Variable::Shock_Height_Shooting(){
    upper_bound = shock_height;
    lower_bound = shock_height;
    Update_Shock_Height((upper_bound+lower_bound)/2);
    if(Get_Landing_Altitude(1./shock_height)<0){
        upper_bound *= 2;
    }
    else{
        lower_bound /=2;
    }
    Bracket_Shock_Height(1.);
    Bracket_Shock_Height(0.1);
    Bracket_Shock_Height(0.01);
    Bracket_Shock_Height(0);
    Update_Shock_Height((upper_bound+lower_bound)/2);
    accretion_column.Integrate(this, 0.25, 1e-4, {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))});
    previous_shock_height = (upper_bound+lower_bound)/2;
}

void Cataclysmic_Variable::Build_Column_Profile(){
    // generate a velocity grid that is has a spacing of ~ kT_grid_spacing
    // interpolate between each point in the RK grid to find a set of velocities to evaluate our integral at
    vector<double> vel_eval = {0.25};
    double kT, dir, kT_new, dv_dkT;
    double kT_left, kT_right;
    int n_possible_vals; // maximum number of grid points that could be in an interval

    kT = accretion_column.y[0][2]*accretion_column.t[0]*pow((non_dim_radius+accretion_column.y[0][0])/(non_dim_radius+1),area_exponent);
    kT *= erg_to_kev*avg_ion_mass*shock_speed*shock_speed/density_const;

    kT_left = kT;

    for(int i=1; i<accretion_column.t.size(); i++){
        kT_right = accretion_column.y[i][2]*accretion_column.t[i]*pow((non_dim_radius+accretion_column.y[i][0])/(non_dim_radius+1),area_exponent);
        kT_right *= erg_to_kev*avg_ion_mass*shock_speed*shock_speed/density_const;

        dir = -1;
        if(kT_right > kT){
            dir = 1;
        }
        n_possible_vals = ceil(abs((kT_right-kT_left)/kT_grid_spacing));
        while(n_possible_vals>0){
            kT_new = kT + dir*kT_grid_spacing;
            if((kT_right-kT_new)*(kT_left-kT_new)<0){ // only true if kT_new is contained in the interval
                dv_dkT = (accretion_column.t[i]-accretion_column.t[i-1])/(kT_right-kT_left);
                kT = kT_new;
                vel_eval.push_back(dv_dkT*(kT-kT_left) + accretion_column.t[i-1]);
                n_possible_vals--;
            }
            else{
                break;
            }
        }
        kT_left = kT_right;
    }
    accretion_column.Integrate(this, 0.25, vel_eval.back(), {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))}, vel_eval);
    int n_points = vel_eval.size();
    velocity.resize(n_points);
    altitude.resize(n_points);
    total_pressure.resize(n_points);
    electron_pressure.resize(n_points);
    electron_density.resize(n_points);
    ion_density.resize(n_points);
    electron_temperature.resize(n_points);
    ion_temperature.resize(n_points);
    volume.resize(n_points);

    double mdot;

    for(int i=0; i<accretion_column.t_eval.size(); i++){
        velocity[i] = accretion_column.t_eval[i]*shock_speed;
        altitude[i] = accretion_column.y_eval[i][0]*shock_height;
        total_pressure[i] = shock_mdot*shock_speed*accretion_column.y_eval[i][1];
        electron_pressure[i] = shock_mdot*shock_speed*accretion_column.y_eval[i][2];
        mdot = (accretion_rate/accretion_area)*pow(1+altitude[i]/radius,-area_exponent);
        electron_density[i] = (mdot/velocity[i])*density_const/avg_ion_mass;
        ion_density[i] = electron_density[i]/avg_atomic_charge;
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*avg_atomic_charge*(total_pressure[i]-electron_pressure[i])/electron_density[i];
    }

    double x0,x1;
    for(int i=1; i<volume.size()-1; i++){
        x0 = (altitude[i+1]+altitude[i])/2;
        x1 = (altitude[i-1]+altitude[i])/2;
        volume[i] = accretion_area*radius*pow(1+x1/radius,area_exponent+1)/(area_exponent+1);
        volume[i] -= accretion_area*radius*pow(1+x0/radius,area_exponent+1)/(area_exponent+1);
    }
    x0 = x1;
    x1 = altitude[volume.size()-1];
    volume[volume.size()-1] = accretion_area*radius*pow(1+x1/radius,area_exponent+1)/(area_exponent+1);
    volume[volume.size()-1] -= accretion_area*radius*pow(1+x0/radius,area_exponent+1)/(area_exponent+1);
    x0 = (altitude[0]+altitude[1])/2;
    x1 = altitude[0];
    volume[0] = accretion_area*radius*pow(1+x1/radius,area_exponent+1)/(area_exponent+1);
    volume[0] -= accretion_area*radius*pow(1+x0/radius,area_exponent+1)/(area_exponent+1);
}

void Cataclysmic_Variable::MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string){
    int n = flux.size();
    double refl_amp, last_refl_amp, sum_refl=0, sum_weights=0;
    valarray<double> apec_flux(n);
    valarray<double> reflected_flux(n);
    valarray<double> flux_error(n);
    valarray<double> apec_parameters(3);
    apec_parameters = {0,metalicity,0};
    valarray<double> refl_parameters(5);
    refl_parameters = {1,0,metalicity,metalicity,incl_angle};

    last_refl_amp = 1-sqrt(1.0-1.0/pow(1+altitude[0]/radius,2));

    for(int i=0; i<altitude.size(); i++){
        apec_parameters[0] = electron_temperature[i];
        if (apec_parameters[0] > 64.0){
            CXX_bremss(energy, apec_parameters, spectrum_num, apec_flux, flux_error, init_string);
        }
        else{
            CXX_apec(energy, apec_parameters, spectrum_num, apec_flux, flux_error, init_string);
        }
        apec_flux *= volume[i]*ion_density[i]*electron_density[i]*1e-14;
        apec_flux /= 4*pi*distance*distance;

        if(refl==1){
            reflected_flux += apec_flux;

            refl_amp = 1-sqrt(1.0-1.0/pow(1+altitude[i]/radius,2));
            sum_refl += refl_amp*ion_density[i]*electron_density[i]*sqrt(electron_temperature[i])*volume[i];
            sum_weights += ion_density[i]*electron_density[i]*sqrt(electron_temperature[i])*volume[i];

            if(abs(refl_amp-last_refl_amp)>reflect_spacing){
                refl_parameters[0] = sum_refl/sum_weights;
                CXX_reflect(energy, refl_parameters, spectrum_num, reflected_flux, flux_error, init_string);
                flux += reflected_flux;
                reflected_flux *= 0;
                sum_refl = 0;
                sum_weights = 0;
                last_refl_amp = refl_amp;
            }
        }
        else{
            flux += apec_flux;
        }
        apec_flux *= 0;
    }
    if(refl==1 && sum_refl>0){
        refl_parameters[0] = sum_refl/sum_weights;
        CXX_reflect(energy, refl_parameters, spectrum_num, reflected_flux, flux_error, init_string);
        flux += reflected_flux;
    }
}

void Cataclysmic_Variable::Print_Properties(){
    cout << "===================================================" << endl;
    cout << "                   mCV Properties                  " << endl;
    cout << "===================================================" << endl;
    cout << " mass:              " << mass/m_sol << " M_solar" << endl;
    cout << " radius:            " << radius/r_sol << " R_solar" << endl;
    cout << " B_field:           " << b_field/1e6 << " MG" << endl;
    if(inverse_mag_radius != 0){
        cout << " R_m/R:             " << (1./inverse_mag_radius)/radius << endl;
    }
    cout << " accretion rate:    " << accretion_rate << " g/s" << endl;
    cout << " accretion rate:    " << shock_mdot << "-->" << accretion_rate/accretion_area << " g/cm2/s" << endl;
    cout << " shock height:      " << shock_height/radius << " (h/R_wd)" << endl;
    cout << " shock temperature: " << electron_temperature[0] << " keV" << endl;
    cout << " cooling ratio:     " << cooling_ratio << endl;
}
