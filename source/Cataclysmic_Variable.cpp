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
using std::floor;

valarray<double> Flow_Equation_Wrapper(double t, valarray<double> y, void* cv_instance){
    return ((Cataclysmic_Variable*)(cv_instance))->Flow_Equation(t, y);
}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double b, double metals, double fractional_area, double theta, double n, double dist, int reflection):
    mass(m), b_field(b), inverse_mag_radius(0), distance(dist), metalicity(metals), pressure_ratio(1.), incl_angle(theta), area_exponent(n),  refl(reflection),
    accretion_column(Flow_Equation_Wrapper, this, 3)
{}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double b, double metals, double luminosity, double fractional_area, double theta, double n, double dist, int reflection):
    mass(m), b_field(b), inverse_mag_radius(0), distance(dist), metalicity(metals), pressure_ratio(1.), incl_angle(theta), area_exponent(n),  refl(reflection),
    accretion_column(Flow_Equation_Wrapper, this, 3)
{
    Set_Radius();
    Set_Accretion_Rate(luminosity);
    accretion_area = fractional_area*4.*pi*radius*radius;
    Set_Abundances(metalicity);
    Set_Shock_Speed(5);
    Set_Cooling_Ratio();
}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double metals, double luminosity, double fractional_area, double theta, double n, double dist, int reflection, double r_m):
    mass(m), inverse_mag_radius(1/r_m), distance(dist), metalicity(metals), pressure_ratio(1.), incl_angle(theta), area_exponent(n),  refl(reflection),
    accretion_column(Flow_Equation_Wrapper, this, 3)
{
    Set_Radius();
    Set_Accretion_Rate(luminosity);
    b_field = sqrt(32*accretion_rate*sqrt(grav_const*mass/pow(inverse_mag_radius,7)))/(radius*radius*radius);
    accretion_area = fractional_area*4.*pi*radius*radius;
    Set_Abundances(metalicity);
    Set_Shock_Speed(5);
    Set_Cooling_Ratio();
}

void Cataclysmic_Variable::Set_Inverse_Mag_Radius(double mag_ratio){
    inverse_mag_radius = 1./(mag_ratio*radius);
    b_field = sqrt(32*accretion_rate*sqrt(grav_const*mass/pow(inverse_mag_radius,7)))/(radius*radius*radius);
    Set_Shock_Speed(5);
    Set_Cooling_Ratio();
}

void Cataclysmic_Variable::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
    avg_ion_mass = (abundances*atomic_mass).sum()*amu_to_g;
    avg_atomic_charge = (abundances*atomic_charge).sum();
    double avg_charge_squared = (abundances*atomic_charge*atomic_charge).sum();
    double avg_charge_sqr_over_mass = (abundances*atomic_charge*atomic_charge/atomic_mass).sum()/amu_to_g;
    density_const = avg_atomic_charge/(1 + m_e*avg_atomic_charge/avg_ion_mass + 1);

    force_const = grav_const*mass/(radius*radius);
    cooling_ratio_const = 8.07e-2*avg_atomic_charge*pow(b_field, 2.85)*pow(avg_ion_mass/density_const,3.85);
    cooling_ratio_const /= gaunt_factor*avg_charge_squared*k_b*k_b*pow(accretion_area,0.425);
    coulomb_log_const = 0.5*log(2*m_e/(pi*alpha*c)) + 1.5*log(avg_ion_mass/(hbar*density_const));
    exchange_const = 4*(alpha*hbar*c)*(alpha*hbar*c)*sqrt(2*pi*m_e*pow((density_const/avg_ion_mass),5))*avg_charge_sqr_over_mass;
    bremss_const = sqrt(512*pi/(27*m_e*m_e*m_e))*alpha*alpha*alpha*hbar*hbar*gaunt_factor;
    bremss_const *= (avg_charge_squared/avg_atomic_charge)*sqrt(pow(density_const/avg_ion_mass,3));
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

void Cataclysmic_Variable::Set_Accretion_Rate(double luminosity){
    accretion_rate = luminosity/(grav_const*mass*((1./radius) - inverse_mag_radius));
}

void Cataclysmic_Variable::Set_Shock_Speed(int num_itters){
    double integral = (39.*sqrt(3.) - 20*pi)/96.; // value of integral from EQ 7a of Wu 1994 DOI: 10.1086/174103
    shock_speed = sqrt(2*grav_const*mass*((1./radius) - inverse_mag_radius));
    shock_height = pow(shock_speed,3.)*integral*sqrt((avg_atomic_charge+1)/avg_atomic_charge)/(2.*bremss_const*accretion_rate);
    for(int i = 0; i < num_itters; i++){
        shock_speed = sqrt(2*grav_const*mass*((1./(radius+shock_height)) - inverse_mag_radius));
        shock_height = pow(shock_speed,3.)*integral*sqrt((avg_atomic_charge+1)/avg_atomic_charge)/(2.*bremss_const*accretion_rate);
    }
}

void Cataclysmic_Variable::Set_Cooling_Ratio(){
    shock_mdot = accretion_rate/(accretion_area*pow(1+shock_height/radius, area_exponent));
    non_dim_radius = radius/shock_height;
    cooling_ratio = cooling_ratio_const*pow(shock_speed,5.85)/pow(shock_mdot, 1.85);
}

void Cataclysmic_Variable::Set_Radius(){
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
    radius = white_dwarf_radius[left_ind] + (delta_r/delta_m)*(mass-white_dwarf_mass[left_ind]);
}

valarray<double> Cataclysmic_Variable::Flow_Equation(double vel, valarray<double> pos_pres_epres){
    double pos = pos_pres_epres[1];
    double press = pos_pres_epres[1];
    double e_press = pos_pres_epres[1];

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

double Cataclysmic_Variable::Get_Landing_Altitude(){
    accretion_column.Integrate(this, 0.25, 1e-4, {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))});
    double slope = Flow_Equation(accretion_column.t.back(),  accretion_column.y.back())[0];
    return accretion_column.y.back()[0] - slope*accretion_column.t.back();
}

void Cataclysmic_Variable::Shock_Height_Shooting(){
    double x_final = Get_Landing_Altitude();
    double upp_bound = shock_height;
    double low_bound = 100;
    if (x_final < 0){
        low_bound = shock_height;
        upp_bound = 10*radius;
        if(upp_bound*inverse_mag_radius > 0.9){
            upp_bound = 0.9/inverse_mag_radius;
        }
    }
    while(abs(x_final) > 1e-8){
        shock_height = (upp_bound+low_bound)/2;
        Set_Cooling_Ratio();
        x_final = Get_Landing_Altitude();
        if(x_final < 0){
            low_bound = shock_height;
        }
        else{
            upp_bound = shock_height;
        }
    }
    shock_height = (upp_bound+low_bound)/2;
    Set_Cooling_Ratio();
    accretion_column.Integrate(this, 0.25, 1e-4, {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))});
}

void Cataclysmic_Variable::Build_Column_Profile(){
    double kT = 0;
    double pressure, pressure_e, n_electron;
    int i = accretion_column.t.size()-1;
    while(kT < 0.5){
        pressure_e = shock_mdot*shock_speed*accretion_column.y[i][2];
        n_electron = pow((non_dim_radius+1)/(non_dim_radius+accretion_column.y[i][0]),area_exponent)/accretion_column.t[i];
        n_electron *= density_const*shock_mdot/(avg_ion_mass*shock_speed);
        kT = erg_to_kev*pressure_e/n_electron;
        i--;
    }
    i++;

    // check profile to make sure delta kT ~< 1 keV everywhere
    pressure_e = shock_mdot*shock_speed*accretion_column.y[0][2];
    n_electron = pow((non_dim_radius+1)/(non_dim_radius+accretion_column.y[0][0]),area_exponent)/accretion_column.t[0];
    n_electron *= density_const*shock_mdot/(avg_ion_mass*shock_speed);
    kT = erg_to_kev*pressure_e/n_electron;
    vector<double> vel_eval;
    double delta_kT, delta_v;
    double n_add_points;
    for(int j=1; j<i; j++){
        pressure_e = shock_mdot*shock_speed*accretion_column.y[j][2];
        n_electron = pow((non_dim_radius+1)/(non_dim_radius+accretion_column.y[j][0]),area_exponent)/accretion_column.t[j];
        n_electron *= density_const*shock_mdot/(avg_ion_mass*shock_speed);
        delta_kT = abs(erg_to_kev*pressure_e/n_electron - kT);
        if(delta_kT > kT_grid_spacing){
            delta_v = accretion_column.t[j] - accretion_column.t[j-1];
            n_add_points = floor(delta_kT/kT_grid_spacing);
            for(int n=1; n<=n_add_points; n++){
                vel_eval.push_back(accretion_column.t[j-1] + n*delta_v/(n_add_points+1));
            }
        }
    }
    accretion_column.Integrate(this, 0.25, accretion_column.t[i], {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))}, vel_eval);
    int n_points = accretion_column.t.size() + vel_eval.size();
    velocity.resize(n_points);
    altitude.resize(n_points);
    electron_temperature.resize(n_points);
    ion_temperature.resize(n_points);
    electron_density.resize(n_points);
    ion_density.resize(n_points);
    electron_pressure.resize(n_points);
    total_pressure.resize(n_points);

    int k=0;
    for(int j = 0; j < accretion_column.t.size(); j++){
        pressure = shock_mdot*shock_speed*accretion_column.y[j][1];
        pressure_e = shock_mdot*shock_speed*accretion_column.y[j][2];
        n_electron = pow((non_dim_radius+1)/(non_dim_radius+accretion_column.y[j][0]),area_exponent)/accretion_column.t[j];
        n_electron *= density_const*shock_mdot/(avg_ion_mass*shock_speed);
        velocity[j+k] = accretion_column.t[j]*shock_speed;
        altitude[j+k] = accretion_column.y[j][0]*shock_height;
        if(accretion_column.t_eval[k] > accretion_column.t[j]){
            pressure = shock_mdot*shock_speed*accretion_column.y_eval[k][1];
            pressure_e = shock_mdot*shock_speed*accretion_column.y_eval[k][2];
            n_electron = pow((non_dim_radius+1)/(non_dim_radius+accretion_column.y_eval[k][0]),area_exponent)/accretion_column.t_eval[k];
            n_electron *= density_const*shock_mdot/(avg_ion_mass*shock_speed);
            velocity[j+k] = accretion_column.t_eval[k]*shock_speed;
            altitude[j+k] = accretion_column.y_eval[k][0]*shock_height;
            k++;
            j--;
        }
        total_pressure[j+k] = pressure;
        electron_pressure[j+k] = pressure_e;
        electron_density[j+k] = n_electron;
        ion_density[j+k] = n_electron/avg_atomic_charge;
        electron_temperature[j+k] = erg_to_kev*pressure_e/n_electron;
        ion_temperature[j+k] = erg_to_kev*avg_atomic_charge*(pressure-pressure_e)/n_electron;
    }
}

void Cataclysmic_Variable::MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string){
    int n = flux.size();
    int m = altitude.size()-1;
    double segment_volume = 0;
    double segment_top, segment_bottom;
    valarray<double> flux_from_layer(n);
    valarray<double> flux_error(n);
    valarray<double> apec_parameters(3);
    valarray<double> refl_parameters(5);


    refl_parameters[0] = 1-sqrt(1.0-1.0/pow(1+altitude[m]/radius,2));
    refl_parameters[1] = 0.;
    refl_parameters[2] = metalicity;
    refl_parameters[3] = metalicity;
    refl_parameters[4] = incl_angle;

    for(int i=0; i<altitude.size(); i++){
        for(int j=0; j<n; j++){
            flux[j] += flux_from_layer[j];
            flux_from_layer[j] = 0;
        }

        apec_parameters[0] = electron_temperature[i];
        apec_parameters[1] = metalicity;
        apec_parameters[2] = 0.0;

        if (apec_parameters[0] > 64.0){
            CXX_bremss(energy, apec_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }
        else{
            CXX_apec(energy, apec_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }

        // normalizes spectrum appropriately
        if(i==0){
            segment_top = altitude[0];
            segment_bottom = (altitude[1]+altitude[0])/2;
        }
        else if(i==altitude.size()-1){
            segment_top = (altitude[i]+altitude[i-1])/2;
            segment_bottom = altitude[i];
        }
        else{
            segment_top = (altitude[i]+altitude[i-1])/2;
            segment_bottom = (altitude[i]+altitude[i+1])/2;
        }
        segment_volume = (accretion_area*radius/(area_exponent+1))*(pow(1+segment_top/radius,area_exponent+1)-pow(1+segment_bottom/radius,area_exponent+1));
        flux_from_layer *= segment_volume*ion_density[i]*electron_density[i]*1e-14;
        flux_from_layer /= 4*pi*distance*distance;

        // apply reflect to each slice
        if (refl == 1){
            refl_parameters[0] = 1-sqrt(1.0-1.0/pow(1+altitude[i]/radius,2));
            CXX_reflect(energy, refl_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }
    }

    for(int j=0; j<n; j++){
        flux[j] += flux_from_layer[j];
        flux_from_layer[j] = 0;
    }

    if(refl == 2){
        CXX_reflect(energy, refl_parameters, spectrum_num, flux, flux_error, init_string);
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
    cout << " accretion rate:    " << shock_mdot << "-->" << accretion_rate/accretion_area << "g/cm2/s" << endl;
    cout << " shock height:      " << shock_height/radius << " (h/R_wd)" << endl;
    cout << " shock temperature: " << electron_temperature[0] << " keV" << endl;
    cout << " cooling ratio:     " << cooling_ratio << endl;
}
