#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "integration.hh"
#include <iostream>
#include <XSFunctions/Utilities/FunctionUtility.h>

using std::copy;
using std::cout;
using std::endl;

valarray<double> Flow_Equation_Wrapper(double t, valarray<double> y, void* cv_instance){
    return ((Cataclysmic_Variable*)(cv_instance))->Flow_Equation(t, y);
}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double b, double metals, double fractional_area, double theta, double dist, int reflection):
    mass(m), b_field(b), inverse_mag_radius(0), distance(dist), metalicity(metals), pressure_ratio(1.), incl_angle(theta), refl(reflection),
    accretion_column(Flow_Equation_Wrapper, this, 3)
{}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double b, double metals, double luminosity, double fractional_area, double theta, double dist, int reflection):
    mass(m), b_field(b), inverse_mag_radius(0), distance(dist), metalicity(metals), pressure_ratio(1.), incl_angle(theta), refl(reflection),
    accretion_column(Flow_Equation_Wrapper, this, 3)
{
    Radius_Shooting(100000);
    Set_Accretion_Rate(luminosity);
    accretion_area = fractional_area*4.*pi*radius*radius;
    accretion_rate /= accretion_area;
    Set_Abundances(metalicity);
    Set_Pre_Shock_Speed(5);
    Set_Cooling_Ratio();
}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double metals, double luminosity, double fractional_area, double theta, double dist, int reflection, double r_m):
    mass(m), inverse_mag_radius(1/r_m), distance(dist), metalicity(metals), pressure_ratio(1.), incl_angle(theta), refl(reflection),
    accretion_column(Flow_Equation_Wrapper, this, 3)
{
    Radius_Shooting(100000);
    Set_Accretion_Rate(luminosity);
    b_field = sqrt(32.*accretion_rate)*pow(grav_const*mass, 1./4.)*pow(r_m,7./4.)*pow(radius,-3);
    accretion_area = fractional_area*4.*pi*radius*radius;
    accretion_rate /= accretion_area;
    Set_Abundances(metalicity);
    Set_Pre_Shock_Speed(5);
    Set_Cooling_Ratio();
}

void Cataclysmic_Variable::Set_Inverse_Mag_Radius(double mag_ratio){
    inverse_mag_radius = 1./(mag_ratio*radius);
    b_field = sqrt(32.*accretion_area*accretion_rate)*pow(grav_const*mass, 1./4.)*pow(1./inverse_mag_radius,7./4.)*pow(radius,-3);
    Set_Pre_Shock_Speed(5);
    Set_Cooling_Ratio();
}

void Cataclysmic_Variable::Set_Cooling_Constants(){
    charge.resize(atomic_charge.size());
    copy(begin(atomic_charge), end(atomic_charge), begin(charge));
    avg_ion_mass = (abundances*atomic_mass).sum()*amu_to_g;
    avg_atomic_charge = (abundances*charge).sum();
    double avg_charge_squared = (abundances*charge*charge).sum();
    thermal_constant = avg_atomic_charge/(avg_atomic_charge + avg_ion_mass/electron_mass);

    bremss_constant = sqrt(32./(27.*pow(pi,3)))*planck_const*planck_const*pow(fine_structure_constant/electron_mass, 3)*gaunt_factor
                      *sqrt(pow(thermal_constant,3))*avg_charge_squared/avg_atomic_charge;

    cooling_ratio_const = 8.07e-2*avg_atomic_charge*pow(electron_mass, 3.85)*pow(shock_ratio-1, 2)/pow(shock_ratio, 5.85)*pow(pressure_ratio/(pressure_ratio+1),2)
                          /(gaunt_factor*avg_charge_squared*boltz_const*boltz_const*pow(accretion_rate, 1.85)*pow(thermal_constant, 3.85));

    exchange_constant = sqrt(2/pow(pi,3))*pow(fine_structure_constant*planck_const*light_speed,2)*sqrt(pow(thermal_constant, 5))
                        /pow(electron_mass,3);
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

void Cataclysmic_Variable::Set_Pre_Shock_Speed(int num_itters){
    double integral = (39.*sqrt(3.) - 20*pi)/96.; // value of integral from EQ 7a of Wu 1994 DOI: 10.1086/174103
    pre_shock_speed = sqrt(2*grav_const*mass*((1./radius) - inverse_mag_radius));
    shock_height = pow(pre_shock_speed,3.)*integral*sqrt((avg_atomic_charge+1)/avg_atomic_charge)/(2.*bremss_constant*accretion_rate);
    for(int i = 0; i < num_itters; i++){
        pre_shock_speed = sqrt(2*grav_const*mass*((1./(radius+shock_height)) - inverse_mag_radius));
        shock_height = pow(pre_shock_speed,3.)*integral*sqrt((avg_atomic_charge+1)/avg_atomic_charge)/(2.*bremss_constant*accretion_rate);
    }
}

void Cataclysmic_Variable::Set_Cooling_Ratio(){
    cooling_ratio = cooling_ratio_const*pow(pre_shock_speed,5.85)*pow(accretion_area, -0.425)*pow(b_field, 2.85);
}

valarray<double> Cataclysmic_Variable::Chandrasekhar_White_Dwarf_Equation(double eta, valarray<double> phi_psi, void* y0){
    double y_0 = *(double *)y0;
    return {phi_psi[1],(-2./eta)*phi_psi[1] - sqrt(pow(phi_psi[0]*phi_psi[0] - 1./(y_0*y_0), 3))};
}

void Cataclysmic_Variable::Radius_Shooting(int max_itter){
    int itters = 0;
    double y_0 = 2;
    double solved_mass = 2*solar_mass;
    const double pressure_constant = pi*pow(electron_mass,4)*pow(light_speed,5)/(3*pow(planck_const,3));
    const double density_constant = (8.*pi/3.)*wd_mol_mass*amu_to_g*pow(electron_mass*light_speed/planck_const,3);
    Integrator white_dwarf(Chandrasekhar_White_Dwarf_Equation, 2);
    double eta0 = 1e-32;
    vector<double> eta;
    valarray<double> bounds;
    valarray<double> start(2);

    while(abs(solved_mass/mass - 1.) > relative_err && itters < max_itter){
        bounds = {1./y_0, 1e3};
        start[0] = 1 - eta0*eta0*sqrt(pow(1.0-1./(y_0*y_0), 3))/2;
        start[1] = 1 - eta0*sqrt(pow(1.0-1./(y_0*y_0), 3));
        white_dwarf.Integrate(&y_0, eta0, 1e5, start, bounds);
        eta = white_dwarf.t;
        solved_mass = (-4*pi/(density_constant*density_constant))*pow((2*pressure_constant/(pi*grav_const)),3./2.)*(eta.back()*eta.back())*white_dwarf.y.back()[1];
        y_0 /= (eta.back()*white_dwarf.y.back()[1]*y_0/2)*(mass/solved_mass - 1) + 1;
        itters++;
    }
    bounds = {1./y_0, 1e3};
    start[0] = 1 - eta0*eta0*sqrt(pow(1.0-1./(y_0*y_0), 3))/2;
    start[1] = 1 - eta0*sqrt(pow(1.0-1./(y_0*y_0), 3));
    white_dwarf.Integrate(&y_0, eta0, 1e5, start, bounds);
    eta = white_dwarf.t;
    radius = sqrt(2*pressure_constant/(pi*grav_const*density_constant*density_constant))*(eta.back()/y_0);
}

valarray<double> Cataclysmic_Variable::Flow_Equation(double vel, valarray<double> pos_pres_epres){
    double position = pos_pres_epres[0];
    double pressure = pos_pres_epres[1];
    double elec_pressure = pos_pres_epres[2];
    double n_e = thermal_constant*accretion_rate/(pre_shock_speed*electron_mass*vel);
    double temp = accretion_rate*pre_shock_speed*elec_pressure/(boltz_const*n_e);
    double coulomb_logarithm = 15.9 + log(temp*1e-8) - 0.5*log(n_e*1e-16);

    double energy_loss = bremss_constant*sqrt(elec_pressure/pow(vel,3))*(1 + cooling_ratio*(pow(shock_ratio,alpha+beta)/pow((shock_ratio-1),alpha))
                                                                           *(pow(((pressure_ratio+1)/pressure_ratio),alpha))*pow(elec_pressure,alpha)*pow(vel,beta));
    double energy_exchange = exchange_constant*coulomb_logarithm*(pressure-((avg_atomic_charge+1)/avg_atomic_charge)*elec_pressure)/sqrt(pow(vel,5))
                             *(abundances*charge*charge*electron_mass/(amu_to_g*atomic_mass*sqrt(pow(elec_pressure + (pressure-elec_pressure)*avg_atomic_charge*electron_mass/(amu_to_g*atomic_mass),3)))).sum();

    double gravity = (grav_const*mass/shock_height)/((radius+position)*(radius+position));

    double dpos_dvel = pre_shock_speed*pre_shock_speed*(5*pressure - 3*vel)/(3*gravity + 2*(shock_height*accretion_rate/pre_shock_speed)*energy_loss);
    double dpress_dvel = (-1.*gravity/(pre_shock_speed*pre_shock_speed*vel))*dpos_dvel - 1.;
    double depress_dvel = (2./3.)*(shock_height*accretion_rate/pow(pre_shock_speed,3))*((energy_loss - energy_exchange/(pre_shock_speed*pre_shock_speed))/vel)
                                                                                      *dpos_dvel - (5./3.)*elec_pressure/vel;
    return {dpos_dvel,dpress_dvel,depress_dvel};
}

void Cataclysmic_Variable::Shock_Height_Shooting(int max_itter){
    double slope, error = 100.;
    int itters = 0;
    while(itters < max_itter && error > absolute_err*1e3){
        accretion_column.Integrate(this, 1./shock_ratio, 1e-4, {1., ((shock_ratio-1)/shock_ratio), ((shock_ratio-1)/shock_ratio)*(pressure_ratio/(pressure_ratio+1))});
        slope = Flow_Equation(accretion_column.t.back(),  accretion_column.y.back())[0];
        error = accretion_column.y.back()[0] - accretion_column.t.back()*slope;
        shock_height *= 1-error;
        error = abs(error);
        pre_shock_speed = sqrt(2*grav_const*mass*((1./(radius+shock_height)) - inverse_mag_radius));
        Set_Cooling_Ratio();
        itters++;
    }

    int j=accretion_column.t.size()-1;
    while((erg_to_kev*electron_mass*pre_shock_speed*pre_shock_speed/thermal_constant)*accretion_column.t[j]*accretion_column.y[j][2] < 1.){
        j--;
    }
    double delta_vel = ((1./shock_ratio) - accretion_column.t[j+1])/(299.);
    vector<double> vel_eval(300);
    for(int i=0; i<300; i++){
        vel_eval[i] = 1./shock_ratio - i*delta_vel;
    }
    accretion_column.Integrate(this, 1./shock_ratio, vel_eval.back()*0.5, {1., ((shock_ratio-1)/shock_ratio), ((shock_ratio-1)/shock_ratio)*(pressure_ratio/(pressure_ratio+1))}, vel_eval);

    vector<valarray<double>> profile;
    double kT,kT_last;
    profile.push_back({accretion_column.t[299], accretion_column.y[299][0],accretion_column.y[299][1],accretion_column.y[299][2]});
    kT_last = erg_to_kev*(electron_mass*pre_shock_speed*pre_shock_speed/thermal_constant)*accretion_column.t[299]*accretion_column.y[299][2];
    for(int i=298; i>=0; i--){
        kT = erg_to_kev*(electron_mass*pre_shock_speed*pre_shock_speed/thermal_constant)*accretion_column.t[i]*accretion_column.y[i][2];
        if(abs(kT-kT_last)>1.){
            profile.push_back({accretion_column.t[i], accretion_column.y[i][0],accretion_column.y[i][1],accretion_column.y[i][2]});
            kT_last = kT;
        }
    }

    velocity.resize(profile.size());
    altitude.resize(profile.size());
    electron_temperature.resize(profile.size());
    ion_temperature.resize(profile.size());
    electron_density.resize(profile.size());
    ion_density.resize(profile.size());
    electron_pressure.resize(profile.size());
    total_pressure.resize(profile.size());
    double vel, alt, pres, epres; // normalzied velocity, altitude, and electron pressure
    int i = 0;
    for(valarray<double> fluid_element:profile){
        vel = fluid_element[0];
        alt = fluid_element[1];
        pres = fluid_element[2];
        epres = fluid_element[3];

        velocity[i] = pre_shock_speed*vel;
        altitude[i] = shock_height*alt;
        total_pressure[i] = accretion_rate*pre_shock_speed*pres;
        electron_pressure[i] = accretion_rate*pre_shock_speed*epres;
        electron_density[i] = (accretion_rate/pre_shock_speed)*(thermal_constant/electron_mass)/vel;
        ion_density[i] = electron_density[i]/avg_atomic_charge;
        electron_temperature[i] = electron_pressure[i]/electron_density[i];
        ion_temperature[i] = (total_pressure[i]-electron_pressure[i])/ion_density[i];
        i++;
    }
}

void Cataclysmic_Variable::MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string){
    int n = flux.size();
    int m = altitude.size()-1;
    double segment_height = 0;
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
            segment_height = abs(altitude[i+1]-altitude[i])/2.;
        }
        else if(i==altitude.size()-1){
            segment_height = abs(altitude[i]-altitude[i-1])/2.;
        }
        else{
            segment_height = 0.5*abs(altitude[i+1]-altitude[i-1]);
        }
        flux_from_layer *= (accretion_area*segment_height*ion_density[i]*electron_density[i]*1e-14)/avg_atomic_charge;
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
    cout << "                   McV Properties                  " << endl;
    cout << "===================================================" << endl;
    cout << " radius: " << radius/solar_radius << " R_solar" << endl;
    if(inverse_mag_radius != 0){
        cout << " R_m/R: " << (1./inverse_mag_radius)/radius << endl;
    }
    cout << " shock height " << shock_height/radius << " (h/R_wd)" << endl;
    cout << " shock temperature " <<  electron_temperature[0] << " keV" << endl;
    cout << " cooling ratio " << cooling_ratio << endl;


}
