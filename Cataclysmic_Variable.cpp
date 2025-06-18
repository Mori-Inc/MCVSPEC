#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "integration.hh"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <XSFunctions/Utilities/FunctionUtility.h>

using std::copy;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::stringstream;

Cataclysmic_Variable::Cataclysmic_Variable(double m, double b, double metals, double luminosity, double fractional_area, double theta, double dist, int reflection):
    mass(m), b_field(b), inverse_mag_radius(0), distance(dist), metalicity(metals), pressure_ratio(1.), incl_angle(theta), refl(reflection),
    accretion_column(Flow_Equation, 2)
{
    Radius_Shooting(100000);
    inverse_mag_radius = 0;
    Set_Accretion_Rate(luminosity);
    accretion_area = fractional_area*4.*pi*radius*radius;
    accretion_rate /= accretion_area;
    Set_Abundances(metalicity);
    Set_Pre_Shock_Speed(5);
    Set_Cooling_Ratio();
}

Cataclysmic_Variable::Cataclysmic_Variable(double m, double metals, double luminosity, double fractional_area, double theta, double dist, int reflection, double r_m):
    mass(m), inverse_mag_radius(1/r_m), distance(dist), metalicity(metals), pressure_ratio(1.), incl_angle(theta), refl(reflection),
    accretion_column(Flow_Equation, 2)
{
    Radius_Shooting(100000);
    Set_Accretion_Rate(luminosity);
    b_field = sqrt(accretion_rate*sqrt(grav_const*mass*pow(r_m,7)))/(2*radius*radius*radius);
    accretion_area = fractional_area*4.*pi*radius*radius;
    accretion_rate /= accretion_area;
    Set_Abundances(metalicity);
    Set_Pre_Shock_Speed(5);
    Set_Cooling_Ratio();
}

void Cataclysmic_Variable::Set_Inverse_Mag_Radius(double mag_ratio){
    inverse_mag_radius = 1./(mag_ratio*radius);
    b_field = sqrt(accretion_rate*sqrt(grav_const*mass*pow(mag_ratio*radius,7)))/(2*radius*radius*radius);
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

    exchange_constant = sqrt(2/pow(pi,3))*pow(fine_structure_constant*planck_const*light_speed,2)*coulomb_logarithm*sqrt(pow(thermal_constant, 5))
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

    kt_const = electron_mass*(avg_atomic_charge + avg_ion_mass/electron_mass)/(1.+avg_atomic_charge);
    electron_density_const = avg_atomic_charge/(avg_atomic_charge + avg_ion_mass/electron_mass);
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
    ifstream file("chandrasekhar.txt");
    string line, radius_string, mass_string;
    double in_mass=0, in_radius=0, old_mass=0, old_radius=0;

    do{
        old_mass = in_mass;
        old_radius = in_radius;
        getline(file, line);
        stringstream stream(line);
        getline(stream, mass_string, ',');
        getline(stream, radius_string, ',');
        in_mass = stod(mass_string);
        in_radius = stod(radius_string);
    }while(mass > old_mass);

    radius = ((in_radius-old_radius)/(in_mass-old_mass))*mass + old_radius - ((in_radius-old_radius)/(in_mass-old_mass))*old_mass;
}

valarray<double> Cataclysmic_Variable::Flow_Equation(double vel, valarray<double> pos_pres, void* my_class_instance){
    Cataclysmic_Variable* cv_instance = (Cataclysmic_Variable*)my_class_instance;
    double cooling_ratio = cv_instance->cooling_ratio;
    double pressure_ratio = cv_instance->pressure_ratio;
    double avg_atomic_charge = cv_instance->avg_atomic_charge;
    valarray<double> abundances = cv_instance->abundances;
    valarray<double> charge = cv_instance->charge;
    double bremss_constant = cv_instance->bremss_constant;
    double exchange_constant = cv_instance->exchange_constant;
    double pre_shock_speed = cv_instance->pre_shock_speed;
    double shock_height = cv_instance->shock_height;
    double accretion_rate = cv_instance->accretion_rate;

    double pressure = pos_pres[1];
    double energy_loss = bremss_constant*sqrt(pos_pres[1]/pow(vel,3))*(1 + cooling_ratio*(pow(shock_ratio,alpha+beta)/pow((shock_ratio-1),alpha))
                                                                           *(pow(((pressure_ratio+1)/pressure_ratio),alpha))*pow(pressure,alpha)*pow(vel,beta));
    double energy_exchange = exchange_constant*(1-vel-((avg_atomic_charge+1)/avg_atomic_charge)*pressure)/sqrt(pow(vel,5))
                             *(abundances*charge*charge*electron_mass/(amu_to_g*atomic_mass*sqrt(pow(pressure + (1-vel-pressure)*avg_atomic_charge*electron_mass/(amu_to_g*atomic_mass),3)))).sum();
    double dpos_dvel = (pow(pre_shock_speed,3)/(shock_height*accretion_rate))*(5. - 8.*vel)/(2*energy_loss);
    double dpress_dvel = ((2.0-8.*vel)/3. - ((5.0-8.*vel)/3.)*(energy_exchange/(pre_shock_speed*pre_shock_speed*energy_loss)))/vel;

    return {dpos_dvel,dpress_dvel};
}

void Cataclysmic_Variable::Shock_Height_Shooting(int max_itter){
    double slope, error = 100.;
    int itters = 0;
    while(itters < max_itter && error > absolute_err*1e3){
        accretion_column.Integrate(this, 1./shock_ratio, 1e-4, {1., ((shock_ratio-1)/shock_ratio)*(pressure_ratio/(pressure_ratio+1))});
        slope = Flow_Equation(accretion_column.t.back(),  accretion_column.y.back(), this)[0];
        error = accretion_column.y.back()[0] - accretion_column.t.back()*slope;
        if(error > 1){
            error = 0.99;
        }
        shock_height *= 1-error;
        error = abs(error);
        pre_shock_speed = sqrt(2*grav_const*mass*((1./(radius+shock_height)) - inverse_mag_radius));
        Set_Cooling_Ratio();
        itters++;
    }

    accretion_column.Integrate(this, 1./shock_ratio, 1e-4, {1., ((shock_ratio-1)/shock_ratio)*(pressure_ratio/(pressure_ratio+1))}, 0.01);

    altitude.resize(accretion_column.t.size());
    electron_temperature.resize(accretion_column.t.size());
    ion_temperature.resize(accretion_column.t.size());
    electron_density.resize(accretion_column.t.size());
    ion_density.resize(accretion_column.t.size());
    double vel, alt, pres, e_temp; // normalzied velocity, altitude, and electron pressure
    int n_points=0;
    altitude[0] = shock_height*accretion_column.y[0][0];
    electron_temperature[0] = erg_to_kev*electron_mass*pre_shock_speed*pre_shock_speed*accretion_column.t[0]*accretion_column.y[0][1]/thermal_constant;
    ion_temperature[0] = erg_to_kev*electron_mass*pre_shock_speed*pre_shock_speed*avg_atomic_charge*vel*(1.0-accretion_column.t[0]-accretion_column.y[0][1])/thermal_constant;
    electron_density[0] = (accretion_rate/pre_shock_speed)*(thermal_constant/electron_mass)/accretion_column.t[0];
    ion_density[0] = electron_density[0]/avg_atomic_charge;

    for(int i = i; i < accretion_column.t.size(); i++){
        vel = accretion_column.t[i];
        alt = accretion_column.y[i][0];
        pres = accretion_column.y[i][1];

        e_temp = erg_to_kev*electron_mass*pre_shock_speed*pre_shock_speed*vel*pres/thermal_constant;
        if(abs(e_temp-electron_temperature[n_points]) < 1 || electron_temperature[n_points]<1){
            continue;
        }

        altitude[i] = shock_height*alt;
        electron_temperature[i] = erg_to_kev*electron_mass*pre_shock_speed*pre_shock_speed*vel*pres/thermal_constant;
        ion_temperature[i] = erg_to_kev*electron_mass*pre_shock_speed*pre_shock_speed*avg_atomic_charge*vel*(1.0-vel-pres)/thermal_constant;
        electron_density[i] = (accretion_rate/pre_shock_speed)*(thermal_constant/electron_mass)/vel;
        ion_density[i] = electron_density[i]/avg_atomic_charge;
        n_points++;
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
        if(electron_temperature[i] < 0.5){
            continue;
        }
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
    cout << " mass: " << mass/solar_mass << " M_solar" << endl;
    cout << " radius: " << radius/solar_radius << " R_solar" << endl;
    cout << " B field: " << b_field/1e6 <<  " MG" << endl;
    if(inverse_mag_radius > 0){
        cout << " R_m/R_wd: " << 1./(radius*inverse_mag_radius) << endl;
    }
    cout << " accretion rate: " << accretion_rate << " g/cm2/s" << endl;
    cout << " shock height " << shock_height/radius << " (h/R_wd)" << endl;
    cout << " shock temperature " <<  electron_temperature[0] << " keV" << endl;
    cout << " cooling ratio " << cooling_ratio << endl;


}
