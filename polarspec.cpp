#include "Cataclysmic_Variable.hh"

extern "C"
void Polarspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    int reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    double b_field = params[1]*1e6; // surface b field [gauss]
    double fractional_area = params[2]; //fractional accretion area
    double luminosity = params[3]*1e33; // luminosity [ergs/s]
    double mass = params[4]*solar_mass; // WD mass [grams]
    double col_abund = params[5]; // accretion column abundance [solar abundances]
    double cos_incl = params[6]; // cos inclination angle
    double source_distance = params[7]*pc_to_cm; // source distnace [cm]

    Cataclysmic_Variable polar(mass, b_field, col_abund, luminosity, fractional_area, cos_incl, source_distance, reflection_sel);
    polar.Shock_Height_Shooting(1000);
    polar.MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}

/*
valarray<double> Flow_Equation(double vel, valarray<double> pos_pres, void* my_class_instance){
    double cooling_ratio = 10.;
    double pressure_ratio = 0.5;
    double avg_atomic_charge = 1.;
    valarray<double> abundances = {9.09672976e-01, 8.88750497e-02, 3.30211290e-04, 1.01883373e-04,
           7.74131702e-04, 1.11889776e-04, 3.45675731e-05, 2.68353528e-06,
           3.22933906e-05, 1.47367022e-05, 3.30211290e-06, 2.08315111e-06,
           4.25726953e-05, 1.61921790e-06};
    valarray<double> charge = {1,2,6,7,8,10,12,13,14,16,18,20,26,28};
    double bremss_constant = cv_instance->bremss_constant;
    double exchange_constant = cv_instance->exchange_constant;
    double pre_shock_speed = cv_instance->pre_shock_speed;
    double shock_height = cv_instance->shock_height;
    double accretion_rate = cv_instance->accretion_rate;


    double pressure = pos_pres[1];
    double energy_loss = bremss_constant*sqrt(pos_pres[1]/pow(vel,3))*(1 + cooling_ratio*(pow(shock_ratio,alpha+beta)/pow((shock_ratio-1),alpha))
                                                                           *(pow(((pressure_ratio+1)/pressure_ratio),alpha))*pow(pressure,alpha)*pow(vel,beta));
    double energy_exchange = exchange_constant*(1-vel-((avg_atomic_charge+1)/avg_atomic_charge)*pressure)/sqrt(pow(vel,5))
                             *(abundances*charge*charge*electron_mass/(atomic_mass*sqrt(pow(pressure + (1-vel-pressure)*avg_atomic_charge*electron_mass/atomic_mass,3)))).sum();
    double dpos_dvel = (pow(pre_shock_speed,3)/(shock_height*accretion_rate))*(5. - 8.*vel)/(2*energy_loss);
    double dpress_dvel = (1./(3*vel))*(2.0-8.*vel - (5.0-8.*vel)*(energy_exchange/(pre_shock_speed*pre_shock_speed*energy_loss)));
    return {dpos_dvel,dpress_dvel};
}
*/
