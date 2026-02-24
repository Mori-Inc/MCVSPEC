#include "Cataclysmic_Variable.hh"
#include "XS_Cataclysmic_Variable.hh"
#include "constants.hh"

XS_Cataclysmic_Variable make_mcv(const RealArray& user_params, const bool is_ip=false, const bool use_lum=false, const bool use_f=false){
    // mass and radius
    double mass = user_params[0]*m_sol;
    double radius = Mass_to_Radius(mass);
    // accretion area in cm2
    int area_ind = is_ip ? 4 : 3;
    double area = user_params[area_ind]*1e15;
    if(use_f){
        double f = user_params[area_ind];
        area = f*4*pi*radius*radius;
    }
    // magnetospheric radius
    double mag_radius = 0;
    double corotation_ratio = 0;
    if(is_ip){
        double p_spin = user_params[1];
        corotation_ratio = user_params[2];
        mag_radius = corotation_ratio*cbrt(grav_const*mass*p_spin*p_spin/(4*pi*pi));
    }
    double inverse_rm = is_ip ? 1./mag_radius : 0;
    // accretion rate
    int mdot_ind = is_ip ? 3 : 2;
    double mdot = area*user_params[mdot_ind];
    if(use_lum){
        mdot = Luminosity_to_Accretion_Rate(user_params[mdot_ind]*1e33, mass, radius, inverse_rm);
    }
    // magnetic field
    double b_field = user_params[1]*1e6;
    if(is_ip){
        b_field = sqrt(32*mdot*sqrt(grav_const*mass*pow(mag_radius,7)))/(radius*radius*radius);
    }
    // abundance, inclination angle, distnace, reflect
    int par_ind = is_ip ? 5 : 4;
    double metallicity = user_params[par_ind];
    double cos_incl = user_params[++par_ind];
    double distance = user_params[++par_ind]*pc_to_cm;
    int refl = user_params[++par_ind];

    return XS_Cataclysmic_Variable(mass, radius, b_field, mdot, area, inverse_rm, corotation_ratio, metallicity, cos_incl, distance, refl);
}

extern "C"
void Polarspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    XS_Cataclysmic_Variable polar = make_mcv(params, false, true, true);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}

extern "C"
void PolarspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    XS_Cataclysmic_Variable polar = make_mcv(params, false, true, false);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}

extern "C"
void PolarspecMdot(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    XS_Cataclysmic_Variable polar = make_mcv(params, false, false, true);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}

extern "C"
void PolarspecMdotArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    XS_Cataclysmic_Variable polar = make_mcv(params, false, false, false);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}

extern "C"
void IPspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    XS_Cataclysmic_Variable ip = make_mcv(params, true, true, true);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string);
    ip.Print_Properties();
}

extern "C"
void IPspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    XS_Cataclysmic_Variable ip = make_mcv(params, true, true, false);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string);
    ip.Print_Properties();
}

extern "C"
void IPspecMdot(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    XS_Cataclysmic_Variable ip = make_mcv(params, true, false, true);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string);
    ip.Print_Properties();
}

extern "C"
void IPspecMdotArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    XS_Cataclysmic_Variable ip = make_mcv(params, true, false, false);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string);
    ip.Print_Properties();
}
