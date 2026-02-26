#include "Cataclysmic_Variable.hh"
#include "XS_Cataclysmic_Variable.hh"
#include "constants.hh"

XS_Cataclysmic_Variable make_mcv(const RealArray& user_params, bool& refl, const bool is_ip=false, const bool use_lum=false, const bool use_f=false){
    White_Dwarf wd;
    Accretion_Column col;
    // invariant pars
    col.shock_pressure_ratio = 0.75;
    // mass and radius
    wd.mass = user_params[0]*m_sol;
    wd.radius = Mass_to_Radius(wd.mass);
    // accretion area in cm2
    int area_ind = is_ip ? 4 : 3;
    col.accretion_area = user_params[area_ind]*1e15;
    if(use_f){
        double f = user_params[area_ind];
        col.accretion_area = f*4*pi*wd.radius*wd.radius;
    }
    // magnetospheric radius
    double mag_radius = 0;
    double corotation_ratio = 0;
    wd.corotation_radius = 1;
    if(is_ip){
        double p_spin = user_params[1];
        corotation_ratio = user_params[2];
        wd.corotation_radius = cbrt(grav_const*wd.mass*p_spin*p_spin/(4*pi*pi));
        mag_radius = corotation_ratio*wd.corotation_radius;
    }
    wd.inverse_mag_radius = is_ip ? 1./mag_radius : 0;
    col.sin_mag_colat = sqrt(wd.inverse_mag_radius);
    // accretion rate
    int mdot_ind = is_ip ? 3 : 2;
    col.accretion_rate = user_params[mdot_ind];
    if(use_lum){
        col.accretion_rate = Luminosity_to_Accretion_Rate(user_params[mdot_ind]*1e33, wd)/col.accretion_area;
    }
    // magnetic field
    wd.b_field = user_params[1]*1e6;
    if(is_ip){
        wd.b_field = sqrt(32*col.accretion_rate*col.accretion_area*sqrt(grav_const*wd.mass*pow(mag_radius,7)))/(wd.radius*wd.radius*wd.radius);
    }
    // abundance, inclination angle, distnace, reflect
    int par_ind = is_ip ? 5 : 4;
    col.metallicity = user_params[par_ind];
    wd.cos_inclination = user_params[++par_ind];
    wd.distance = user_params[++par_ind]*pc_to_cm;
    refl = user_params[++par_ind];

    return XS_Cataclysmic_Variable(wd,col);
}

extern "C"
void Polarspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable polar = make_mcv(params, refl, false, true, true);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    polar.Print_Properties();
}

extern "C"
void PolarspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable polar = make_mcv(params, refl, false, true, false);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    polar.Print_Properties();
}

extern "C"
void PolarspecMdot(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable polar = make_mcv(params, refl, false, false, true);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    polar.Print_Properties();
}

extern "C"
void PolarspecMdotArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable polar = make_mcv(params, refl, false, false, false);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    polar.Print_Properties();
}

extern "C"
void IPspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable ip = make_mcv(params, refl, true, true, true);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    ip.Print_Properties();
}

extern "C"
void IPspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable ip = make_mcv(params, refl, true, true, false);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    ip.Print_Properties();
}

extern "C"
void IPspecMdot(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable ip = make_mcv(params, refl, true, false, true);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    ip.Print_Properties();
}

extern "C"
void IPspecMdotArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable ip = make_mcv(params, refl, true, false, false);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    ip.Print_Properties();
}
