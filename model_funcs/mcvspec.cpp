#include "Cataclysmic_Variable.hh"
#include "XS_Cataclysmic_Variable.hh"
#include "constants.hh"

XS_Cataclysmic_Variable make_mcv(const RealArray& user_params, bool& refl, const bool is_ip=false, const bool use_lum=false, const bool use_f=false){
    White_Dwarf wd;
    Accretion_Column col;
    Tolerance tol;
    tol.absolute_error = 1e-8;
    tol.relative_error = 1e-6;
    tol.kT_grid_spacing = 0.5; // keV
    tol.altitude_grid_spacing = 0.1; // fraction of shock height
    // invariant pars
    col.shock_pressure_ratio = 0.75;

    size_t par_ind = 0;
    // mass and radius
    wd.mass = user_params[par_ind]*m_sol;
    wd.radius = Mass_to_Radius(wd.mass);
    // b_field
    if(is_ip){
        double p_spin = user_params[++par_ind];
        double corotation_ratio = user_params[++par_ind];
        wd.corotation_radius = cbrt(grav_const*wd.mass*p_spin*p_spin/(4*pi*pi));
        wd.inverse_mag_radius = 1./(corotation_ratio*wd.corotation_radius);
    }
    else{
        wd.b_field = user_params[++par_ind]*1e6;
        wd.corotation_radius = 1;
        wd.inverse_mag_radius = 0;
    }
    col.sin_mag_colat = sqrt(wd.inverse_mag_radius);
    // accretion rate and area
    col.accretion_rate = user_params[++par_ind];
    if(use_f){
        double f = user_params[++par_ind];
        col.accretion_area = f*4*pi*wd.radius*wd.radius;
    }
    else{
        col.accretion_area = user_params[++par_ind]*1e15;
    }
    if(use_lum){
        double lum = col.accretion_rate*1e33;
        col.accretion_rate = Luminosity_to_Accretion_Rate(lum, wd)/col.accretion_area;
    }

    if(is_ip){
        wd.b_field = sqrt(32*col.accretion_rate*col.accretion_area*sqrt(grav_const*wd.mass/pow(wd.inverse_mag_radius,7)))/(wd.radius*wd.radius*wd.radius);
    }

    col.metallicity = user_params[++par_ind];
    wd.cos_inclination = user_params[++par_ind];
    wd.distance = user_params[++par_ind]*pc_to_cm;
    refl = user_params[++par_ind];

    return XS_Cataclysmic_Variable(wd,col,tol);
}

extern "C"
void Polarspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable polar = make_mcv(params, refl, false, true, true);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    polar.Set_TCL();
}

extern "C"
void PolarspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable polar = make_mcv(params, refl, false, true, false);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    polar.Set_TCL();
}

extern "C"
void PolarspecMdot(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable polar = make_mcv(params, refl, false, false, true);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    polar.Set_TCL();
}

extern "C"
void PolarspecMdotArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable polar = make_mcv(params, refl, false, false, false);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    polar.Set_TCL();
}

extern "C"
void IPspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable ip = make_mcv(params, refl, true, true, true);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    ip.Set_TCL();
}

extern "C"
void IPspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable ip = make_mcv(params, refl, true, true, false);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    ip.Set_TCL();
}

extern "C"
void IPspecMdot(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable ip = make_mcv(params, refl, true, false, true);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    ip.Set_TCL();
}

extern "C"
void IPspecMdotArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    bool refl;
    XS_Cataclysmic_Variable ip = make_mcv(params, refl, true, false, false);
    ip.XS_Spectrum(energy, spectrum_num, flux, init_string, refl);
    ip.Set_TCL();
}
