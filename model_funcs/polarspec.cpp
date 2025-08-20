#include "XS_Cataclysmic_Variable.hh"

extern "C"
void Polarspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    double mass = params[0]*m_sol; // WD mass [grams]
    double b_field = params[1]*1e6; // surface b field [gauss]
    double luminosity = params[2]*1e33; // luminosity [ergs/s]
    double fractional_area = params[3]; //fractional accretion area
    double col_abund = params[4]; // accretion column abundance [solar abundances]
    double cos_incl = params[5]; // cos inclination angle
    double area_exponent = params[6];
    double source_distance = params[7]*pc_to_cm; // source distnace [cm]
    int reflection_sel = params[8]; // how to apply reflection 0 = off, 1 = on

    double radius = Cataclysmic_Variable::Get_Radius(mass);
    double area = fractional_area*4*pi*radius*radius;
    double mdot = Cataclysmic_Variable::Get_Accretion_Rate(luminosity, mass, radius, 0);
    XS_Cataclysmic_Variable polar(mass, radius, b_field, mdot, 0, 0, col_abund, area, cos_incl, area_exponent, source_distance, reflection_sel);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}

extern "C"
void PolarspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    double mass = params[0]*m_sol; // WD mass [grams]
    double b_field = params[1]*1e6; // surface b field [gauss]
    double luminosity = params[2]*1e33; // luminosity [ergs/s]
    double area = params[3]*1e15; //fractional accretion area
    double col_abund = params[4]; // accretion column abundance [solar abundances]
    double cos_incl = params[5]; // cos inclination angle
    double area_exponent = params[6];
    double source_distance = params[7]*pc_to_cm; // source distnace [cm]
    int reflection_sel = params[8]; // how to apply reflection 0 = off, 1 = on

    double radius = Cataclysmic_Variable::Get_Radius(mass);
    double mdot = Cataclysmic_Variable::Get_Accretion_Rate(luminosity, mass, radius, 0);
    XS_Cataclysmic_Variable polar(mass, radius, b_field, mdot, 0, 0, col_abund, area, cos_incl, area_exponent, source_distance, reflection_sel);
    polar.XS_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}
