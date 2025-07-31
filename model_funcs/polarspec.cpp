#include "Cataclysmic_Variable.hh"

extern "C"
void Polarspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    int reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on
    double b_field = params[1]*1e6; // surface b field [gauss]
    double fractional_area = params[2]; //fractional accretion area
    double luminosity = params[3]*1e33; // luminosity [ergs/s]
    double mass = params[4]*m_sol; // WD mass [grams]
    double col_abund = params[5]; // accretion column abundance [solar abundances]
    double cos_incl = params[6]; // cos inclination angle
    double area_exponent = params[7];
    double source_distance = params[8]*pc_to_cm; // source distnace [cm]

    double radius = Cataclysmic_Variable::Get_Radius(mass);
    double area = fractional_area*4*pi*radius*radius;
    double mdot = Cataclysmic_Variable::Get_Accretion_Rate(luminosity, mass, radius, 0);
    Cataclysmic_Variable polar(mass, radius, b_field, mdot, 0, col_abund, area, cos_incl, area_exponent, source_distance, reflection_sel);
    polar.Shock_Height_Shooting();
    polar.Build_Column_Profile();
    polar.MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}

extern "C"
void PolarspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    int reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on
    double b_field = params[1]*1e6; // surface b field [gauss]
    double area = params[2]*1e15; // accretion area [cm^2]
    double luminosity = params[3]*1e33; // luminosity [ergs/s]
    double mass = params[4]*m_sol; // WD mass [grams]
    double col_abund = params[5]; // accretion column abundance [solar abundances]
    double cos_incl = params[6]; // cos inclination angle
    double area_exponent = params[7];
    double source_distance = params[8]*pc_to_cm; // source distnace [cm]

    double radius = Cataclysmic_Variable::Get_Radius(mass);
    double mdot = Cataclysmic_Variable::Get_Accretion_Rate(luminosity, mass, radius, 0);
    Cataclysmic_Variable polar(mass, radius, b_field, mdot, 0, col_abund, area, cos_incl, area_exponent, source_distance, reflection_sel);
    polar.Shock_Height_Shooting();
    polar.Build_Column_Profile();
    polar.MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}
