#include "XS_Cataclysmic_Variable.hh"
#include "constants.hh"

extern "C"
void IPspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    double mass = params[0]*m_sol; // WD mass [grams]
    double p_spin = params[1]; // spin period [s]
    double r_m_ratio = params[2]; // ratio between corotation radius and sping period
    double luminosity = params[3]*1e33; // luminosity [ergs/s]
    double fractional_area = params[3]; //fractional accretion area
    double col_abund = params[5]; // accretion column abundance [solar abundances]
    double cos_incl = params[6]; // cos inclination angle
    double area_exponent = params[7];
    double source_distance = params[8]*pc_to_cm; // source distnace [cm]
    int reflection_sel = params[9]; // how to apply reflection 0 = off, 1 = on

    double mag_radius = r_m_ratio*cbrt(grav_const*mass*p_spin*p_spin/(4*pi*pi)); // spin equilibrium
    double inverse_mag_radius = 1./mag_radius;
    double radius = Cataclysmic_Variable::Get_Radius(mass);
    double area = fractional_area*4*pi*radius*radius;
    double mdot = Cataclysmic_Variable::Get_Accretion_Rate(luminosity, mass, radius, inverse_mag_radius);
    double b_field = sqrt(32*mdot*sqrt(grav_const*mass/pow(inverse_mag_radius,7)))/(radius*radius*radius);

    XS_Cataclysmic_Variable intermediate_polar(mass, radius, b_field, mdot, inverse_mag_radius, col_abund, area, cos_incl, area_exponent, source_distance, reflection_sel);
    intermediate_polar.Shock_Height_Shooting();
    intermediate_polar.Build_Column_Profile();
    intermediate_polar.XS_Spectrum(energy, spectrum_num, flux, init_string);
    intermediate_polar.Print_Properties();
}

extern "C"
void IPspecArea(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    double mass = params[0]*m_sol; // WD mass [grams]
    double p_spin = params[1]; // spin period [s]
    double r_m_ratio = params[2]; // ratio between corotation radius and sping period
    double luminosity = params[3]*1e33; // luminosity [ergs/s]
    double area = params[3]*1e15; //accretion area [cm^2]
    double col_abund = params[5]; // accretion column abundance [solar abundances]
    double cos_incl = params[6]; // cos inclination angle
    double area_exponent = params[7];
    double source_distance = params[8]*pc_to_cm; // source distnace [cm]
    int reflection_sel = params[9]; // how to apply reflection 0 = off, 1 = on

    double inverse_mag_radius = 1./(r_m_ratio*cbrt(4*pi*pi*grav_const*mass*p_spin*p_spin));
    double radius = Cataclysmic_Variable::Get_Radius(mass);
    double mdot = Cataclysmic_Variable::Get_Accretion_Rate(luminosity, mass, radius, inverse_mag_radius);
    double b_field = sqrt(32*mdot*sqrt(grav_const*mass/pow(inverse_mag_radius,7)))/(radius*radius*radius);

    XS_Cataclysmic_Variable intermediate_polar(mass, radius, b_field, mdot, inverse_mag_radius, col_abund, area, cos_incl, area_exponent, source_distance, reflection_sel);
    intermediate_polar.Shock_Height_Shooting();
    intermediate_polar.Build_Column_Profile();
    intermediate_polar.XS_Spectrum(energy, spectrum_num, flux, init_string);
    intermediate_polar.Print_Properties();
}
