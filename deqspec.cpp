#include "Cataclysmic_Variable.hh"

extern "C"
void Deqspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    int reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    double mag_ratio = params[1]; // magnetospheric to wd radius ratio
    double fractional_area = params[2]; //fractional accretion area
    double luminosity = params[3]*1e33; // luminosity [ergs/s]
    double mass = params[4]*solar_mass; // WD mass [grams]
    double col_abund = params[5]; // accretion column abundance [solar abundances]
    double cos_incl = params[6]; // cos inclination angle
    double source_distance = params[7]*pc_to_cm; // source distnace [cm]

    Cataclysmic_Variable intermediate_polar(mass, col_abund, luminosity, fractional_area, cos_incl, source_distance, reflection_sel, 1e10);
    intermediate_polar.Set_Inverse_Mag_Radius(mag_ratio);
    intermediate_polar.Shock_Height_Shooting(1000);
    intermediate_polar.MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    intermediate_polar.Print_Properties();
}
