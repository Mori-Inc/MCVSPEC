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
    double mass = params[4]*m_sol; // WD mass [grams]
    double col_abund = params[5]; // accretion column abundance [solar abundances]
    double cos_incl = params[6]; // cos inclination angle
    double source_distance = params[7]*pc_to_cm; // source distnace [cm]

    Cataclysmic_Variable polar(mass, b_field, col_abund, luminosity, fractional_area, cos_incl, source_distance, reflection_sel);
    polar.Shock_Height_Shooting(1000);
    polar.MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    polar.Print_Properties();
}
