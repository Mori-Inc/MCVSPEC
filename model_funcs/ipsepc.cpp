#include "Cataclysmic_Variable.hh"

#include <chrono>
using namespace std::chrono;

extern "C"
void IPspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    int reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    double p_spin = params[1]; // spin period [s]; assumes spin equilibrium
    double r_m_ratio = params[2]; // ratio between corotation radius and sping period
    double fractional_area = params[3]; //fractional accretion area
    double luminosity = params[4]*1e33; // luminosity [ergs/s]
    double mass = params[5]*m_sol; // WD mass [grams]
    double col_abund = params[6]; // accretion column abundance [solar abundances]
    double cos_incl = params[7]; // cos inclination angle
    double area_exponent = params[8];
    double source_distance = params[9]*pc_to_cm; // source distnace [cm]

    double mag_radius = r_m_ratio*cbrt(4*pi*pi*grav_const*mass*p_spin*p_spin);

    auto start = high_resolution_clock::now();
    Cataclysmic_Variable intermediate_polar(mass, col_abund, luminosity, fractional_area, cos_incl, area_exponent, source_distance, reflection_sel, mag_radius);
    std::cout << "Object Creation: " << duration_cast<milliseconds>(high_resolution_clock::now()-start).count() << " ms" << std::endl;
    start = high_resolution_clock::now();
    intermediate_polar.Shock_Height_Shooting();
    std::cout << "Solution Found:  " << duration_cast<milliseconds>(high_resolution_clock::now()-start).count() << " ms" << std::endl;
    start = high_resolution_clock::now();
    intermediate_polar.Build_Column_Profile();
    std::cout << "Profile Built:   " << duration_cast<milliseconds>(high_resolution_clock::now()-start).count() << " ms" << std::endl;
    start = high_resolution_clock::now();
    intermediate_polar.MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    std::cout << "Spectrum Produced:  " << duration_cast<milliseconds>(high_resolution_clock::now()-start).count() << " ms" << std::endl;
    intermediate_polar.Print_Properties();
}
