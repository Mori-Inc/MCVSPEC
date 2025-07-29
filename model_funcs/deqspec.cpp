#include "Cataclysmic_Variable.hh"

#include <chrono>
using namespace std::chrono;

extern "C"
void Deqspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    int reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    double mag_ratio = params[1]; // magnetospheric to wd radius ratio
    double fractional_area = params[2]; //fractional accretion area
    double luminosity = params[3]*1e33; // luminosity [ergs/s]
    double mass = params[4]*m_sol; // WD mass [grams]
    double col_abund = params[5]; // accretion column abundance [solar abundances]
    double cos_incl = params[6]; // cos inclination angle
    double area_exponent = params[7];
    double source_distance = params[8]*pc_to_cm; // source distnace [cm]

    auto start = high_resolution_clock::now();
    Cataclysmic_Variable intermediate_polar(mass, col_abund, luminosity, fractional_area, cos_incl, area_exponent, source_distance, reflection_sel, 1e10);
    intermediate_polar.Set_Inverse_Mag_Radius(mag_ratio);
    std::cout << "Object Creation: " << duration_cast<microseconds>(high_resolution_clock::now()-start).count() << " ms" << std::endl;
    start = high_resolution_clock::now();
    intermediate_polar.Shock_Height_Shooting();
    std::cout << "Solution Found:  " << duration_cast<microseconds>(high_resolution_clock::now()-start).count() << " ms" << std::endl;
    start = high_resolution_clock::now();
    intermediate_polar.Build_Column_Profile();
    std::cout << "Profile Built:   " << duration_cast<microseconds>(high_resolution_clock::now()-start).count() << " ms" << std::endl;
    start = high_resolution_clock::now();
    intermediate_polar.MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    std::cout << "Spectrum Produced:  " << duration_cast<microseconds>(high_resolution_clock::now()-start).count() << " ms" << std::endl;
    intermediate_polar.Print_Properties();
}
