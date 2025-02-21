#include "mcvspec.hh"

extern "C"
void Deqspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    mag_ratio = params[1]; // magnetospheric to wd radius ratio
    fractional_area = params[2]; //fractional accretion area
    luminosity = params[3]*1e33; // luminosity [ergs/s]
    mass = params[4]*solar_mass; // WD mass [grams]
    col_abund = params[5]; // accretion column abundance [solar abundances]
    wd_abund = params[6]; // wd surface abundance [solar abundances]
    cos_incl = params[7]; // cos inclination angle

    wd_radius = Calculate_White_Dwarf_Radius(mass);
    accretion_area = fractional_area*4*pi*wd_radius*wd_radius;
    mag_radius = wd_radius*mag_ratio;
    accretion_rate = Calculate_Accretion_Rate(mass, luminosity, wd_radius, mag_radius);
    specific_accretion = accretion_rate/accretion_area;
    b_field = Calculate_Magnetic_Field(mass, accretion_rate, wd_radius, mag_radius);
    free_fall_velocity = Calculate_Free_Fall_Velocity(mass, wd_radius, 0., mag_radius);
    shock_height = Calculate_B_Free_Shock_Height(free_fall_velocity, specific_accretion);

    free_fall_velocity = Calculate_Free_Fall_Velocity(mass, wd_radius, shock_height, mag_radius);
    shock_height = Calculate_B_Free_Shock_Height(free_fall_velocity, specific_accretion);
    shock_temperature = Calculate_Shock_Temperature(free_fall_velocity);
    epsilon_shock = Calculate_Epsilon(b_field, shock_temperature, specific_accretion, free_fall_velocity, shock_height);

    Shock_Height_Shooting(0.005, 100, true);
    MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    std::cout << "h/R_wd = " << shock_height/wd_radius << std::endl;
}
