#include "mcvspec.hh"

extern "C"
void IPspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    p_spin = params[1]; // spin period [s]; assumes spin equilibrium
    fractional_area = params[2]; //fractional accretion area
    luminosity = params[3]*1e33; // luminosity [ergs/s]
    mass = params[4]*solar_mass; // WD mass [grams]
    col_abund = params[5]; // accretion column abundance [solar abundances]
    cos_incl = params[7]; // cos inclination angle
    source_distance = params[8]*pc_to_cm; // source distnace [cm]

    Set_Abundances(col_abund);
    wd_radius = Calculate_White_Dwarf_Radius(mass);
    accretion_area = fractional_area*4*pi*wd_radius*wd_radius;
    mag_radius = cbrt(4*pi*pi*grav_const*mass*p_spin*p_spin);
    accretion_rate = Calculate_Accretion_Rate(mass, luminosity, wd_radius, mag_radius);
    specific_accretion = accretion_rate/accretion_area;
    b_field = Calculate_Magnetic_Field(mass, accretion_rate, wd_radius, mag_radius);
    free_fall_velocity = Calculate_Free_Fall_Velocity(mass, wd_radius, 0., mag_radius);
    shock_height = Calculate_B_Free_Shock_Height(free_fall_velocity, specific_accretion);

    free_fall_velocity = Calculate_Free_Fall_Velocity(mass, wd_radius, shock_height, mag_radius);
    shock_height = Calculate_B_Free_Shock_Height(free_fall_velocity, specific_accretion);
    epsilon_shock = Calculate_Epsilon(free_fall_velocity);

    Shock_Height_Shooting(absolute_err, 100, true);
    MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
    std::cout << "h/R_wd = " << shock_height/wd_radius << std::endl;
}
