#include "XS_Cataclysmic_Variable.hh"
#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "xsTypes.h"
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <funcWrappers.h>
#include <cmath>
#include <iostream>

XS_Cataclysmic_Variable::XS_Cataclysmic_Variable(White_Dwarf wd, Accretion_Column col, Tolerance tol):
    Cataclysmic_Variable(wd, col, tol)
{
    Set_Abundances();
    Find_Shock_Position();
    Build_Column_Profile();
}

void XS_Cataclysmic_Variable::Set_Abundances(){
    abundances.resize(n_elements);
    abundances[0] = FunctionUtility::getAbundance(atomic_charge[0]);

    abundances[1] = FunctionUtility::getAbundance(atomic_charge[1]);
    double abund_sum=abundances[0]+abundances[1];
    for(size_t i = 2; i<n_elements; i++){
        abundances[i] = accretion_column.metallicity*FunctionUtility::getAbundance(atomic_charge[i]);
        abund_sum += abundances[i];
    }
    for(size_t i=0; i<abundances.size(); i++){
        abundances[i] /= abund_sum;
    }
    Set_Cooling_Constants();
}

const void XS_Cataclysmic_Variable::XS_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string, const bool do_refl){
    if(!valid_solution || position.size()==1){
        std::cout << "No valid solution found!" << std::endl;
        std::cout << "Column does not reach WD surface for any shock height" << std::endl;
        flux = std::nan("");
        return;
    }
    int n = flux.size();
    double alt, ion_density, cosi=white_dwarf.cos_inclination;
    double refl_amp;
    RealArray flux_integrand[2] = {RealArray(n), RealArray(n)};
    RealArray direct_flux(n);
    RealArray reflected_flux(n);
    RealArray flux_error(n);
    RealArray apec_parameters = {0,0,accretion_column.metallicity,0};
    RealArray brem_parameters = {0};
    RealArray refl_parameters = {-1,0,accretion_column.metallicity,accretion_column.metallicity,white_dwarf.cos_inclination};
    // refl_amp = -1 means only return reflected spectrum, this ensures that reflection can be done separately to apec

    if(electron_temperature[0] > 86.0 || ion_temperature[0] > 86.0){
        brem_parameters[0] = electron_temperature[0];
        CXX_bremss(energy, brem_parameters, spectrum_num, flux_integrand[1], flux_error, init_string);
        flux_integrand[1] *= 3.02e-15/1e-14;
    }
    else{
        apec_parameters[0] = electron_temperature[0];
        apec_parameters[1] = ion_temperature[0];
        CXX_tapec(energy, apec_parameters, spectrum_num, flux_integrand[1], flux_error, init_string);
    }
    ion_density = electron_density[0]/avg_atomic_charge;
    flux_integrand[1] *= volume_element[0]*electron_density[0]*ion_density*1e-14;
    flux_integrand[1] /= 4*pi*white_dwarf.distance*white_dwarf.distance;

    for(size_t i=1; i<position.size(); i++){
        flux_integrand[0] = flux_integrand[1];
        flux_integrand[1] = 0.;

        if (electron_temperature[i] > 86.0 || ion_temperature[i] > 86.0){
            brem_parameters[0] = electron_temperature[i];
            CXX_bremss(energy, brem_parameters, spectrum_num, flux_integrand[1], flux_error, init_string);
            flux_integrand[1] *= 3.02e-15/1e-14;
        }
        else{
            apec_parameters[0] = electron_temperature[i];
            apec_parameters[1] = ion_temperature[i];
            CXX_tapec(energy, apec_parameters, spectrum_num, flux_integrand[1], flux_error, init_string);
        }
        ion_density = electron_density[i]/avg_atomic_charge;
        flux_integrand[1] *= volume_element[i]*electron_density[i]*ion_density*1e-14;
        flux_integrand[1] /= 4*pi*white_dwarf.distance*white_dwarf.distance;
        direct_flux = 0.5*(position[i]-position[i-1])*(flux_integrand[1]+flux_integrand[0]);

        if(do_refl){
            alt = 0.5*(altitude[i]+altitude[i-1])/white_dwarf.radius;
            if(alt<0){
                alt=0.;
            }
            refl_amp = 1 - (sqrt(alt*(alt+2))/(1+alt))*(1 - (3*cosi*cosi - 1)/(8*(1+alt)*(1+alt)));
            reflected_flux += refl_amp*direct_flux;
        }
        flux += direct_flux;
    }
    if(do_refl){
        CXX_reflect(energy, refl_parameters, spectrum_num, reflected_flux, flux_error, init_string);
        flux += reflected_flux;
    }
}

void XS_Cataclysmic_Variable::Set_TCL() const {
    if(altitude.size() < 1){
        return;
    }
    FunctionUtility::loadDbValue("R_wd", white_dwarf.radius); // cm
    FunctionUtility::loadDbValue("B_0", white_dwarf.b_field/1e6); // MG
    if(white_dwarf.inverse_mag_radius != 0){
        FunctionUtility::loadDbValue("R_m", 1./white_dwarf.inverse_mag_radius); // cm
    }
    FunctionUtility::loadDbValue("Mdot", accretion_column.accretion_rate*accretion_column.accretion_area); // g/s
    FunctionUtility::loadDbValue("mdot", accretion_column.accretion_rate); // g/cm2/s
    FunctionUtility::loadDbValue("h_s", altitude[0]); // cm
    FunctionUtility::loadDbValue("kT_s", ion_temperature[0]); // keV
}
