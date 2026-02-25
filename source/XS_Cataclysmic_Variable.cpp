#include "XS_Cataclysmic_Variable.hh"
#include "constants.hh"
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <funcWrappers.h>
#include <cmath>
#include <iostream>

XS_Cataclysmic_Variable::XS_Cataclysmic_Variable(double m, double r, double b, double mdot, double area, double inv_r_m, double r_m_ratio, double metals, double theta, double dist, int reflection):
    Cataclysmic_Variable(m,r,b,mdot,area,inv_r_m,r_m_ratio,metals,theta,0.75,1e-8,dist,reflection)
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
        abundances[i] = metallicity*FunctionUtility::getAbundance(atomic_charge[i]);
        abund_sum += abundances[i];
    }
    for(size_t i=0; i<abundances.size(); i++){
        abundances[i] /= abund_sum;
    }
    Set_Cooling_Constants();
}

const void XS_Cataclysmic_Variable::XS_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string){
    if(!valid_solution){
        std::cout << "No valid solution found!" << std::endl;
        std::cout << "Column does not reach WD surface for any shock height" << std::endl;
        flux = std::nan("");
        return;
    }
    int n = flux.size();
    double alt;
    double refl_amp;
    RealArray apec_flux(n);
    RealArray reflected_flux(n);
    RealArray flux_error(n);
    RealArray apec_parameters = {0,0,metallicity,0};
    RealArray refl_parameters = {-1,0,metallicity,metallicity,incl_angle};
    // refl_amp = -1 means only return reflected spectrum, this ensures that reflection can be done separately to apec

    for(size_t i=0; i<altitude.size(); i++){
        apec_parameters[0] = electron_temperature[i];
        apec_parameters[1] = ion_temperature[i];
        if (electron_temperature[i] > 64.0 || ion_temperature[i] > 64.0){
            CXX_bremss(energy, apec_parameters, spectrum_num, apec_flux, flux_error, init_string);
        }
        else{
            CXX_tapec(energy, apec_parameters, spectrum_num, apec_flux, flux_error, init_string);
        }
        apec_flux *= volume[i]*electron_density[i]*(electron_density[i]/avg_atomic_charge)*1e-14;
        apec_flux /= 4*pi*distance*distance;
        flux += apec_flux;

        if(refl==1){
            alt = altitude[i]<0 ? 0. : altitude[i];
            refl_amp = 1-sqrt(1.0-1.0/pow(1+alt/radius,2));
            reflected_flux += refl_amp*apec_flux;
        }
        apec_flux *= 0;
    }
    if(refl==1){
        CXX_reflect(energy, refl_parameters, spectrum_num, reflected_flux, flux_error, init_string);
        flux += reflected_flux;
    }
}
