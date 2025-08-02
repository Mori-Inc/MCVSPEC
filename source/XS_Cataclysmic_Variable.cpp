#include "XS_Cataclysmic_Variable.hh"

XS_Cataclysmic_Variable::XS_Cataclysmic_Variable(double m, double r, double b, double mdot, double inv_r_m, double metals, double area, double theta, double n, double dist, int reflection):
    Cataclysmic_Variable(m,r,b,mdot,inv_r_m,area,theta,n,dist,reflection)
{
    metalicity = metals;
    Set_Abundances(metals);
    Guess_Shock_Height();
    Shock_Height_Shooting();
    Build_Column_Profile();
}

void XS_Cataclysmic_Variable::Set_Abundances(double metalicity){
    abundances.resize(atomic_charge.size());
    abundances[0] = FunctionUtility::getAbundance(atomic_charge[0]);

    abundances[1] = FunctionUtility::getAbundance(atomic_charge[1]);
    for(int i = 2; i < 14; i++){
        abundances[i] = metalicity*FunctionUtility::getAbundance(atomic_charge[i]);
    }
    abundances = abundances/abundances.sum();
    Set_Cooling_Constants();
}

void XS_Cataclysmic_Variable::XS_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string){
    int n = flux.size();
    double refl_amp;
    valarray<double> apec_flux(n);
    valarray<double> reflected_flux(n);
    valarray<double> flux_error(n);
    valarray<double> apec_parameters = {0,metalicity,0};
    valarray<double> refl_parameters = {-1,0,metalicity,metalicity,incl_angle};
    // refl_amp = -1 means only return reflected spectrum, this ensures that reflection can be done separately to apec

    for(int i=0; i<altitude.size(); i++){
        apec_parameters[0] = electron_temperature[i];
        if (apec_parameters[0] > 64.0){
            CXX_bremss(energy, apec_parameters, spectrum_num, apec_flux, flux_error, init_string);
        }
        else{
            CXX_apec(energy, apec_parameters, spectrum_num, apec_flux, flux_error, init_string);
        }
        apec_flux *= volume[i]*ion_density[i]*electron_density[i]*1e-14;
        apec_flux /= 4*pi*distance*distance;
        flux += apec_flux;

        if(refl==1){
            refl_amp = 1-sqrt(1.0-1.0/pow(1+altitude[i]/radius,2));
            reflected_flux += refl_amp*apec_flux;
        }
        apec_flux *= 0;
    }
    if(refl==1){
        CXX_reflect(energy, refl_parameters, spectrum_num, reflected_flux, flux_error, init_string);
        flux += reflected_flux;
    }
}
