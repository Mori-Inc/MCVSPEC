#include "constants.hh"
#include "gaunt.hh"

double gaunt::gaunt_factor(double kT){
    double gamma_sqr = ryd_to_erg/kT;
    double log_gam = log10(gamma_sqr);
    double p,q;
    if(kT > 1./erg_to_kev){
        int i = (log_gam-gamma_interp[0])/0.1;
        return ((gaunt_interp[i+1]-gaunt_interp[i])/(gamma_interp[i+1]-gamma_interp[i]))*(log_gam-gamma_interp[i]) + gaunt_interp[i];
    }
    else if(log_gam<=0.8){
        p = a_high[0] + log_gam*(a_high[1] + log_gam*(a_high[2] + log_gam*(a_high[3] + log_gam*a_high[4])));
        q = b_high[0] + log_gam*(b_high[1] + log_gam*(b_high[2] + log_gam*(b_high[3] + log_gam*b_high[4])));
    }
    else{
        p = a_low[0] + log_gam*(a_low[1] + log_gam*(a_low[2] + log_gam*(a_low[3] + log_gam*a_low[4])));
        q = b_low[0] + log_gam*(b_low[1] + log_gam*(b_low[2] + log_gam*(b_low[3] + log_gam*b_low[4])));
    }
    return p/q;
}
