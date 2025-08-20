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

/*
    Copyright (c) 2014, Peter A.M. van Hoof.

    This program is provided 'as-is', without any express or implied warranty. In
    no event will the author be held liable for any damages arising from the use
    of this program.

    Permission is granted to anyone to use this program for any purpose, including
    commercial applications, and to alter it and redistribute it freely, subject
    to the following restrictions:

    1. The origin of this program must not be misrepresented; you must not claim
       that you created the original program. If you use this program in a product,
       an acknowledgment in the product documentation would be appreciated but
       is not required.
    2. Altered program versions must be plainly marked as such, and must not be
       misrepresented as being the original program.
    3. This notice may not be removed or altered from any further distribution.

    Peter A.M. van Hoof
    Royal Observatory of Belgium
    Ringlaan 3
    B-1180 Brussels
    Belgium
    p.vanhoof@oma.be
*/

/*
 * This program was retrieved from https://data.nublado.org/gauntff/ on Aug 18, 2025
 * It was modified heavily from its original form
 * The non-relativistic formula are the same with simplification for use in the relevant temperature range
 * For high temperature plasmas, linear interpolation is used over the relativsitc quantities from the same source
 */
