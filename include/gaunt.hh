#pragma once

namespace gaunt{
    // coefficients for exact non-relativistic gaunt factor
    const double a_high[5] = {1.43251926625281e+00, 3.50626935257777e-01, 4.36183448595035e-01, 6.03536387105599e-02, 3.66626405363100e-02};
    const double b_high[5] = {1.00000000000000e+00, 2.92525161994346e-01, 4.05566949766954e-01, 5.62573012783879e-02, 3.33019373823972e-02};
    const double a_low[5] = {1.45481634667278e+00, -9.55399384620923e-02, 1.46327814151538e-01, -1.41489406498468e-02, 2.76891413242655e-03};
    const double b_low[5] = {1.00000000000000e+00, 3.31149751183539e-02, 1.31127367293310e-01, -1.32658217746618e-02, 2.74809263365693e-03};
    // interpolation grid for relativistic gaunt factor
    const double gamma_interp[26] = {-4.3, -4.2, -4.1, -4., -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1,
                                    -3., -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2. , -1.9, -1.8};
    const double gaunt_interp[26] = {1.97078, 1.74022, 1.5677 , 1.44298, 1.35169, 1.2879, 1.24188, 1.21068, 1.18857, 1.1742 , 1.16452,
                                     1.15894, 1.15595, 1.15531, 1.15629, 1.1588 , 1.16247, 1.16731, 1.17318, 1.18009, 1.188 ,1.19693,
                                     1.20686, 1.21782, 1.22984, 1.24287};

    double gaunt_factor(double);
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
