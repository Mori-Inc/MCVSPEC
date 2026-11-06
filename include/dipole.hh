#pragma once
#include <cmath>

struct Dipole{
    double u, w_0, a_0;

    Dipole(double uu):u(uu), w_0(sqrt(1-uu)){
        double r, proj, conv, metric[3];
        update_coordinates(w_0, r, proj, conv, metric);
        a_0 = metric[0];
    }

    double get_radial_coordinate(double w) const{
        const double w2 = w*w;

        const double _w = -u*u/(64*w2*w2);
        const double _x = (u*u*u*u + 256*w2/27)/(16384*w2*w2*w2*w2);
        const double _sp = cbrt(-0.5*_w + sqrt(_x));
        const double _sm = cbrt(-0.5*_w - sqrt(_x));
        const double ssum = _sp+_sm;
        const double sdif = _sp-_sm;
        const double _y = sqrt(ssum*ssum + 3*sdif*sdif);

        return sqrt(_y-_sp-_sm) - u/(4*w2*_y);
    }

    // solves for the relevant coordinate transforms for a dipole geometry with r=1 -> stellar surface
    void update_coordinates(double w, double& r, double& proj_r_w, double& convergance, double scale_factors[3]) const{

        r = get_radial_coordinate(w);

        const double r2 = r*r;
        const double r3 = r2*r;
        const double r4 = r3*r;

        const double costheta = w*r2;

        const double psi = 1./sqrt(1+3*costheta*costheta);

        /*  For this dipole coordiante system the scale factor h_w is equal to
        the product h_u*h_phi. We make use of this fact to make our computation
        safer by only computing r3*psi, the product. This is advantageous since
        h_u ~ 1/sin(theta) which = nan for theta=0. Since we only use the u and
        phi scale factors to compute the cross sectional area A ~ h_u*h_phi we
        simply store the area in h_u and set h_phi=1.

        The complete scale factors are:
        scale_factors[0] = r3*psi;
        scale_factors[1] = r2*psi/sintheta;
        scale_factors[2] = r*sintheta;
        */

        scale_factors[0] = r3*psi;
        scale_factors[1] = scale_factors[0];
        scale_factors[2] = 1;

        convergance = -3*w*r4*psi*psi*psi*psi*(8 - 5*u*r);
        proj_r_w = -2*costheta*psi;
    }
};
