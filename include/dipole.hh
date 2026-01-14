#pragma once
#include<valarray>
#include <cmath>

using std::valarray;

struct Dipole{
    double u, w_0, a_0;

    Dipole(double uu):u(uu), w_0(sqrt(1-uu)){
        double r, proj, conv, metric[3];
        update_coordinates(w_0, r, proj, conv, metric);
        a_0 = metric[0]*metric[2];
    }

    // solves for the relevant coordinate transforms for a dipole geometry with r=1 -> stellar surface
    void update_coordinates(double w, double& r, double& proj_r_w, double& convergance, double metric[3]) const{
        const double w2 = w*w;

        const double canalle_w = -u*u/(64*w2*w2);
        const double canalle_x = (u*u*u*u + 256*w2/27)/(16384*w2*w2*w2*w2);
        const double canalle_sp = cbrt(-0.5*canalle_w + sqrt(canalle_x));
        const double canalle_sm = cbrt(-0.5*canalle_w - sqrt(canalle_x));
        const double ssum = canalle_sp+canalle_sm;
        const double sdif = canalle_sp-canalle_sm;
        const double canalle_y = sqrt(ssum*ssum + 3*sdif*sdif);

        r = sqrt(canalle_y-canalle_sp-canalle_sm) - u/(4*w2*canalle_y);

        const double r2 = r*r;
        const double r3 = r2*r;
        const double r4 = r3*r;

        const double denom = 1 + 3*w2*r4;
        metric[0] = r2/sqrt(r*u*denom);
        metric[2] = sqrt(u*r3);
        metric[1] = metric[0]*metric[2];

        convergance = -(3*w*r4/(denom*denom))*(5*w2*r4 + 3);

        const double dr_dw = -2*w*r4/(4*w2*r3 + u);
        proj_r_w = (r/(2*metric[1]))*(2*w*r4 + (3*u + 6*w2*r3)*dr_dw);
    }
};
