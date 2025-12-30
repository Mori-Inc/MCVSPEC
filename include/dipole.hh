#pragma once
#include<valarray>
#include <cmath>

using std::valarray;

struct Dipole{
    double u, w_0, a_0, dr_dw_0;
    Dipole(double uu):u(uu), w_0(sqrt(1-uu)){
        double r, proj, conv, metric[3];
        update_coordinates(w_0, r, dr_dw_0, proj, conv, metric);
        a_0 = metric[0]*metric[1];
    }

    void set_bounds(double& upper_bound, double& lower_bound){
        upper_bound = w_0 + 1e-6/dr_dw_0;
        lower_bound = w_0/4;
    }

    // solves for the relevant coordinate transforms for a dipole geometry with r=1 -> stellar surface
    void update_coordinates(double w, double& r, double& dr_dw, double& proj_r_w, double& convergance, double metric[3]) const{
        const double w2 = w*w;
        const double w4 =  w2*w2;
        const double w8 = w4*w4;

        const double canalle_w = -u*u/(64*w4);
        const double canalle_x = (u*u*u*u + 256*w2/27)/(16384*w8);
        const double canalle_sp = cbrt(-0.5*canalle_w + sqrt(canalle_x));
        const double canalle_sm = cbrt(-0.5*canalle_w - sqrt(canalle_x));
        const double ssum = canalle_sp+canalle_sm;
        const double sdif = canalle_sp-canalle_sm;
        const double canalle_y = sqrt(ssum*ssum + 3*sdif*sdif);

        r = sqrt(canalle_y-canalle_sp-canalle_sm) - u/(4*w2*canalle_y);
        const double r2 = r*r;
        const double r3 = r2*r;
        const double r4 = r2*r2;

        const double du = 4*w2*r3 + u;
        const double dr_du = -r/du;
        dr_dw = -2*w*r4/du;
        const double dr_du_dw = -(dr_dw + dr_du*(8*w*r3 + 12*w2*r2*dr_dw))/du;

        const double a = r+3*u*dr_du;
        const double b = 3*w*r2*dr_du;
        const double c = r+3*w*dr_dw;
        metric[0] = sqrt((r/(4*u))*a*a + b*b);
        metric[1] = sqrt(r4*c*c + 2.25*u*r*dr_dw*dr_dw);
        metric[2] = sqrt(u*r3);

        proj_r_w = dr_dw/metric[1];

        const double ts = 1./u + 2*dr_du/r + dr_du_dw/dr_dw;
        const double te = 1./w + 3.5*dr_dw/r + dr_du_dw/dr_du;
        convergance = (3*u*r*r4/(metric[0]*metric[0]*metric[2]*metric[2]));
        convergance *= (dr_dw/4)*(1 + 3*dr_du*u/r)*ts + 3*(w2*r2*dr_du*dr_du)*te;
    }
};
