#pragma once
#include<valarray>
#include <cmath>

using std::valarray;

struct Dipole{
    double u;
    Dipole(double uu):u(uu){}

    // solves for the relevant coordinate transforms for a dipole geometry with r=1 -> stellar surface
    void update_coordinates(double w, double& r, double& proj_r_w, double& convergance, double metric[3]) const{
        const double canalle_x = (u*u*u*u + 256*w*w/27)/(16384*pow(w,8));
        const double canalle_w = -u*u/(64*w*w*w*w);
        const double canalle_s_p = cbrt(-0.5*canalle_w + sqrt(canalle_x));
        const double canalle_s_m = cbrt(-0.5*canalle_w - sqrt(canalle_x));
        const double ssum = canalle_s_p+canalle_s_m;
        const double sdif = canalle_s_p-canalle_s_m;
        const double canalle_y = sqrt(ssum*ssum + 3*sdif*sdif);

        r = sqrt(canalle_y-canalle_s_p-canalle_s_m) - u/(4*w*w*canalle_y);

        const double du = (4*w*w*r*r*r + u);
        const double dr_du = -r/du;
        const double r4 = r*r*r*r;
        const double dr_dw = -2*w*r4/du;
        const double dr_du_dw = -(dr_dw + dr_du*(8*w*r*r*r + 12*w*w*r*r*dr_dw))/(4*w*w*r*r*r + u);
        const double ts = 1./u + 2*dr_du/r + dr_du_dw/dr_dw;
        const double te = 1./w + 7*dr_dw/(2*r) + dr_du_dw/dr_du;

        const double a = r+3+u*dr_du;
        const double b = w*r*dr_du;
        metric[0] = sqrt(r*a*a/(4*u) + 9*r*r*b*b);
        const double c = r+3*w*dr_dw;
        metric[1] = sqrt(r4*c*c + 2.25*u*r*dr_dw*dr_dw);
        metric[2] = sqrt(u*r*r*r);
        proj_r_w = (3*u*r - 2)*(w*r4 - r)/(metric[1]*(4*w*w*r*r*r + u));
        convergance = 3*u*r*r4*(0.25*ts*dr_dw*(1. + 3*dr_du*u/r) + 3*te*b*b)/(metric[0]*metric[0]*metric[2]*metric[2]);
    }
};
