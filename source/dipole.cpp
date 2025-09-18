#include "dipole.hh"
#include <cmath>

void Dipole::Solve_Coordinates(double w){
    const double canalle_x = (u*u*u*u + 256*w*w/27)/(16384*pow(w,8));
    const double canalle_w = -u*u/(64*w*w*w*w);
    const double canalle_s_p = cbrt(-0.5*canalle_w + sqrt(canalle_x));
    const double canalle_s_m = cbrt(-0.5*canalle_w - sqrt(canalle_x));
    const double canalle_y = sqrt(pow(canalle_s_p+canalle_s_m,2) + 3*pow(canalle_s_p-canalle_s_m, 2));
    const double du = (4*w*w*r*r*r + u);
    const double dr_du = -r/du;
    const double dr_dw = -2*w*pow(r,4)/du;
    const double dr_du_dw = -(dr_dw + dr_du*(8*w*r*r*r + 12*w*w*r*r*dr_dw))/(4*w*w*r*r*r + u);
    const double ts = 1./u + 2*dr_du/r + dr_du_dw/dr_dw;
    const double te = 1./w + 7*dr_dw/(2*r) + dr_du_dw/dr_du;

    r = sqrt(canalle_y-canalle_s_p-canalle_s_m) - u/(4*w*w*canalle_y);
    metric[0] = sqrt(r*pow(r+3+u*dr_du,2)/(4*u) + pow(3*w*r*r*dr_du,2));
    metric[1] = sqrt(pow(r,4)*pow(r+3*w*dr_dw,2) + 2.25*u*r*dr_dw*dr_dw);
    metric[2] = sqrt(u*r*r*r);
    proj_r_w = (3*u*r - 2)*(w*pow(r,4) - r)/(metric[1]*(4*w*w*r*r*r + u));
    convergance = 3*u*pow(r,5)*(0.25*ts*dr_dw*(1. + 3*dr_du*u/r) + 3*te*pow(w*r*dr_du,2))/pow(metric[0]*metric[2],2);
}

double Dipole::Get_Radial_Distance(){
    return r;
}
double Dipole::Get_Projection(){
    return proj_r_w;
}
double Dipole::Get_Metric(int i){
    return metric[i];
}
double Dipole   ::Get_Convergence(){
    return convergance;
}
