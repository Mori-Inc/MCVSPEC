#include "mcvspec.hh"

valarray<double> Normalized_Position_Derivative(double vel, valarray<double> pos, void* eps_s){
    double epsilon_s = *(double *)eps_s;
    double prefactor = pow(free_fall_velocity,3)/(2.*shock_height*bremss_const*specific_accretion);
    double cyclotron = 1/(1+epsilon_s*pow(4.,alpha+beta)*pow(3.,-alpha)*pow(1-vel,alpha)*pow(vel,beta));
    return {prefactor*vel*vel*(5.0-8.0*vel)/sqrt(vel*(1-vel))*cyclotron};
}
