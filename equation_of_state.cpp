#include "equation_of_state.hh"

valarray<double> Chandrasekhar_White_Dwarf_Equation(double eta, valarray<double> phi_psi, void* y0){
    double y_0 = *(double *)y0;
    return {phi_psi[1],(-2./eta)*phi_psi[1] - sqrt(pow(phi_psi[0]*phi_psi[0] - 1./(y_0*y_0), 3))};
}
