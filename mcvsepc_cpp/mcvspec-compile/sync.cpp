#include <cmath>
#include <iostream>
#include <vector>

#include "xsTypes.h"

//extern "C"
void synchrotron(const RealArray& E, const RealArray& params, int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	  const string& initString)
{
    int switchCase = params[0];
    Real E0;
    Real E1;

    int NE = E.size();
    fluxArray.resize(NE-1);
    fluxErrArray.resize(0);
    switch (switchCase)
    {
    // k = ratio of upstream and downstream B strengths
    // single-valued B field
    case 1:
        /* k = 1 */
        E0 = params[1];
        for(int i=0;i<NE;++i){
            fluxArray[i] = std::pow(E[i], -2)*std::pow(1+0.46*std::pow(E[i]/E0,0.6),11.0/4.8)*exp(-sqrt(E[i]/E0));
        }
        break;

    case 2:
        /* k = 1/√11 */
        E0 = params[1];
        for(int i=0;i<NE;++i){
            fluxArray[i] = std::pow(1+0.38*std::pow(E[i]/E0,0.5),11.0/4.0)*exp(-sqrt(E[i]/E0))/std::pow(E[i],2);
        }
        break;

    // distributed B field
    case 3:
        /* k = 1 */
        E1 = params[1];
        for(int i=0;i<NE;++i){
            fluxArray[i] = std::pow(E[i], -2)*std::pow(1+0.172*std::pow(E[i]/E1,0.46),25./(12.*0.46))*exp(-std::pow(E[i]/E1,1.0/3.0));
        }
        break;

    case 4:
        /* k = 1/√11 */
        E1 = params[1];
        for(int i=0;i<NE;i++){
            fluxArray[i] = std::pow(E[i], -2)*std::pow(1+0.185*std::pow(E[i]/E1,0.4),25./(12.*0.4))*exp(-std::pow(E[i]/E1,1.0/3.0));
        }
        break;

    default:
        std::cerr << "Invalid formula selection. Please choose an integer from 1-4." << std::endl;
        break;
    }

}