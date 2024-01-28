#include <cmath>
#include <iostream>
#include "xsTypes.h"

using std::pow;

void synchrotron(const RealArray& E, const RealArray& params, int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	  const string& initString)
{
    Real E0 = params[0];
	fluxArray.resize(E.size()-1);
	for(int i=0; i<fluxArray.size(); i++){
		fluxArray[i] = pow(E[i], -2)*pow(1+0.46*pow(E[i]/E0,0.6),11.0/4.8)*exp(-sqrt(E[i]/E0));
	}

}
