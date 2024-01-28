// Voigt profile emission line
//  Parameters are :
//    0:   ecenter, line center (keV)
//    1:   sigma from Gaussian  (keV)
//    2:   gamma from Lorentzian  (keV)

// Calculation via the Faddeeva function w(z) = exp(-z^2) erfc(-iz)
// in XSUtil/Numerics/Faddeeva.{cxx,h} from http://ab-initio.mit.edu/Faddeeva
//    Voigt(x,sigma,gamma) = Re[w(z)] / sigma*sqrt(2pi)
// where z = (x + i*gamma)/sigma*sqrt(2) and x = E-ecenter.

#include <XSFunctions/functionMap.h>
#include <xsTypes.h>

// uses calcLine which is in calcLines.cxx
void calcLine(const RealArray& energyArray, const Real ecenter,
	      const RealArray& lineParams, const Real lineflux, const Real crtLevel,
	      const int lineShape, const bool qspeedy, RealArray& fluxArray);

void voigtLine(const RealArray& energyArray, const RealArray& params, 
	       int spectrumNumber, RealArray& fluxArray, 
	       RealArray& fluxErrArray, const string& initString)
{
  const Real crtLevel = 1.0e-8;
  RealArray lineParams(2);
  lineParams[0] = params[1];
  lineParams[1] = params[2];

  fluxArray.resize(energyArray.size()-1);
  fluxArray = 0.0;

  calcLine(energyArray, params[0], lineParams, (Real)1.0, crtLevel,
	   2, false, fluxArray);
  fluxErrArray.resize(0);

  return;
}
