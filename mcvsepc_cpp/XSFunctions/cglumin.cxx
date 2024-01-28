// Convolution model that can be used to create a parameter which is the 
// luminosity in a particular energy range for the input redshift
// Parameters are     energ_lo     Source frame low energy over which to 
//                                 calculate luminosity
//                    energ_hi     Source frame high energy over which to 
//                                 calculate luminosity
//                    dist         Distance to source in kpc
//                    lumin        log (base 10) Luminosity in units of erg/s

#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

#include <cmath>

void cglumin (
         const RealArray& energyArray, 
         const RealArray& params, 
         int spectrumNumber,
         RealArray& fluxArray, 
         RealArray& fluxErrArray,
         const string& initString)
{

  using namespace std;
  using namespace Numerics;

  // constant to calculate the luminosity in the correct units
  static const Real LUMCON (9.5234e42);

  Real eMin (params[0]);
  Real eMax (params[1]);
  Real dist (params[2]);
  Real lumin (pow(10.0, params[3]));

  // Integrate the input array between elow and ehi into fluxsum.

  pair<Real,Real> fluxsum (integrationKernel(energyArray,fluxArray,eMin,eMax));

  // lumin is in erg/s
  // so luminsum needs to be in erg/s
  // fluxsum.second is the erg flux, in cgs (erg/s/cm^2)
  // Luminosity is 4*pi*d^2
  // LUMCON is the conversion factor between kpc^2 and cm^2
  
  Real luminsum = fluxsum.second * 4 * M_PI * pow(dist,2) * LUMCON;

  // Scale by lumin/luminsum

  fluxArray *= lumin/luminsum;
  fluxErrArray *= lumin/luminsum;

  return;

}
