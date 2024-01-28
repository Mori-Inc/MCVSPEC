#include <xsTypes.h>
#include <functionMap.h>

// prototype for call to vwDem model in vwDem.cxx.

void vwDem(const RealArray& energyArray, const RealArray& params,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString);

// XSPEC model subroutine to calculate collisional plasma with a power-law DEM
//
//  dY/dT = 0                   if   T <= beta*Tmax
//          Ytot*T^alpha        if   beta*Tmax < T < Tmax
//          0                   if   Tmax <= T
// where Y is the emission measure and Ytot is the total integrated over
// beta*Tmax to Tmax.
// 
// Parameters:
//    params[0] = Maximum temperature (Tmax)
//    params[1] = Ratio of minimum to maximum temperature (beta = Tmin/Tmax)
//    params[2] = Inverse slope (p = 1/alpha)
//    params[3] = nH (cm^-3)  Fixed at 1 for most applications
//    params[4] = abundance
//    params[5] = redshift
//    params[6] = switch(0=calculate MEKAL model, 1=interpolate MEKAL model,
//                       2=AtomDB model)

void wDem(const RealArray& energyArray, const RealArray& params,
	  int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	  const string& initString)
{

  RealArray pparams(20);

  for (int i=0; i<4; i++) pparams[i] = params[i];
  pparams[4] = 1.0;
  for (int i=5; i<18; i++) pparams[i] = params[4];
  pparams[18] = params[5];
  pparams[19] = params[6];

  vwDem(energyArray, pparams, spectrumNumber, flux, fluxErr, 
	initString);

  return;

}

