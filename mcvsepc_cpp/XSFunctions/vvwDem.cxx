#include <xsTypes.h>
#include <functionMap.h>
#include <cmath>


// function from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
                        const IntegerVector& Zarray, const RealArray& abun, 
                        const Real dens, const Real z, const RealArray& Tarr, 
                        const RealArray& DEMarr, const int ifl, const bool qtherm, 
                        const Real velocity, RealArray& fluxArray, 
                        RealArray& fluxErrArray);


// function from vwDem.cxx
void calculatePowerLawDEM(const Real Tmin, const Real Tmax, 
			  const Real invSlope, RealArray& Tarray, 
			  RealArray& DEMarray);



// XSPEC model subroutine to calculate collisional plasma with a power-law DEM
//
//  dY/dT = 0                   if   T <= beta*Tmax
//          Ytot*T^alpha        if   beta*Tmax < T < Tmax
//          0                   if   Tmax <= T
// where Y is the emission measure and Ytot is the total integrated over
// beta*Tmax to Tmax.
// 
// Parameters:
//       0 = Maximum temperature (Tmax)
//       1 = Ratio of minimum to maximum temperature (beta = Tmin/Tmax)
//       2 = Inverse slope (p = 1/alpha)
//       3 = nH (cm^-3)  Fixed at 1 for most applications
//   4..33 = abundances
//      34= redshift
//      35 = switch(0=calculate MEKAL model, 1=interpolate MEKAL model,
//                       2=AtomDB model)


void vvwDem(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	    const string& initString)
{
   const Real Tmax = params[0];
   Real Tmin = params[1]*params[0];
   if ( Tmin < 0.0 ) Tmin = 0.0;
   const Real invSlope = params[2];

   // set up the T and DEM arrays

   RealArray Tarray, DEMarray;
   calculatePowerLawDEM(Tmin, Tmax, invSlope, Tarray, DEMarray);

   // set up all the variables to pass to calcMultiTempPlasma

   int swtch = static_cast<int>(params[35]);
   int plasmaType(6);
   if ( swtch == 0 ) {
     plasmaType = 3;
   } else if ( swtch == 1 ) {
     plasmaType = 4;
   } else if ( swtch == 2 ) {
     plasmaType = 6;
   }

   const Real density = params[3];
   const Real redshift = params[34];

   RealArray abun(30);
   IntegerVector Zarray(30);
   for (size_t i=0; i<abun.size(); i++) {
     abun[i] = params[i+4];
     Zarray[i] = i+1;
   }

   const bool qtherm = false;
   const Real velocity = 0.0;

   calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, density,
		       redshift, Tarray, DEMarray, spectrumNumber, 
		       qtherm, velocity, flux, fluxErr);

}
