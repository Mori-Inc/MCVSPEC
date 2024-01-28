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
//    params(0) = Maximum temperature (Tmax)
//    params(1) = Ratio of minimum to maximum temperature (beta = Tmin/Tmax)
//    params(2) = Inverse slope (p = 1/alpha)
//    params(3) = nH (cm^-3)  Fixed at 1 for most applications
//    params(4) = He abundance
//    params(5) = C   "
//    params(6) = N   "
//    params(7) = O   "
//    params(8) = Ne  "
//    params(9) = Na  "
//    params(10)= Mg  "
//    params(11)= Al  "
//    params(12)= Si  "
//    params(13)= S   "
//    params(14)= Ar  "
//    params(15)= Ca  "
//    params(16)= Fe  "
//    params(17)= Ni  " 
//    params(18)= redshift
//    params(19) = switch(0=calculate MEKAL model, 1=interpolate MEKAL model,
//                       2=AtomDB model)


void vwDem(const RealArray& energyArray, const RealArray& params,
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

   int swtch = static_cast<int>(params[19]);
   int plasmaType(6);
   if ( swtch == 0 ) {
     plasmaType = 3;
   } else if ( swtch == 1 ) {
     plasmaType = 4;
   } else if ( swtch == 2 ) {
     plasmaType = 6;
   }

   const Real density = params[3];
   const Real redshift = params[18];

   RealArray abun(14);
   for (size_t i=0; i<abun.size(); i++) abun[i] = params[i+4];
   const int elements[] = {2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};
   IntegerVector Zarray(14);
   for (size_t i=0; i<14; i++) Zarray[i] = elements[i];

   const bool qtherm = false;
   const Real velocity = 0.0;

   calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, density,
		       redshift, Tarray, DEMarray, spectrumNumber, 
		       qtherm, velocity, flux, fluxErr);

}


void calculatePowerLawDEM(const Real Tmin, const Real Tmax, 
			  const Real invSlope, RealArray& Tarray, 
			  RealArray& DEMarray)
{

  // *******************************************************************
  // set up arrays of temperature and DEM values
  // use nT temperature steps running linearly from Tmin to Tmax

  const int nT(51);

  Tarray.resize(nT);
  DEMarray.resize(nT);

  Real Tdelta = (Tmax - Tmin) / (nT-1);
  for (int i=0; i<nT; i++) Tarray[i] = Tmin + Tdelta * (i + 0.5);

   // special case of p = 0 which places all the EM in the highest temperature
   // bin

   if ( invSlope == 0 ) {

     for (int i=0; i<nT-1; i++) DEMarray[i] = 0.0;
     DEMarray[nT-1] = 1.0;

   } else {

     Real slope = 1.0/invSlope;
     
     for (int i=0; i<nT; i++) {
       Real T1 = Tarray[i] - Tdelta/2.0;
       Real T2 = Tarray[i] + Tdelta/2.0;
       if ( slope == -1.0 ) {
	 DEMarray[i] = log(T2) - log(T1);
       } else {
	 DEMarray[i] = (pow(T2,slope+1.0)-pow(T1,slope+1.0))/(slope+1.0);
       }
     }

     // normalize to total emissivity
     Real demsum(0.0);
     for (int i=0; i<nT; i++) demsum += DEMarray[i];
     DEMarray /= demsum;

   }

   // end of set up arrays of temperature and DEM values
   // *******************************************************************

   return;
}
