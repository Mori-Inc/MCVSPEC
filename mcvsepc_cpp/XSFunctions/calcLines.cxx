// Generic code for calculating lines with various shapes

#define GAUSS 0
#define LORENTZ 1
#define VOIGT 2

#define SQRT2 1.4142135623730950488

#include <xsTypes.h>
#include <XSstreams.h>
#include <XSUtil/Numerics/AdaptiveIntegrate.h>
#include <XSUtil/Numerics/BinarySearch.h>
#include <XSUtil/Numerics/Faddeeva.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <cmath>
#include <sstream>
#include <iostream>

using std::atan;
using namespace Numerics;


void calcManyLines(const RealArray& energyArray, const RealArray& ecenterArray,
		   const std::vector<RealArray>& lineParamsArray,
		   const RealArray& linefluxArray, const Real crtLevel,
		   const int lineShape, const bool qspeedy, RealArray& fluxArray);
void calcLine(const RealArray& energyArray, const Real ecenter, 
	      const RealArray& lineParams, const Real lineflux, const Real crtLevel, 
	      const int lineShape, const bool qspeedy, RealArray& fluxArray);

Real lineFraction(const int lineShape, const Real energy, const Real ecenter, 
		  const RealArray& lineParams, const bool qspeedy);
Real gaussFraction(const Real deltasigma, const bool qspeedy);
Real lorentzFraction(const Real deltasigma, const bool qspeedy);
Real voigtFraction(const Real energy, const Real ecenter, const Real sigma, const Real gamma);
RealArray voigtLineIntegrand(const RealArray& x, void *p);


// to calculate line shapes. If qspeedy is true then uses a tabulation
// otherwise does a proper calculation. Lines are calculated down to
// crtlevel*(flux in center bin)
// Arguments:
//      energyArray                 model energy bins
//      ecenterArray                array of line center energies
//      lineParamsArray             vector of line parameter arrays
//      linefluxArray               array of line fluxes
//      crtLevelIn                  critical level down which to calculate lines
//                                    can be overridden by xset LINECRITLEVEL
//      lineShape                   0==Gauss, 1==Lorentz, 2==Voigt
//      qspeedy                     if true then do fast but less accurate calculation
//      fluxArray                   input/output fluxes - note that this routine
//                                    does not initialize input flux to zero but uses +=

void calcManyLines(const RealArray& energyArray, const RealArray& ecenterArray,
		   const std::vector<RealArray>& lineParamsArray,
		   const RealArray& linefluxArray, const Real crtLevelIn, 
		   const int inLineShape, const bool qspeedy, RealArray& fluxArray)
{
  // Find out whether we want to override the crtLevel

  Real crtLevel(crtLevelIn);
  string pname = "LINECRITLEVEL";
  string pvalue = FunctionUtility::getModelString(pname);
  if ( pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    std::istringstream iss(pvalue);
    Real tmpVal(crtLevel);
    if (!(iss >> tmpVal) || !iss.eof()) {
      std::ostringstream oss;
      oss << "Failed to read value from LINECRITLEVEL - assuming "
	  << crtLevel << "\n";
      FunctionUtility::xsWrite(oss.str(),10);
    } else {
      crtLevel = tmpVal;
    }
  }
  
  // Find the bins containing the line centers. This assumes that the ecenter
  // array is in increasing order of energy. If the center energy falls outside
  // the energy array then lineCenterBin will be set to < 0. Otherwise
  // energyArray[lineCenterBin] < ecenterArray <= energyArray[lineCenterBin+1]

  IntegerVector lineCenterBin(ecenterArray.size());
  lineCenterBin = BinarySearch(energyArray, ecenterArray);

  Real efirst = energyArray[0];
  int nE = energyArray.size();
  Real elast = energyArray[nE-1];

  // useful to have a non-const lineShape variable so we can change it
  // for Voigt special cases
  int lineShape = inLineShape;

  // Loop around the lines ignoring those with zero line flux

  for (size_t iline=0; iline<ecenterArray.size(); iline++) {

    Real lineflux = linefluxArray[iline];
    if ( lineflux == 0.0 ) continue;

    int icen = lineCenterBin[iline];
    Real ecenter = ecenterArray[iline];

    // handle the special cases for Voigt where we can use Gauss or Lorentz

    RealArray lineParams;
    if ( (lineShape == VOIGT) && (lineParamsArray[iline][0] <= 0.0) ) {
      lineShape = LORENTZ;
      lineParams.resize(1);
      lineParams[0] = lineParamsArray[iline][1];
    } else if ( (lineShape == VOIGT) && (lineParamsArray[iline][1] <= 0.0) ) {
      lineShape = GAUSS;
      lineParams.resize(1);
      lineParams[0] = lineParamsArray[iline][0];
    } else {
      lineParams.resize(lineParamsArray[iline].size());
      lineParams = lineParamsArray[iline];
    }

    // first do case of zero width line. assume for the moment that this occurs
    // if all the lineParams are zero

    bool deltaFunction(true);
    for (size_t i=0; i<lineParams.size(); i++) {
      if (lineParams[i] != 0.0) {
	deltaFunction = false;
	break;
      }
    }

    if ( deltaFunction ) {
      if ( icen >= 0 ) fluxArray[icen] += lineflux;
      continue;
    }

    // If the line center is below first bin then don't calculate the 
    // lower part of the line. If line center is above the first bin
    // then just calculate part of line within the energy range.

    if ( ecenter >= efirst ) {
      int ielow = icen;
      Real alow = 0.0;
      Real fractionInsideRange = lineFraction(lineShape, efirst, ecenter, lineParams, qspeedy);
      if ( ecenter > elast ) {
	ielow = nE-2;
	alow = lineFraction(lineShape, elast, ecenter, lineParams, qspeedy);
	fractionInsideRange -= alow;
      }

      // Do the low energy part of the line

      Real lineSum = 0.0;
      Real ahi;
      while ( ielow >= 0 ) {
	ahi = lineFraction(lineShape, energyArray[ielow], ecenter, lineParams, qspeedy);
	Real fract = ahi-alow;

	fluxArray[ielow] += fract*lineflux;
	lineSum += fract;
	if ( (fractionInsideRange-lineSum) < crtLevel ) {
	  // Too many sigma away so stop now and add the rest of the line into this
	  // bin. Not strictly correct but shouldn't matter and ensures that the total
	  // flux is preserved.
	  fluxArray[ielow] += (fractionInsideRange-lineSum)*lineflux;
	  ielow = 0;
	}
	alow = ahi;
	ielow -= 1;
      }
    }

    // If the line center is above the last bin then don't calculate 
    // the upper part of the line. If line center is below first bin then 
    // just calculate the part of the line within energy range.

    if ( ecenter <= elast ) {
      Real alow = 0.0;
      int ielow = icen;
      Real fractionInsideRange = lineFraction(lineShape, elast, ecenter, lineParams, qspeedy);
      if ( ecenter <  efirst ) {
	ielow = 0;
	alow = lineFraction(lineShape, energyArray[ielow], ecenter, lineParams, qspeedy);
	fractionInsideRange -= alow;
      }

      // Do the high energy part of the line

      Real lineSum = 0.0;
      Real ahi;
      while ( ielow < nE-1 ) {
	ahi = lineFraction(lineShape, energyArray[ielow+1], ecenter, lineParams, qspeedy);
	Real fract = ahi-alow;
	fluxArray[ielow] += fract*lineflux;
	lineSum += fract;
	if ( (fractionInsideRange-lineSum) < crtLevel ) {
	  // Too many sigma away so stop now and add the rest of the Gaussian into this
	  // bin. Not strictly correct but shouldn't matter and ensures that the total
	  // flux is preserved.
	  fluxArray[ielow] += (fractionInsideRange-lineSum)*lineflux;
	  ielow = nE;
	}
	alow = ahi;
	ielow += 1;
      }
    }
	  
  }

  return;
}

void calcLine(const RealArray& energyArray, const Real ecenter, 
	      const RealArray& lineParams, const Real lineflux, const Real crtLevel, 
	      const int lineShape, const bool qspeedy, RealArray& fluxArray)
{
  RealArray ecenterArray(ecenter,1);
  std::vector<RealArray> lineParamsArray(1,lineParams);
  RealArray linefluxArray(lineflux,1);
  calcManyLines(energyArray, ecenterArray, lineParamsArray, linefluxArray, 
		crtLevel, lineShape, qspeedy, fluxArray);
  return;
}

// Calculates the fraction of the line between energy and ecenter.
// If lineShape==GAUSS then lineParams[0] is sigma
// If lineShape==LORENTZ then lineParams[1] is width
// If lineShape==VOIGT then lineParams[0] is sigma and lineParams[1] is gamma

Real lineFraction(const int lineShape, const Real energy, const Real ecenter,
		  const RealArray& lineParams, const bool qspeedy)
{
  static bool first(true);
  static Real saveEcenter, saveWidth, lnorm;
  if ( lineShape == GAUSS ) {
    return gaussFraction(fabs(energy-ecenter)/lineParams[0], qspeedy)/2.0;
  } else if ( lineShape == LORENTZ ) {
    if ( first || ecenter != saveEcenter || saveWidth != lineParams[0] ) {
      saveEcenter = ecenter;
      saveWidth = lineParams[0];
      lnorm = 1.0/(M_PI/2.0 - atan(-2.0*ecenter/lineParams[0]));
      first = false;
    }
    Real lfrac = lorentzFraction(fabs(energy-ecenter)/lineParams[0], qspeedy)/2.0;
    return lfrac * lnorm;
  } else if ( lineShape == VOIGT ) {
    return voigtFraction(energy, ecenter, lineParams[0], lineParams[1]);
  } else {
    return 0.0;
  }
}

Real gaussFraction(const Real deltasigma, const bool qspeedy)
{

  // Function to return the integral of a Gaussian(0,1) distribution from 
  // -deltasigma to +deltasigma
  // If qspeedy=true interpolates on a previously calculated grid of erf calculations
  // while if qspeedy=false then calls erf each time.

  static const size_t GaussMaxStep(1200);
  static const size_t GaussMax(6);
  static const Real GaussStep = (Real)GaussMax / GaussMaxStep;

  static RealArray tabErf(GaussMaxStep+1);
  static bool first(true);

  // if the first time through then calculate the erf grid on which we
  // will interpolate

  if ( qspeedy && first ) {
    for (size_t i=0; i<=GaussMaxStep; i++) {
      tabErf[i] = erf(i*GaussStep/SQRT2);
    }
    first = false;
  }
  
  // how many sigmas away we are
  Real x = fabs( deltasigma );

  if ( qspeedy ) {

    // Now interpolate from the table
    size_t index = (size_t)(x/GaussStep);

    // If we're past the edge of the tabulated data return 1.0
    if ( index >= GaussMaxStep ) return 1.0;

    Real remainder = (x - (Real)index*GaussStep) * (1.0/GaussStep);

    // Do the interpolation
    return (1.0-remainder)*tabErf[index] +
      remainder*tabErf[index+1];

  } else {

    return erf(x/SQRT2);

  }

}


Real lorentzFraction(const Real deltasigma, const bool qspeedy)
{
  // Function to return the integral of a Lorentzian(0,1) distribution from 
  // -deltasigma to +deltasigma. There is an additional normalization factor
  // which depends on the line energy and width of 1/(pi/2 - arctan(-2*E0/W))
  // and should be applied by the routine calling lorentzFraction.
  // If qspeedy=true interpolates on a previously calculated grid of
  // calculations for low values, but pi-1/x for large values
  // while if qspeedy=false then does calculation each time. The difference
  // in values for the fast and slow methods is less than 1e-3.

  static const size_t LorentzMaxStep(250);
  static const size_t LorentzMax(5);
  static const Real LorentzStep = (Real)LorentzMax / LorentzMaxStep;
  
  static RealArray tabLor(LorentzMaxStep+1);
  static bool first(true);

  // how many sigmas away we are
  Real x = fabs( deltasigma );

  if ( qspeedy ) {

    // if the first time through then calculate the grid on which we
    // will interpolate
    if ( first ) {
      for (size_t i=0; i<=LorentzMaxStep; i++) {
        tabLor[i] = 2*atan(2*i*LorentzStep);
      }
      first = false;
    }
    
    // Now interpolate from the table
    Real xstep = x*(1.0/LorentzStep);
    size_t index = (size_t)(xstep);

    if ( index >= LorentzMaxStep ) {
      // If we're past the edge of the tabulated data return
      // asymptotic function for the line wings
      return M_PI - 1/x;
    } else {
      Real remainder = xstep - index;

      // Do the interpolation
      return (1.0-remainder)*tabLor[index] +
        remainder*tabLor[index+1];
    }

  } else {

    return 2*atan(2*x);

  }
}

Real voigtFraction(const Real energy, const Real ecenter, const Real sigma, const Real gamma)
{
  // Function to return the fraction of the Voigt line profile between energy and ecenter
  // using adaptive integration of the Faddeeva function

  RealArray p(3);
  p[0] = ecenter;
  p[1] = sigma;
  p[2] = gamma;

  Real emin, emax;
  if ( energy > ecenter ) {
    emin = ecenter;
    emax = energy;
  } else {
    emin = energy;
    emax = ecenter;
  }

  Real fraction, fractionError;
  Real precision(1.0e-6);
  AdaptiveIntegrate<voigtLineIntegrand>(emin, emax, &p, precision, fraction, fractionError);

  return fraction;
}

// evaluate the Voigt line integrand with p used to pass in the parameter values

RealArray voigtLineIntegrand(const RealArray& x, void *p)
{
  RealArray OutArray(x.size());

  RealArray params(3);
  params = *(static_cast<RealArray*>(p));

  Real center = params[0];
  Real sigma = params[1];
  Real gamma = params[2];

  Real vnorm = 1.0/sigma/sqrt(2*M_PI);

  Real znorm = 1.0/sigma/sqrt(2.0);
  for (size_t i = 0; i<x.size(); i++) {
    std::complex<double> z(x[i]-center,gamma);
    z *= znorm;
    std::complex<double> wofz = Faddeeva::w(z);
    OutArray[i] = vnorm * wofz.real();
  }

  return OutArray;
}
