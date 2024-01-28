
#ifndef XSFORTRAN_H
#define XSFORTRAN_H
#include <sys/types.h>
#include <stdbool.h>

/* functions that the fortran library needs that are implemented in C++ now. */


#ifdef __cplusplus
extern "C"
{
#endif

  char* FGMODF(void);
  char* FGDATD(void);
  float FGABND(const char* element);
  float FGABNZ(const int Z);
  float FGTABN(const char* table, const char* element);
  float FGTABZ(const char* table, const int Z);
  char* FGELTI(const int index);
  int   FGNELT();
  char* FGABFL();
  void  FPABFL(const char* fname);
  char* FGAPTH();
  void  FPAPTH(const char* abunDir);
  char* FGXSCT(void); 
  char* FGSOLR(void);
  char* FGMSTR(const char* dname);
  void  FPDATD(const char* dataDir);
  void  FPSOLR(const char* table, int* ierr);
  void  FPXSCT(const char* csection, int* ierr);
  void  FPMSTR(const char* value1, const char* value2);
  void  FPSLFL(float rvalue[], int nvalue, int *ierr);
  void  RFLABD(const char* fname, int *ierr);
  void  FNINIT(void);
  int   xs_write(char*, int);
  int   xs_read(const char*, char*, int*);
  float csmgq0(void);
  float csmgh0(void);
  float csmgl0(void);
  void  csmpq0(const float q0);
  void  csmph0(const float H0);
  void  csmpl0(const float lambda0);
  void  csmpall(const float H0, const float q0, const float lambda0);
  float fzsq(const float z, const float q0, const float lambda);
  int   DGNFLT(int);
  float DGFILT(int, const char*);
  bool  DGQFLT(int, const char*);
  void  DPFILT(int, const int, const char*, const float);
  void  DCLFLT();
  double GDBVAL(const char*);
  void  PDBVAL(const char*, double);
  void  CDBASE();
  char* FGATDV();
  void  FPATDV(const char*);
  char* FGNEIV();
  void  FPNEIV(const char*);

  /* These are not xspec's chatter levels.  They are intended for
     use by other programs which link to XSFunctions. */
  int   FGCHAT(void);
  void  FPCHAT(int chat);
  /* This retrieves xspec's chatter levels, or 0 if not calling
     from xspec. */
  void  xs_getChat(int* cons, int* log);

  /* Copies xspec's version string into buffer. 
     Will copy at most buffSize-1 characters, and pads with
     a trailing NULL.  Returns -1 if the copy was truncated, 
     else returns 0. */
  int   xs_getVersion(char* buffer, int buffSize);

  /* Interface for user C and Fortran subroutines requiring use of
     Gamma functions in the Numerics subdirectory of the XSUtil library. */
  float xs_erf( float x);
  float xs_erfc(float x);
  float gammap(float a, float x);
  float gammq(float a, float x);

  /* Interface for tableInterpolate with spectrum number */
  void tabint(const float* ear, const int ne, const float* param,
	      const int npar, const char* filenm, int ifl,
	      const char* tabtyp, float* photar, float* photer);
  /* Interface for tableInterpolate passing xflt information instead of
     spectrum number */
  void tabintxflt(const float* ear, const int ne, const float* param,
		  const int npar, const char* filenm,
		  const char **xfltname, const float *xfltvalue,
		  const int nxflt,
		  const char* tabtyp, float* photar, float* photer);

  /* Interface for user C and Fortran subroutines requiring the functions
     in LinearInterp.h in the XSUtil/Numerics library */
  bool findFirstBins(const int nCurr, const float* currBins, const int nTarg,
		     const float* targBins, const float fFUZZY, int* currStart,
		     int* targStart);
  void initBins(const int nCurr, const float* currBins, const int nTarg,
		const float* targBins, const float fFUZZY,
		int* currStart, int* targStart, int* sBin, int* eBin,
		float* sWeight, float* eWeight);
  void rebinBins(const int nIn, const float* inArray, const int* sBin,
		 const int* eBin, const float* sWeight,
		 const float* eWeight, const int nOut, float* outArray,
		 const float lowVal, const float highVal);
  void interpBins(const int nIn, const float* inArray, const int* sBin,
		  const int* eBin, const float* sWeight,
		  const float* eWeight, const int nOut, float* outArray,
		  const float lowVal, const float highVal);
  void gainRebin(const int nIn, const float* inArray, const int* sBin,
		 const int* eBin, const float* sWeight,
		 const float* eWeight, const int nOut, float* outArray);
  void linInterpInteg(const int nCurr, const float* currPoints,
		      const float* indValues, const int nTarg,
		      const float* targBins, float* outValues,
		      const float lowVal, const float highVal);
  /* double precision versions */
  bool dfindFirstBins(const int nCurr, const double* currBins, const int nTarg,
		     const double* targBins, const double fFUZZY, int* currStart,
		     int* targStart);
  void dinitBins(const int nCurr, const double* currBins, const int nTarg,
		const double* targBins, const double fFUZZY,
		int* currStart, int* targStart, int* sBin, int* eBin,
		double* sWeight, double* eWeight);
  void drebinBins(const int nIn, const double* inArray, const int* sBin,
		 const int* eBin, const double* sWeight,
		 const double* eWeight, const int nOut, double* outArray,
		 const double lowVal, const double highVal);
  void dinterpBins(const int nIn, const double* inArray, const int* sBin,
		  const int* eBin, const double* sWeight,
		  const double* eWeight, const int nOut, double* outArray,
		  const double lowVal, const double highVal);
  void dgainRebin(const int nIn, const double* inArray, const int* sBin,
		 const int* eBin, const double* sWeight,
		 const double* eWeight, const int nOut, double* outArray);
  void dlinInterpInteg(const int nCurr, const double* currPoints,
		      const double* indValues, const int nTarg,
		      const double* targBins, double* outValues,
		      const double lowVal, const double highVal);

  /* handy routines to get constants defined in Numerics.h */
  double getkeVtoA();
  double getkeVtoHz();
  double getkeVtoErg();
  double getkeVtoJy();
  double getdegtorad();
  double getlightspeed();
  double getamu();
  
#ifdef __cplusplus                 
}
#endif

#endif
