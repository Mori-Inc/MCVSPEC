
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/IosHolder.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSUtil/Numerics/CosmologyFunction.h>
#include <XSUtil/Numerics/Gamma.h>
#include <XSUtil/Numerics/LinearInterp.h>
#include <XSUtil/Numerics/Numerics.h>
#include <cstdlib>
#include <string>
#include <cstring>
#include <unistd.h>

char* FGDATD()
{
        return const_cast<char*>(FunctionUtility::managerPath().c_str());

} 

char* FGMODF()
{
        return const_cast<char*>(FunctionUtility::modelDataPath().c_str());       
} 

char* FGXSCT()
{
        return const_cast<char*>(FunctionUtility::XSECT().c_str());

}

char* FGSOLR()
{
        return const_cast<char*>(FunctionUtility::ABUND().c_str());

}

float FGABND(const char* element)
{
  return FunctionUtility::getAbundance(string(element));      
}

float FGABNZ(const int Z)
{
  return FunctionUtility::getAbundance((size_t)Z);      
}

float FGTABN(const char* table, const char* element)
{
  return FunctionUtility::getAbundance(string(table), string(element));      
}

float FGTABZ(const char* table, const int Z)
{
  return FunctionUtility::getAbundance(string(table), (size_t)Z);      
}

char* FGELTI(const int index)
{
  return const_cast<char*>(FunctionUtility::elements((size_t)index).c_str());      
}

int FGNELT()
{
  return (int)FunctionUtility::NELEMS();      
}

char* FGABFL()
{
  return const_cast<char*>(FunctionUtility::abundanceFile().c_str());      
}

void FPABFL(const char* fname)
{
  FunctionUtility::abundanceFile(string(fname));      
  return;
}

char* FGAPTH()
{
  return const_cast<char*>(FunctionUtility::abundPath().c_str());      
}

void FPAPTH(const char* abunDir)
{
  FunctionUtility::abundPath(string(abunDir));      
  return;
}


char* FGMSTR(const char* dname)
{
        static string value;
        value = FunctionUtility::getModelString(string(dname));
        if (value == FunctionUtility::NOT_A_KEY())
        {
           value.erase();
        }
        return const_cast<char*>(value.c_str());      
}

void FPDATD(const char* dataDir)
{
   FunctionUtility::managerPath(string(dataDir));
}

void FPSOLR(const char* table, int* ierr)
{
   string tableName = string(table);
   tableName = XSutility::lowerCase(tableName);
   if (tableName == "file")
   {
      FunctionUtility::ABUND(tableName);
      *ierr = 0;
   }
   else
   {
      if (FunctionUtility::checkAbund(tableName))
      {
         FunctionUtility::ABUND(tableName);
         *ierr = 0;
      }
      else
      {
         *ierr = 1;
      }
   }
}

void FPXSCT(const char* csection, int* ierr)
{
   string tableName = string(csection);
   tableName = XSutility::lowerCase(tableName);
   if (FunctionUtility::checkXsect(tableName))
   {
      FunctionUtility::XSECT(tableName);
      *ierr = 0;
   }
   else
   {
      *ierr = 1;
   }
}

void FPMSTR(const char* value1, const char* value2)
{
   string svalue1 = XSutility::upperCase(string(value1));
   if (svalue1 == "INITIALIZE")
   {
      FunctionUtility::eraseModelStringDataBase();
   }
   else
   {
      FunctionUtility::setModelString(svalue1, string(value2));
   }
}

void FPSLFL(float rvalue[], int nvalue, int *ierr)
{
   // Load the values of the file solar abundance table.
   *ierr = 0;
   size_t nelem = static_cast<size_t>(nvalue);
   if (nelem != FunctionUtility::NELEMS())
   {
      *ierr = 1;
   } 
   else
   {
      std::vector<float> rvect(nelem);
      for (size_t i=0; i<nelem; ++i)
      {
         rvect[i] = rvalue[i];
      }
      FunctionUtility::abundanceVectors("file",rvect);
   }  
}

void RFLABD(const char* fname, int *ierr)
{
   string fNameStr(fname);
   *ierr = 0;
   try
   {
      FunctionUtility::readNewAbundances(fNameStr);
   }
   catch (...)
   {
      *ierr = 1;
   }
}

void FNINIT()
{
   // get directory paths
   const char* cHeadasLoc = getenv("HEADAS");
   const char* cHomeLoc = getenv("HOME");

   if (!cHeadasLoc) {
     *IosHolder::outHolder() << "\n***Error: HEADAS environment variable not set."
			     << std::endl;
     return;
   }
   string HOME;
   if (cHomeLoc) {
      HOME = string(cHomeLoc);
   }
   const string HEADAS(cHeadasLoc);

   string spectralLoc = HEADAS + "/../spectral/";
   if (access(spectralLoc.c_str(),R_OK)) {
     *IosHolder::outHolder() << "\n***Error: Unable to find spectral directory."
			     << std::endl;
     return;          
   }

   // Initialize the directory for the mode data files. Use the
   // XSPEC_MDATA_DIR environment variable if it is set, otherwise
   // use the standard directory.
   string modelDataLoc;
   const char* cMdataLoc = getenv("XSPEC_MDATA_DIR");
   if (cMdataLoc)
   {
      modelDataLoc = string(cMdataLoc);
   }
   else
   {
      modelDataLoc = spectralLoc + "modelData/";
   }

   string managerLoc = spectralLoc + "manager";

   FunctionUtility::modelDataPath(modelDataLoc);
   FunctionUtility::managerPath(managerLoc);

   // read the default settings
   string defaultSettingsFileName(managerLoc+"/Xspec.init");
   std::map<string,string> defaultSettings = XSutility::readSettingsFile(defaultSettingsFileName);

   // read the user's settings
   string userSettingsFileName(HOME+"/.xspec/Xspec.init");
   std::map<string,string> userSettings = XSutility::readSettingsFile(userSettingsFileName);
   
   // Read in the abundance and crosssections.dat files.  
   string dummy1, dummy2;
   FunctionUtility::readInitializers(dummy1, dummy2);

   int ierr=1;

   std::map<string,string>::iterator it;
   string value;

   // Initialize the solar abundance table in use.
   it = userSettings.find("ABUND");
   if ( it != userSettings.end() ) {
     value = it->second;
   } else {
     it = defaultSettings.find("ABUND");
     if ( it != defaultSettings.end() ) {
       value = it->second;
     } else {
       value = "angr";
     }
   }
   
   FPSOLR(value.c_str(), &ierr);
   if (ierr) {
     *IosHolder::outHolder() << "\n***Error: Failed to set Solar abundance table."<<std::endl;
     return;
   }

   // Initialize the photoelectric cross-sections in use.
   it = userSettings.find("XSECT");
   if ( it != userSettings.end() ) {
     value = it->second;
   } else {
     it = defaultSettings.find("XSECT");
     if ( it != defaultSettings.end() ) {
       value = it->second;
     } else {
       value = "bcmc";
     }
   }

   FPXSCT(value.c_str(), &ierr);
   if (ierr) {
     *IosHolder::outHolder() << "\n***Error: Failed to set photoelectric cross-section"<<std::endl;
     return;
   }

   // Initialize the list of model string parameters
   FPMSTR("INITIALIZE"," ");

   // Initialize the ATOMDB_VERSION and NEI version
   it = userSettings.find("ATOMDB_VERSION");
   if ( it != userSettings.end() ) {
     value = it->second;
   } else {
     it = defaultSettings.find("ATOMDB_VERSION");
     if ( it != defaultSettings.end() ) {
       value = it->second;
     } else {
       value = "3.0.9";
     }
   }
   FunctionUtility::atomdbVersion(value);

   it = userSettings.find("NEI_VERSION");
   if ( it != userSettings.end() ) {
     value = it->second;
   } else {
     it = defaultSettings.find("NEI_VERSION");
     if ( it != defaultSettings.end() ) {
       value = it->second;
     } else {
       value = "3.0.4";
     }
   }
   FunctionUtility::neiVersion(value);
}

// Emulates the xanlib routine FXWRITE.
int xs_write(char *wrtstr,int  idest)
{
  FunctionUtility::xsWrite(string(wrtstr), idest);
  return 0;
}

// Emulates XCREAD, etc, but assumes interactive input.
int xs_read(const char *prompt,char * buffer,int* ierr)
{
  XSstream* xsin = dynamic_cast<XSstream*>(IosHolder::inHolder());
  char  bk[] = "/*";
  if (xsin)
  {      
     char* temp_prompt =  new char [strlen(IosHolder::xsPrompt())+1]; 
     XSstream::setPrompter(*xsin,string(prompt));
     *xsin >> buffer;
     XSstream::setPrompter(*xsin,string(temp_prompt));
     delete [] temp_prompt;
  }
  else
  {
     *IosHolder::inHolder() >> buffer;
  }

  if( !strncmp(buffer, bk, strlen(bk)) )
  {
        *ierr = -1;
        return 0;
  }

  *ierr = 0;

  return 0;
}

float csmgq0()
{
   return FunctionUtility::getq0();
}

float csmgh0()
{
   return FunctionUtility::getH0();
}

float csmgl0()
{
   return FunctionUtility::getlambda0();
}

void csmpq0(const float q0)
{
  FunctionUtility::setq0(q0);
}

void csmph0(const float H0)
{
  FunctionUtility::setH0(H0);
}

void csmpl0(const float lambda0)
{
  FunctionUtility::setlambda0(lambda0);
}

void csmpall(const float H0, const float q0, const float lambda0)
{
  FunctionUtility::setFunctionCosmoParams((Real)H0, (Real)q0, (Real)lambda0);
}

float fzsq(const float z, const float q0, const float lambda)
{
   Numerics::FZSQ val;
   return (float)val(z,q0,lambda);
}

int DGNFLT(int ifl)
{
  return FunctionUtility::getNumberXFLT(ifl);
}

float DGFILT(int ifl, const char* key)
{
  return FunctionUtility::getXFLT(ifl, string(key));
}

bool DGQFLT(int ifl, const char* key)
{
  return FunctionUtility::inXFLT(ifl, string(key));
}

void DPFILT(int ifl, int nfilt, const char* key, const float keyval)
{
  std::map<string, Real> filtmap;
  filtmap.insert(std::pair<string,Real>(string(key),keyval));

  FunctionUtility::loadXFLT(ifl, filtmap);
  return;
}

void DCLFLT()
{
  FunctionUtility::clearXFLT();
  return;
}

double GDBVAL(const char* keyword)
{
  return FunctionUtility::getDbValue(string(keyword));
}

void PDBVAL(const char* keyword, double value)
{
  FunctionUtility::loadDbValue(string(keyword), value);
  return;
}

void CDBASE()
{
  FunctionUtility::clearDb();
  return;
}

char* FGATDV()
{
  return const_cast<char*>(FunctionUtility::atomdbVersion().c_str());
}

void FPATDV(const char* version)
{
  FunctionUtility::atomdbVersion(string(version));
  return;
}

char* FGNEIV()
{
  return const_cast<char*>(FunctionUtility::neiVersion().c_str());
}

void FPNEIV(const char* version)
{
  FunctionUtility::neiVersion(string(version));
  return;
}

int FGCHAT()
{
   return FunctionUtility::xwriteChatter();
}

void FPCHAT(int chat)
{
   FunctionUtility::xwriteChatter(chat);
}

// Emulates xanlib/xparse's xgtcht, but only when this library is
// linked into xspec's executable.  For external programs this
// simply returns levels = 0.
void xs_getChat(int* cons, int* log)
{
  XSstream* xsout = dynamic_cast<XSstream*>(IosHolder::outHolder());
  if (xsout)
  {
     *cons = xsout->consoleChatterLevel();
     *log = xsout->logChatterLevel();
  }
  else
  {
     *cons = 0;
     *log = 0;
  }   
}

int xs_getVersion(char* buffer, int buffSize)
{
   const string& versStr(XSutility::xs_version());
   const int len = static_cast<int>(versStr.length());
   int isOK = -1;
   if (buffSize > len)
   {
      versStr.copy(buffer, len);
      buffer[len] = 0;
      isOK = 0;
   }
   else if (buffSize > 0)
   {
      versStr.copy(buffer, buffSize-1);
      buffer[buffSize-1] = 0;
   }
   return isOK;
}


float xs_erf(float x)
{
        Numerics::Erf E;
        Real y(x);
        return float(E(y));
}

float xs_erfc(float x)
{
        Numerics::Erfc E;
        Real y(x);
        return float(E(y));
}

float gammap(float a, float x)
{
        Numerics::GammaP G;
        Real y(x);
        Real b(a);
        return float(G(b,y));
}



float gammq(float a, float x)
{
        Numerics::GammaQ G;
        Real y(x);
        Real b(a);
        return float(G(b,y));
}

void tabint(const float* ear, const int ne, const float* param,
	    const int npar, const char* filenm, int ifl,
	    const char* tabtyp, float* photar, float* photer)
{
  string fileName(filenm);
  string initString("");
  string tableType(tabtyp);

  RealArray energyArray(.0,ne+1);
  RealArray parameters(npar);
  RealArray fluxArray, fluxErrArray;

  for (size_t i=0; i<(size_t)ne+1; i++) energyArray[i] = ear[i];
  for (size_t i=0; i<(size_t)npar; i++) parameters[i] = param[i];

  FunctionUtility::tableInterpolate(energyArray, parameters, fileName, ifl,
				    fluxArray, fluxErrArray, initString,
				    tableType, true);

  for (size_t i=0; i<(size_t)ne; i++) photar[i] = fluxArray[i];
  if ( fluxErrArray.size() > 0 ) {
    for (size_t i=0; i<(size_t)ne; i++) photer[i] = fluxErrArray[i];
  } else {
    for (size_t i=0; i<(size_t)ne; i++) photer[i] = 0.0;
  }

  return;
}

void tabintxflt(const float* ear, const int ne, const float* param,
		const int npar, const char* filenm,
		const char **xfltname, const float *xfltvalue,
		const int nxflt,
		const char* tabtyp, float* photar, float* photer)
{
  string fileName(filenm);
  string initString("");
  string tableType(tabtyp);

  RealArray energyArray(.0,ne+1);
  RealArray parameters(npar);
  RealArray fluxArray, fluxErrArray;

  for (size_t i=0; i<(size_t)ne+1; i++) energyArray[i] = ear[i];
  for (size_t i=0; i<(size_t)npar; i++) parameters[i] = param[i];

  std::map<string,Real> spectrumXFLT;
  for (size_t i=0; i<(size_t)nxflt; i++) {
    spectrumXFLT[string(xfltname[i])] = (Real)xfltvalue[i];
  }

  FunctionUtility::tableInterpolate(energyArray, parameters, fileName,
				    fluxArray, fluxErrArray, initString,
				    tableType, true, spectrumXFLT);

  for (size_t i=0; i<(size_t)ne; i++) photar[i] = fluxArray[i];
  if ( fluxErrArray.size() > 0 ) {
    for (size_t i=0; i<(size_t)ne; i++) photer[i] = fluxErrArray[i];
  } else {
    for (size_t i=0; i<(size_t)ne; i++) photer[i] = 0.0;
  }

  return;
}

/* Interface for user C and Fortran subroutines requiring the functions
   in LinearInterp.h in the XSUtil/Numerics library */
bool findFirstBins(const int nCurr, const float* currBins, const int nTarg,
		   const float* targBins, const float fFUZZY, int* currStart,
		   int* targStart)
{
  RealArray currentBins(nCurr);
  RealArray targetBins(nTarg);
  const Real FUZZY((Real)fFUZZY);
  size_t currentStart((size_t)*currStart);
  size_t targetStart((size_t)*targStart);

  for (size_t i=0; i<currentBins.size(); i++) currentBins[i] = currBins[i];
  for (size_t i=0; i<targetBins.size(); i++) targetBins[i] = targBins[i];
  
  bool found = Numerics::Rebin::findFirstBins(currentBins, targetBins, FUZZY,
					      currentStart, targetStart);

  currStart = (int*)&currentStart;
  targStart = (int*)&targetStart;

  return found;
}

void initBins(const int nCurr, const float* currBins, const int nTarg,
	      const float* targBins, const float fFUZZY,
	      int* currStart, int* targStart, int* sBin, int* eBin,
	      float* sWeight, float* eWeight)
{
  RealArray currentBins(nCurr);
  RealArray targetBins(nTarg);
  const Real FUZZY((Real)fFUZZY);
  size_t currentStart((size_t)*currStart);
  size_t targetStart((size_t)*targStart);

  for (size_t i=0; i<currentBins.size(); i++) currentBins[i] = currBins[i];
  for (size_t i=0; i<targetBins.size(); i++) targetBins[i] = targBins[i];

  IntegerVector startBin(nTarg-1), endBin(nTarg-1);
  RealArray startWeight(nTarg-1), endWeight(nTarg-1);
  
  Numerics::Rebin::initializeBins(targetBins, currentBins, FUZZY, targetStart,
				  currentStart, startBin, endBin, startWeight,
				  endWeight);

  currStart = (int*)&currentStart;
  targStart = (int*)&targetStart;

  // NB this assumes that memory was grabbed for sBin, eBin, sWeight,
  // and eWeight in the calling routine
  
  for (size_t i=0; i<startBin.size(); i++) sBin[i] = (int)startBin[i];
  for (size_t i=0; i<endBin.size(); i++) eBin[i] = (int)endBin[i];
  for (size_t i=0; i<startWeight.size(); i++) sWeight[i] = (float)startWeight[i];
  for (size_t i=0; i<endWeight.size(); i++) eWeight[i] = (float)endWeight[i];

  return;
}

void rebinBins(const int nIn, const float* inArray, const int* sBin,
	       const int* eBin, const float* sWeight,
	       const float* eWeight, const int nOut, float* outArray,
	       const float lowVal, const float highVal)
{
  RealArray inputArray(nIn), outputArray(nOut);
  IntegerVector startBin(nOut), endBin(nOut);
  RealArray startWeight(nOut), endWeight(nOut);
  Real lowValue((Real)lowVal);
  Real highValue((Real)highVal);

  for (size_t i=0; i<inputArray.size(); i++) inputArray[i] = (Real)inArray[i];
  for (size_t i=0; i<startBin.size(); i++) startBin[i] = (int)sBin[i];
  for (size_t i=0; i<endBin.size(); i++) endBin[i] = (int)eBin[i];
  for (size_t i=0; i<startWeight.size(); i++) startWeight[i] = (Real)sWeight[i];
  for (size_t i=0; i<endWeight.size(); i++) endWeight[i] = (Real)eWeight[i];

  Numerics::Rebin::rebin(inputArray, startBin, endBin, startWeight, endWeight,
			 outputArray, lowValue, highValue);

  for (size_t i=0; i<outputArray.size(); i++) outArray[i] = (float)outputArray[i];
  
  return;  
}

void interpBins(const int nIn, const float* inArray, const int* sBin,
		const int* eBin, const float* sWeight,
		const float* eWeight, const int nOut, float* outArray,
		const float lowVal, const float highVal)
{
  RealArray inputArray(nIn), outputArray(nOut);
  IntegerVector startBin(nOut), endBin(nOut);
  RealArray startWeight(nOut), endWeight(nOut);
  Real lowValue((Real)lowVal);
  Real highValue((Real)highVal);

  for (size_t i=0; i<inputArray.size(); i++) inputArray[i] = (Real)inArray[i];
  for (size_t i=0; i<startBin.size(); i++) startBin[i] = (int)sBin[i];
  for (size_t i=0; i<endBin.size(); i++) endBin[i] = (int)eBin[i];
  for (size_t i=0; i<startWeight.size(); i++) startWeight[i] = (Real)sWeight[i];
  for (size_t i=0; i<endWeight.size(); i++) endWeight[i] = (Real)eWeight[i];

  Numerics::Rebin::interpolate(inputArray, startBin, endBin, startWeight,
			       endWeight, outputArray, lowValue, highValue);

  for (size_t i=0; i<outputArray.size(); i++) outArray[i] = (float)outputArray[i];
  
  return;  
}

void gainRebin(const int nIn, const float* inArray, const int* sBin,
	       const int* eBin, const float* sWeight,
	       const float* eWeight, const int nOut, float* outArray)
{
  RealArray inputArray(nIn), outputArray(nOut);
  IntegerVector startBin(nOut), endBin(nOut);
  RealArray startWeight(nOut), endWeight(nOut);

  for (size_t i=0; i<inputArray.size(); i++) inputArray[i] = (Real)inArray[i];
  for (size_t i=0; i<startBin.size(); i++) startBin[i] = (int)sBin[i];
  for (size_t i=0; i<endBin.size(); i++) endBin[i] = (int)eBin[i];
  for (size_t i=0; i<startWeight.size(); i++) startWeight[i] = (Real)sWeight[i];
  for (size_t i=0; i<endWeight.size(); i++) endWeight[i] = (Real)eWeight[i];

  Numerics::Rebin::gainRebin(inputArray, startBin, endBin, startWeight,
			     endWeight, outputArray);

  for (size_t i=0; i<outputArray.size(); i++) outArray[i] = (float)outputArray[i];
  
  return;  

}

void linInterpInteg(const int nCurr, const float* currPoints,
		    const float* indValues, const int nTarg,
		    const float* targBins, float* outValues,
		    const float lowVal, const float highVal)
{
  RealArray currentPoints(nCurr);
  RealArray targetBins(nTarg);
  RealArray inputdValues(nCurr);
  RealArray outputValues(nTarg-1);
  Real lowValue((Real)lowVal);
  Real highValue((Real)highVal);

  for (size_t i=0; i<currentPoints.size(); i++) currentPoints[i] = (Real)currPoints[i];
  for (size_t i=0; i<inputdValues.size(); i++) inputdValues[i] = (Real)indValues[i];
  for (size_t i=0; i<targetBins.size(); i++) targetBins[i] = (Real)targBins[i];
  
  Numerics::Rebin::linInterpInteg(currentPoints, inputdValues, targetBins,
				  outputValues, lowValue, highValue);

  // Note that this assumes that the memory for outValues has been acquired
  // by the calling routine
  
  for (size_t i=0; i<outputValues.size(); i++) outValues[i] = (float)outputValues[i];

  return;
}
/* double precision versions of above */
bool dfindFirstBins(const int nCurr, const double* currBins, const int nTarg,
		   const double* targBins, const double fFUZZY, int* currStart,
		   int* targStart)
{
  RealArray currentBins(nCurr);
  RealArray targetBins(nTarg);
  const Real FUZZY((Real)fFUZZY);
  size_t currentStart((size_t)*currStart);
  size_t targetStart((size_t)*targStart);

  for (size_t i=0; i<currentBins.size(); i++) currentBins[i] = currBins[i];
  for (size_t i=0; i<targetBins.size(); i++) targetBins[i] = targBins[i];
  
  bool found = Numerics::Rebin::findFirstBins(currentBins, targetBins, FUZZY,
					      currentStart, targetStart);

  currStart = (int*)&currentStart;
  targStart = (int*)&targetStart;

  return found;
}

void dinitBins(const int nCurr, const double* currBins, const int nTarg,
	      const double* targBins, const double fFUZZY,
	      int* currStart, int* targStart, int* sBin, int* eBin,
	      double* sWeight, double* eWeight)
{
  RealArray currentBins(nCurr);
  RealArray targetBins(nTarg);
  const Real FUZZY((Real)fFUZZY);
  size_t currentStart((size_t)*currStart);
  size_t targetStart((size_t)*targStart);

  for (size_t i=0; i<currentBins.size(); i++) currentBins[i] = currBins[i];
  for (size_t i=0; i<targetBins.size(); i++) targetBins[i] = targBins[i];

  IntegerVector startBin(nTarg-1), endBin(nTarg-1);
  RealArray startWeight(nTarg-1), endWeight(nTarg-1);
  
  Numerics::Rebin::initializeBins(targetBins, currentBins, FUZZY, targetStart,
				  currentStart, startBin, endBin, startWeight,
				  endWeight);

  currStart = (int*)&currentStart;
  targStart = (int*)&targetStart;

  // NB this assumes that memory was grabbed for sBin, eBin, sWeight,
  // and eWeight in the calling routine
  
  for (size_t i=0; i<startBin.size(); i++) sBin[i] = (int)startBin[i];
  for (size_t i=0; i<endBin.size(); i++) eBin[i] = (int)endBin[i];
  for (size_t i=0; i<startWeight.size(); i++) sWeight[i] = (double)startWeight[i];
  for (size_t i=0; i<endWeight.size(); i++) eWeight[i] = (double)endWeight[i];

  return;
}

void drebinBins(const int nIn, const double* inArray, const int* sBin,
	       const int* eBin, const double* sWeight,
	       const double* eWeight, const int nOut, double* outArray,
	       const double lowVal, const double highVal)
{
  RealArray inputArray(nIn), outputArray(nOut);
  IntegerVector startBin(nOut), endBin(nOut);
  RealArray startWeight(nOut), endWeight(nOut);
  Real lowValue((Real)lowVal);
  Real highValue((Real)highVal);

  for (size_t i=0; i<inputArray.size(); i++) inputArray[i] = (Real)inArray[i];
  for (size_t i=0; i<startBin.size(); i++) startBin[i] = (int)sBin[i];
  for (size_t i=0; i<endBin.size(); i++) endBin[i] = (int)eBin[i];
  for (size_t i=0; i<startWeight.size(); i++) startWeight[i] = (Real)sWeight[i];
  for (size_t i=0; i<endWeight.size(); i++) endWeight[i] = (Real)eWeight[i];

  Numerics::Rebin::rebin(inputArray, startBin, endBin, startWeight, endWeight,
			 outputArray, lowValue, highValue);

  for (size_t i=0; i<outputArray.size(); i++) outArray[i] = (double)outputArray[i];
  
  return;  
}

void dinterpBins(const int nIn, const double* inArray, const int* sBin,
		const int* eBin, const double* sWeight,
		const double* eWeight, const int nOut, double* outArray,
		const double lowVal, const double highVal)
{
  RealArray inputArray(nIn), outputArray(nOut);
  IntegerVector startBin(nOut), endBin(nOut);
  RealArray startWeight(nOut), endWeight(nOut);
  Real lowValue((Real)lowVal);
  Real highValue((Real)highVal);

  for (size_t i=0; i<inputArray.size(); i++) inputArray[i] = (Real)inArray[i];
  for (size_t i=0; i<startBin.size(); i++) startBin[i] = (int)sBin[i];
  for (size_t i=0; i<endBin.size(); i++) endBin[i] = (int)eBin[i];
  for (size_t i=0; i<startWeight.size(); i++) startWeight[i] = (Real)sWeight[i];
  for (size_t i=0; i<endWeight.size(); i++) endWeight[i] = (Real)eWeight[i];

  Numerics::Rebin::interpolate(inputArray, startBin, endBin, startWeight,
			       endWeight, outputArray, lowValue, highValue);

  for (size_t i=0; i<outputArray.size(); i++) outArray[i] = (double)outputArray[i];
  
  return;  
}

void dgainRebin(const int nIn, const double* inArray, const int* sBin,
	       const int* eBin, const double* sWeight,
	       const double* eWeight, const int nOut, double* outArray)
{
  RealArray inputArray(nIn), outputArray(nOut);
  IntegerVector startBin(nOut), endBin(nOut);
  RealArray startWeight(nOut), endWeight(nOut);

  for (size_t i=0; i<inputArray.size(); i++) inputArray[i] = (Real)inArray[i];
  for (size_t i=0; i<startBin.size(); i++) startBin[i] = (int)sBin[i];
  for (size_t i=0; i<endBin.size(); i++) endBin[i] = (int)eBin[i];
  for (size_t i=0; i<startWeight.size(); i++) startWeight[i] = (Real)sWeight[i];
  for (size_t i=0; i<endWeight.size(); i++) endWeight[i] = (Real)eWeight[i];

  Numerics::Rebin::gainRebin(inputArray, startBin, endBin, startWeight,
			     endWeight, outputArray);

  for (size_t i=0; i<outputArray.size(); i++) outArray[i] = (double)outputArray[i];
  
  return;  

}

void dlinInterpInteg(const int nCurr, const double* currPoints,
		    const double* indValues, const int nTarg,
		    const double* targBins, double* outValues,
		    const double lowVal, const double highVal)
{
  RealArray currentPoints(nCurr);
  RealArray targetBins(nTarg);
  RealArray inputdValues(nCurr);
  RealArray outputValues(nTarg-1);
  Real lowValue((Real)lowVal);
  Real highValue((Real)highVal);

  for (size_t i=0; i<currentPoints.size(); i++) currentPoints[i] = (Real)currPoints[i];
  for (size_t i=0; i<inputdValues.size(); i++) inputdValues[i] = (Real)indValues[i];
  for (size_t i=0; i<targetBins.size(); i++) targetBins[i] = (Real)targBins[i];
  
  Numerics::Rebin::linInterpInteg(currentPoints, inputdValues, targetBins,
				  outputValues, lowValue, highValue);

  // Note that this assumes that the memory for outValues has been acquired
  // by the calling routine
  
  for (size_t i=0; i<outputValues.size(); i++) outValues[i] = (double)outputValues[i];

  return;
}

// Handy routine to get constants defined in Numerics.h
double getkeVtoA() { return Numerics::KEVTOA; }
double getkeVtoHz() { return Numerics::KEVTOHZ; }
double getkeVtoErg() { return Numerics::KEVTOERG; }
double getkeVtoJy() { return Numerics::KEVTOJY; }
double getdegtorad() { return Numerics::DEGTORAD; }
double getlightspeed() { return Numerics::LIGHTSPEED; }
double getamu() { return Numerics::AMU; }
