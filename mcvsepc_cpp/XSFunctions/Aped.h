// Classes for handling the APED data. ApedTemperatureRecord is the complete data
// for a given temperature while ApedElementRecord is the complete data for a given
// element within ApedTemperatureRecord and ApedIonRecord is the complete data for
// a given ion within ApedElementRecord.


#include <xsTypes.h>

using namespace std;



class ApedIonRecord{
 public:

  int m_Ion;
  RealArray m_ContinuumEnergy;
  RealArray m_ContinuumFlux;
  RealArray m_ContinuumFluxError;
  RealArray m_PseudoContinuumEnergy;
  RealArray m_PseudoContinuumFlux;
  RealArray m_PseudoContinuumFluxError;
  RealArray m_LineEnergy;
  RealArray m_LineEnergyError;
  RealArray m_LineEmissivity;
  RealArray m_LineEmissivityError;
  IntegerVector m_ElementDriver;
  IntegerVector m_IonDriver;
  IntegerVector m_UpperLevel;
  IntegerVector m_LowerLevel;

  ApedIonRecord();    // default constructor
  ~ApedIonRecord();   // destructor

  ApedIonRecord& operator=(const ApedIonRecord&);    // deep copy

};

class ApedElementRecord{
 public:

  int m_AtomicNumber;
  vector<ApedIonRecord> m_IonRecord;

  ApedElementRecord();    // default constructor
  ~ApedElementRecord();   // destructor

  void LoadIonRecord(ApedIonRecord input);

  ApedElementRecord& operator=(const ApedElementRecord&);    // deep copy

};

class ApedTemperatureRecord{
 public:

  Real m_Temperature;
  vector<ApedElementRecord> m_ElementRecord;

  ApedTemperatureRecord();     //default constructor
  ~ApedTemperatureRecord();    //destructor

  void LoadElementRecord(ApedElementRecord input);

  ApedTemperatureRecord& operator=(const ApedTemperatureRecord&);    // deep copy


};

class Aped{
 public:

  vector<ApedTemperatureRecord> m_TemperatureRecord;

  // These store the information in the initial PARAMETERS extension

  RealArray m_Temperatures;
  IntegerVector m_NelementsLine;
  IntegerVector m_NelementsCoco;
  IntegerVector m_Nline;
  IntegerVector m_Ncont;

  string m_coconame;
  string m_linename;

  bool m_noLines;
  bool m_noResonanceLines;
  bool m_thermalBroadening;
  Real m_velocityBroadening;
  Real m_minimumLinefluxForBroadening;
  bool m_broadenPseudoContinuum;
  bool m_logTempInterpolation;
  bool m_multiThread;

  Aped();     //default constructor
  ~Aped();    //destructor

  void SetNoLines(const bool qno);
  void SetNoResonanceLines(const bool qnores);
  void SetThermalBroadening(const bool qtherm);
  void SetVelocityBroadening(const Real velocity);
  void SetMinimumLinefluxForBroadening(const Real flux);
  void SetBroadenPseudoContinuum(const bool qpbroad);
  void SetLogTempInterpolation(const bool qltinterp);
  void SetMultiThread(const bool qmulti);

  void SetLineSpecs(const bool qno, const bool qnores, const bool qtherm, 
		    const Real velocity, const Real flux, const bool qpbroad,
		    const bool qltinterp, const bool qmulti);
  
  // Reads the continuum and line files and stores the names and
  // temperatures. Does not read the actual continuum and line data.

  int Read(string cocofilename, string linefilename);

  int ReadTemperature(const int TemperatureIndex);
  int ReadTemperature(const vector<int>& TemperatureIndex);

  int NumberTemperatures();     // return number of tabulated temperatures
  int NumberElements();         // return number of tabulated elements
  int NumberIons(int Z);        // return number of ions across all temperatures
                                // for element Z.

  // return true if the data has been loaded for the requested temperature index.

  bool IsTemperatureLoaded(const int TemperatureIndex);

  Aped& operator=(const Aped&);    // deep copy


  // SumEqSpectra sums spectra for the temperatures and associated DEMs. It
  // assumes ionization equilibrium so the object must have been loaded from
  // the ionization equilibrium input files. Versions with single and multiple
  // temperatures for the following methods.

  void SumEqSpectra(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Dem, RealArray& fluxArray, RealArray& fluxErrArray);

  void SumEqSpectra(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Dem, RealArray& fluxArray, RealArray& fluxErrArray);

  // case where the temperature used for the thermal broadening differs from that
  // for the ionization

  void SumEqSpectra(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput, 
		    const Real& Tbinput, const Real& Dem,
		    RealArray& fluxArray, RealArray& fluxErrArray);

  void SumEqSpectra(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput, 
		    const RealArray& Tbinput, const RealArray& Dem,
		    RealArray& fluxArray, RealArray& fluxErrArray);

  // SumNeqSpectra is the equivalent case for NEI where the input ionization 
  // fractions are given as well as the temperatures. Each temperature requires
  // an ionization fraction for all the relevant levels for each element.

  void SumNeqSpectra(const RealArray& energyArray, 
		     const IntegerVector& Zinput, const RealArray& abundance,
		     const Real Redshift, const Real& Tinput,
		     const vector<RealArray>& IonFrac, 
		     RealArray& fluxArray, RealArray& fluxErrArray);

  void SumNeqSpectra(const RealArray& energyArray, 
		     const IntegerVector& Zinput, const RealArray& abundance,
		     const Real Redshift, const RealArray& Tinput,
		     const vector<vector<RealArray> >& IonFrac, 
		     RealArray& fluxArray, RealArray& fluxErrArray);

  // case where the temperature used for the thermal broadening differs from that
  // for the ionization

  void SumNeqSpectra(const RealArray& energyArray, 
		     const IntegerVector& Zinput, const RealArray& abundance,
		     const Real Redshift, const Real& Tinput,
		     const Real& Tbinput, 
		     const vector<RealArray>& IonFrac, 
		     RealArray& fluxArray, RealArray& fluxErrArray);

  void SumNeqSpectra(const RealArray& energyArray, 
		     const IntegerVector& Zinput, const RealArray& abundance,
		     const Real Redshift, const RealArray& Tinput,
		     const RealArray& Tbinput, 
		     const vector<vector<RealArray> >& IonFrac, 
		     RealArray& fluxArray, RealArray& fluxErrArray);

  // SumSpectra underlies both SumEqSpectra and SumNeqSpectra

  void SumSpectra(const RealArray& energyArray, 
		  const IntegerVector& Zinput, const RealArray& abundance,
		  const Real Redshift, const RealArray& Tinput,
		  const RealArray& Tbinput, const vector<vector<RealArray> >& IonFrac, 
		  const bool isCEI, RealArray& fluxArray, RealArray& fluxErrArray);

  int calcContinuumSpectrumForTemperature(const size_t TRecordIndex, const IntegerVector& Zinput,
					  const RealArray& abunZ, const vector<RealArray>& IonFrac,
					  const RealArray& sourceFrameEnergy, const bool isCEI,
					  RealArray& fluxArray);

  int calcLineSpectrumForTemperature(const size_t TRecordIndex, const IntegerVector& Zinput,
				     const RealArray& abunZ, const vector<RealArray>& IonFrac,
				     const RealArray& sourceFrameEnergy, const bool isCEI,
				     const Real minLineflux, RealArray& fluxArray, Real& maxLineflux);

};

// Handy routines

Real getAtomicMass(const int& AtomicNumber);

bool getApedFileNames(string& version, string& continuumFile, string& lineFile);

bool getNEIApedFileNames(string& version, string& eigenVersion, string& continuumFile, string& lineFile);

void apedInterpFlux(const RealArray& inputEnergy, const RealArray& inputFlux, 
		    const Real& coeff, const RealArray& energyArray, 
		    RealArray& fluxArray);

// Wrap-up routines that can be called by outside routines. Creates an internal static Aped 
// object then calls SumSpectra. For the CEI case checks the APECROOT variable to
// find the files to read. The RS case is for the historical Raymond-Smith CEI.
// For the NEI case checks the NEIVERS and NEIAPECROOT variables. Note versions 
// with both single and multiple temperatures. For the CEI case there is the option of
// not calculating any lines

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Tbinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Tbinput,
		    const RealArray& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    const bool noLines, RealArray& fluxArray, RealArray& fluxErrArray);

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Dem, const bool qtherm, const Real velocity,
		    const bool noLines, RealArray& fluxArray, RealArray& fluxErrArray);

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Tbinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    const bool noLines, RealArray& fluxArray, RealArray& fluxErrArray);

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Tbinput,
		    const RealArray& Dem, const bool qtherm, const Real velocity,
		    const bool noLines, RealArray& fluxArray, RealArray& fluxErrArray);

int calcRSSpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcRSSpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcRSSpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Tbinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcRSSpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Tbinput,
		    const RealArray& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcNEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const vector<RealArray>& IonFrac, 
		    const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcNEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const vector<vector<RealArray> >& IonFrac, 
		    const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcNEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Tbinput,
		    const vector<RealArray>& IonFrac, 
		    const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

int calcNEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Tbinput,
		    const vector<vector<RealArray> >& IonFrac, 
		    const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray);

