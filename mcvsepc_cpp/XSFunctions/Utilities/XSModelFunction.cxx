//   Read the documentation to learn more about C++ code generator
//   versioning.
//	  %X% %Q% %Z% %W%

// XSModelFunction
#include <XSFunctions/Utilities/XSModelFunction.h>
#include <XSFunctions/Utilities/XSCall.h>
#include <XSUtil/Utils/IosHolder.h>
// sstream
#include <sstream>
#include <fstream>
#include <queue>
#include <memory>
#include <utility>

// Class XSModelFunction::NoSuchComponent 

XSModelFunction::NoSuchComponent::NoSuchComponent()
   : YellowAlert()  
{
  // probably a gratuitous constructor...
}


// Class XSModelFunction 
CompInfoCacheType XSModelFunction::s_compInfoCache;
const string XSModelFunction::s_NOT_FOUND = string("hN#@xnNOT_FOUNDp!749w9q");

void XSModelFunction::updateComponentList (const string& modelDatFile, bool isStandard)
{
  // This function sets up a cache of model component attributes for 
  // use with the model command. It reads information from modelDatFile and
  // saves ComponentInfo records as follows:

  // ComponentInfo:  name, modelDatFile, type, error flag, additional info string,
  //                 is standard flag, minimum energy, maximum energy

  // initial use of this function is to initialize the entries in the model.dat
  // file. It should also be possible to treat user models as load module addins
  // using this technique [this is the reason for the otherwise redundant model.dat
  // file argument. Cache lookup is provided by compInfoCache and fullName, which
  // are complicated by the fact that models sometimes begin with the same
  // few letters. I chose to resolve this by using a multimap and storing entries
  // against two letter keys, as the least messy way of resolving abbreviations
  // that users use all the time.


  using namespace std;

  ifstream modfile;
  modfile.open(modelDatFile.c_str());

  if (!modfile) 
  {
        string message = "Cannot open model definition file: " + modelDatFile;
        throw RedAlert(message);
  }
  string line("");
  string prevline("");
  size_t linesRead(0);

  /*    the model.dat format consists of stanzas such as 
        bbody          1   1.e-20     1.e20          xsblbd    add  0
        kT      keV     3.0   1.e-4   1.e-2   100.      200.      0.01

        which are delimited by blank lines. The easiest way of telling we have
        a new entry is that the previous line was blank.
  */
  size_t nPar(0);
  Real upper(0);
  Real lower(0);
  string functionName("");
  string type("");
  string name("");
  int errorFlag(0);
  string infoString("");
  const string WS(" \t\n\r");

  while (!modfile.eof()) {

    // prevline && line are blank initially. after some reading has taken place
    // we look for the case where prevline is blank and line is not.
    getline(modfile,line);
    ++linesRead;
    if ((prevline.find_first_not_of(WS) == string::npos) && line.length() != 0) {

      istringstream s(line);  
      s >>  name >> nPar >> lower >> upper >> functionName >> type >> errorFlag; 
      if (!s)
      {
         string errMsg("Format error in first line of ");
         errMsg += name + string(" entry in model description file.\n");
	 errMsg += string(" line is: ") + line + "\n";
         throw YellowAlert(errMsg);
      }
      // Remaining parameters are optional.  They may be:
      // <dependency flag> followed by <info string>,
      // <dependency flag> OR <info string>, or nothing at all. 
      string testString;
      bool dependencyFlag = false;
      s >> testString;
      if (testString.length())
      {
         // At least 1 optional param.  Is it a bool?
         istringstream issTest(testString);
         bool testBool = false;
         if (!(issTest >> testBool) || !issTest.eof())
         {
            // Not a bool, assume info string.
            infoString = testString;
         }
         else
         {
            dependencyFlag = testBool;
            // Still need to test for info string.
            testString.erase();
            if (s >> testString)
            {
               infoString = testString;
            } 
         }
         testString.erase();
         if (s >> testString)
         {
            // This should catch case of dependency flag 
            // following info string.
            *IosHolder::outHolder() << "***Warning: In first line of " << name << " entry of model description file,"
               << "\n     the string \"" << testString << "\" following \"" << infoString
               << "\" will be ignored." << std::endl;
         }
      }
      

      // likely to be empty, but that case ought to be harmless.
      bool error(errorFlag != 0);
      // If an already existing component record has the same 
      // full name, then remove it.
      CompInfoCacheType::iterator match = compIteratorExactMatchName(name);
      if (match != s_compInfoCache.end()) s_compInfoCache.erase(match);

      // read the parameter lines and put them in parStrings
      vector<string> parStrings(nPar);
      for (size_t ipar=0; ipar<nPar; ipar++) {
	getline(modfile,line);
	parStrings[ipar] = line;
	++linesRead;
      }

      ComponentInfo compInfo(name,type,error,infoString,isStandard,lower,upper,parStrings);
      compInfo.isPythonModel(false);
      compInfo.isSpecDependent(dependencyFlag);
      addCompInfoToCache(XSutility::lowerCase(name.substr(0,2)), compInfo);

    }

    prevline = line;       
  }
}

bool XSModelFunction::addMdefToComponentList(const string& name, const string& type,
					     const string& expression,
					     const Real& minEnergy, const Real& maxEnergy)
{
  // Don't let in a name the same as not_found indicator.  It's extremely
  // unlikely anyone could match the not_found name accidentally.
  if (name == s_NOT_FOUND) {
    string err(name);
    err += " is an internally reserved name and is unavailable for mdefine.\n";
    throw YellowAlert(name);
  }

  // If this name already exists then allow a current mdefine model with the same name to
  // be overridden otherwise throw

  try {
    ComponentInfo testCompInfo = XSModelFunction::compMatchName(name);
    if ( !testCompInfo.isMdefineModel() ) {
      string err(name);
      err += " is already in use and cannot be redefined.\n";
      throw YellowAlert(name);
    } else {
      delFromComponentList(name);
    }
  } catch(XSModelFunction::NoSuchComponent&) {
  }

  // we have to go through a convoluted process to get the number of parameters and construct
  // the default parameter strings

  std::unique_ptr<MdefExpression> expr(new MdefExpression(std::make_pair(minEnergy,maxEnergy),type,name));
  expr->init(expression);
  const std::vector<string>& parNames = expr->distinctParNames();
  std::vector<string> parStrings(parNames.size());
  for (size_t i=0; i<parStrings.size(); i++) {
    parStrings[i] = parNames[i] + " \" \" 1.0 -1.0e22 -1.0e22 1.0e22 1.0e22 0.01";
  }

  // we also need the information to set isSpecDependent

  bool isSpecDependent = expr->callsSpecDependentFunctions();

  // set up the function pointer

  XSCallBase* callExpr = (XSCallBase*)(new XSCall<MdefExpression>(expr.get()));
  expr.release();

  // store the component info in the cache. note that the mdef expression is stored in the infoString member
  
  ComponentInfo inCompInfo(name, type, false, expression, false, minEnergy, maxEnergy, parStrings);
  inCompInfo.isPythonModel(false);
  inCompInfo.isSpecDependent(isSpecDependent);
  inCompInfo.isMdefineModel(true);
  inCompInfo.functionPtr(callExpr);
  addCompInfoToCache(XSutility::lowerCase(name.substr(0,2)), inCompInfo);

  return true;
}

string XSModelFunction::delFromComponentList(const string& name)
{
  string compType(s_NOT_FOUND);
  // better be an exact match if we are going to remove a component
  if ( !isExactMatchName(name) ) return compType;

  ComponentInfo compInfo = compIteratorExactMatchName(name)->second;
  compType = compInfo.type();
  // reset the ComponentInfo object to clean up memory then erase
  // entry from the cache.
  compInfo.reset();
  s_compInfoCache.erase(compIteratorExactMatchName(name));
  return compType;
}


void XSModelFunction::printComponentList (std::ostream& s)
{
   using namespace std;
   // this code is not state of the art regarding flexibility for
   // future expansion but it is expected that adding model types will
   // be a rare occurrence.    
   static vector<string> addMods;
   static vector<string> mulMods;
   static vector<string> conMods;
   static vector<string> mixMods;
   static vector<string> acnMods;     
   static vector<string> amxMods;     
   static size_t lastCount = 0;
   size_t nCon = 0;
   size_t nMix = 0;
   size_t nAcn = 0;  
   size_t nAmx = 0;  

   size_t cacheSize = s_compInfoCache.size();

   if (cacheSize != lastCount)
   {
        // then iterate through the map and build the lists of different model types.
        CompInfoCacheType::const_iterator  comp = s_compInfoCache.begin();
        CompInfoCacheType::const_iterator  compEnd = s_compInfoCache.end();
        addMods.clear();
        mulMods.clear();
        conMods.clear();
        mixMods.clear();
        acnMods.clear();
        amxMods.clear();
        while ( comp != compEnd)
        {
                const ComponentInfo& current = comp->second;
                string name = current.name();
                if (!current.isStandardXspec())
                {
                   // Append local user components with a '*' or a '#' in printout.
		   // The former means a local model and the latter an mdefine model.
		  if (current.isMdefineModel()) {
		    name += "#";
		  } else {
		    name += "*";
		  }
                }
                if ( current.type() == "add")
                {
                        addMods.push_back(name);
                }
                else if ( current.type() == "mul" )
                {
                        mulMods.push_back(name);
                } 
                else if (current.type() == "con")
                {
                        conMods.push_back(name);
                }       
                else if (current.type() == "mix")
                {
                        mixMods.push_back(name);
                } 
                else if (current.type() == "acn")
                {
                        acnMods.push_back(name);
                }
                else if (current.type() == "amx")
                {
                        amxMods.push_back(name);
                } 
                ++comp;       
        }

        nCon = conMods.size();
        nMix = mixMods.size();
        nAcn = acnMods.size();             
        nAmx = amxMods.size();             
        sort(addMods.begin(),addMods.end());
        sort(mulMods.begin(),mulMods.end());
        if (nMix) sort(mixMods.begin(),mixMods.end());
        if (nCon) sort(conMods.begin(),conMods.end());
        if (nAcn) sort(acnMods.begin(),acnMods.end());
        if (nAmx) sort(amxMods.begin(),amxMods.end());
        lastCount = s_compInfoCache.size();
   }        


   const int LINE = 6;
   const int WIDTH = 12;

   s.setf(std::ios_base::left);
   s << " Additive Models: \n";

   XSutility::printStrings(addMods,s,LINE,WIDTH);

   s << "\n Multiplicative Models: \n";

   XSutility::printStrings(mulMods,s,LINE,WIDTH);

   if (conMods.size())
   {
        s << "\n Convolution Models: \n";
        XSutility::printStrings(conMods,s,LINE,WIDTH);
   }      

   if (mixMods.size())
   {
        s << "\n Mixing Models: \n";
        XSutility::printStrings(mixMods,s,LINE,WIDTH);
   }

   if (acnMods.size())
   {
        s << "\n Pile-up Models: \n";
        XSutility::printStrings(acnMods,s,LINE,WIDTH);
   }

   if (amxMods.size())
   {
        s << "\n Mixing pile-up Models: \n";
        XSutility::printStrings(amxMods,s,LINE,WIDTH);
   }

   s.setf(std::ios_base::right);

   s << "\n Table models may be used with the commands atable/mtable/etable";
   s << "\n\t atable{</path/to/tablemodel.mod>}";
   s << "\n and are described at:";
   s << "\n\t heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html";
   s << "\n";

   s << "\n Additional models are available at:";
   s << "\n\t heasarc.gsfc.nasa.gov/docs/xanadu/xspec/newmodels.html";
   s << endl; 
}

void XSModelFunction::printDebugComponentInfo (std::ostream& s)
{
  using namespace std;
  CompInfoCacheType::const_iterator  comp = s_compInfoCache.begin();
  CompInfoCacheType::const_iterator  compEnd = s_compInfoCache.end();
  while ( comp != compEnd ) {
    s << "-----" << comp->first << "-----" << "\n";
    s << comp->second.printInfo() << "\n";
    comp++;
  }
}

std::vector<ComponentInfo> XSModelFunction::compsMatchAbbrevName (const string& name)
{
  // entries in the name cache have two letter keys, so there is multiplicity.
  if (name.size() < 2)
  {
        throw XSparse::InputError("Model name must be at least two characters");         
  }
  const string abbrev = XSutility::lowerCase(name.substr(0,2));
  size_t n (s_compInfoCache.count(abbrev)); 
  if (n == 0) throw NoSuchComponent();
  std::vector<ComponentInfo> matches(n);
  std::pair<CompInfoCacheType::const_iterator,CompInfoCacheType::const_iterator> ff(s_compInfoCache.equal_range(abbrev));

  CompInfoCacheType::const_iterator m(ff.first);
  size_t index(0);
  while (m != ff.second)
  {
        matches[index] = m->second;
        ++m;
        ++index;       
  }
  return matches;  
}

ComponentInfo XSModelFunction::compMatchName (const string& fullName)
{
  // so far we know all of the matches start with the same two letters as 
  // fullname.

  std::vector<ComponentInfo> matches(compsMatchAbbrevName(fullName));
  size_t m(fullName.size());
  size_t n(matches.size());
  size_t i(0);
  while ( i < n )
  {
        if (m <= matches[i].name().size()) 
        {
	  if (matches[i].isMdefineModel())
	  {

	    // Special case for mdefine models.  Since mdefine models
	    // are not defined in models.dat in any proper order, we
	    // have to ensure an exact match rather than a lazy
	    // abreviation-allowed match.  Also, historically, mdefine
	    // models require an exact match and not an abbrev-match.
	    if (XSutility::lowerCase(fullName) ==
		XSutility::lowerCase(matches[i].name()))
	      return matches[i];
	  } else {
	        // Internal model which may be abbreviated
                // example: take bexriv as the intended model.
                // if the users type anything less than 'bexri' they
                // will get bexra, because both will be in the matched list.
                // but it should return bexrav correctly.
                // the test should guarantee a match of "bbody" against anything
                // shorter than "bbodyr".

                if (XSutility::lowerCase(fullName) == 
		    XSutility::lowerCase(matches[i].name().substr(0,m))) 
                        return matches[i];
	  }
        }
        ++i;
  }
  throw NoSuchComponent();    
}

std::vector<ComponentInfo> XSModelFunction::allComponentsInList()
{
  std::vector<ComponentInfo> allComps;

  // iterate through the cache
  CompInfoCacheType::const_iterator comp = s_compInfoCache.begin();
  CompInfoCacheType::const_iterator compEnd = s_compInfoCache.end();
  while ( comp != compEnd ) {
    const ComponentInfo& current = comp->second;
    allComps.push_back(current);
    comp++;
  }
  return allComps;
}

std::vector<ComponentInfo> XSModelFunction::allMdefineComponentsInList()
{
  std::vector<ComponentInfo> allComps;

  // iterate through the cache
  CompInfoCacheType::const_iterator comp = s_compInfoCache.begin();
  CompInfoCacheType::const_iterator compEnd = s_compInfoCache.end();
  while ( comp != compEnd ) {
    const ComponentInfo& current = comp->second;
    if ( current.isMdefineModel() ) allComps.push_back(current);
    comp++;
  }
  return allComps;
}

CompInfoCacheType::iterator XSModelFunction::compIteratorExactMatchName (const string& fullName)
{
  const string abbrev = XSutility::lowerCase(fullName.substr(0,2));
  std::pair<CompInfoCacheType::iterator,CompInfoCacheType::iterator> ff(s_compInfoCache.equal_range(abbrev));
  CompInfoCacheType::iterator match = s_compInfoCache.end();
  CompInfoCacheType::iterator it = ff.first;
  while (it != ff.second)
  {
     if (XSutility::lowerCase(fullName) == XSutility::lowerCase(it->second.name()))
     {
        match = it;
        break;
     }
     ++it;
  }
  return match;
}

bool XSModelFunction::isExactMatchName(const string& fullName)
{
  CompInfoCacheType::const_iterator itComp = compIteratorExactMatchName(fullName);
  if ( itComp == s_compInfoCache.end() ) return false;
  return true;
}

void XSModelFunction::saveMdefComponents(std::ostream& os)
{
  std::ostringstream mdefInfo;
  bool activelyFindingComponents = true;
  std::set<string> componentsAlreadyOutput;

  std::vector<ComponentInfo> mdefComponents(allMdefineComponentsInList());

  // Keep iterating until we stop finding mdefine components to output
  // This is also an infinite-loop-preventer.
  while (activelyFindingComponents)
  {
    activelyFindingComponents = false;

    // Loop through all Mdefine components and find components that either
    // have no dependencies, or have all dependencies written out already.
    for (size_t i=0; i<mdefComponents.size(); i++)
    {
      ComponentInfo current = mdefComponents[i];

      // Already been written once, don't write again
      if (componentsAlreadyOutput.count(current.name()) > 0)
	continue;

      XSCall<MdefExpression>* xsMdefExpr = (XSCall<MdefExpression>*)functionPointer(current.name());
      if (! xsMdefExpr ) continue; // shouldn't happen!!
      MdefExpression* mdefExpr = xsMdefExpr->generator();

      bool doOutput = true;  // should we output this entry?
      std::set<string>::const_iterator itDep = mdefExpr->usingOtherMdefs().begin();
      while (itDep != mdefExpr->usingOtherMdefs().end())
      {

	if ( componentsAlreadyOutput.count(*itDep) <= 0 )
	  doOutput = false; // dependency has not been output yet

	itDep ++;
      }

      if (doOutput)
      {
	mdefInfo << "mdefine " << current.name() << " "
		 << current.infoString() << " : " << current.type() <<"\n";
	// We are still actively finding new components to output
	activelyFindingComponents = true;
	// Add to list of components that have been output
	componentsAlreadyOutput.insert(current.name());
      }
    }
  }

  if (mdefInfo.str().length())
    os << mdefInfo.str() << std::flush;
}

void XSModelFunction::getMdefDeletionDependencies(const string& toDelete, 
        std::set<string>& affected)
{
    // Mdef 'toDelete' component should already have been checked to exist.
    XSCall<MdefExpression>* deleteComp = (XSCall<MdefExpression>*)functionPointer(toDelete);
    if (!deleteComp)
    {
       throw RedAlert("Programmer error in getMdefDeletionDependencies");
    }
    std::queue<string> namesToCheck;
    namesToCheck.push(toDelete);
    affected.insert(toDelete);   
        
    while (!namesToCheck.empty())
    {
       const string& nameToDelete = namesToCheck.front();
       string lcDelete(XSutility::lowerCase(nameToDelete));
       // loop over all mdefine components now defined
       std::vector<ComponentInfo> mdefComponents = allMdefineComponentsInList();
       for (size_t i=0; i<mdefComponents.size(); i++) {
	 string compName(mdefComponents[i].name());
	 // No need to check the nameToDelete component for its
	 //  dependencies.  We KNOW it should be deleted, and
	 //  it's already been inserted into 'affected' set.
	 if (lcDelete != XSutility::lowerCase(compName)) {
	   XSCall<MdefExpression>* testMdefComp = (XSCall<MdefExpression>*)functionPointer(compName);
	   MdefExpression* expr = testMdefComp->generator();
	   if (expr->usingOtherMdefs().find(lcDelete) != expr->usingOtherMdefs().end()) {
	     affected.insert(compName);
	     namesToCheck.push(compName);
	   }
	 }
       }
       namesToCheck.pop();
    }
}

size_t XSModelFunction::numberParameters(const string& fullName)
{
   CompInfoCacheType::const_iterator itComp = compIteratorExactMatchName(fullName);
   const ComponentInfo& current = itComp->second;
   const std::vector<string>& pInfo = current.parInfo();
   return pInfo.size();
}

void XSModelFunction::deleteAllFunctionPointers()
{
  CompInfoCacheType::iterator itFm = s_compInfoCache.begin();
  CompInfoCacheType::iterator itFmEnd = s_compInfoCache.end();
  while (itFm != itFmEnd) {
    itFm->second.deleteFunctionPtr();
    ++itFm;
  }
}

bool XSModelFunction::hasFunctionPointer(const string name)
{
  bool found = false;
  CompInfoCacheType::const_iterator itFm = compIteratorExactMatchName(name);
  if ( itFm != s_compInfoCache.end() ) {
    if ( itFm->second.functionPtr() != 0 ) found = true;
  }
  return found;
}

XSCallBase* XSModelFunction::functionPointer(const string name)
{
  CompInfoCacheType::const_iterator itFm = compIteratorExactMatchName(name);
  if ( itFm != s_compInfoCache.end() ) {
    ComponentInfo compInfo = itFm->second;
    return compInfo.functionPtr();
  } else {
    throw NoSuchComponent();
  }
}

void XSModelFunction::addFunctionPointer(const string name, XSCallBase* funcPointer)
{
  CompInfoCacheType::iterator itFm = compIteratorExactMatchName(name);
  if ( itFm != s_compInfoCache.end() ) {
    itFm->second.functionPtr(funcPointer);
  } else {
    throw NoSuchComponent();
  }
}

void XSModelFunction::deleteFunctionPointer(const string name)
{
  CompInfoCacheType::iterator itFm = compIteratorExactMatchName(name);
  if ( itFm != s_compInfoCache.end() ) {
    itFm->second.deleteFunctionPtr();
  } else {
    throw NoSuchComponent();
  }
}

