//   Read the documentation to learn more about C++ code generator
//   versioning.
//	  %X% %Q% %Z% %W%

#ifndef XSMODELFUNCTION_H
#define XSMODELFUNCTION_H 1

// xsTypes
#include <xsTypes.h>
// Error
#include <XSUtil/Error/Error.h>
// ComponentInfo
#include <XSFunctions/Utilities/ComponentInfo.h>
//#ifdef INITPACKAGE
//#       include <funcType.h>
//#else
#       include <XSFunctions/Utilities/funcType.h>
//#endif

#include <XSUtil/Parse/XSparse.h>  // invalid exception
#include <XSFunctions/Utilities/MdefExpression.h>
#include <map>

class MixUtility;

// This class consists entirely of static functions and data members.
//  No instances of this are generated.

class XSModelFunction 
{
  public:

      XSModelFunction() = delete;
      XSModelFunction(const XSModelFunction &right) = delete;
      XSModelFunction & operator=(const XSModelFunction &right) = delete;

    class NoSuchComponent : public YellowAlert  //## Inherits: <unnamed>%3FBA6F620269
    {
      public:
          NoSuchComponent();

      protected:
      private:
      private: //## implementation
    };

      // add to the component cache the models listed in modelDatFile. this should be followed by separate
      // calls to addFunctionPointer for each model (see eg functionMap.cxx).
      static void updateComponentList (const string& modelDatFile, bool isStandard = false);

      // add a ComponentInfo object to the cache. shortName is the two character abbreviation.
      static void addCompInfoToCache (const string& shortName, const ComponentInfo& value);

      // add a new component from an mdefine command
      static bool addMdefToComponentList(const string& name, const string& type,
				     const string& expression, const Real& minEnergy,
				     const Real& maxEnergy);

      // delete a component from the cache. the name must be an exact match
      static string delFromComponentList(const string& name);

      // print a user-friendly list of the available components
      static void printComponentList (std::ostream& s);

      // print lots of debug output on the available components
      static void printDebugComponentInfo (std::ostream& s);

      // name is the name of the component typed by the user. The vector of ComponentInfo objects
      // returned are all those with names whose first two characters match those in name.
      static std::vector<ComponentInfo> compsMatchAbbrevName (const string& name);

      // compMatchName returns the single ComponentInfo which matches to fullName by the abbreviation rules
      static ComponentInfo compMatchName (const string& fullName);

      // compIteratorExactMatchName returns an iterator to the single ComponentInfo entry in the cache
      // which exactly matches to fullName.
      static CompInfoCacheType::iterator compIteratorExactMatchName (const string& fullName);

      // return true if there is an exact match to a component in the cache
      static bool isExactMatchName(const string& fullName);

      // A vector of all the ComponentInfo currently defined
      static std::vector<ComponentInfo> allComponentsInList();

      // A vector of all the ComponentInfo from mdefine currently defined
      static std::vector<ComponentInfo> allMdefineComponentsInList();

      // The number of parameters of the component which matches to fullName by the abbreviation rules
      static size_t numberParameters(const string& fullName);

      // writes the commands to make the mdefine'd components for use in the save command. Note that
      // this can be complicated if mdefine'd components depend on others because the commands must
      // then appear in the correct order
      static void saveMdefComponents(std::ostream& os);

      // if an mdefined component is to be deleted then check for any other mdefined components which
      // depend on it
      static void getMdefDeletionDependencies(const string& toDelete,
		  std::set<string>& affected);

      // delete the function pointers for all the components
      static void deleteAllFunctionPointers();

      // check whether the component has a non-zero function pointer
      static bool hasFunctionPointer(const string name);

      // return the function pointer for the component
      static XSCallBase* functionPointer(const string name);

      // set the function pointer for the named component, which must exist
      static void addFunctionPointer(const string name, XSCallBase* funcPointer);

      // delete the function pointer for the named component and set it to zero
      static void deleteFunctionPointer(const string name);

      // return the string used to indicate that a component was not found.
      static const string& NOT_FOUND();

  public:
  protected:

  private: //## implementation
    // Data Members for Associations
      // This contains information on all the model components defined, those read
      // from model.dat and mixmodel.dat, user models, and mdefined models.
      static CompInfoCacheType s_compInfoCache;
      
      static const string s_NOT_FOUND;

};
// Class XSModelFunction::NoSuchComponent 

// Class XSModelFunction 

inline void XSModelFunction::addCompInfoToCache (const string& shortName, const ComponentInfo& value)
{
  s_compInfoCache.insert(CompInfoCacheType::value_type(shortName,value));
}

inline const string& XSModelFunction::NOT_FOUND()
{
  return s_NOT_FOUND;
}


#endif
