//   Read the documentation to learn more about C++ code generator
//   versioning.
//	  %X% %Q% %Z% %W%

#ifndef COMPONENTINFO_H
#define COMPONENTINFO_H 1
#include "xsTypes.h"
#include "XSCall.h"

struct ComponentInfo 
{
  public:
  ComponentInfo();
  ComponentInfo (const string& componentName, const string& componentType, bool errorFlag, string addString, bool isStandard, Real eMin, Real eMax, std::vector<string>& pInfo);

  const string printInfo() const;
  void reset ();

  void name(const string componentName);
  const string& name() const;
  void type(const string componentType);
  const string& type() const;
  void error(const bool hasError);
  const bool& error() const;
  void infoString(const string componentInfoString);
  const string& infoString() const;
  void isStandardXspec(const bool standard);
  const bool& isStandardXspec() const;
  void isPythonModel(const bool python);
  const bool& isPythonModel() const;
  void isSpecDependent(const bool dependent);
  const bool& isSpecDependent() const;
  void isMdefineModel(const bool mdefine);
  const bool& isMdefineModel() const;
  void minEnergy(const Real eMin);
  const Real& minEnergy() const;
  void maxEnergy(const Real eMax);
  const Real& maxEnergy() const;
  void parInfo(const std::vector<string>& pInfo);
  const std::vector<string>& parInfo() const;
  void functionPtr(XSCallBase* fPtr);
  XSCallBase* functionPtr() const;
  void deleteFunctionPtr();


  protected:
  private:
  private: //## implementation
  // Data Members for Class Attributes
  string m_name;
  string m_type;
  //	true if the model calculates errors.
  bool m_error;
  //	An additional information string that can be read from
  //	the
  //	"model.dat" file to supply anything that new model
  //	implementations might need.
  //
  //	Initial use might be a generic 'component that needs
  //	to read from a fits file' filename and directory
  //    If this is an mdefine component then contains the
  //    expression string
  string m_infoString;
  bool m_isStandardXspec;
  // Currently only Python-coded models are using the following 2 flags.
  //   For other models, isSpecDependent is determined when reading
  //   the model.dat file in Component::read().
  bool m_isPythonModel;
  bool m_isSpecDependent;
  // Flag for whether the component is from mdefine
  bool m_isMdefineModel;
  Real m_minEnergy;
  Real m_maxEnergy;
  // Vector of strings describing the parameters
  std::vector<string> m_parInfo;
  // Pointer to the template class to calculate the component
  XSCallBase* m_functionPtr;
};

typedef std::multimap<string,ComponentInfo> CompInfoCacheType;


// Class ComponentInfo 

inline void ComponentInfo::name(const string componentName)
{
  m_name = componentName;
}
inline const string& ComponentInfo::name() const
{
  return m_name;
}
inline void ComponentInfo::type(const string componentType)
{
  m_type = componentType;
}
inline const string& ComponentInfo::type() const
{
  return m_type;
}
inline void ComponentInfo::error(const bool hasError)
{
  m_error = hasError;
}
inline const bool& ComponentInfo::error() const
{
  return m_error;
}
inline void ComponentInfo::infoString(const string componentInfoString)
{
  m_infoString = componentInfoString;
}
inline const string& ComponentInfo::infoString() const
{
  return m_infoString;
}
inline void ComponentInfo::isStandardXspec(const bool standard)
{
  m_isStandardXspec = standard;
}
inline const bool& ComponentInfo::isStandardXspec() const
{
  return m_isStandardXspec;
}
inline void ComponentInfo::isPythonModel(const bool python)
{
  m_isPythonModel = python;
}
inline const bool& ComponentInfo::isPythonModel() const
{
  return m_isPythonModel;
}
inline void ComponentInfo::isSpecDependent(const bool dependent)
{
  m_isSpecDependent = dependent;
}
inline const bool& ComponentInfo::isSpecDependent() const
{
  return m_isSpecDependent;
}
inline void ComponentInfo::isMdefineModel(const bool mdefine)
{
  m_isMdefineModel = mdefine;
}
inline const bool& ComponentInfo::isMdefineModel() const
{
  return m_isMdefineModel;
}
inline void ComponentInfo::minEnergy(const Real eMin)
{
  m_minEnergy = eMin;
}
inline const Real& ComponentInfo::minEnergy() const
{
  return m_minEnergy;
}
inline void ComponentInfo::maxEnergy(const Real eMax)
{
  m_maxEnergy = eMax;
}
inline const Real& ComponentInfo::maxEnergy() const
{
  return m_maxEnergy;
}
inline void ComponentInfo::parInfo(const std::vector<string>& pInfo)
{
  m_parInfo = pInfo;
}
inline const std::vector<string>& ComponentInfo::parInfo() const
{
  return m_parInfo;
}
inline void ComponentInfo::functionPtr(XSCallBase* fPtr)
{
  m_functionPtr = fPtr;
}
inline XSCallBase* ComponentInfo::functionPtr() const
{
  return m_functionPtr;
}
inline void ComponentInfo::deleteFunctionPtr()
{
  delete m_functionPtr;
  m_functionPtr=0;
}


#endif
