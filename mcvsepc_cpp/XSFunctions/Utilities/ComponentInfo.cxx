//   Read the documentation to learn more about C++ code generator
//   versioning.
//	  %X% %Q% %Z% %W%

#include <sstream>

// ComponentInfo
#include <XSFunctions/Utilities/ComponentInfo.h>


// Class ComponentInfo 


ComponentInfo::ComponentInfo()
      : m_name(""),
        m_type("null"),
        m_error(false),
        m_infoString(""),
        m_isStandardXspec(false),
        m_isPythonModel(false),
        m_isSpecDependent(false),
        m_isMdefineModel(false),
	m_minEnergy(0.0),
	m_maxEnergy(1.0e10),
	m_parInfo(std::vector<string>()),
	m_functionPtr(0)
{
}

ComponentInfo::ComponentInfo (const string& componentName, const string& componentType, bool errorFlag, string addString, bool isStandard, Real eMin, Real eMax, std::vector<string>& pInfo)
      : m_name(componentName),
        m_type(componentType),
        m_error(errorFlag),
        m_infoString(addString),
        m_isStandardXspec(isStandard),
        m_isPythonModel(false),
        m_isSpecDependent(false),
	m_isMdefineModel(false),
	m_minEnergy(eMin),
	m_maxEnergy(eMax),
	m_parInfo(pInfo)
{
}

const string ComponentInfo::printInfo() const
{
  std::stringstream outStr;
  outStr << m_name << ": " << m_type << ", " << m_error << ", " << m_infoString
	 << ", " << m_isStandardXspec << ", " << m_isPythonModel << ", "
	 << m_isSpecDependent << ", " << m_isMdefineModel << ", " << m_minEnergy
	 << ", " << m_maxEnergy << "\n";
  for (size_t i=0; i<m_parInfo.size(); i++) outStr << m_parInfo[i] << "\n";
  outStr << "pointer = " << m_functionPtr << "\n";
  return outStr.str();
}

void ComponentInfo::reset ()
{
  m_name = "";
  m_type = "nul";    
  m_error = false;
  m_infoString = "";
  m_isStandardXspec = false;
  m_isPythonModel = false;
  m_isSpecDependent = false;
  m_isMdefineModel = false;
  m_minEnergy = 0.0;
  m_maxEnergy = 1.e10;
  m_parInfo = std::vector<string>();
  delete m_functionPtr;
  m_functionPtr = 0;
}
