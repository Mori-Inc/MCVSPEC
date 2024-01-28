//	Interface for calling model calculation functions

#ifndef XSCALL_H
#define XSCALL_H 1

#include <xsTypes.h>
#include <XSUtil/Error/Error.h>
#include <XSFunctions/Utilities/MdefExpression.h>

// Virtual base class for the template class to inherit from.

class XSCallBase
{
public:

  virtual ~XSCallBase() = 0;

  virtual void operator () (const RealArray& energyArray, const RealArray& params, int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, const string& initString = string()) const = 0;
  virtual void operator () (const EnergyPointer& energyArray, const std::vector<Real>& parameterValues, GroupFluxContainer& fluxArrays, GroupFluxContainer& fluxErrArrays, MixUtility* mixGenerator = 0, const string& modelName = string()) const = 0;
  virtual MixUtility* getUtilityObject() const = 0;

protected:
  XSCallBase();
  XSCallBase(const XSCallBase &right);

private:
  XSCallBase& operator=(const XSCallBase &right);
};

inline XSCallBase::XSCallBase()
{
}

inline XSCallBase::XSCallBase(const XSCallBase &right)
{
}

inline XSCallBase::~XSCallBase()
{
}

template <typename T>
class XSCall : public XSCallBase  //## Inherits: <unnamed>%3E319B420391
{
  public:
      XSCall(const XSCall< T > &right);
      XSCall (T* generator);
      virtual ~XSCall();

      virtual void operator () (const RealArray& energyArray, const RealArray& params, int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, const string& initString = string()) const;
      virtual void operator () (const EnergyPointer& energyArray, const std::vector<Real>& parameterValues, GroupFluxContainer& fluxArrays, GroupFluxContainer& fluxErrArrays, MixUtility* mixGenerator = 0, const string& modelName = string()) const;
      // More complicated models (ie. mix) may need to create MixUtility objects.  
      // For everything else, the default implementation of this function returns NULL.
      // Note that the MixUtility type isn't defined until a higher level library.
      virtual MixUtility* getUtilityObject() const;
      T* generator () const;
      void generator (T* value);

  protected:
  private:
      XSCall();
      XSCall< T > & operator=(const XSCall< T > &right);

  private: //## implementation
    // Data Members for Class Attributes
      T* m_generator;

};

// Parameterized Class XSCall 

template <typename T>
MixUtility* XSCall<T>::getUtilityObject() const
{
   return 0;
}

template <typename T>
inline T* XSCall<T>::generator () const
{
  return m_generator;
}

template <typename T>
inline void XSCall<T>::generator (T* value)
{
  m_generator = value;
}

// Parameterized Class XSCall 

template <typename T>
XSCall<T>::XSCall(const XSCall<T> &right)
  : XSCallBase(right),
    m_generator(right.m_generator) // default non-owning - shallow copy
                                   // (assume it's a function pointer)
{
}

template <typename T>
XSCall<T>::XSCall (T* generator)
  : m_generator(generator)
{
}


template <typename T>
XSCall<T>::~XSCall()
{
}


template <typename T>
void XSCall<T>::operator () (const RealArray& energyArray, const RealArray& params, int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, const string& initString) const
{
   throw RedAlert("Model function is missing due to specialized template build error (1).");
}

template <typename T>
void XSCall<T>::operator () (const EnergyPointer& energyArray, const std::vector<Real>& parameterValues, GroupFluxContainer& fluxArrays, GroupFluxContainer& fluxErrArrays, MixUtility* mixGenerator, const string& modelName) const
{
   throw RedAlert("Model function is missing due to specialized template build error (2).");
}

template  <> 
void XSCall<xsf77Call>::operator() (const RealArray& energyArray, const RealArray& params,
                        int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, 
                        const string& initString) const;

template  <> 
void XSCall<xsF77Call>::operator() (const RealArray& energyArray, const RealArray& params,
                        int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, 
                        const string& initString) const;

template  <> 
void XSCall<xsccCall>::operator() (const RealArray& energyArray, const RealArray& params,
                        int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, 
                        const string& initString) const;


template  <> 
void XSCall<XSCCall>::operator() (const RealArray& energyArray, const RealArray& params,
                        int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, 
                        const string& initString) const;

template  <> 
void XSCall<XSMixCCall>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator, 
                        const string& modelName) const;

template  <> 
void XSCall<xsmixcall>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator, 
                        const string& modelName) const;

template <>
void XSCall<MdefExpression>::operator() (const RealArray& energyArray, const RealArray& params,
                        int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, 
                        const string& initString) const;

template <>
XSCall<MdefExpression>::~XSCall();

template <>
XSCall<MdefExpression>::XSCall(const XSCall<MdefExpression> &right);


#endif
