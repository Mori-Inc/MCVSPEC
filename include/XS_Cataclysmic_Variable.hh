#pragma once

#include "Cataclysmic_Variable.hh"
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <xsTypes.h>
#include <funcWrappers.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

class XS_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        XS_Cataclysmic_Variable(double,double,double,double,double,double,double,double,double,double,double,int);
        void XS_Spectrum(const RealArray&, const int, RealArray&, const string&);
    protected:
        void Set_Abundances(double) override;
};
