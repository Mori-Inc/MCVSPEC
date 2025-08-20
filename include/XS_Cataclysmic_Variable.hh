#pragma once

#include "Cataclysmic_Variable.hh"
#include <xsTypes.h>

class XS_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        XS_Cataclysmic_Variable(double,double,double,double,double,double,double,double,double,double,double,int);
        void XS_Spectrum(const RealArray&, const int, RealArray&, const string&);
    protected:
        void Set_Abundances(double) override;
};
