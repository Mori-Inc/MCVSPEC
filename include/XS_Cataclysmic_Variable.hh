#pragma once

#include "Cataclysmic_Variable.hh"
#include <xsTypes.h>

class XS_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        XS_Cataclysmic_Variable(White_Dwarf, Accretion_Column);
        const void XS_Spectrum(const RealArray&, const int, RealArray&, const string&, const bool);
    protected:
        void Set_Abundances() override;
};
