#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Cataclysmic_Variable.hh"
#include <iostream>
using std::cout;
using std::endl;

namespace py = pybind11;

class Py_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        Py_Cataclysmic_Variable(double m, double b, double metals, double luminosity, double fractional_area, double theta, double dist, int reflection):
            Cataclysmic_Variable(m,b,metals,fractional_area,theta,dist,reflection)
        {
            Radius_Shooting(100000);
            inverse_mag_radius = 0;
            Set_Accretion_Rate(luminosity);
            accretion_area = fractional_area*4.*pi*radius*radius;
            accretion_rate /= accretion_area;
            Set_Abundances(metalicity);
            Set_Pre_Shock_Speed(5);
            Set_Cooling_Ratio();
        }
        void Set_Abundances(double) override{
            abundances.resize(atomic_charge.size());
            abundances = {1.00e+00, 9.77e-02, 1.45e-11, 1.41e-11, 3.98e-10, 3.63e-04,
                1.12e-04, 8.51e-04, 3.63e-08, 1.23e-04, 2.14e-06, 3.80e-05,
                2.95e-06, 3.55e-05, 2.82e-07, 1.62e-05, 3.16e-07, 3.63e-06,
                1.32e-07, 2.29e-06, 1.26e-09, 9.77e-08, 1.00e-08, 4.68e-07,
                2.45e-07, 4.68e-05, 8.32e-08, 1.78e-06, 1.62e-08, 3.98e-08};
            abundances = abundances/abundances.sum();
            Set_Cooling_Constants();
        }
        void MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string) override{
            std::cout << "Spectrum creation is not yet implented in python" << std::endl;
        }
};

PYBIND11_MODULE(mcvspec, module) {
    py::class_<Py_Cataclysmic_Variable>(module, "Polar")
        .def(py::init<double,double,double,double,double,double,double,int>())
        .def("execute", &Py_Cataclysmic_Variable::Shock_Height_Shooting);
}
