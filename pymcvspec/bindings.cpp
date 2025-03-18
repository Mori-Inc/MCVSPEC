#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
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
            Set_Cooling_Ratio();;
        }
        void Set_Abundances(const double m) override{
            abundances.resize(atomic_charge.size());
            abundances = {1.00e+00, 9.77e-02, m*3.63e-04, m*1.12e-04, m*8.51e-04, m*1.23e-04,
                          m*3.80e-05, m*2.95e-06, m*3.55e-05, m*1.62e-05, m*3.63e-06, m*2.29e-06,
                          m*4.68e-05, m*1.78e-06};
            abundances = abundances/abundances.sum();
            Set_Cooling_Constants();
        }
        void MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string) override{
            std::cout << "Spectrum creation is not yet implented in python" << std::endl;
        }
        py::array_t<double> Get_Altitude(){
            py::array_t<double> array(altitude.size());
            py::detail::unchecked_mutable_reference<double, 1> np_array = array.mutable_unchecked<1>();
            for(int i = 0; i < altitude.size(); i++){
                np_array(i) = altitude[i];
            }
            return array;
        }
        py::array_t<double> Get_Electron_Temperature(){
            py::array_t<double> array(electron_temperature.size());
            py::detail::unchecked_mutable_reference<double, 1> np_array = array.mutable_unchecked<1>();
            for(int i = 0; i < electron_temperature.size(); i++){
                np_array(i) = electron_temperature[i];
            }
            return array;
        }
        py::array_t<double> Get_Ion_Temperature(){
            py::array_t<double> array(ion_temperature.size());
            py::detail::unchecked_mutable_reference<double, 1> np_array = array.mutable_unchecked<1>();
            for(int i = 0; i < ion_temperature.size(); i++){
                np_array(i) = ion_temperature[i];
            }
            return array;
        }
        py::array_t<double> Get_Electron_Density(){
            py::array_t<double> array(electron_density.size());
            py::detail::unchecked_mutable_reference<double, 1> np_array = array.mutable_unchecked<1>();
            for(int i = 0; i < electron_density.size(); i++){
                np_array(i) = electron_density[i];
            }
            return array;
        }
        py::array_t<double> Get_Ion_Density(){
            py::array_t<double> array(ion_density.size());
            py::detail::unchecked_mutable_reference<double, 1> np_array = array.mutable_unchecked<1>();
            for(int i = 0; i < ion_density.size(); i++){
                np_array(i) = ion_density[i];
            }
            return array;
        }
};

PYBIND11_MODULE(mcvspec, module) {
    py::class_<Py_Cataclysmic_Variable>(module, "Polar")
        .def(py::init<double,double,double,double,double,double,double,int>())
        .def("execute", &Py_Cataclysmic_Variable::Shock_Height_Shooting)
        .def("altitude", &Py_Cataclysmic_Variable::Get_Altitude)
        .def("electron_temperature", &Py_Cataclysmic_Variable::Get_Electron_Temperature)
        .def("ion_temperature", &Py_Cataclysmic_Variable::Get_Ion_Temperature)
        .def("electron_density", &Py_Cataclysmic_Variable::Get_Electron_Density)
        .def("ion_density", &Py_Cataclysmic_Variable::Get_Ion_Density);
}
