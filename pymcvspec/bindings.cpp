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
            Set_Radius();
            inverse_mag_radius = 0;
            Set_Accretion_Rate(luminosity);
            accretion_area = fractional_area*4.*pi*radius*radius;
            accretion_rate /= accretion_area;
            Set_Abundances(metalicity);
            Set_Pre_Shock_Speed(5);
            Set_Cooling_Ratio();
            cout << bremss_constant << endl;
            cout << cooling_ratio << endl;
        }
        Py_Cataclysmic_Variable(double m, double metals, double luminosity, double fractional_area, double theta, double dist, int reflection, double r_m):
            Cataclysmic_Variable(m,r_m,metals,fractional_area,theta,dist,reflection)
        {
            inverse_mag_radius = 1/r_m;
            Set_Radius();
            Set_Accretion_Rate(luminosity);
            b_field = sqrt(32.*accretion_rate)*pow(grav_const*mass, 1./4.)*pow(r_m,7./4.)*pow(radius,-3);
            accretion_area = fractional_area*4.*pi*radius*radius;
            accretion_rate /= accretion_area;
            Set_Abundances(metalicity);
            Set_Pre_Shock_Speed(5);
            Set_Cooling_Ratio();
        }
        void Set_Abundances(const double m) override{
            abundances.resize(atomic_charge.size());
            abundances = {1.00e+00, 9.77e-02, m*3.63e-04, m*1.12e-04, m*8.51e-04, m*1.23e-04,
                          m*3.80e-05, m*2.95e-06, m*3.55e-05, m*1.62e-05, m*3.63e-06, m*2.29e-06,
                          m*4.68e-05, m*1.78e-06};
            abundances = abundances/abundances.sum();
            Set_Cooling_Constants();
        }
        void Set_Pressure_Ratio(double sigma_s){
            pressure_ratio = sigma_s;
        }
        void MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string) override{
            std::cout << "Spectrum creation is not yet implented in python" << std::endl;
        }
        py::array_t<double> Valarray_to_Numpy(valarray<double>* arr){
            py::array_t<double> array(arr->size());
            py::detail::unchecked_mutable_reference<double, 1> np_array = array.mutable_unchecked<1>();
            for(int i = 0; i < arr->size(); i++){
                np_array(i) = (*arr)[i];
            }
            return array;
        }
        py::array_t<double> Get_Altitude(){
            return Valarray_to_Numpy(&altitude);
        }
        py::array_t<double> Get_Electron_Temperature(){
            return Valarray_to_Numpy(&electron_temperature);
        }
        py::array_t<double> Get_Ion_Temperature(){
            return Valarray_to_Numpy(&ion_temperature);
        }
        py::array_t<double> Get_Electron_Density(){
            return Valarray_to_Numpy(&electron_density);
        }
        py::array_t<double> Get_Ion_Density(){
            return Valarray_to_Numpy(&ion_density);
        }
        py::array_t<double> Get_Total_Pressure(){
            return Valarray_to_Numpy(&total_pressure);
        }
        py::array_t<double> Get_Electron_Pressure(){
            return Valarray_to_Numpy(&electron_pressure);
        }
        double Get_Radius(){
            return radius;
        }
        double Get_Accretion_Rate(){
            return accretion_rate;
        }
        double Get_Shock_Height(){
            return shock_height;
        }
};

PYBIND11_MODULE(mcvspec, module) {
    py::class_<Py_Cataclysmic_Variable>(module, "Polar", py::module_local())
        .def(py::init<double,double,double,double,double,double,double,int>(),
            py::arg("mass") = 0.7*solar_mass, py::arg("b_field") = 1e7,
            py::arg("metalicity") = 1., py::arg("luminosity") = 1e33,
            py::arg("area_frac") = 1e-4, py::arg("cos_incl_angle") = 0.5,
            py::arg("src_distance") = 200*pc_to_cm, py::arg("refl_on") = 1)
        .def("flow_eq", &Py_Cataclysmic_Variable::Flow_Equation)
        .def("execute", &Py_Cataclysmic_Variable::Shock_Height_Shooting)
        .def("set_pressure_ratio", &Py_Cataclysmic_Variable::Set_Pressure_Ratio)
        .def("altitude", &Py_Cataclysmic_Variable::Get_Altitude)
        .def("electron_temperature", &Py_Cataclysmic_Variable::Get_Electron_Temperature)
        .def("ion_temperature", &Py_Cataclysmic_Variable::Get_Ion_Temperature)
        .def("electron_density", &Py_Cataclysmic_Variable::Get_Electron_Density)
        .def("ion_density", &Py_Cataclysmic_Variable::Get_Ion_Density)
        .def("electron_pressure", &Py_Cataclysmic_Variable::Get_Electron_Pressure)
        .def("total_pressure", &Py_Cataclysmic_Variable::Get_Total_Pressure)
        .def("radius", &Py_Cataclysmic_Variable::Get_Radius)
        .def("m_dot", &Py_Cataclysmic_Variable::Get_Accretion_Rate)
        .def("shock_height", &Py_Cataclysmic_Variable::Get_Shock_Height);
    py::class_<Py_Cataclysmic_Variable>(module, "Intermediate_Polar")
        .def(py::init<double,double,double,double,double,double,int,double>(),
            py::arg("mass") = 0.7*solar_mass, py::arg("metalicity") = 1.,
            py::arg("luminosity") = 1e33, py::arg("area_frac") = 1e-4,
            py::arg("cos_incl_angle") = 0.5, py::arg("src_distance") = 200*pc_to_cm,
            py::arg("refl_on") = 1, py::arg("mag_radius") = 2*solar_radius)
        .def("execute", &Py_Cataclysmic_Variable::Shock_Height_Shooting)
        .def("set_pressure_ratio", &Py_Cataclysmic_Variable::Set_Pressure_Ratio)
        .def("altitude", &Py_Cataclysmic_Variable::Get_Altitude)
        .def("electron_temperature", &Py_Cataclysmic_Variable::Get_Electron_Temperature)
        .def("ion_temperature", &Py_Cataclysmic_Variable::Get_Ion_Temperature)
        .def("electron_density", &Py_Cataclysmic_Variable::Get_Electron_Density)
        .def("ion_density", &Py_Cataclysmic_Variable::Get_Ion_Density)
        .def("electron_pressure", &Py_Cataclysmic_Variable::Get_Electron_Pressure)
        .def("total_pressure", &Py_Cataclysmic_Variable::Get_Total_Pressure)
        .def("radius", &Py_Cataclysmic_Variable::Get_Radius)
        .def("m_dot", &Py_Cataclysmic_Variable::Get_Accretion_Rate)
        .def("shock_height", &Py_Cataclysmic_Variable::Get_Shock_Height);
}
