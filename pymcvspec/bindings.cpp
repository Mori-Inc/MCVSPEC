#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include <iostream>

namespace py = pybind11;

class Py_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        Py_Cataclysmic_Variable(double m, double r, double b, double mdot, double inv_r_m, double metals, double area, double theta, double n, double dist, int reflection):
            Cataclysmic_Variable(m,r,b,mdot,inv_r_m,area,theta,n,dist,reflection)
        {
            metalicity = metals;
            Set_Abundances(metals);
            Guess_Shock_Height();
            Shock_Height_Shooting();
            Build_Column_Profile();
        }
        void Set_Abundances(const double m) override{
            abundances.resize(atomic_charge.size());
            abundances = {1.00e+00, 9.77e-02, m*3.63e-04, m*1.12e-04, m*8.51e-04, m*1.23e-04,
                          m*3.80e-05, m*2.95e-06, m*3.55e-05, m*1.62e-05, m*3.63e-06, m*2.29e-06,
                          m*4.68e-05, m*1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
            abundances = abundances/abundances.sum();
            Set_Cooling_Constants();
        }
        void Set_Pressure_Ratio(double sigma_s){
            pressure_ratio = sigma_s;
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

PYBIND11_MODULE(_pymcvspec, module) {
    module.def("_mass_to_radius", &Cataclysmic_Variable::Get_Radius, "Returns the radius (cm) for a corresponding WD mass (g)");
    module.def("_luminosity_to_mdot", &Cataclysmic_Variable::Get_Accretion_Rate, "Returns the radius (cm) for a corresponding WD mass (g)");
    py::class_<Py_Cataclysmic_Variable>(module, "_cataclysmic_variable", py::module_local())
        .def(py::init<double,double,double,double,double,double,double,double,double,double,int>(),
            py::arg("mass") = 0.7*m_sol, py::arg("radius") = 0.01*r_sol, py::arg("b_field") = 1e7,
            py::arg("mdot") = 1e15, py::arg("inv_r_m") = 0., py::arg("metalicity") = 1.,
            py::arg("area") = 1e15, py::arg("cos_incl_angle") = 0.5, py::arg("area_exp") = 0,
            py::arg("src_distance") = 200*pc_to_cm, py::arg("refl_on") = 1)
        .def("flow_eq", &Py_Cataclysmic_Variable::Flow_Equation)
        .def("set_pressure_ratio", &Py_Cataclysmic_Variable::Set_Pressure_Ratio)
        .def("get_altitude", &Py_Cataclysmic_Variable::Get_Altitude)
        .def("get_electron_temperature", &Py_Cataclysmic_Variable::Get_Electron_Temperature)
        .def("get_ion_temperature", &Py_Cataclysmic_Variable::Get_Ion_Temperature)
        .def("get_electron_density", &Py_Cataclysmic_Variable::Get_Electron_Density)
        .def("get_ion_density", &Py_Cataclysmic_Variable::Get_Ion_Density)
        .def("get_electron_pressure", &Py_Cataclysmic_Variable::Get_Electron_Pressure)
        .def("get_total_pressure", &Py_Cataclysmic_Variable::Get_Total_Pressure)
        .def("get_radius", &Py_Cataclysmic_Variable::Get_Radius)
        .def("get_m_dot", &Py_Cataclysmic_Variable::Get_Accretion_Rate)
        .def("get_shock_height", &Py_Cataclysmic_Variable::Get_Shock_Height)
        .def("print", &Py_Cataclysmic_Variable::Print_Properties);
}
