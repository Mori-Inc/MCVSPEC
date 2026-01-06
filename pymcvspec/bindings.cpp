#include <cstddef>
#include <pybind11/buffer_info.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "dipole.hh"

#include <iostream>

namespace py = pybind11;

class Py_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        Py_Cataclysmic_Variable(double m, double r, double b, double mdot, double inv_r_m, double r_m_ratio, double area, double metals, double theta, double dist, int reflection):
            Cataclysmic_Variable(m,r,b,mdot,inv_r_m,r_m_ratio,area,metals,theta,dist,reflection)
        {
            Set_Abundances();
            Guess_Shock_Height();
            Shock_Height_Shooting();
            Build_Column_Profile();
        }
        void Set_Abundances() override{
            abundances.resize(atomic_charge.size());
            abundances = {1.00e+00, 9.77e-02, 3.63e-04, 1.12e-04, 8.51e-04, 1.23e-04,
                          3.80e-05, 2.95e-06, 3.55e-05, 1.62e-05, 3.63e-06, 2.29e-06,
                          4.68e-05, 1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
            for(uint i=2; i<abundances.size(); i++){
                abundances[i] *= metalicity;
            }
            abundances = abundances/abundances.sum();
            Set_Cooling_Constants();
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
        py::array_t<double> Get_Velocity(){
            return Valarray_to_Numpy(&velocity);
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
        py::array_t<double> Get_Density(){
            return Valarray_to_Numpy(&density);
        }
        py::array_t<double> Get_Total_Pressure(){
            return Valarray_to_Numpy(&total_pressure);
        }
        py::array_t<double> Get_Electron_Pressure(){
            return Valarray_to_Numpy(&electron_pressure);
        }
        py::array_t<double> Get_Volume(){
            return Valarray_to_Numpy(&volume);
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
        double Get_Avg_Atomic_Charge(){
            return avg_atomic_charge;
        }
};

PYBIND11_MODULE(_pymcvspec, module) {
    module.def("_mass_to_radius", &Cataclysmic_Variable::Get_Radius, "Returns the radius (cm) for a corresponding WD mass (g)");
    module.def("_luminosity_to_mdot", &Cataclysmic_Variable::Get_Accretion_Rate, "Returns the radius (cm) for a corresponding WD mass (g)");
    py::class_<Py_Cataclysmic_Variable>(module, "_cataclysmic_variable", py::module_local())
        .def(py::init<double,double,double,double,double,double,double,double,double,double,int>(),
            py::arg("mass") = 0.7*m_sol, py::arg("radius") = 0.01*r_sol, py::arg("b_field") = 1e7,
            py::arg("mdot") = 1e15, py::arg("inv_r_m") = 0., py::arg("r_m_ratio") = 1., py::arg("metalicity") = 1.,
            py::arg("area") = 1e15, py::arg("cos_incl_angle") = 0.5,
            py::arg("src_distance") = 200*pc_to_cm, py::arg("refl_on") = 1)
        .def("set_shock_height", &Py_Cataclysmic_Variable::Update_Shock_Height)
        .def("get_landing", &Py_Cataclysmic_Variable::Get_Landing_Altitude)
        .def("get_altitude", &Py_Cataclysmic_Variable::Get_Altitude)
        .def("get_velocity", &Py_Cataclysmic_Variable::Get_Velocity)
        .def("get_electron_temperature", &Py_Cataclysmic_Variable::Get_Electron_Temperature)
        .def("get_ion_temperature", &Py_Cataclysmic_Variable::Get_Ion_Temperature)
        .def("get_electron_density", &Py_Cataclysmic_Variable::Get_Electron_Density)
        .def("get_density", &Py_Cataclysmic_Variable::Get_Density)
        .def("get_electron_pressure", &Py_Cataclysmic_Variable::Get_Electron_Pressure)
        .def("get_total_pressure", &Py_Cataclysmic_Variable::Get_Total_Pressure)
        .def("get_volume", &Py_Cataclysmic_Variable::Get_Volume)
        .def("get_radius", &Py_Cataclysmic_Variable::Get_Radius)
        .def("get_m_dot", &Py_Cataclysmic_Variable::Get_Accretion_Rate)
        .def("get_shock_height", &Py_Cataclysmic_Variable::Get_Shock_Height)
        .def("get_avg_charge", &Py_Cataclysmic_Variable::Get_Avg_Atomic_Charge)
        .def("print", &Py_Cataclysmic_Variable::Print_Properties);

    py::class_<Dipole>(module, "_dipole", py::module_local())
        .def(py::init<double>(),py::arg("u")=0.03)
        .def("update_geo", [](Dipole& self, double w, double& r, double& dr_dw, double& proj_r_w, double& convergance, py::array_t<double> metric){
            py::buffer_info buf = metric.request();
            double* ptr = static_cast<double*>(buf.ptr);
            self.update_coordinates(w, r, dr_dw, proj_r_w, convergance, ptr);
            return py::make_tuple(r, dr_dw, proj_r_w, convergance);
        });
}
