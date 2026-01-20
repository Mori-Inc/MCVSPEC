#include <pybind11/buffer_info.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "dipole.hh"

#include <vector>

using std::memcpy;
using std::transform;

namespace py = pybind11;

template <typename T>
vector<T> Numpy_to_Vector(py::array_t<T, py::array::c_style | py::array::forcecast> np_array){
    py::buffer_info array_info = np_array.request();
    if (array_info.ndim != 1){
        throw std::runtime_error("Only Handles 1D Numpy Arrays");
    }
    const size_t arr_len =  static_cast<size_t>(array_info.shape[0]);
    const T* source = static_cast<const T*>(array_info.ptr);
    vector<T> cpp_vec(arr_len);
    memcpy(cpp_vec.data(), source, arr_len*sizeof(T));
    return cpp_vec;
}

template <typename T>
py::array_t<T> Vector_to_Numpy(const std::vector<T>& cpp_vec) {
    py::array_t<T> np_array(static_cast<py::size_t>(cpp_vec.size()));
    py::buffer_info array_info = np_array.request();
    T* array_pointer = static_cast<T*>(array_info.ptr);
    memcpy(array_pointer, cpp_vec.data(), cpp_vec.size()*sizeof(T));
    return np_array;
}

class Py_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        Py_Cataclysmic_Variable(double m, double r, double b, double mdot, double area, double inv_r_m, double r_m_ratio, double metals, double theta, double dist, int reflection):
            Cataclysmic_Variable(m,r,b,mdot,area,inv_r_m,r_m_ratio,metals,theta,dist,reflection)
        {
            Set_Abundances();
        }
        void Set_Abundances() override{
            abundances.resize(atomic_charge.size());
            abundances = {1.00e+00, 9.77e-02, 3.63e-04, 1.12e-04, 8.51e-04, 1.23e-04,
                          3.80e-05, 2.95e-06, 3.55e-05, 1.62e-05, 3.63e-06, 2.29e-06,
                          4.68e-05, 1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
            double total=0;
            for(uint i=2; i<abundances.size(); i++){
                abundances[i] *= metalicity;
                total += abundances[i];
            }
            transform(abundances.begin(),abundances.end(),abundances.begin(),[total](double x) {return x/total;});
            Set_Cooling_Constants();
        }
        py::array_t<double> Get_Altitude(){
            return Vector_to_Numpy(altitude);
        }
        py::array_t<double> Get_Velocity(){
            return Vector_to_Numpy(velocity);
        }
        py::array_t<double> Get_Electron_Temperature(){
            return Vector_to_Numpy(electron_temperature);
        }
        py::array_t<double> Get_Ion_Temperature(){
            return Vector_to_Numpy(ion_temperature);
        }
        py::array_t<double> Get_Electron_Density(){
            return Vector_to_Numpy(electron_density);
        }
        py::array_t<double> Get_Density(){
            return Vector_to_Numpy(density);
        }
        py::array_t<double> Get_Total_Pressure(){
            return Vector_to_Numpy(total_pressure);
        }
        py::array_t<double> Get_Electron_Pressure(){
            return Vector_to_Numpy(electron_pressure);
        }
        py::array_t<double> Get_Volume(){
            return Vector_to_Numpy(volume);
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
            py::arg("mdot") = 1e15, py::arg("area") = 1e15, py::arg("inv_r_m") = 0., py::arg("r_m_ratio") = 1.,
            py::arg("metalicity") = 1., py::arg("cos_incl_angle") = 0.5,
            py::arg("src_distance") = 200*pc_to_cm, py::arg("refl_on") = 1)
        .def("set_shock_height", &Py_Cataclysmic_Variable::Update_Shock_Position)
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
        .def("print", &Py_Cataclysmic_Variable::Print_Properties)
        .def("flow_equation", [](Py_Cataclysmic_Variable& self, double s, py::array_t<double, py::array::c_style | py::array::forcecast> state_py){
            vector<double> state = Numpy_to_Vector(state_py);
            vector<double> deriv(state.size());
            self.Flow_Equation(s, state, deriv);
            return Vector_to_Numpy(deriv);
        });

    py::class_<Dipole>(module, "_dipole", py::module_local())
        .def(py::init<double>(),py::arg("u")=0.03)
        .def("update_geo", [](Dipole& self, double w){
            double r, proj_r_w, convergance, metric[3];
            self.update_coordinates(w, r, proj_r_w, convergance, metric);
            py::array_t<double> array(3);
            py::detail::unchecked_mutable_reference<double, 1> np_array = array.mutable_unchecked<1>();
            for(int i = 0; i < 3; i++){
                np_array(i) = metric[i];
            }
            return py::make_tuple(r, proj_r_w, convergance, array);
        });
}
