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
static py::array_t<T> Vector_to_Numpy(const std::vector<T>& cpp_vec) {
  py::array_t<T> np_array(cpp_vec.size());
  memcpy(np_array.mutable_data(), cpp_vec.data(), cpp_vec.size()*sizeof(T));
  return np_array;
}

static std::vector<double> Numpy_to_Vector(py::handle obj) {
    py::array array = py::array::ensure(obj);
    if (!array) {
        throw py::type_error("Input must be array-like.");
    }
    py::array_t<double, py::array::c_style | py::array::forcecast> np_array(array);

    if (np_array.ndim() != 1) {
        throw py::value_error("Input must be 1-d");
    }

    vector<double> cpp_vec((size_t)np_array.shape(0));
    memcpy(cpp_vec.data(), np_array.data(), cpp_vec.size()*sizeof(double));
    return cpp_vec;
}

class Py_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        Py_Cataclysmic_Variable(double m, double r, double b, double mdot, double area, double inv_r_m, double r_m_ratio, double metals, double theta, double pressure_ratio, double u, double dist, int reflection):
            Cataclysmic_Variable(m,r,b,mdot,area,inv_r_m,r_m_ratio,metals,theta,pressure_ratio,u,dist,reflection)
        {
            Set_Abundances();
        }

        void Set_Abundances() override{
            abundances.resize(atomic_charge.size());
            abundances = {1.00e+00, 9.77e-02, 3.63e-04, 1.12e-04, 8.51e-04, 1.23e-04,
                          3.80e-05, 2.95e-06, 3.55e-05, 1.62e-05, 3.63e-06, 2.29e-06,
                          4.68e-05, 1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
            double total= abundances[0]+abundances[1];
            for(uint i=2; i<abundances.size(); i++){
                abundances[i] *= metalicity;
                total += abundances[i];
            }
            transform(abundances.begin(),abundances.end(),abundances.begin(),[total](double x) {return x/total;});
            Set_Cooling_Constants();
        }

        void Solve_Profile(){
            Determine_Shock_Position();
            Build_Column_Profile();
        }

        double Get_Mass(){return mass;}
        double Get_B_Field(){return b_field;}
        double Get_inv_Mag_Radius(){return inverse_mag_radius;}
        double Get_Corotation_Ratio(){return corotation_ratio;}
        double Get_Distance(){return distance;}
        double Get_Mdot(){return accretion_rate;}
        double Get_Area(){return accretion_area;}
        double Get_Abund(){return metalicity;}
        double Get_Shock_Ratio(){return pressure_ratio;}
        double Get_Inclination_Angle(){return incl_angle;}
        double Get_Shock_Height(){return shock_height;}
        double Get_mBar(){return avg_ion_mass;}
        double Get_ZBar(){return avg_atomic_charge;}
        double Get_Dens_Const(){return density_const;}
        double Get_Exch_Const(){return exchange_const;}
        double Get_Bremss_Const(){return bremss_const;}
        double Get_Cycl_Const(){return cyclotron_const;}
        double Get_Length_Conv(){return length_conv;}
        double Get_Mass_Conv(){return mass_conv;}
        double Get_Time_Conv(){return time_conv;}
        double Get_Vel_Conv(){return vel_conv;}
        double Get_Vol_Conv(){return volume_conv;}
        double Get_Energy_Conv(){return energy_conv;}
        double Get_Density_Conv(){return density_conv;}
        vector<double>& Get_Altitude(){return altitude;}
        vector<double>& Get_Volume(){return volume;}
        vector<double>& Get_Velocity(){return velocity;}
        vector<double>& Get_Density(){return density;}
        vector<double>& Get_Pressure(){return total_pressure;}
        vector<double>& Get_Electron_Pressure(){return electron_pressure;}
        vector<double>& Get_Electron_Density(){return electron_density;}
        vector<double>& Get_Electron_Temperature(){return electron_temperature;}
        vector<double>& Get_Ion_Temperature(){return ion_temperature;}
        vector<double>& Get_Abundance(){return abundances;}

};

PYBIND11_MODULE(_pymcvspec, module) {
    module.def("_mass_to_radius", &Cataclysmic_Variable::Get_Radius, "Returns the radius (cm) for a corresponding WD mass (g)");
    module.def("_luminosity_to_mdot", &Cataclysmic_Variable::Get_Accretion_Rate, "Returns the radius (cm) for a corresponding WD mass (g)");
    py::class_<Py_Cataclysmic_Variable>(module, "_cataclysmic_variable", py::module_local())
        .def(py::init<double,double,double,double,double,double,double,double,double,double,double,double,int>(),
            py::arg("mass") = 0.7*m_sol, py::arg("radius") = 0.01*r_sol, py::arg("b_field") = 1e7,
            py::arg("mdot") = 1e15, py::arg("area") = 1e15, py::arg("inv_r_m") = 0., py::arg("r_m_ratio") = 1.,
            py::arg("metalicity") = 1., py::arg("cos_incl_angle") = 0.5, py::arg("shock_ratio") = 0.75,
            py::arg("shock_coord") = 1e-8, py::arg("src_distance") = 200*pc_to_cm, py::arg("refl_on") = 1)
        .def_property_readonly("mass", &Py_Cataclysmic_Variable::Get_Mass)
        .def_property_readonly("b_field", &Py_Cataclysmic_Variable::Get_B_Field)
        .def_property_readonly("inv_r_m", &Py_Cataclysmic_Variable::Get_inv_Mag_Radius)
        .def_property_readonly("corotation_ratio", &Py_Cataclysmic_Variable::Get_Corotation_Ratio)
        .def_property_readonly("distance", &Py_Cataclysmic_Variable::Get_Distance)
        .def_property_readonly("accretion_rate", &Py_Cataclysmic_Variable::Get_Mdot)
        .def_property_readonly("accretion_area", &Py_Cataclysmic_Variable::Get_Area)
        .def_property_readonly("metalicity", &Py_Cataclysmic_Variable::Get_Abund)
        .def_property_readonly("shock_ratio", &Py_Cataclysmic_Variable::Get_Shock_Ratio)
        .def_property_readonly("inclination_angle", &Py_Cataclysmic_Variable::Get_Inclination_Angle)
        .def_property_readonly("shock_height", &Py_Cataclysmic_Variable::Get_Shock_Height)
        .def_property_readonly("average_ion_mass", &Py_Cataclysmic_Variable::Get_mBar)
        .def_property_readonly("average_ion_charge", &Py_Cataclysmic_Variable::Get_ZBar)
        .def_property_readonly("density_const", &Py_Cataclysmic_Variable::Get_Dens_Const)
        .def_property_readonly("exchange_const", &Py_Cataclysmic_Variable::Get_Exch_Const)
        .def_property_readonly("bremss_const", &Py_Cataclysmic_Variable::Get_Bremss_Const)
        .def_property_readonly("cyclotron_const", &Py_Cataclysmic_Variable::Get_Cycl_Const)
        .def_property_readonly("length_converter", &Py_Cataclysmic_Variable::Get_Length_Conv)
        .def_property_readonly("mass_converter", &Py_Cataclysmic_Variable::Get_Mass_Conv)
        .def_property_readonly("time_converter", &Py_Cataclysmic_Variable::Get_Time_Conv)
        .def_property_readonly("velocity_converter", &Py_Cataclysmic_Variable::Get_Vel_Conv)
        .def_property_readonly("volume_converter", &Py_Cataclysmic_Variable::Get_Vol_Conv)
        .def_property_readonly("energy_converter", &Py_Cataclysmic_Variable::Get_Energy_Conv)
        .def_property_readonly("density_converter", &Py_Cataclysmic_Variable::Get_Density_Conv)
        .def_property_readonly("altitude", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Altitude());})
        .def_property_readonly("volume", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Volume());})
        .def_property_readonly("velocity", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Velocity());})
        .def_property_readonly("density", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Density());})
        .def_property_readonly("total_pressure", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Pressure());})
        .def_property_readonly("electron_pressure", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Electron_Pressure());})
        .def_property_readonly("electron_density", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Electron_Density());})
        .def_property_readonly("electron_temperature", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Electron_Temperature());})
        .def_property_readonly("ion_temperature", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Ion_Temperature());})
        .def_property_readonly("abundance", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Abundance());})
        .def("solve", &Py_Cataclysmic_Variable::Solve_Profile)
        .def("print", &Py_Cataclysmic_Variable::Print_Properties)
        .def("flow_equation", [](Py_Cataclysmic_Variable& self, double s, py::array_t<double, py::array::c_style | py::array::forcecast> state_py){
            vector<double> state = Numpy_to_Vector(state_py);
            vector<double> deriv(state.size());
            self.Flow_Equation(s, state, deriv);
            return Vector_to_Numpy(deriv);
        });

    py::class_<Dipole>(module, "_dipole", py::module_local())
        .def(py::init<double>(),py::arg("u")=0.03)
        .def_readonly("u", &Dipole::u)
        .def_readonly("w0", &Dipole::w_0)
        .def_readonly("a0", &Dipole::a_0)
        .def("update_coordinates", [](Dipole& self, double w){
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
