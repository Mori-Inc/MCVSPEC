#include <pybind11/attr.h>
#include <pybind11/buffer_info.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "dipole.hh"
#include <vector>

namespace py = pybind11;

using std::vector;

template <typename T>
static py::array_t<T> Vector_to_Numpy(const std::vector<T>& cpp_vec) {
  py::array_t<T> np_array(cpp_vec.size());
  std::memcpy(np_array.mutable_data(), cpp_vec.data(), cpp_vec.size()*sizeof(T));
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
    std::memcpy(cpp_vec.data(), np_array.data(), cpp_vec.size()*sizeof(double));
    return cpp_vec;
}

template <size_t n_dim>
inline void Numpy_to_State(const double* np_ptr, State<n_dim>& state) {
    std::memcpy(state.data(), np_ptr, n_dim*sizeof(double));
}

template <size_t n_dim>
inline void State_to_Numpy(const State<n_dim>& state, double* np_ptr) {
    std::memcpy(np_ptr, state.data(), n_dim*sizeof(double));
}

class Py_Cataclysmic_Variable : public Cataclysmic_Variable {
    public:
        static constexpr size_t n_dim = Cataclysmic_Variable::n_dim;
        Py_Cataclysmic_Variable(White_Dwarf wd, Accretion_Column col, Tolerance tol):
            Cataclysmic_Variable(wd, col, tol)
        {
            Set_Abundances();
        }

        void Set_Abundances() override{
            abundances.resize(n_elements);
            abundances = {1.00e+00, 9.77e-02, 3.63e-04, 1.12e-04, 8.51e-04, 1.23e-04,
                          3.80e-05, 2.95e-06, 3.55e-05, 1.62e-05, 3.63e-06, 2.29e-06,
                          4.68e-05, 1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
            double total= abundances[0]+abundances[1];
            for(size_t i=2; i<abundances.size(); i++){
                abundances[i] *= accretion_column.metallicity;
                total += abundances[i];
            }
            std::transform(abundances.begin(),abundances.end(),abundances.begin(),[total](double x) {return x/total;});
            Set_Cooling_Constants();
        }

        int Solve_Profile(){
            Find_Shock_Position();
            Build_Column_Profile();
            if(!valid_solution){
                return -1;
            }
            return 0;
        }

        const double Get_Mass() const {return white_dwarf.mass;}
        const double Get_B_Field() const {return white_dwarf.b_field;}
        const double Get_inv_Mag_Radius() const {return white_dwarf.inverse_mag_radius;}
        const double Get_Corotation_Radius() const {return white_dwarf.corotation_radius;}
        const double Get_Distance() const {return white_dwarf.distance;}
        const double Get_Mdot() const {return accretion_column.accretion_rate;}
        const double Get_Area() const {return accretion_column.accretion_area;}
        const double Get_Abund() const {return accretion_column.metallicity;}
        const double Get_Shock_Ratio() const {return accretion_column.shock_pressure_ratio;}
        const double Get_Cos_Inclination_Angle() const {return white_dwarf.cos_inclination;}
        const double Get_mBar() const {return avg_ion_mass;}
        const double Get_ZBar() const {return avg_atomic_charge;}
        const double Get_Dens_to_ne() const {return mass_to_number_density;}
        const double Get_Exch_Const() const {return exchange_const;}
        const double Get_Bremss_Const() const {return bremss_const;}
        const double Get_Cycl_Const() const {return cyclotron_const;}
        const double Get_Length_Conv() const {return length_conv;}
        const double Get_Mass_Conv() const {return mass_conv;}
        const double Get_Time_Conv() const {return time_conv;}
        const double Get_Vel_Conv() const {return vel_conv;}
        const double Get_Vol_Conv() const {return volume_conv;}
        const double Get_Energy_Conv() const {return energy_conv;}
        const double Get_Density_Conv() const {return density_conv;}
        const double Get_Column_Coord() const {return geometry.u;}
        const vector<double>& Get_Position() const {return position;}
        const vector<double>& Get_Altitude() const {return altitude;}
        const vector<double>& Get_Volume_Element() const {return volume_element;}
        const vector<double>& Get_Velocity() const {return velocity;}
        const vector<double>& Get_Density() const {return density;}
        const vector<double>& Get_Pressure() const {return total_pressure;}
        const vector<double>& Get_Electron_Pressure() const {return electron_pressure;}
        const vector<double>& Get_Electron_Density() const {return electron_density;}
        const vector<double>& Get_Electron_Temperature() const {return electron_temperature;}
        const vector<double>& Get_Ion_Temperature() const {return ion_temperature;}
        const vector<double>& Get_Abundance() const {return abundances;}
        State<n_dim> state{}, deriv{};
};

PYBIND11_MODULE(_pymcvspec, module) {
    module.attr("_atomic_charges") = py::cast(atomic_charge);
    module.attr("_atomic_masses") = py::cast(atomic_charge);
    module.def("_mass_to_radius", Mass_to_Radius, "Returns the radius (cm) for a corresponding WD mass (g)");
    module.def("_luminosity_to_mdot", Luminosity_to_Accretion_Rate, "Returns the accretion rate (g/s) for a corresponding luminosity (erg/s), mass (g), and radius (cm)");
    py::class_<White_Dwarf>(module, "_white_dwarf", py::module_local())
        .def(py::init<double, double, double, double, double, double, double>(),
            py::arg("mass") = 0.7*m_sol, py::arg("radius") = 0.01*r_sol, py::arg("b_field") = 1e7,
            py::arg("cos_incl") = 0.5, py::arg("inv_mag_rad") = 0, py::arg("corot_rad") = 1,
            py::arg("distance") = 200*pc_to_cm)
        .def_readwrite("mass", &White_Dwarf::mass)
        .def_readwrite("radius", &White_Dwarf::radius)
        .def_readwrite("b_field", &White_Dwarf::b_field)
        .def_readwrite("cos_inclination", &White_Dwarf::cos_inclination)
        .def_readwrite("inverse_mag_radius", &White_Dwarf::inverse_mag_radius)
        .def_readwrite("corotation_radius", &White_Dwarf::corotation_radius)
        .def_readwrite("distance", &White_Dwarf::distance);
    py::class_<Accretion_Column>(module, "_accretion_column", py::module_local())
        .def(py::init<double, double, double, double, double>(),
            py::arg("mdot")=1, py::arg("area")=1e15, py::arg("metallicity")=1,
            py::arg("shock_ratio")=0.75, py::arg("sin_colat")=0)
        .def_readwrite("mdot", &Accretion_Column::accretion_rate)
        .def_readwrite("area", &Accretion_Column::accretion_area)
        .def_readwrite("metallicity", &Accretion_Column::metallicity)
        .def_readwrite("shock_ratio", &Accretion_Column::shock_pressure_ratio)
        .def_readwrite("sin_colatitude", &Accretion_Column::sin_mag_colat);
    py::class_<Tolerance>(module, "_tolerance", py::module_local())
        .def(py::init<double, double, double, double>(),
            py::arg("abserr")=1e-8, py::arg("relerr")=1e-6, py::arg("dkT")=0.5, py::arg("dz")=0.1)
        .def_readwrite("abs_err", &Tolerance::absolute_error)
        .def_readwrite("rel_err", &Tolerance::relative_error)
        .def_readwrite("dkT", &Tolerance::kT_grid_spacing)
        .def_readwrite("dz", &Tolerance::altitude_grid_spacing);
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
            for(size_t i = 0; i < 3; i++){
                np_array(i) = metric[i];
            }
            return py::make_tuple(r, proj_r_w, convergance, array);
        });
    py::class_<Py_Cataclysmic_Variable>(module, "_cataclysmic_variable", py::module_local())
        .def(py::init<White_Dwarf, Accretion_Column, Tolerance>(),
                py::arg("white_dwarf"), py::arg("accretion_column"), py::arg("tol"))
        .def_property_readonly("mass", &Py_Cataclysmic_Variable::Get_Mass)
        .def_property_readonly("b_field", &Py_Cataclysmic_Variable::Get_B_Field)
        .def_property_readonly("inv_r_m", &Py_Cataclysmic_Variable::Get_inv_Mag_Radius)
        .def_property_readonly("corotation_radius", &Py_Cataclysmic_Variable::Get_Corotation_Radius)
        .def_property_readonly("distance", &Py_Cataclysmic_Variable::Get_Distance)
        .def_property_readonly("accretion_rate", &Py_Cataclysmic_Variable::Get_Mdot)
        .def_property_readonly("accretion_area", &Py_Cataclysmic_Variable::Get_Area)
        .def_property_readonly("metallicity", &Py_Cataclysmic_Variable::Get_Abund)
        .def_property_readonly("shock_ratio", &Py_Cataclysmic_Variable::Get_Shock_Ratio)
        .def_property_readonly("cos_inclination_angle", &Py_Cataclysmic_Variable::Get_Cos_Inclination_Angle)
        .def_property_readonly("average_ion_mass", &Py_Cataclysmic_Variable::Get_mBar)
        .def_property_readonly("average_ion_charge", &Py_Cataclysmic_Variable::Get_ZBar)
        .def_property_readonly("density_to_ne", &Py_Cataclysmic_Variable::Get_Dens_to_ne)
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
        .def_property_readonly("column_coord", &Py_Cataclysmic_Variable::Get_Column_Coord)
        .def_property_readonly("position", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Position());})
        .def_property_readonly("altitude", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Altitude());})
        .def_property_readonly("volume_element", [](Py_Cataclysmic_Variable& self) { return Vector_to_Numpy(self.Get_Volume_Element());})
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
            Numpy_to_State(state_py.data(), self.state);
            self.Flow_Equation(s, self.state, self.deriv);
            py::array_t<double> np_array(self.n_dim);
            State_to_Numpy(self.deriv, np_array.mutable_data());
            return np_array;
        });
}
