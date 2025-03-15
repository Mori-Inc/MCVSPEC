#include <pybind11/pybind11.h>
#include "../include/Cataclysmic_Variable.hh"

namespace py = pybind11;

PYBIND11_PLUGIN(mcvspec) {
    py::module m("mcvspec", "MCVSPEC");
    py::class_<Cataclysmic_Variable>(m, "Polar")
        .def(py::init<double,double,double,double,double,double,double,int>())
        .def("execute", &Cataclysmic_Variable::Shock_Height_Shooting);
}
