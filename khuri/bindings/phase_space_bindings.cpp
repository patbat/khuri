#include "phase_space.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"

namespace py = pybind11;
using phase_space::Complex;

PYBIND11_MODULE(_khuri_phase_space, m) {
    m.doc() = "Two-particle phase spaces with different analytic structure";

    m.def("rho",
        py::vectorize(py::overload_cast<double, double>(phase_space::rho)),
        py::arg("mass"),
        py::arg("s"));

    m.def("rho",
        py::vectorize(py::overload_cast<double, Complex>(phase_space::rho)),
        py::arg("mass"),
        py::arg("s"));

    m.def("sigma",
        py::vectorize(py::overload_cast<double, double>(phase_space::sigma)),
        py::arg("mass"),
        py::arg("s"));

    m.def("sigma",
        py::vectorize(py::overload_cast<double, Complex>(phase_space::sigma)),
        py::arg("mass"),
        py::arg("s"));
}
