#include "phase_space.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_MODULE(_khuri_phase_space, m) {
    m.doc() = "Two-particle phase spaces with different analytic structure";

    m.def("rho", py::vectorize(phase_space::rho),
        py::arg("mass"),
        py::arg("s"));

    m.def("sigma", py::vectorize(phase_space::sigma),
        py::arg("mass"),
        py::arg("s"));
}
