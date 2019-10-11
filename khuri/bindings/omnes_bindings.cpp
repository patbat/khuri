#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"
#include "pybind11/functional.h"
#include "gsl_interface.h"
#include "omnes.h"

namespace py = pybind11;
using omnes::Omnes;
using gsl::Function;
using gsl::Settings;

PYBIND11_MODULE(_khuri_omnes, m) {
    m.doc() = "The Omnes function.";

    py::class_<Omnes<gsl::Cquad>>(m, "OmnesCquad")
        .def(py::init<const Function&, double, double, Settings>(),
             py::arg("phase"),
             py::arg("threshold"),
             py::arg("minimal_distance") = 1e-10,
             py::arg("config") = Settings{})
        .def(py::init<const Function&, double, double, double, double,
                      Settings>(),
             py::arg("phase"),
             py::arg("threshold"),
             py::arg("constant"),
             py::arg("cut"),
             py::arg("minimal_distance") = 1e-10,
             py::arg("config") = Settings{})
        .def("__call__", py::vectorize(&Omnes<gsl::Cquad>::operator()),
                py::arg("s"));

    py::class_<Omnes<gsl::Qag>>(m, "OmnesQag")
        .def(py::init<const Function&, double, double, Settings>(),
             py::arg("phase"),
             py::arg("threshold"),
             py::arg("minimal_distance") = 1e-10,
             py::arg("config") = Settings{})
        .def(py::init<const Function&, double, double, double, double,
                      Settings>(),
             py::arg("phase"),
             py::arg("threshold"),
             py::arg("constant"),
             py::arg("cut"),
             py::arg("minimal_distance") = 1e-10,
             py::arg("config") = Settings{})
        .def("__call__", py::vectorize(&Omnes<gsl::Qag>::operator()),
                py::arg("s"));
}
