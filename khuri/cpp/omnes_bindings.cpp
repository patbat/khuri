#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"
#include "gsl_interface.h"
#include "omnes.h"

namespace py = pybind11;
using omnes::Omnes;
using gsl::Function;
using gsl::Settings;


PYBIND11_MODULE(omnes, m) {
    m.doc() = "The Omnes function.";

    py::class_<Omnes>(m, "Omnes")
        .def(py::init<const Function&, double, double, Settings>())
        .def(py::init<const Function&, double, double, double, double,
                      Settings>())
        .def("__call__", py::vectorize(&Omnes::operator()),
                py::arg("s"));
}
