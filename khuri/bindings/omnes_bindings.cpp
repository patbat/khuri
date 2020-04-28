#include "omnes.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"
#include "pybind11/functional.h"
#include "gsl_interface.h"

namespace py = pybind11;
using omnes::Omnes;
using gsl::Function;
using gsl::Settings;

template<typename T>
void create_binding(py::module& m, const std::string& name)
{
    py::class_<Omnes<T>>(m, name.c_str())
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
        .def("__call__", py::vectorize(&Omnes<T>::operator()),
                py::arg("s"));
}

template<typename T>
void second_sheet_binding(py::module& m, const std::string& name)
{
    m.def(name.c_str(),
        // py::vectorize(omnes::second_sheet<T>),
        omnes::second_sheet<T>,
        "Evaluate the Omnes function on the second Riemann sheet.",
        py::arg("omnes_function"),
        py::arg("amplitude"),
        py::arg("s"));
}

PYBIND11_MODULE(_khuri_omnes, m) {
    m.doc() = "The Omnes function.";

    create_binding<gsl::Cquad>(m, "OmnesCquad");
    create_binding<gsl::Qag>(m, "OmnesQag");

    second_sheet_binding<gsl::Cquad>(m, "second_sheet_cquad");
    second_sheet_binding<gsl::Qag>(m, "second_sheet_quag");
}
