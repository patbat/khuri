#include "curved_omnes.h"
#include "piecewise.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"

namespace py = pybind11;
using omnes::OmnesF;
using curved_omnes::CFunction;
using curved_omnes::CurvedOmnes;

template<typename Curve>
void create_binding(py::module&m, const std::string& name)
{
    py::class_<CurvedOmnes>(m, ("CurvedOmnes" + name).c_str())
        .def(py::init<OmnesF, CFunction, Curve>(),
             py::arg("omnes"),
             py::arg("amplitude"),
             py::arg("curve"))
        .def("__call__", py::vectorize(&CurvedOmnes::operator()),
             py::arg("mandelstam_s"))
        .def("original", &CurvedOmnes::original);
}

PYBIND11_MODULE(_khuri_curved_omnes, m) {
    m.doc() = "The Omnes function with cut along a somewhat general curve.";

    create_binding<piecewise::Real>(m, "Real");
    create_binding<piecewise::Vector_decay>(m, "VectorDecay");
    create_binding<piecewise::Adaptive>(m, "Adaptive");
}
