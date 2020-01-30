#include "khuri_treiman.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

namespace py = pybind11;
using khuri_treiman::Complex;
using khuri_treiman::Curve;
using khuri_treiman::Piecewise;

PYBIND11_MODULE(_khuri_khuri_treiman, m) {
    m.doc() = "Khuri Treiman equations.";

    py::class_<Curve>(m, "Curve", "A curve in the complex plane.")
        .def("curve_func", py::vectorize(&Curve::curve_func),
             "Evaluate the curve at `x`.",
             py::arg("x"))
        .def("derivative_func", py::vectorize(&Curve::derivative_func),
             "Evaluate the derivative of the curve at `x`.",
             py::arg("x"))
        .def("hits", &Curve::hits,
             "Determine whether `s` hits the curve.\n\n"
             "If `s` lies on the curve, return the parameter values marking"
             " the beginning and the end of the segment that is hit.",
             py::arg("x"))
        .def("boundaries", &Curve::boundaries);

    py::class_<Piecewise, Curve> piecewise(m, "Piecewise");
    piecewise
        .def(py::init<const std::vector<Complex>&,
                      const std::vector<Piecewise::Para>&>())
        .def("lower", &Piecewise::lower,
             "Return the parameter value corresponding to the start of the"
             " curve.")
        .def("upper", &Piecewise::upper,
             "Return the parameter value corresponding to the end of the"
             " curve.")
        .def("piece_index", py::vectorize(&Piecewise::piece_index),
                "Return the number of the segment corresponding to the"
                " parameter value `x`.",
             py::arg("x"));

    py::enum_<Piecewise::Para>(piecewise, "Para")
        .value("linear", Piecewise::Para::linear)
        .value("quadratic", Piecewise::Para::quadratic)
        .export_values();
}
