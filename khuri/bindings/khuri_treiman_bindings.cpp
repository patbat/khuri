#include "khuri_treiman.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

namespace py = pybind11;
using khuri_treiman::Complex;
using khuri_treiman::Piecewise;

PYBIND11_MODULE(_khuri_khuri_treiman, m) {
    m.doc() = "Khuri Treiman equations.";

    py::class_<Piecewise> piecewise(m, "Piecewise");
    piecewise
        .def(py::init<const std::vector<Complex>&,
                      const std::vector<Piecewise::Para>&>())
        .def("curve_func", py::vectorize(&Piecewise::curve_func),
             "Evalute the piecewise linear curve at `x`.",
             py::arg("x"))
        .def("derivative_func", py::vectorize(&Piecewise::derivative_func),
             "Evalute the derivative of the piecewise linear curve at `x`.",
             py::arg("x"))
        .def("hits", &Piecewise::hits,
             py::arg("x"))
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
