#include "khuri_treiman.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;
using khuri_treiman::Complex;
using khuri_treiman::Piecewise;

PYBIND11_MODULE(_khuri_khuri_treiman, m) {
    m.doc() = "Khuri Treiman equations.";

    py::class_<Piecewise> piecewise(m, "Piecewise");
    piecewise.def(py::init<const std::vector<Complex>&,
                           const std::vector<Piecewise::Para>&>());

    py::enum_<Piecewise::Para>(piecewise, "Para")
        .value("linear", Piecewise::Para::linear)
        .value("quadratic", Piecewise::Para::quadratic)
        .export_values();
}
