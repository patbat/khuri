#include "khuri_treiman.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

#include <vector>
#include <type_traits>

namespace py = pybind11;
using khuri_treiman::Complex;
using khuri_treiman::Curve;
using khuri_treiman::Grid;
using khuri_treiman::Piecewise;

template<typename T>
void create_grid_binding(py::module& m, const std::string& type_name)
{
    using V = std::add_lvalue_reference_t<std::add_const_t<T>>;
    using G = Grid<T>;
    const std::string name = "Grid" + type_name;
    const std::string init_docstring =
        "Parameters\n"
        "----------\n"
        "t: The continuous curve in the x-plane.\n"
        "x_sizes: The number of knots along the (different segements of the)"
        " curve in the x-plane.\n"
        "z_size: The number of knots along the line in the z-plane.\n";
    py::class_<G, T>(m, name.c_str())
        .def(py::init<V, std::vector<std::size_t>, std::size_t>(),
             init_docstring.c_str())
        .def("__call__", py::vectorize(&G::operator()),
             py::arg("x_index"),
             py::arg("z_index"))
        .def("x_parameter_values", &G::x_parameter_values)
        .def("x", &G::x,
             py::arg("x_index"))
        .def("derivative", &G::derivative,
             py::arg("x_index"))
        .def("z", &G::z,
             py::arg("z_index"))
        .def("x_size", &G::x_size)
        .def("z_size", &G::z_size)
        .def("x_parameter_lower", &G::x_parameter_lower)
        .def("x_parameter_upper", &G::x_parameter_upper);
}

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

    py::class_<Piecewise, Curve> piecewise(m, "Piecewise",
                                           "a piecewise linear path in the"
                                           " complex plane");
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

    py::class_<khuri_treiman::Real, Piecewise>(m, "Real",
                                            "linear curve along the real axis")
        .def(py::init<double, double>());

    py::class_<khuri_treiman::Vector_decay, Piecewise>(m, "Vector_decay",
                                            "curve for the decay as described"
                                            " in the paper by Gasser and"
                                            " Rusetsky")
        .def(py::init<double, double, double>());

    py::class_<khuri_treiman::Adaptive, Piecewise>(m, "Adaptive",
                                            "curve for arbitrary virtualities"
                                            " above the three-pion threshold"
                                            " and arbitrary pion masses")
        .def(py::init<double, double, double>());

    create_grid_binding<khuri_treiman::Real>(m, "Real");
    create_grid_binding<khuri_treiman::Vector_decay>(m, "Vector_decay");
    create_grid_binding<khuri_treiman::Adaptive>(m, "Adaptive");
}
