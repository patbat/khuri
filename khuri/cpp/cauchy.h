#ifndef CAUCHY_HEADER_H
#define CAUCHY_HEADER_H

#include "facilities.h"
#include "gsl_interface.h"

#include <algorithm>
#include <complex>
#include <functional>
#include <initializer_list>
#include <tuple>
#include <vector>

/// Facilities for dealing with complex valued functions.
namespace cauchy {
using Complex = std::complex<double>;
using Curve = std::function<Complex(double)>;
using Complex_function = std::function<Complex(const Complex&)>;
using gsl::Interval;

// -- Basic facilities --------------------------------------------------------

std::vector<double> real(const std::vector<Complex>& vec);
    ///< Return the real parts of the elements of `vec`.

std::vector<double> imag(const std::vector<Complex>& vec);
    ///< Return the imaginary parts of the elements of `vec`.

// -- Piecewise defined functions ---------------------------------------------

/// Piecewise defined function.

/// @tparam Argument needs to provide `operator<=`.
template<class Return, class Argument>
class Piecewise_function {
public:
    using Function = std::function<Return(Argument)>;
        // A `Piecewise_function` is a `Function` itself.

    Piecewise_function(std::initializer_list<Function> pieces,
            std::initializer_list<Argument> boundaries);
        ///< Invariant: `boundaries` needs to be sorted in strictly ascending
        ///< order. `pieces` needs to contain at least one element and
        ///< `boundaries` needs to contain one more element than `pieces`.
        ///<
        ///< The domain of definition of `pieces[i]` is given as
        ///< (`boundaries[i]`,`boundaries[i+1]`].

    Piecewise_function(const Function& f, const Argument& left,
            const Argument& right);
        ///< The domain of definition of `f` is given as [`left`,`right`].

    Return operator()(const Argument& x) const;
        ///< Evaluate curve at `x`.

    template<class F>
    F for_each_piece(F f);
        ///< Call `f` for each element of `pieces`.

        ///< @tparam F needs to meet the requirements of `FunctionObject`.
        ///< Its signature needs to equal `f(Function)`, `f(const Function&)`
        ///< or `f(Function&)`.
private:
    std::vector<Function> pieces;
    std::vector<Argument> boundaries_;
    typename std::vector<Argument>::const_iterator second;
        // points to the second element of `boundaries_`
};

using Piecewise_curve = Piecewise_function<Complex,double>;

template<class Return, class Argument>
Piecewise_function<Return,Argument>::Piecewise_function(
        std::initializer_list<Function> pieces,
        std::initializer_list<Argument> boundaries)
    : pieces{pieces}, boundaries_{boundaries}
{
    if (!pieces.size())
        throw std::invalid_argument{"Piecewise_function needs to contain at \
least one curve"};
    if (boundaries_.size() != pieces.size()+1)
        throw std::invalid_argument{"Piecewise_function needs to contain one \
curve less than boundaries."};
    const bool sorted{std::is_sorted(boundaries_.begin(),boundaries_.end(),
            std::less_equal<Argument>{})};
    if (!sorted)
        throw std::invalid_argument{"Piecewise_function's boundaries need to \
be sorted in strictly ascending order"};

    second = std::next(boundaries_.begin());
}

template<class Return, class Argument>
Piecewise_function<Return,Argument>::Piecewise_function(const Function& f,
        const Argument& left,
        const Argument& right)
    : pieces{f}, boundaries_{left,right}
{
    if (left>=right)
        throw std::invalid_argument{"Piecewise_function's boundaries need to \
be sorted in strictly ascending order"};
    second = std::next(boundaries_.begin());
}

template<class Return, class Argument>
Return Piecewise_function<Return,Argument>::operator()(const Argument& x) const
{
    if (x<boundaries_.front() || boundaries_.back()<x)
        throw std::invalid_argument{"Piecewise_function cannot be evaluated \
outside domain of definition"};

    const auto it{std::find_if(second,boundaries_.end(),
            [x](const Argument& y)
            {
                return x<=y;
            })};
    const auto dist{std::distance(second,it)};
    return pieces[dist](x);
}

template<class Return, class Argument>
template<class F>
F Piecewise_function<Return,Argument>::for_each_piece(F f)
{
    return std::for_each(pieces.begin(),pieces.end(),f);
}

// -- Integration -------------------------------------------------------------

std::tuple<Complex,double,double> c_integrate(const Curve& c,
        double lower, double upper, const gsl::Integration& integrate);
    ///< Integrate `c` in the interval [`lower`,`upper`] using `integrate`.
    ///< Return the value of the integral, the error of the real part and the
    ///< error of the imaginary part.
    ///<
    ///< Note: There is no class for this task, since in many situations one
    ///< needs to integrate both real valued and complex valued functions.
    ///< Using `c_integrate`, the same instance of `Integration` can be used for
    ///< both.

std::tuple<Complex,double,double> c_integrate(const Complex_function& f,
        const Curve& c, const Curve& c_derivative, double lower, double upper,
        const gsl::Integration& integrate);
    ///< Integrate `f` along c` in the interval [`lower`,`upper`] using
    ///< `integrate`. Return the value of the integral, the error of the real
    ///< part and the error of the imaginary part.

// -- Interpolation -----------------------------------------------------------

/// @brief Interpolate data provided as pairs \f$(x_i,y_i)\f$, here \f$y_i\f$
/// is a complex number.
class Interpolate {
public:
    Interpolate(const Interval& x, const std::vector<Complex>& y,
            gsl::Interpolation_method m);
        ///< The sizes of `x` and `y` need to be the same.

    Complex operator()(double x) const;
        ///< Return the value of the (interpolated) data at point `x`.

        ///< If evaluated outside the interval (`front()`,`back()`), the
        ///< boundary values are returned.

    double front() const noexcept {return real_part.front();}
    double back() const noexcept {return real_part.back();}
private:
    gsl::Interpolate real_part;
    gsl::Interpolate imaginary_part;
};

Interpolate sample(const Complex_function& f, const Curve& c,
        const Interval& i, gsl::Interpolation_method m);
    ///< Interpolate `f` along `c`.

    ///< Return interpolator for data pairs
    ///< ( i[k], f(c(i[k])) ). It is assumed that `i` is sorted (in ascending
    ///< order) and contains at least 2 elements. It is guaranteed that the
    ///< returned `Interpolate` instance works at boundaries of `i`.

Interpolate sample(const Curve& c, const Interval& i,
        gsl::Interpolation_method m);
    ///< Interpolate `c` along `i`.

    ///< Return interpolator for data pairs
    ///< ( i[k], c(i[k]) ). It is assumed that `i` is sorted (in ascending
    ///< order) and contains at least 2 elements. It is guaranteed that the
    ///< returned `Interpolate` instance works at boundaries of `i`.
} // cauchy

#endif // CAUCHY_HEADER_H
