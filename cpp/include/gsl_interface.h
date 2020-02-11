#ifndef GSL_INTERFACE_H
#define GSL_INTERFACE_H

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <utility>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

/// @brief Interface to some of the gsl facilities (including more modern error
/// handling).
namespace gsl {
using Value = std::pair<double,double>; // a value and its error
using Function = std::function<double(double)>;

/// By representing intervals used e.g. for interpolation in this way, the user
/// can choose regions with dense sampling and ones with coarse sampling.
/// An `Interval` needs to be sorted in ascending order.
using Interval = std::vector<double>;

// -- Error handling ----------------------------------------------------------

struct Turn_off_gsl_errors_helper {
    Turn_off_gsl_errors_helper() {gsl_set_error_handler_off();}
        // Turn off default GSL error handling, sucht that error handling with
        // exceptions can be invoked.
};

class Turn_off_gsl_errors {
    // Turn off the default GSL error handling as soon as this header file is
    // included.
    const static Turn_off_gsl_errors_helper t;
};

void check(int status);
    ///< Check whether a `status` returned by a GSL routine indicates failure.
    ///< If so, throw an appropriate exception.

template<class Function, class ...Arguments>
void call(Function f, Arguments... a)
    /// @brief Call GSL routine `f` with arguments `a` and throw an exception
    /// if an error occurs.
    ///
    /// Note: All GSL routines that return an int to indicate the status
    /// (e.g. success, errors) should be invoked only by using `call`.
{
    int status{f(a...)};
    check(status);
}

class Error : public std::exception {
public:
    Error(const std::string& message) : message{message} {}
    const char* what() const noexcept {return message.data();}
private:
    std::string message;
};

struct Allocation_error : Error {
    using Error::Error;
};

struct Divergence_error : Error {
    using Error::Error;
};

struct Subdivision_error : Error {
    using Error::Error;
};

struct Roundoff_error : Error {
    using Error::Error;
};

struct Bad_integrand_error : Error {
    using Error::Error;
};

struct Domain_error : Error {
    using Error::Error;
};


// -- Integration: Gauss-Legendre  --------------------------------------------

/// Deleter needed for Gauss-Legendre integration
struct Glfixed_deleter {
    void operator()(gsl_integration_glfixed_table* p)
    {
        gsl_integration_glfixed_table_free(p);
    }
};

/// Integration via Gauss-Legendre quadrature.
class Gauss_Legendre {
public:
    explicit Gauss_Legendre(std::size_t n);
        ///< Allocate resources for a `n`-point integration scheme.
    Gauss_Legendre(const Gauss_Legendre& other);
    Gauss_Legendre(Gauss_Legendre&& other) noexcept;
    Gauss_Legendre& operator=(const Gauss_Legendre& other);
    Gauss_Legendre& operator=(Gauss_Legendre&& other) noexcept;

    std::pair<double,double> point(double lower, double upper, std::size_t i)
        const;
        ///< @brief Return the `i`th (point,weight) pair for an integration in
        ///< the interval [`lower`,`upper`].
    double operator()(Function f, double lower, double upper) const;
        ///< Integrate `f` from `lower` to `upper`.
    void resize(std::size_t s);
        ///< Adjust the number of points of the integration scheme.
    std::size_t size() const noexcept;
        ///< Return the number of points of the integration scheme.
    ~Gauss_Legendre() noexcept;
private:
    std::size_t _size;
    std::unique_ptr<gsl_integration_glfixed_table,Glfixed_deleter> table;
};


// -- Integration: adaptive routines ------------------------------------------

struct Integration {
    virtual Value operator()(Function f, double lower, double upper) const=0;
        ///< Integrate the function `f` in the interval [`lower`,`upper`].

        ///< Both `lower` and `upper` are allowed to be infinity
        ///< (use e.g. `std::numeric_limits<double>::infinity()`).

    virtual ~Integration() {}
};

template<typename Space, std::enable_if_t<std::is_same<gsl_integration_workspace, Space>::value, int> = 0>
Space* allocate(std::size_t size)
    ///< Allocate the resources needed by adaptive integration routines.

    ///< @tparam Space needs to equal the type of a GSL workspace.
{
    return gsl_integration_workspace_alloc(size);
}

template<typename Space, std::enable_if_t<std::is_same<gsl_integration_cquad_workspace, Space>::value, int> = 0>
Space* allocate(std::size_t size)
    ///< Allocate the resources needed by adaptive integration routines.

    ///< @tparam Space needs to equal the type of a GSL workspace.
{
    return gsl_integration_cquad_workspace_alloc(size);
}

inline void deallocate(gsl_integration_workspace* workspace)
{
    gsl_integration_workspace_free(workspace);
}

inline void deallocate(gsl_integration_cquad_workspace* workspace)
{
        gsl_integration_cquad_workspace_free(workspace);
}

/// Manage the resources needed by adaptive GSL integration routines.
template<class Space>
class Workspace {
public:
    explicit Workspace(std::size_t space=0);
        ///< Create a workspace of size `space`.
    Workspace(const Workspace& other);
    Workspace(Workspace&& other);

    Workspace& operator=(const Workspace& other);
    Workspace& operator=(Workspace&& other);

    void resize(std::size_t space);
        ///< Adjust the size of the workspace to equal space.
    void reserve(std::size_t space);
        ///< Change the size of the workspace (never decrease size).
    std::size_t size() const noexcept {return available_space;}

    Space* data() const noexcept {return workspace;}

    ~Workspace() noexcept;
private:
    std::size_t available_space;
    Space* workspace;
};

template<class Space>
Workspace<Space>::Workspace(std::size_t space)
    : available_space{space}, workspace{allocate<Space>(space)}
{
}

template<class Space>
Workspace<Space>::Workspace(const Workspace<Space>& other)
    : available_space{other.size()},
    workspace{allocate<Space>(available_space)}
{
}

template<class Space>
Workspace<Space>::Workspace(Workspace<Space>&& other)
    : available_space{other.size()},
    workspace{other.workspace}
{
    other.available_space = 0;
    other.workspace = nullptr;
}

template<class Space>
Workspace<Space>& Workspace<Space>::operator=(const Workspace<Space>& other)
{
    resize(other.size());
    return *this;
}

template<class Space>
Workspace<Space>& Workspace<Space>::operator=(Workspace<Space>&& other)
{
    // update workspace
    resize(0);
    available_space = other.size();
    workspace = other.workspace;

    // empty old workspace
    other.available_space = 0;
    other.workspace = nullptr;

    return *this;
}

template<class Space>
void Workspace<Space>::resize(std::size_t space)
{
    if (space!=size()) {
        deallocate(workspace);
        available_space = space;
        workspace = allocate<Space>(size());
    }
}

template<class Space>
void Workspace<Space>::reserve(std::size_t space)
{
    if (space>size())
        resize(space);
}

template<class Space>
Workspace<Space>::~Workspace() noexcept
{
    deallocate(workspace);
}

using Qag_workspace = Workspace<gsl_integration_workspace>;
using Cquad_workspace = Workspace<gsl_integration_cquad_workspace>;

// Note:
// Up to now, there are two classes derived from `Integration`, namely `Qag`
// and `Cquad`. Both have a small overlap, i.e. the trivial functions
// (`set_absolute`, `size` etc.) are the same. However, there is no class in
// between `Integration` and `Qag`/`Cquad`. This is due to the fact that the
// other GSL integration routines do not make use of this overlap (have a
// different overlap). Hence, an additional class would lead to a more
// complicated class hierarchy without beeing of use in future extentions of
// this interface.

struct Settings {
    double absolute_precision{0.0};
    double relative_precision{1e-7};
    std::size_t space{1000};
};

/// @brief Integration of one function or multiple functions using GSL CQUAD
/// adaptive routines. This routine is able to handle more difficult
/// integrands compared to `Qag`.
class Cquad : public Integration {
public:
    Cquad(const Settings& set=Settings{});
        // If `absolute_precision` is set to zero, `relative_precision` is used
        // and vice versa. `space` denotes the size of the workspace used by
        // the gsl integration routine.
    Cquad(const Cquad&)=default;
    Cquad(Cquad&&)=default;
    Cquad& operator=(const Cquad&)=default;
    Cquad& operator=(Cquad&&)=default;
    ~Cquad() noexcept {}

    Value operator()(Function f, double lower, double upper) const override;

    void reserve(std::size_t space);
        ///< Change the size of the workspace used by the gsl integration
        ///< routine.
    void set_absolute(double abs) noexcept {absolute_precision = abs;}
    void set_relative(double rel) noexcept {relative_precision = rel;}

    double absolute() const noexcept {return absolute_precision;}
    double relative() const noexcept {return relative_precision;}
    std::size_t size() const noexcept {return workspace.size();}
private:
    double absolute_precision;
    double relative_precision;
    Cquad_workspace workspace;
};

/// @brief Integration of one function or multiple functions using GSL QAG
/// adaptive routines.
class Qag : public Integration {
public:
    Qag(const Settings& set=Settings{});
        ///< If `absolute_precision` is set to zero, `relative_precision` is
        ///< used and vice versa. `space` denotes the size of the workspace
        ///< used by the gsl integration routine.
    Qag(const Qag&)=default;
    Qag(Qag&&)=default;
    Qag& operator=(const Qag&)=default;
    Qag& operator=(Qag&&)=default;
    ~Qag() noexcept {}

    Value operator()(Function f, double lower, double upper) const override;

    void reserve(std::size_t space);
        // Change the size of the workspace used by the gsl integration
        // routine.
    void set_absolute(double abs) noexcept {absolute_precision = abs;}
    void set_relative(double rel) noexcept {relative_precision = rel;}

    double absolute() const noexcept {return absolute_precision;}
    double relative() const noexcept {return relative_precision;}
    std::size_t size() const noexcept {return workspace.size();}
private:
    double absolute_precision;
    double relative_precision;
    std::size_t limit;
    Qag_workspace workspace;
};

// -- Interpolation -----------------------------------------------------------

/// @brief These methods can be used by the interpolation routine accessed via
/// `Interpolate`.
class Interpolation_method {
public:
    enum Method {
        linear,
        polynomial,
        cubic,
        cubic_periodic,
        akima,
        akima_periodic,
        steffen
    };

    Interpolation_method(Method m);
    Interpolation_method(const gsl_interp_type* p) : method{p} {}

    const gsl_interp_type* get_method() const noexcept {return method;}
    unsigned int min_size() const {return gsl_interp_type_min_size(method);}
        ///< Return the minimal number of required points, such that the method
        ///< can be used for interpolation.
private:
    // The default copy constructor etc. are fine, since this pointer does not
    // handle resources.
    const gsl_interp_type* method;
};

/// Interpolation of 1 dimensional data provided as pairs (x_i,y_i).
class Interpolate {
public:
    Interpolate(const Interval& x, const std::vector<double>& y,
            Interpolation_method m, bool tolerant=true);
        ///< The sizes of `x` and `y` need to be the same. `tolerant` influences
        ///< the behavour of `operator()` at the boundaries and beyond, see
        ///< below.
    Interpolate(const Interpolate& other);
    Interpolate(Interpolate&& other);
    Interpolate& operator=(const Interpolate& other);
    Interpolate& operator=(Interpolate&& other);
    ~Interpolate() noexcept;

    // No default constructor, since an empty data set cannot be interpolated
    // in a meaningful way.

    double operator()(double x) const;
        ///< Return the value of the (interpolated) data at point `x`.
        ///< If `tolerant==true`, the interpolator will return the boundary
        ///< values if `x` is outside the interpolated interval. Otherwise
        ///< the interpolator works only for values in the interval
        ///< [`front()`,`back()`]. In this case, evaluation outside the interval
        ///< will throw.

    double front() const noexcept {return x_data.front();}
    double back() const noexcept {return x_data.back();}

    bool is_tolerant() noexcept {return tolerant;}
    void be_tolerant() noexcept {tolerant = true;}
    void be_strict() noexcept {tolerant = false;}

private:
    // The GSL routine needs a copy of the data. This copy is stored
    // explicitly (instead of implicitly by the GSL 1d higher-level interface)
    // to allow for copy and move constructors.
    std::vector<double> x_data;
    std::vector<double> y_data;

    Interpolation_method method;
    bool tolerant;
    gsl_interp_accel* acc{gsl_interp_accel_alloc()};
    gsl_interp* spline;
};

Interpolate sample(const Function& f, const Interval& i,
        Interpolation_method m, bool tolerant=true);
    ///< Interpolate `f` along `i`.

template<class InIterA, class InIterB>
Interpolate make_interpolate(InIterA first1, InIterA last1, InIterB first2,
        InIterB last2, Interpolation_method m)
    /// Make an interpolator for the data [first1,last1) and [first2,last2).

    /// Both data sets need to be of the same size and contain doubles.
    /// `InIterA` and `InIterB` need to meet the requirements of InputIterator.
    /// `InIterA` and `InIterB` need to point to doubles.
{
    constexpr bool check1{std::is_same<double,
        typename std::iterator_traits<InIterA>::value_type>::value};
    constexpr bool check2{std::is_same<double,
        typename std::iterator_traits<InIterB>::value_type>::value};
    static_assert(check1 && check2, "value_types need to be the same");

    // The data is copied into a vector since the gsl routine works for data
    // stored in a plain array only.
    std::vector<double> x_values{first1,last1};
    std::vector<double> y_values{first2,last2};

    return Interpolate{x_values,y_values,m};
}

template<class ContainerA, class ContainerB>
Interpolate make_interpolate(ContainerA a, ContainerB b,
        Interpolation_method m)
    /// Make an Interpolator for the data stored in `a` and `b`.

    /// `ContainerA` and `ContainerB` need to meet the requirements of
    /// `Container`. Both need to have double as value_type.
{
    return make_interpolate(a.begin(),a.end(),b.begin(),b.end(),m);
}

// -- Helper-functions --------------------------------------------------------

template<class Number>
bool signed_interval(Number& a, Number& b)
    /// @brief Check if [a,b] is a valid integral (i.e. a<=b).
    /// Swap a and b if it is not.

    /// @tparam Number needs to meet the requirements of LessThanComparable.
{
    if (a>b) {
        std::swap(a,b);
        return false;
    }
    return true;
}
} // gsl

#endif // GSL_INTERFACE_H
