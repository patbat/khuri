#include "gsl_interface.h"

namespace gsl {
// -- Error handling ----------------------------------------------------------
const Turn_off_gsl_errors_helper Turn_off_gsl_errors::t{};

void check(int status)
{
    std::string message{gsl_strerror(status)};
    switch (status) {
        case 0: // everything worked fine
            return;
        case GSL_ENOMEM:  // could not allocate enough space
            throw Allocation_error{message};
        case GSL_EDIVERGE: // integral is divergent or too slowly convergent
            throw Divergence_error{message};
        case GSL_EMAXITER: // maximum number of subdivisions in integration 
                           // exceeded
            throw Subdivision_error{message};
        case GSL_EROUND: // roundoff error
            throw Roundoff_error{message};
        case GSL_ESING: // bad integrand behaviour (e.g. non-integrable
                        // singularity)
            throw Bad_integrand_error{message};
        case GSL_EDOM: // domain error
            throw Domain_error{message};
        default:
            throw Error{message};
    }
}


// -- Integration: general ----------------------------------------------------

double unwrap(double x, void* p)
    // Call *p with x, where *p should be a std::function<double(double)>.
    // `unwrap` provides the signature needed by the gsl integration routine.
{
    auto fp = static_cast<std::function<double(double)>*>(p);
    return (*fp)(x);
}


// -- Integration: Gauss-Legendre  --------------------------------------------

Gauss_Legendre::Gauss_Legendre(std::size_t s)
    : _size{s}, table{gsl_integration_glfixed_table_alloc(s)}
{
}

Gauss_Legendre::Gauss_Legendre(const Gauss_Legendre& other)
    : _size{other.size()}, table{gsl_integration_glfixed_table_alloc(_size)}
{
}

Gauss_Legendre::Gauss_Legendre(Gauss_Legendre&& other) noexcept
    : _size{other.size()}, table{std::move(other.table)}
{
}

Gauss_Legendre& Gauss_Legendre::operator=(const Gauss_Legendre& other)
{
    resize(other.size());
    return *this;
}

Gauss_Legendre& Gauss_Legendre::operator=(Gauss_Legendre&& other) noexcept
{
    _size = other.size();
    table = std::move(other.table);
    return *this;
}

std::pair<double,double> Gauss_Legendre::point(double lower, double upper,
        std::size_t i) const
{
    if (i>=_size)
        throw std::out_of_range{"requested value exceeds number of knots"};
    double point{-1}, weight{-1}; // useless values to indicate failure
    gsl_integration_glfixed_point(lower,upper,i,&point,&weight,table.get());
    return std::make_pair(point,weight);
}

double Gauss_Legendre::operator()(Function f, double lower,
        double upper) const
{
    gsl_function wrapper;
    wrapper.function = unwrap;
    // Use the `void*`, which points to parameters usually, to pass a
    // `Function` to the GSL routine. `unwrap` takes care of everything, i.e.
    // applies this function. This way, the GSL routine can be used to
    // integrate both function objects and ordinary functions.
    wrapper.params = &f; 
    return gsl_integration_glfixed(&wrapper,lower,upper,table.get());
}

void Gauss_Legendre::resize(std::size_t s)
{
    _size = s;
    table = std::unique_ptr<gsl_integration_glfixed_table,Glfixed_deleter>{
            gsl_integration_glfixed_table_alloc(_size)};
}

std::size_t Gauss_Legendre::size() const noexcept
{
    return _size;
}

Gauss_Legendre::~Gauss_Legendre() noexcept
{
}


// -- Integration: adaptive routines ------------------------------------------

Qag::Qag(const Settings& set)
: absolute_precision{set.absolute_precision},
    relative_precision{set.relative_precision},
    limit{set.space},
    workspace{set.space}
{
}

Value Qag::operator()(Function f, double lower, double upper) const
{
    int sign{signed_interval(lower,upper) ? 1 : -1};

    gsl_function wrapper;
    wrapper.function = unwrap;
    // Use the `void*`, which points to parameters usually, to pass a
    // `Function` to the GSL routine. `unwrap` takes care of everything, i.e.
    // applies this function. This way, the GSL routine can be used to
    // integrate both function objects and ordinary functions.
    wrapper.params = &f; 

    double result{0.0};
    double error{0.0};

    bool lower_inf{std::isinf(lower)};
    bool upper_inf{std::isinf(upper)};

    if (lower_inf && upper_inf)
        call(gsl_integration_qagi,&wrapper,absolute_precision,
                relative_precision,limit,workspace.data(),&result,&error);
    else if (lower_inf)
        call(gsl_integration_qagil,&wrapper,upper,absolute_precision,
                relative_precision,limit,workspace.data(),&result,&error);
    else if (upper_inf)
        call(gsl_integration_qagiu,&wrapper,lower,absolute_precision,
                relative_precision,limit,workspace.data(),&result,&error);
    else
        call(gsl_integration_qags,&wrapper,lower,upper,absolute_precision,
                relative_precision,limit,workspace.data(),&result,&error);

    return Value{sign*result,error};
}

void Qag::reserve(std::size_t space)
{
    workspace.reserve(space);
    limit = space;
}

Cquad::Cquad(const Settings& set)
: absolute_precision{set.absolute_precision},
    relative_precision{set.relative_precision},
    workspace{set.space}
{
}

Value Cquad::operator()(Function f, double lower, double upper) const
{
    int sign{signed_interval(lower,upper) ? 1 : -1};

    bool lower_inf{std::isinf(lower)};
    bool upper_inf{std::isinf(upper)};

    // `gsl_integration_cquad` does not provide functions for the integration
    // of infinite intervals. Hence, the required change of variables is
    // performed explicitly.
    Function integrand{f};
    if (lower_inf && upper_inf) {
        integrand = [&f](double x){return (f((1-x)/x) + f((x-1)/x)) / (x*x);};
        lower = 0.0;
        upper = 1.0;
    }
    else if (lower_inf) {
        integrand = [&f,upper](double x){return f(upper+(x-1)/x) / (x*x);};
        lower = 0.0;
        upper = 1.0;
    }
    else if (upper_inf) {
        integrand = [&f,lower](double x){return f(lower+(1-x)/x) / (x*x);};
        lower = 0.0;
        upper = 1.0;
    }

    gsl_function wrapper;
    wrapper.function = unwrap;
    wrapper.params = &integrand;
    double result{0.0};
    double error{0.0};

    std::size_t evaluations{0}; // dummy variable for function call
    call(gsl_integration_cquad,&wrapper,lower,upper,absolute_precision,
            relative_precision,workspace.data(),&result,&error,&evaluations);
    return Value{sign*result,error};
}

void Cquad::reserve(std::size_t space)
{
    workspace.reserve(space);
}

// -- Interpolation -----------------------------------------------------------

Interpolate::Interpolate(const std::vector<double>& x,
        const std::vector<double>& y, Interpolation_method m, bool tolerant)
    : x_data{x}, y_data{y}, method{m}, tolerant{tolerant}
{
    if (x_data.size()!=y_data.size())
        throw std::invalid_argument("x and y need to have the same size");
    if (x_data.size()<method.min_size())
        throw std::invalid_argument("not enough data points for the choosen \
interpolation method");
    spline = gsl_interp_alloc(method.get_method(),x_data.size());
    call(gsl_interp_init,spline,x_data.data(),y_data.data(),x.size());
}

Interpolate::Interpolate(const Interpolate& other)
    : x_data{other.x_data}, y_data{other.y_data}, method{other.method},
    tolerant{other.tolerant}
{
    spline = gsl_interp_alloc(method.get_method(),x_data.size());
    call(gsl_interp_init,spline,x_data.data(),y_data.data(),x_data.size());
}

Interpolate::Interpolate(Interpolate&& other)
    : x_data{std::move(other.x_data)}, y_data{std::move(other.y_data)},
    method{other.method}, tolerant{other.tolerant},
    acc{other.acc}, spline{other.spline}
{
    other.acc = nullptr;
    other.spline = nullptr;
}

Interpolate& Interpolate::operator=(const Interpolate& other)
{
    // free old data
    gsl_interp_free(spline);

    // copy
    x_data = other.x_data;
    y_data = other.y_data;
    method = other.method;
    tolerant = other.tolerant;
    spline = gsl_interp_alloc(method.get_method(),x_data.size());
    call(gsl_interp_init,spline,x_data.data(),y_data.data(),x_data.size());

    return *this;
}

Interpolate& Interpolate::operator=(Interpolate&& other)
{
    // free old data
    gsl_interp_free(spline);

    // copy
    x_data = std::move(other.x_data);
    y_data = std::move(other.y_data);
    method = other.method;
    spline = other.spline;
    tolerant = other.tolerant;

    // empty old object
    other.acc = nullptr;
    other.spline = nullptr;

    return *this;
}

double Interpolate::operator()(double x) const
{
    double result{};
    if(tolerant) {
        if (x<front())
            x = front();
        else if (x>back())
            x = back();
    }
    call(gsl_interp_eval_e,spline,x_data.data(),y_data.data(),x,acc,&result);
    return result;
}

Interpolate::~Interpolate() noexcept
{
    gsl_interp_free(spline);
    gsl_interp_accel_free(acc);
}
 
Interpolation_method::Interpolation_method(Interpolation_method::Method m)
{
    switch (m) {
        case Interpolation_method::linear:
            method =  gsl_interp_linear;
            break;
        case Interpolation_method::polynomial:
            method =  gsl_interp_polynomial;
            break;
        case Interpolation_method::cubic:
            method =  gsl_interp_cspline;
            break;
        case Interpolation_method::cubic_periodic:
            method =  gsl_interp_cspline_periodic;
            break;
        case Interpolation_method::akima:
            method =  gsl_interp_akima;
            break;
        case Interpolation_method::akima_periodic:
            method =  gsl_interp_akima_periodic;
            break;
        case Interpolation_method::steffen:
            method =  gsl_interp_steffen;
            break;
    }
}

Interpolate sample(const Function& f, const Interval& i,
        Interpolation_method m, bool tolerant)
{
    std::vector<double> y_values(i.size());
    std::transform(i.cbegin(),i.cend(),y_values.begin(),f);

    return Interpolate{i,y_values,m,tolerant};
}
} // gsl
