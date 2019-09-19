#include "gsl_interface.h"
#include "gtest/gtest.h"
#include <cmath>
#include <vector>

using gsl::Gauss_Legendre;
using gsl::Interpolate;

void test_integration(const Gauss_Legendre& g)
{
    double calculated{g([](double x){return 1.0;},0.0,1.0)};
    constexpr double tolerance{1e-6};
    EXPECT_NEAR(calculated,1.0,tolerance);
}

TEST(GaussLegendre, CopyConstructor)
{
    constexpr std::size_t size{100};
    Gauss_Legendre g{size};
    Gauss_Legendre g2{g};
    EXPECT_EQ(g.size(),size);
    EXPECT_EQ(g2.size(),size);

    test_integration(g);
    test_integration(g2);
}

TEST(GaussLegendre, MoveConstructor)
{
    constexpr std::size_t size{100};
    Gauss_Legendre g{size};
    Gauss_Legendre g2{std::move(g)};
    EXPECT_EQ(g2.size(),size);
    test_integration(g2);
}

TEST(GaussLegendre, CopyAsignment)
{
    constexpr std::size_t size{100};
    Gauss_Legendre g{size};
    Gauss_Legendre g2{size*2};
    EXPECT_EQ(g.size(),size);
    EXPECT_EQ(g2.size(),size*2);
    g2 = g;
    EXPECT_EQ(g2.size(),g.size());

    test_integration(g);
    test_integration(g2);

}

TEST(GaussLegendre, MoveAsignment)
{
    constexpr std::size_t size{100};
    Gauss_Legendre g{size};
    Gauss_Legendre g2{size*2};
    EXPECT_EQ(g.size(),size);
    EXPECT_EQ(g2.size(),size*2);
    g2 = std::move(g);
    EXPECT_EQ(g2.size(),size);

    test_integration(g2);
}

TEST(GaussLegendre, Size)
{
    constexpr std::size_t size{100};
    Gauss_Legendre g{size};
    EXPECT_EQ(g.size(),size);
}

TEST(GaussLegendre, Integrate)
{
    constexpr std::size_t size{3};
    constexpr double tolerance{1e-2};
    Gauss_Legendre g{size};
    const auto f{[](double x){return 2.0*std::pow(x,5) - x*x + 3.5*x - 1.0;}};
    double calculated{g(f,-2.0,5.0)};
    double expected{5172.42};
    EXPECT_NEAR(calculated,expected,tolerance);

    calculated = g(f,5.0,-2.0);
    EXPECT_NEAR(calculated,-expected,tolerance);
}

TEST(GaussLegendre, Point)
{
    constexpr std::size_t size{2};
    Gauss_Legendre g{size};
    constexpr double lower{-1.0};
    constexpr double upper{1.0};
    std::vector<std::pair<double,double>> points_weights{
        {-1.0/std::sqrt(3.0), 1.0},
        {1.0/std::sqrt(3.0), 1.0},
    };
    for (std::size_t i{0}; i<size; ++i) {
        auto point{g.point(lower,upper,i)};
        auto expected{points_weights[i]};
        EXPECT_DOUBLE_EQ(point.first,expected.first);
        EXPECT_DOUBLE_EQ(point.second,expected.second);
    }
}

TEST(GaussLegendre, Resize)
{
    constexpr std::size_t size{3};
    constexpr std::size_t new_size{70};

    Gauss_Legendre g{size};
    g.resize(new_size);

    EXPECT_EQ(g.size(),new_size);

    test_integration(g);
}

TEST(Interpolate, Sample)
{
    std::vector<double> knots{1,2,3,4,5};
    auto f{[](double x){return 2.0*x;}};
    Interpolate i{gsl::sample(f,knots,gsl::Interpolation_method::linear)};
    double x{2.5};
    EXPECT_DOUBLE_EQ(f(x),i(x));
}

TEST(Interpolate, Throw)
{
    std::vector<double> knots{1};
    auto f{[](double x){return 2.0*x;}};
    ASSERT_THROW(gsl::sample(f,knots,gsl::Interpolation_method::linear),
            std::invalid_argument);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
