#include "gtest/gtest.h"
#include "omnes.h"
#include "test_common.h"

double phase(double s)
{
    return 1.0 + 2.0 / s;
}

TEST(Omnes, Call)
{
    const auto omnes_function{omnes::Omnes(phase, 4.0, 1e-10)};
    expect_near(omnes_function(0.0), {1.0, 0.0});
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
