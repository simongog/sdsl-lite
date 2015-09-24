#include "gtest/gtest.h"
@SDSL_INCLUDE_ALL@

namespace
{
// The fixture for testing the compilation of all header files.
class CompileTest : public ::testing::Test { };

//! Test constructors
TEST_F(CompileTest, Compile) { }

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
