#include <stdexcept>

#include "crpropa/electricField/ElectricField.h"

#include "gtest/gtest.h"

using namespace crpropa;

// simple test for uniform electric fields copied from magnetic fields
TEST(testUniformElectricField, SimpleTest) {
	UniformElectricField E(Vector3d(-1, 5, 3));
	Vector3d e = E.getField(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(e.getX(), -1);
	EXPECT_DOUBLE_EQ(e.getY(), 5);
	EXPECT_DOUBLE_EQ(e.getZ(), 3);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}