#include <stdexcept>

#include "crpropa/electricField/ElectricField.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(testUniformElectricField, SimpleTest) {
	UniformElectricField E(Vector3d(-1, 5, 3));
	Vector3d e = E.getField(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(e.getX(), -1);
	EXPECT_DOUBLE_EQ(e.getY(), 5);
	EXPECT_DOUBLE_EQ(e.getZ(), 3);
}

TEST(testElectricFieldList, SimpleTest) {
	ElectricFieldList E;
	// Test a list of three electric fields
	E.addField(new UniformElectricField(Vector3d(1, 0, 0)));
	E.addField(new UniformElectricField(Vector3d(0, 2, 0)));
	E.addField(new UniformElectricField(Vector3d(0, 0, 3)));
	Vector3d e = E.getField(Vector3d(0.));
	EXPECT_DOUBLE_EQ(e.x, 1);
	EXPECT_DOUBLE_EQ(e.y, 2);
	EXPECT_DOUBLE_EQ(e.z, 3);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}