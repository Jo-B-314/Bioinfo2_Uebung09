#include <gtest/gtest.h>

#include <Test.hpp>

TEST(TestSuite1, add)
{
	Uebung01::Test test;
	EXPECT_NEAR(test.add(2.4, 4.2), 6.6, 0.000000001);
}
