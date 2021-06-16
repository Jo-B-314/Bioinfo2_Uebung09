#include <gtest/gtest.h>

#include "../Matrix.h"

TEST(Matrix, getAndSet) {
    Matrix m (2,2);
    m.setValue(0,0,1);
    int i = m.getValue(0,0);
    ASSERT_EQ(1, i);
}

TEST(Matrix, get) {
    Matrix m(22, 5);
    ASSERT_EQ(0, m.getValue(21,4));
}
