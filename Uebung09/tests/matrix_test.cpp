#include <gtest/gtest.h>

#include "../Matrix.h"

TEST(Matrix, getAndSet) {
    Matrix::Matrix m (2,2);
    m.setValue(0,0,1);
    int i = m.getValue(0,0);
    ASSERT_EQ(1, i);
}
