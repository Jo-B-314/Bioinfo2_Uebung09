#ifndef BIOINFO_MATRIX_H
#define BIOINFO_MATRIX_H

#include <vector>

class Matrix {
  
public: 
    
    using Array = std::vector<int>; //row
    using Matrixtype = std::vector<Array>; //column
    
    /**
     * default constructor 
     */
    Matrix(int size_s,int size_p);
    
    /**
     * sets the given value at pos Dij
     */
    void setValue(int i, int j, int value);
    
    /**
     * returns the value at pos Dij
     */
    int getValue(int i, int j);
    
    /**
     * returns the number of rows
     */
    size_t getRowNum();
    
    /**
     * returns the number of columns
     */
    size_t getColumnNum();
    
private:
    
    Matrixtype matrix_;
};

Matrix::Matrix(int size_s, int size_p) {
    matrix_.resize(size_p);
    //with resize we have something in our arrays so we have no undefined behavior in getValue
    for (int i = 0; i < size_p; i++) {
        Matrix::Array a;
        a.resize(size_s);
        matrix_[i]=a;
    }
}

void Matrix::setValue(int i, int j, int value) {
    matrix_[i][j] = value;
}

int Matrix::getValue(int i, int j) {
    return matrix_[i][j];
}

#endif
