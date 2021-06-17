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
     * resizes the matrix
     */
    void resize(int new_p);

    
    
private:
    
    Matrixtype matrix_;
};

Matrix::Matrix(int size_s, int size_p) {
    matrix_.resize(size_p);
    //we initialize our Matrix with zeros
    for (int i = 0; i < size_p; i++) {
        Matrix::Array a;
        a.resize(size_s);
        for (int j = 0; j < size_s; j++) {
            a[j] = 0;
        }
        matrix_[i]=a;
    }
}

void Matrix::resize(int new_p) {
    int old_p = (int) matrix_.size();
    matrix_.resize(new_p);
    for (int i = old_p; i < new_p; i++) {
        Matrix::Array a;
        int s = (int) matrix_[0].size();
        a.resize(s);
        for (int j = 0; j < s; j++) {
            a[j] = 0;
        }
        matrix_[i] = a;
    }
}

void Matrix::setValue(int i, int j, int value) {
    //has to be like this because otherwise we get a segfault
    matrix_[j][i] = value;
}

int Matrix::getValue(int i, int j) {
    return matrix_[j][i];
}



#endif
