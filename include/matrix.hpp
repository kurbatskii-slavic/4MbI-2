#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <limits>
#ifndef MATRIX_HPP
#define MATRIX_HPP


struct Matrix;
struct HouseholderMatrix;
struct ScalarMatrix;

double norm(const std::vector<double> &x, const size_t &shift=0);
std::vector<double> operator*(const Matrix &A, const std::vector<double> &x);
const Matrix operator*(const Matrix &A, const Matrix &B);
std::vector<double> operator-(const std::vector<double> &self, const std::vector<double> &other);
void operator-=(std::vector<double> &self, const std::vector<double> &other);
std::vector<double> operator+(const std::vector<double> &self, const std::vector<double> &other);
void operator+=(std::vector<double> &self, const std::vector<double> &other);
std::vector<double> operator*(const std::vector<double> x, const double &a);
void operator*=(std::vector<double> &self, const double &a);
std::vector<double> operator/(const std::vector<double> x, const double &a);

struct HouseholderMatrix // struct for Householder matrices
{
    std::vector<double> w; // normal vector
    size_t shift = 0;
    int coeff = 1;
    HouseholderMatrix(){};
    HouseholderMatrix(const std::vector<double> &v): w(v / norm(v)) {} // constructor
};

struct Matrix // struct for matrices
{
    std::vector<std::vector<double>> arr; // array of columns (convenient for calculations)
    size_t rows, cols; // dimensions
    bool transposed = false;
    Matrix(): rows(0), cols(0) {}
    Matrix(size_t m, size_t n) : rows(m), cols(n), arr(n, std::vector<double>(m, 0)) {} // constructors
    double &operator()(size_t i, size_t j) { return transposed ? arr[i][j] : arr[j][i]; } // element access
    double operator()(size_t i, size_t j) const { return transposed ? arr[i][j] : arr[j][i]; } // element access
    Matrix operator-(const Matrix &B) const {
        Matrix C = *this;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                C(i, j) -= B(i, j);
            }    
        }
        return C;
    }
    void transpose() { transposed = !transposed; }
    std::vector<double> &operator[](size_t i) { return arr[i]; } // get i-column
    std::vector<double> operator[](size_t i) const { return arr[i]; } // get i-column
};

struct ScalarMatrix
{
    size_t n;
    double a;
    ScalarMatrix(): n(0), a(0){}
    ScalarMatrix(size_t n, double a): n(n), a(a) {}
    double get_a() const { return a; };
    Matrix operator+(const Matrix &A) const {
        Matrix C = A;
        for (size_t i = 0; i < A.cols; i++) {
            C(i, i) += a;
        }
        return C;
    }
};

double dot_product(const std::vector<double> &x, const std::vector<double> &y, const size_t &shift=0);
std::ostream &operator<<(std::ostream &os, const Matrix& A);
std::ostream &operator<<(std::ostream &os, const std::vector<double>& v);
std::ostream &operator<<(std::ostream &os, const HouseholderMatrix& T);
std::istream &operator>>(std::istream &is, Matrix& A);
double maximum_norm(const std::vector<double> &x);
double matrix_maximum_norm(const Matrix &A);


std::vector<double> solve_triangular_system(const Matrix &R, const std::vector<double> &f);
void matvec(const HouseholderMatrix &H, std::vector<double> &v);


void make_reflection(HouseholderMatrix &H, const std::vector<double> &x, const std::vector<double> &y, const size_t &shift);
void get_reflection(std::vector<double> &ref, const std::vector<double> &x, const size_t &shift);


#endif
