#include "matrix.hpp"

double 
dot_product(const std::vector<double> &x, const std::vector<double> &y, const size_t &shift) // vector dot_product
{
    double result = 0;
    auto j = y.begin() + shift;
    for (auto i = x.begin() + shift; i < x.end(); i++, j++) {
        result += (*i) * (*j);
    }
    return result;
}

double
norm(const std::vector<double> &x, const size_t &shift) // second Holder norm (||.||_2)
{
    return std::sqrt(dot_product(x, x, shift));
}

double 
mean_square_norm(const std::vector<double> &x) // mean square vector norm
{
    return norm(x) / std::sqrt(x.size());
}

std::ostream 
&operator<<(std::ostream &os, const Matrix& A) // cout overloading
{
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            os << A(i, j) << ' ';
        }
        os << std::endl;
    }
    return os;
}

std::ostream 
&operator<<(std::ostream &os, const std::vector<double>& v)
{
    for (auto i: v) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
    return os;
}

std::ostream 
&operator<<(std::ostream &os, const std::vector<size_t>& v)
{
    for (auto i: v) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
    return os;
}

std::istream 
&operator>>(std::istream &is, Matrix& A) // cin overloading
{
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            is >> A(i, j);
        }
    }
    return is;
}

std::vector<double> 
operator-(const std::vector<double> &self, const std::vector<double> &other) // "x - y" overloading
{
    std::vector<double> res = self;
    res -= other;
    return res;
}


void
operator-=(std::vector<double> &self, const std::vector<double> &other) // "x -= y" overloading
{
    if (self.size() == other.size()) {
        for (size_t i = 0; i < other.size(); i++) {
            self[i] -= other[i];
        }
    }
}

std::vector<double> 
operator+(const std::vector<double> &self, const std::vector<double> &other) // "x + y" overloading
{
    std::vector<double> res = self;
    res += other;
    return res;
}


void
operator+=(std::vector<double> &self, const std::vector<double> &other) // "x += y" overloading
{
    if (self.size() == other.size()) {
        for (size_t i = 0; i < other.size(); i++) {
            self[i] += other[i];
        }
    }
}

std::vector<double>
operator*(const Matrix &A, const std::vector<double> &x) // "A * v" overloading
{
    std::vector<double> result(x.size(), 0);
    for (size_t i = 0; i < A.rows; i++) {
        double sum = 0;
        for (size_t j = 0; j < A.cols; j++) {
            sum += x[j] * A(i, j);
        }    
        result[i] = sum;
    }
    return result;
}

const Matrix 
operator*(const Matrix &A, const Matrix &B)
{
    Matrix C(A.rows, A.cols);
    for (size_t i = 0; i < B.cols; i++) {
        C[i] = A * B[i];
    }
    return C;
}

std::vector<double> 
operator*(const std::vector<double> x, const double &a) // "v * const" overloading
{
    std::vector<double> result = x;
    result *= a;
    return result;
}

void
operator*=(std::vector<double> &self, const double &a)
{
    for (size_t i = 0; i < self.size(); i++) {
        self[i] *= a;
    }   
}

std::vector<double> 
operator/(const std::vector<double> x, const double &a) // "v / const" overloading
{
    std::vector<double> result = x;
    for (size_t i = 0; i < x.size(); i++) {
        result[i] /= a;
    }
    return result;
}

std::vector<double> 
solve_triangular_system(const Matrix &R, const std::vector<double> &f) // Gauss method
{
    size_t n = f.size() - 1;
    std::vector<double> x(f.size());
    for (int m = n; m >= 0; m--) {
        double s = 0;
        for (size_t j = m + 1; j <= n; j++) {
            s += R(m, j) * x[j];
        }
        x[m] = (f[m] - s) / R(m, m);
    }
    return x;
}
