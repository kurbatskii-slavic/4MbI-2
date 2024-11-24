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
&operator<<(std::ostream &os, const HouseholderMatrix& T)
{
    std::cout << T.w;
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

double
maximum_norm(const std::vector<double> &x) // vector maximum_norm
{
    double m = 0;    
    for (auto x_i: x) {
        double abs = x_i > 0 ? x_i : -x_i;
        m = abs > m ? abs : m;
    }
    return m;
}

double
matrix_maximum_norm(const Matrix &A) // matrix maximum norm
{
    double m = 0;
    for (size_t i = 0; i < A.rows; i++) {
        double sum = 0;
        for (size_t j = 0; j < A.cols; j++) {
            double abs = A(i, j) > 0 ? A(i, j) : -A(i, j);
            sum += abs;
        }
        m = sum > m ? sum : m;
    } 
    return m;
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

void 
matvec(const HouseholderMatrix &H, std::vector<double> &v) // reflection operation (H * v)
{
    v -= H.w * (2 * dot_product(H.w, v, H.shift));
    if (H.coeff == -1) {
        for (auto i = v.begin() + H.shift; i < v.end(); i++) {
            (*i) *= H.coeff;
        }
    }
}

void 
make_reflection(HouseholderMatrix &H, const std::vector<double> &x, const std::vector<double> &y, const size_t &shift) // Make H: x -> y
{
    int sign = dot_product(x, y, shift) > 0 ? 1 : -1;
    std::vector<double> v = x + y * sign;
    for (size_t i = 0; i < shift; i++) {
       v[i] = 0;
    }
    H.w = v / norm(v);
    H.coeff = -sign;
    H.shift = shift;
}

void
get_reflection(std::vector<double> &ref, const std::vector<double> &x, const size_t &shift)
{
    ref = std::vector<double>(ref.size(), 0);
    double sub_vec_norm = norm(x, shift);
    for (size_t i = 0; i < shift; i++) {
        ref[i] = x[i];
    } 
    ref[shift] = sub_vec_norm;
}