#include "matrix.hpp"
#include <vector>
#include <cmath>

double
tau(const double &t0, const double &r0, const size_t m, const size_t n)
{
    const double pi = 3.14159265;
    return t0 / (1 + r0 * std::cos(pi * (n - 1) / (2 * m)));
}

std::vector<double>
Chebyshev_solve(const Matrix& A, const std::vector<double>& f,
const std::vector<double>& y0, double a, double b, size_t m, 
const std::vector<size_t> &mask) // Chebyshev method
{
    double t0 = 2 / (a + b); // approximation for optimal tau from simple-iteration method
    double cond = b / a; // approximation for cond(A)
    double r0 = (1 - 1 / cond) / (1 + 1 / cond); // approximation for ro from simple-iteration method
    double t_k;
    std::vector<double> y(y0);
    for (size_t k = 0; k < m; k++) { // iterations
        t_k = tau(t0, r0, m, mask[k]); // count parameter
        y = (f - A * y) * t_k + y; // update y
    }
    return y;
}

