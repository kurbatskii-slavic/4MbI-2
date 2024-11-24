#include "matrix.hpp"
#include <vector>

void
QRGram_Shmidt(Matrix &Q, Matrix &R, const Matrix &A) // QR decomposition using Gramm-Shmidt orthogonalization algorithm
{
    std::vector<double> b_k;
    for (size_t k = 0; k < Q.cols; k++) {
        b_k = A[k];
        for (size_t s = 0; s < k; s++) {
            R(s, k) = dot_product(A[k], Q[s]);
            b_k -= Q[s] * R(s, k);
        }
        double b_k_norm = norm(b_k);
        R(k, k) = b_k_norm;
        Q[k] = b_k / b_k_norm;
    }
}