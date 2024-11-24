#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <limits>
#include <chrono>
#include <random>


#include "matrix.hpp"
#include "Gram-Shmidt.hpp"
#include "generate_rand.hpp"
#include "spectrum.hpp"

enum
{
    SCALE = 1000,
    BASE = 10
};

int
main(int argc, char** argv)
{
    const size_t SIZE = std::strtol(argv[1], nullptr, BASE);
    Matrix tmp(SIZE, SIZE); 
    std::cin >> tmp; // read matrix from file
    ScalarMatrix I(SIZE, 1);
    Matrix A = I + tmp, Q(SIZE, SIZE), R(SIZE, SIZE); // matrices
    std::vector<double> x = generate_random(-1, 1, SIZE); // generate random vector
    std::vector<double> f = A * x; // right side of Ax = f system
    QRGram_Shmidt(Q, R, A); // Gram-Shmidt QR-decomposition
    Q.transpose();
    std::vector<double> f_n = Q * f;
    std::vector<double> x_f = solve_triangular_system(R, f_n);
    double f_error = norm(f - A * x_f); // calculate f error 
    double x_error = norm(x - x_f); // calculate x error
    std::cout << "||x - x_f|| = " <<  x_error << std::endl;
    std::cout << "||f - Ax_f|| = " << f_error << std::endl;
    std::cout << '\n' << Matrix_circles(A).get_range() << '\n';
    return 0;
}