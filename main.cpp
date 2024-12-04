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
#include "iterative_methods.hpp"

enum
{
    SCALE = 1000,
    BASE = 10
};

int
main(int argc, char** argv)
{
    const size_t SIZE = std::strtol(argv[1], nullptr, BASE); // Matrix size
    Matrix tmp(SIZE, SIZE); 
    std::cin >> tmp; // read matrix from file
    ScalarMatrix I(SIZE, 1);
    Matrix A = I + tmp, Q(SIZE, SIZE), R(SIZE, SIZE); // matrices

    std::vector<double> x = generate_random(-1, 1, SIZE); // generate random vector
    std::vector<double> f = A * x; // right side of Ax = f system
    QRGram_Shmidt(Q, R, A); // Gram-Shmidt QR-decomposition

    Q.transpose();
    std::vector<double> f_n = Q * f;
    std::vector<double> x_f = solve_triangular_system(R, f_n); // solve system
    
    double f_error = norm(f - A * x_f); // calculate f error 
    double x_error = norm(x - x_f); // calculate x error
    std::cout << "Gram-Shmidt: \n";
    std::cout << "||x - x_f|| = " <<  x_error << std::endl;
    std::cout << "||f - Ax_f|| = " << f_error << std::endl;
    std::cout << "||x - x_f|| / ||x|| = " << x_error / norm(x) << std::endl;

    Matrix_circles circles(A);
    std::cout << "\nSpectrum range:\n";
    std::cout << circles.get_range() << "\n\n"; // Matrix spectrum range
    double a = circles.get_range().first, b = circles.get_range().second; // Chebyshev method parameters

    std::vector<double> y0(SIZE, 0), y(SIZE, 0);
    size_t count = 1;
    std::vector<size_t> mask = {1};
    while (norm(x - y) >= x_error) {
        count *= 2;
        std::vector<size_t> new_mask(count);
        for (size_t k = 0; k < count; k++) {
            if (k % 2 == 0) {
                new_mask[k] = mask[k / 2];
            } else {
                new_mask[k] = 2 * count - mask[k / 2]; 
            }
        }
        mask = new_mask;
        y = Chebyshev_solve(A, f, y0, a, b, count, mask);
    }
    f_error = norm(f - A * y); // calculate f error 
    x_error = norm(x - y); // calculate x error
    std::cout << "Chebyshev: \n";
    std::cout << "||x - y|| = " <<  x_error << std::endl;
    std::cout << "||f - Ay|| = " << f_error << std::endl;
    std::cout << "||x - y|| / ||x|| = " << x_error / norm(x) << std::endl;
    std::cout << "Iterations count = " << count << std::endl;

    std::vector<double> errors = Chebyshev_anylyze(A, f, y0, a, b, count, mask, x); // iteration errors
    #ifdef DISPLAY_ERRORS
        std::cout << errors;
    #endif
    return 0;
}