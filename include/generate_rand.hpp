#include <random>
#include <vector>

std::vector<double>
generate_random(double left, double right, size_t size)
{
    std::random_device rd; // random numbers generation
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> d(left, right); // numbers in [-1, 1]
    std::vector<double> x(size); // analytical solution
    for (size_t i = 0; i < size; i++) {
        x[i] = d(gen); // create random vector
    }
    return x;
}