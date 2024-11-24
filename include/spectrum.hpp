#include "matrix.hpp"
#include <vector>

class Gershgorin_circle
{
    double radius, center;
public:
    Gershgorin_circle(double r, double c): radius(r), center(c) {}    
    friend std::ostream &operator<<(std::ostream &os, const Gershgorin_circle& B);
};

class Matrix_circles
{
    std::vector<Gershgorin_circle> arr;
    std::pair<double, double> range;
public:
    friend std::ostream &operator<<(std::ostream &os, const Matrix_circles& v);
    Gershgorin_circle operator[](size_t i) const { return arr[i]; }
    Gershgorin_circle &operator[](size_t i) { return arr[i]; }
    size_t size() const { return arr.size(); }
    Matrix_circles(const Matrix& A)
    {
        range.first = 0, range.second = 0;
        double center, radius;
        for (size_t i = 0; i < A.rows; i++) {
            center = A(i, i);
            radius = 0;
            for (size_t j = 0; j < A.cols; j++) {
                if (i != j) [[likely]] {
                    radius += std::abs(A(i, j));
                }
            }
            range.first = std::min(range.first, center - radius);
            range.second = std::max(range.second, center + radius);
            arr.push_back(Gershgorin_circle(radius, center));
        }
    }

};

std::ostream 
&operator<<(std::ostream &os, const Gershgorin_circle& B)
{
    os << "Center = " << B.center <<  ", Radius = " << B.radius;
    return os;
}


std::ostream 
&operator<<(std::ostream &os, const Matrix_circles& v) // cout overloading
{
    for (size_t i = 0; i < v.size(); i++) {
        os << "Circle " << i + 1 << ": " << v[i] << std::endl;  
    }
    os << "Range: [" << v.range.first << ", " << v.range.second << ']';
    return os;
}