#include "Integral.hpp"
namespace Maths {


    double integral(std::function<double(double)> f,
        double a, double b, int n)
    {
        double h = (b - a) / n;
        double sum = 0.0;

        for (int i = 0; i < n; i++)
        {
            double x1 = a + i * h;
            double x2 = a + (i + 1) * h;

            sum += (f(x1) + f(x2)) * h / 2.0;
        }

        return sum;
    }
}