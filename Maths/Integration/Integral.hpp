#pragma once
#include <functional>
namespace Maths {


    double integral(std::function<double(double)> f,
        double a, double b, int n);

}