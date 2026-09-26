#pragma once

#include <vector>
namespace Maths {

    struct SplineCoefficients
    {
        double a, b, c, d;
    };

    std::vector<SplineCoefficients> BuildSpline(
        const std::vector<double>& x,
        const std::vector<double>& y);


    double CubicSplineInterpolation(
        const std::vector<double>& x,
        const std::vector<SplineCoefficients>& coeffs,
        double X);

    double LinearInterpolation(double x0, double y0,
        double x1, double y1,
        double x);


    double BilinearInterpolation(
        double x0, double x1,
        double y0, double y1,
        double Q00, double Q10,
        double Q01, double Q11,
        double x, double y);

}



