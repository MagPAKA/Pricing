#include "Interpolation.hpp"

#include <vector>
#include <cmath>
#include <stdexcept>

namespace Maths
{

    double LinearInterpolation(double x0, double y0,
        double x1, double y1, double x)
    {
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }


    double BilinearInterpolation(
        double x0, double x1,
        double y0, double y1,
        double Q00, double Q10,
        double Q01, double Q11,
        double x, double y)
    {
        double R1 = Q00 + (Q10 - Q00) * (x - x0) / (x1 - x0);
        double R2 = Q01 + (Q11 - Q01) * (x - x0) / (x1 - x0);

        double P = R1 + (R2 - R1) * (y - y0) / (y1 - y0);

        return P;
    }


    std::vector<SplineCoefficients> BuildSpline(
        const std::vector<double>& x,
        const std::vector<double>& y)
    {
        int n = x.size();

        if (n < 2 || y.size() != x.size())
            throw std::invalid_argument("invalid data");

        
        std::vector<double> h(n - 1);

        for (int i = 0; i < n - 1; i++)
        {
            h[i] = x[i + 1] - x[i];

            if (h[i] <= 0)
                throw std::invalid_argument("x must be increasing");
        }

       
        std::vector<double> alpha(n, 0.0);
        std::vector<double> l(n, 0.0);
        std::vector<double> mu(n, 0.0);
        std::vector<double> z(n, 0.0);

        for (int i = 1; i < n - 1; i++)
        {
            alpha[i] =
                (3.0 / h[i]) * (y[i + 1] - y[i])
                - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
        }

        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;

        for (int i = 1; i < n - 1; i++)
        {
            l[i] = 2.0 * (x[i + 1] - x[i - 1])
                - h[i - 1] * mu[i - 1];

            mu[i] = h[i] / l[i];

            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n - 1] = 1.0;
        z[n - 1] = 0.0;

        std::vector<double> c(n, 0.0);
        std::vector<double> b(n - 1);
        std::vector<double> d(n - 1);

        for (int j = n - 2; j >= 0; j--)
        {
            c[j] = z[j] - mu[j] * c[j + 1];

            b[j] = (y[j + 1] - y[j]) / h[j]
                - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;

            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }

        std::vector<SplineCoefficients> coeffs(n - 1);

        for (int i = 0; i < n - 1; i++)
        {
            coeffs[i] = {
                y[i],  // a
                b[i],  // b
                c[i],  // c
                d[i]   // d
            };
        }

        return coeffs;
    }

    double CubicSplineInterpolation(
        const std::vector<double>& x,
        const std::vector<SplineCoefficients>& coeffs,
        double X)
    {
        int n = x.size();

        if (n < 2 || coeffs.size() != n - 1)
            throw std::invalid_argument("Donnees invalides");

        // Recherche de l'intervalle [x[i], x[i+1]]
        int i;

        if (X <= x[0])
        {
            i = 0;
        }
        else if (X >= x[n - 1])
        {
            i = n - 2;
        }
        else
        {
            i = 0;

            while (i < n - 2 && X > x[i + 1])
                i++;
        }

        double dx = X - x[i];

        // S_i(x) = a + b*dx + c*dx² + d*dx³
        return coeffs[i].a
            + coeffs[i].b * dx
            + coeffs[i].c * dx * dx
            + coeffs[i].d * dx * dx * dx;
    }

}
