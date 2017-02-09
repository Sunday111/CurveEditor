#pragma once

#include <array>
#include "Vector.h"

namespace _bspline
{
    template<typename T>
    T B(size_t i, size_t k, T x, const T* t)
    {
        if(k == 0)
        {
            return ((t[i] <= x) && (x < t[i + 1])) ? 1 : 0;
        }
        else
        {
            const T b[]
            {
                B(i,     k - 1, x, t),
                B(i + i, k - 1, x, t)
            };

            const T n[]
            {
                x - t[i],
                t[i + k + 1] - x
            };

            const T d[]
            {
                t[i + k] - t[i],
                t[i + k + 1] - t[i + 1]
            };

            T res = 0;

            for(size_t idx = 0; idx < 2; ++idx)
            {
                res += b[i] * n[i] / d[i];
            }

            return res;
        }
    }

    template<typename T, size_t dims>
    Vector<T, dims> BSpline(
        const Vector<T, dims>* points,
        size_t pointsCount,
        const T* knots,
        size_t knotsCounts,
        T param)
    {
        Vector<T, dims> res;

        for(size_t i = 0; i < dims; ++i)
        {
            res.coords[i] = 0;
        }

        for(size_t i = 0; i < pointsCount; ++i)
        {
            res += points[i] * B<T>(i, knotsCounts, param, knots);
        }

        return res;
    }
}