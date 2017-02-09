#pragma once

#include "TunableParametricFunction.h"

enum class HyperbolaDataIndex
{
    StartPoint,
    FirstAxis,
    SecondAxis,
    __Last
};

template<typename T, size_t dimensions>
class Hyperbola :
    public TunableParametricFunctionT<T, dimensions, HyperbolaDataIndex>
{
public:
    explicit Hyperbola(
        T a, T b,
        const Vector& start,
        const Vector& ix,
        const Vector& iy)
    {
        m_pts[Index<HyperbolaDataIndex::StartPoint>()] = start;
        m_pts[Index<HyperbolaDataIndex::FirstAxis>()] = start + ix * a;
        m_pts[Index<HyperbolaDataIndex::SecondAxis>()] = start + iy * b;
    }

    Vector Evaluate(T t) override
    {
        t = t * 2 - 1;

        const auto s = m_pts[Index<HyperbolaDataIndex::StartPoint>()];

        Vector res = s;

        res += (m_pts[Index<HyperbolaDataIndex::FirstAxis>()] - s) * std::cosh(t);
        res += (m_pts[Index<HyperbolaDataIndex::SecondAxis>()] - s) * std::sinh(t);

        return res;
    }
};