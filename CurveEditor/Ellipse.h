#pragma once

#include <cmath>

#include "TunableParametricFunction.h"

enum class EllipseDataIndex
{
    Center,
    FirstAxis,
    SecondAxis,
    __Last
};

template<typename T, size_t dimensions>
class Ellipse :
    public TunableParametricFunctionT<T, dimensions, EllipseDataIndex>
{
public:
    explicit Ellipse(
        T r1, T r2,
        const Vector& center,
        const Vector& ix,
        const Vector& iy)
    {
        m_pts[Index<EllipseDataIndex::Center>()] = center;
        m_pts[Index<EllipseDataIndex::FirstAxis>()] = center + ix * r1;
        m_pts[Index<EllipseDataIndex::SecondAxis>()] = center + iy * r2;
    }

    Vector Evaluate(T t) override
    {
        constexpr const T pi = 3.14159265358979323846;
        t *= pi * 2;

        const auto c = m_pts[Index<EllipseDataIndex::Center>()];

        Vector res = c;

        res += (m_pts[Index<EllipseDataIndex::FirstAxis>()] - c) * std::cos(t);
        res += (m_pts[Index<EllipseDataIndex::SecondAxis>()] - c) * std::sin(t);

        return res;
    }
};