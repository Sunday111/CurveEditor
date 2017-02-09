#pragma once

#include "TunableParametricFunction.h"

enum class ParabolaDataIndex
{
    StartPoint,
    FirstAxis,
    SecondAxis,
    __Last
};

template<typename T, size_t dimensions>
class Parabola :
    public TunableParametricFunctionT<T, dimensions, ParabolaDataIndex>
{
public:
    explicit Parabola(
        T p,
        const Vector& start,
        const Vector& ix,
        const Vector& iy)
    {
        m_pts[Index<ParabolaDataIndex::StartPoint>()] = start;
        m_pts[Index<ParabolaDataIndex::FirstAxis>()] = start + ix * p;
        m_pts[Index<ParabolaDataIndex::SecondAxis>()] = start + iy * p;
    }

    Vector Evaluate(T t) override
    {
        t = t * 2 - 1;

        return
            m_pts[Index<ParabolaDataIndex::StartPoint>()] +
            m_pts[Index<ParabolaDataIndex::FirstAxis>()] * t * t * 0.5 +
            m_pts[Index<ParabolaDataIndex::SecondAxis>()] * t;
    }
};