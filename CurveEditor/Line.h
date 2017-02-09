#pragma once

#include "TunableParametricFunction.h"

enum class LineDataIndex
{
    StartPoint,
    EndPoint,
    __Last
};

template<typename T, size_t dimensions>
class Line :
    public TunableParametricFunctionT<T, dimensions, LineDataIndex>
{
public:
    explicit Line(
        const Vector& startPoint,
        const Vector& endPoint)
    {
        m_pts[Index<LineDataIndex::StartPoint>()] = startPoint;
        m_pts[Index<LineDataIndex::EndPoint>()] = endPoint;
    }

    Vector Evaluate(T t) override
    {
        return
            m_pts[Index<LineDataIndex::StartPoint>()] * (1 - t) +
            m_pts[Index<LineDataIndex::EndPoint>()] * t;
    }
};