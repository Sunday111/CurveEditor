#pragma once

#include "Vector.h"

template<typename T, size_t dimensions>
class TunableFunction
{
public:
    using Vector = Vector<T, dimensions>;

    virtual size_t GetControlPointsCount() const = 0;
    virtual bool GetControlPoint(size_t index, const Vector** pp) const = 0;
    virtual bool SetControlPoint(size_t index, const Vector& val) = 0;
    virtual ~TunableFunction() = default;
};