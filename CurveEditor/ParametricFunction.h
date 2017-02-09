#pragma once

#include "Vector.h"

template<typename T, size_t dimensions>
class ParametricFunction
{
public:
    using Vector = Vector<T, dimensions>;
    virtual Vector Evaluate(T t) = 0;
    virtual ~ParametricFunction() = default;
};