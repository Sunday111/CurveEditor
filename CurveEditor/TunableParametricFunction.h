#pragma once

#include "ParametricFunction.h"
#include "TunableFunction.h"

template<typename T, size_t dimensions>
class TunableParametricFunction :
    public ParametricFunction<T, dimensions>,
    public TunableFunction<T, dimensions>
{
public:
    using Vector = typename ParametricFunction<T, dimensions>::Vector;
};

template<typename T, size_t dimensions, typename _enum>
class TunableParametricFunctionT :
    public TunableParametricFunction<T, dimensions>
{
protected:
    template<_enum e>
    constexpr size_t Index() const { return static_cast<size_t>(e); }

public:
    size_t GetControlPointsCount() const override
    {
        return Index<_enum::__Last>();
    }

    bool GetControlPoint(size_t index, const Vector** pp) const override
    {
        if(index >= Index<_enum::__Last>())
        {
            return false;
        }

        *pp = m_pts + index;

        return true;
    }

    bool SetControlPoint(size_t index, const Vector& val) override
    {
        if(index >= Index<_enum::__Last>())
        {
            return false;
        }

        m_pts[index] = val;

        return true;
    }

protected:
    Vector m_pts[static_cast<size_t>(_enum::__Last)];
};