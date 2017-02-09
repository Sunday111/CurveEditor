#pragma once


template<class T, size_t dimensions>
class Vector
{
public:

    explicit Vector(const T(&_coords)[dimensions])
    {
        for(int i = 0; i < dimensions; ++i)
        {
            coords[i] = _coords[i];
        }
    }

    explicit Vector()
    {
    }

    Vector Cross(const Vector& v)
    {

    }

    T DistanceTo(const Vector& p) const
    {
        return (*this - p).Length();
    }

    T Length() const
    {
        T result = 0;

        for(int i = 0; i < dimensions; ++i)
        {
            result += coords[i] * coords[i];
        }

        return result;
    }

    Vector operator* (const T& val) const
    {
        Vector result = *this;

        for(int i = 0; i < dimensions; ++i)
        {
            result.coords[i] *= val;
        }

        return result;
    }

    Vector& operator*= (const T& val)
    {
        for(int i = 0; i < dimensions; ++i)
        {
            coords[i] *= val;
        }

        return *this;
    }

    Vector operator+ (const Vector& p) const
    {
        Vector result = *this;

        for(int i = 0; i < dimensions; ++i)
        {
            result.coords[i] += p.coords[i];
        }

        return result;
    }

    Vector& operator+= (const Vector& p)
    {
        for(int i = 0; i < dimensions; ++i)
        {
            coords[i] += p.coords[i];
        }

        return *this;
    }

    Vector operator-() const
    {
        Vector result;

        for(int i = 0; i < dimensions; ++i)
        {
            result.coords[i] = -coords[i];
        }

        return result;
    }

    Vector operator- (const Vector& p) const
    {
        Vector result(*this);

        for(int i = 0; i < dimensions; ++i)
        {
            result.coords[i] -= p.coords[i];
        }

        return result;
    }

    Vector operator/ (const T& val) const
    {
        Vector result(*this);

        for(int i = 0; i < dimensions; ++i)
        {
            result.coords[i] /= val;
        }

        return result;
    }

    T coords[dimensions];
};