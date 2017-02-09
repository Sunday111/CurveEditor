#pragma once

#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include <sstream>
#include <unordered_map>

#include "Vector.h"

long long unsigned gcd(long long unsigned a, long long unsigned b)
{
    int c;
    while (a != 0) {
        c = a; a = b%a;  b = c;
    }
    return b;
}

/* The number of combinations from n by i */
inline long double Combinations(size_t n, size_t i)
{
    using lt = long long unsigned;

    std::unordered_map<lt, lt> numerator, denominator;

    auto factorial = [](std::unordered_map<lt, lt>& map, size_t n)
    {
        for (size_t i = 2; i <= n; ++i)
        {
            auto iVal = map.find(i);

            if (iVal == map.end())
            {
                map[i] = 1;
            }
            else
            {
                ++iVal->second;
            }
        }
    };

    factorial(numerator, n);
    factorial(denominator, i);
    factorial(denominator, n - i);

    long double result = 1;

    for (auto iNumerator = numerator.begin(); iNumerator != numerator.end(); ++iNumerator)
    {
        auto iDenominator = denominator.find(iNumerator->first);

        bool updateNumerator = false;
        bool updateDenominator = false;

        if (iDenominator == denominator.end())
        {
            updateNumerator = true;
        }
        else if (iNumerator->second < iDenominator->second)
        {
            iDenominator->second -= iNumerator->second;
            iNumerator->second = 0;

            updateDenominator = true;
        }
        else if (iDenominator->second < iNumerator->second)
        {
            iNumerator->second -= iDenominator->second;
            iDenominator->second = 0;

            updateNumerator = true;
        }
        else
        {
            iNumerator->second = 0;
            iDenominator->second = 0;
        }

        if (updateNumerator)
        {
            result *= std::pow(iNumerator->first, iNumerator->second);
        }

        if (updateDenominator)
        {
            result /= std::pow(iDenominator->first, iDenominator->second);
        }
    }

    return result;

    // n! / (i!(n-i)!)
    //return Factorial(n) / (Factorial(i) * Factorial(n - i));
}

template<class T>
inline T BernsteinPolynomial(size_t n, size_t i, const T& t)
{
    return
        Combinations(n, i) *
        std::pow(t, i) *
        std::pow(1 - t, n - i);
}

class CombinationsCache
{
public:

    const std::vector<size_t>& GetValues(size_t n) const
    {
        auto iCache = m_d.find(n);

        if (iCache == m_d.end())
        {
            auto& values = m_d[n];

            values.reserve(n + 1);

            for (size_t i = 0; i < n + 1; ++i)
            {
                values.push_back(Combinations(n, i));
            }

            return values;
        }

        return iCache->second;
    }

private:
    mutable std::unordered_map<size_t, std::vector<size_t>> m_d;
};

template<class T>
inline T BernsteinPolynomial(size_t n, size_t i, const T& t, CombinationsCache* cache)
{
    const T combinations = cache == nullptr ? Combinations(n, i) : cache->GetValue(n, i);

    return
        combinations *
        std::pow(t, i) *
        std::pow(1 - t, n - i);
}

enum class BezierCurvePointType
{
    Strong_Normal,
    Strong_Smooth,
    Weak
};

enum class BezierCurveStrongPointType
{
    Normal = static_cast<int>(BezierCurvePointType::Strong_Normal),
    Smooth = static_cast<int>(BezierCurvePointType::Strong_Smooth)
};

template<class T, size_t dimensions>
class BezierStrongPoint :
    public Vector<T, dimensions>
{
public:
    using Base = Vector<T, dimensions>;

    explicit BezierStrongPoint() :
        m_pointType(BezierCurveStrongPointType::Normal)
    {}

    explicit BezierStrongPoint(BezierCurveStrongPointType pointType) :
        m_pointType(pointType)
    {}

    explicit BezierStrongPoint(const T(&_coords)[dimensions], BezierCurveStrongPointType pointType) :
        Base(_coords),
        m_pointType(pointType)
    {}

    BezierCurveStrongPointType GetPointType() const { return m_pointType; }

private:
    BezierCurveStrongPointType m_pointType;
};

template<class T, size_t dimensions>
class BezierPoints
{
private:
    enum { First = 0, Last = 1 };

public:
    using Vector = Vector<T, dimensions>;
    using StrongPoint = BezierStrongPoint<T, dimensions>;
    using WalkThroughThePointsFn = bool(*)(size_t pointIndex, BezierCurvePointType pointType, Vector* point, void* userContext);
    using WalkThroughThePointsConstFn = bool(*)(size_t pointIndex, BezierCurvePointType pointType, const Vector* point, void* userContext);

    /* Strong point can be referenced by another segment */
    using StrongPointPtr = std::shared_ptr<StrongPoint>;

    explicit BezierPoints(
            const StrongPointPtr& sp, const Vector& ep,
            size_t order, BezierCurveStrongPointType pointType) :
        m_order(order)
    {
        assert(sp != nullptr);

        m_strongPoints[First] = sp;
        m_strongPoints[Last] = std::make_shared<StrongPoint>(ep.coords, pointType);

        if (order < 2)
        {
            return;
        }

        const size_t innerPoints = order - 1;

        m_weakPoints.resize(innerPoints);
        auto dir = ep - *sp;

        for (size_t i = 0; i < innerPoints; ++i)
        {
            m_weakPoints[i] = (dir * (i + 1)) / (innerPoints + 1) + *sp;
            m_weakPoints[i].coords[1] += dir.coords[1] / 5;
        }
    }
    BezierPoints(const BezierPoints&) = delete;
    BezierPoints(BezierPoints&& uref) :
        m_order(uref.m_order),
        m_weakPoints(std::move(uref.m_weakPoints))
    {
        m_strongPoints[First] = (std::move(uref.m_strongPoints[First]));
        m_strongPoints[Last] = (std::move(uref.m_strongPoints[Last]));
    }

    const size_t PointsCount() const { return m_order + 1; }
    const size_t InnerPoints() const { return m_order - 1; }
    const size_t Order() const { return m_order; }

    T GetParameter(const T& coord) const
    {
        const T range[]
        {
            m_strongPoints[First]->coords[0],
            m_strongPoints[Last]->coords[0]
        };

        if (coord > range[1])
        {
            return -1;
        }

        return (coord - range[0]) / (range[1] - range[0]);
    }

    StrongPointPtr& RFirstPoint() { return m_strongPoints[First]; }
    const StrongPointPtr& RFirstPoint() const { return m_strongPoints[First]; }
    StrongPointPtr& RLastPoint() { return m_strongPoints[Last]; }
    const StrongPointPtr& RLastPoint() const { return m_strongPoints[Last]; }

    bool InsertWeakPoint(const Vector& p)
    {
        for(size_t pointIndex = 0; pointIndex < m_weakPoints.size(); ++pointIndex)
        {
            const Vector& point = m_weakPoints[pointIndex];

            if (p.coords[0] < point.coords[0])
            {
                m_weakPoints.insert(m_weakPoints.begin() + pointIndex, p);
                ++m_order;
                return true;
            }
        }

        if (p.coords[0] < m_strongPoints[Last]->coords[0])
        {
            m_weakPoints.push_back(p);
            ++m_order;
            return true;
        }

        return false;
    }

    int FindPoint(const Vector& p, const T& tol, BezierCurvePointType* pointType = nullptr) const
    {
        assert(m_strongPoints[First] != nullptr && m_strongPoints[Last] != nullptr);

        const size_t pointsCount = PointsCount();

        size_t resIdx = 0;
        T minDist = std::numeric_limits<T>::max();

        for (size_t pointIndex = 0; pointIndex < pointsCount; ++pointIndex)
        {
            const T dist = (*this)[pointIndex].DistanceTo(p);
            if (dist < minDist)
            {
                minDist = dist;
                resIdx = pointIndex;
            }
        }

        const Vector& resPt = (*this)[resIdx];

        for (int coordIdx = 0; coordIdx < dimensions; ++coordIdx)
        {
            if (std::abs(resPt.coords[coordIdx] - p.coords[coordIdx]) > tol)
            {
                return -1;
            }
        }

        if (pointType != nullptr)
        {
            if (resIdx == 0)
            {
                *pointType = static_cast<BezierCurvePointType>(m_strongPoints[First]->GetPointType());
            }
            else if (resIdx == m_order)
            {
                *pointType = static_cast<BezierCurvePointType>(m_strongPoints[Last]->GetPointType());
            }
            else
            {
                *pointType = BezierCurvePointType::Weak;
            }
        }

        return static_cast<int>(resIdx);
    }

    bool GetPointInfo(size_t pointIndex, Vector** point = nullptr, BezierCurvePointType* pointType = nullptr)
    {
        return const_cast<const BezierPoints*>(this)->GetPointInfo(pointIndex, point, pointType);
    }

    bool GetPointInfo(size_t pointIndex, const Vector** point = nullptr, BezierCurvePointType* pointType = nullptr) const
    {
        if (pointIndex > m_order)
        {
            return false;
        }

        if (point != nullptr)
        {
            *point = &operator[](pointIndex);
        }

        if (pointType != nullptr)
        {
            if (pointIndex == 0)
            {
                *pointType = static_cast<BezierCurvePointType>(m_strongPoints[First]->GetPointType());
            }
            else if (pointIndex == m_order)
            {
                *pointType = static_cast<BezierCurvePointType>(m_strongPoints[Last]->GetPointType());
            }
            else
            {
                *pointType = BezierCurvePointType::Weak;
            }
        }

        return true;
    }

    BezierPoints SplitInPoint(const Vector& p, BezierCurveStrongPointType pointType)
    {
        BezierPoints retVal(m_strongPoints[First], p, 1, pointType);
        m_strongPoints[First] = retVal.m_strongPoints[Last];

        size_t pointIndex = 0;

        for (; pointIndex < m_weakPoints.size(); ++pointIndex)
        {
            const Vector& point = m_weakPoints[pointIndex];

            if (point.coords[0] > p.coords[0])
            {
                break;
            }

            retVal.m_weakPoints.push_back(point);
        }

        retVal.m_order = 1 + retVal.m_weakPoints.size();

        m_weakPoints.erase(m_weakPoints.begin(), m_weakPoints.begin() + pointIndex);
        m_order = 1 + m_weakPoints.size();

        return retVal;
    }

    bool RemoveWeakPoint(size_t index)
    {
        if (index == 0 || index == m_order)
        {
            return false;
        }

        --m_order;
        m_weakPoints.erase(m_weakPoints.begin() + index - 1);

        return true;
    }

    bool Merge(BezierPoints& points)
    {
        // Merge with previous segment
        if (points.m_strongPoints[Last] == m_strongPoints[First])
        {
            m_strongPoints[First] = points.m_strongPoints[First];

            points.m_weakPoints.reserve(points.m_weakPoints.size() + m_weakPoints.size());

            for (auto iPoint = m_weakPoints.rbegin();
                iPoint != m_weakPoints.rend(); ++iPoint)
            {
                points.m_weakPoints.push_back(*iPoint);
            }

            m_weakPoints = std::move(points.m_weakPoints);

            m_order = 1 + m_weakPoints.size();

            return true;
        }
        // Merge with next segment
        else if (m_strongPoints[Last] == points.m_strongPoints[First])
        {
            m_strongPoints[Last] = std::move(points.m_strongPoints[Last]);
            m_weakPoints.insert(m_weakPoints.end(),
                points.m_weakPoints.begin(), points.m_weakPoints.end());

            m_order = 1 + m_weakPoints.size();

            return true;
        }
        // Bad segment
        else
        {
            return false;
        }
    }

    BezierPoints& operator=(const BezierPoints&) = delete;

    BezierPoints& operator=(BezierPoints&& uref)
    {
        m_order = uref.m_order;
        m_strongPoints[First] = std::move(uref.m_strongPoints[First]);
        m_weakPoints = std::move(uref.m_weakPoints);
        m_strongPoints[Last] = std::move(uref.m_strongPoints[Last]);

        return *this;
    }

    Vector& operator[](size_t idx)
    {
        const size_t pointsCount = PointsCount();
        assert(idx < pointsCount);

        if (idx == 0)
        {
            return *m_strongPoints[First];
        }

        if (idx == pointsCount - 1)
        {
            return *m_strongPoints[Last];
        }

        return m_weakPoints[idx - 1];
    }

    const Vector& operator[](size_t idx) const
    {
        return const_cast<BezierPoints*>(this)->operator[](idx);
    }

    bool WalkThroughThePoints(size_t startIndex, WalkThroughThePointsFn f, void* userContext)
    {
        size_t weakIdx;

        if (startIndex == 0)
        {
            StrongPointPtr& first = m_strongPoints[First];

            if (f(0, static_cast<BezierCurvePointType>(first->GetPointType()), first.get(), userContext))
            {
                return true;
            }

            weakIdx = 0;
        }
        else
        {
            weakIdx = startIndex - 1;
        }

        for (; weakIdx < m_weakPoints.size(); ++weakIdx)
        {
            if (f(weakIdx + 1, BezierCurvePointType::Weak, &m_weakPoints[weakIdx], userContext))
            {
                return true;
            }
        }

        if (startIndex <= m_order)
        {
            StrongPointPtr& last = m_strongPoints[Last];

            if (f(m_order, static_cast<BezierCurvePointType>(last->GetPointType()), last.get(), userContext))
            {
                return true;
            }
        }

        return false;
    }

    bool WalkThroughThePoints(size_t startIndex, WalkThroughThePointsConstFn f, void* userContext) const
    {
        return const_cast<BezierPoints*>(this)->WalkThroughThePoints(startIndex, 
            reinterpret_cast<WalkThroughThePointsFn>(f), userContext);
    }

    size_t GetOrder() const
    {
        return m_order;
    }

protected:
    BezierPoints() :
        m_strongPoints[Last](std::make_unique<Vector>())
    {}

    size_t m_order;
    //First and last strong point in segment
    StrongPointPtr m_strongPoints[2];
    std::vector<Vector> m_weakPoints;
};

/* This is the class that represents single BezierCurve
 */
template<class T, size_t dimensions>
class BezierSegment
{
public:
    using Points = BezierPoints<T, dimensions>;
    using Vector = typename Points::Vector;
    using StrongPoint = typename Points::StrongPoint;
    using StrongPointPtr = typename Points::StrongPointPtr;

    explicit BezierSegment(const std::shared_ptr<StrongPoint>& start, const Vector& end, size_t order, BezierCurveStrongPointType endPointType, CombinationsCache* cache) :
        m_points(start, end, order, endPointType),
        m_cache(cache)
    {}
    BezierSegment(const BezierSegment&) = delete;
    BezierSegment(BezierSegment&& uref) :
        m_points(std::move(uref.m_points)),
        m_cache(uref.m_cache)
    {}

    Points& RPoints() { return m_points; }
    const Points& RPoints() const { return m_points; }

    /* Split this bezier curve.
     * Produces the new segment which must be placed before this segment
     */
    BezierSegment SplitInPoint(const Vector& p, BezierCurveStrongPointType pointType)
    {
        return BezierSegment(m_points.SplitInPoint(p, pointType), m_cache);
    }

    Vector operator() (const T& t) const
    {
        Vector result;
        result *= 0;

        const size_t pointsCount = m_points.PointsCount();
        const size_t order = m_points.Order();

        auto& values = m_cache->GetValues(order);

        for (size_t i = 0; i < pointsCount; ++i)
        {
            result +=
                m_points[i] *
                values[i] *
                std::pow(t, i) *
                std::pow(1 - t, order - i);
        }

        return result;
    }

    BezierSegment& operator=(const BezierSegment&) = delete;
    BezierSegment& operator=(BezierSegment&& uref)
    {
        m_points = std::move(uref.m_points);
        m_cache = uref.m_cache;
        return *this;
    }

private:
    BezierSegment(Points&& points, CombinationsCache* cache) :
        m_points(std::move(points)),
        m_cache(cache)
    {}

    Points m_points;
    CombinationsCache* m_cache;
};

/* This class is array of bezier curves */
template<class T, size_t dimensions>
class BezierCurve
{
public:
    using Segment = BezierSegment<T, dimensions>;
    using Vector = typename Segment::Vector;
    using StrongPoint = typename Segment::StrongPoint;
    using StrongPointPtr = typename Segment::StrongPointPtr;
    using WalkThroughThePointsFn = bool(*)(
        size_t indexInSegment, size_t indexInCurve, Segment* segment,
        BezierCurvePointType pointType, Vector* point, void* userContext);
    using WalkThroughThePointsConstFn = bool(*)(
        size_t indexInSegment, size_t indexInCurve, const Segment* segment,
        BezierCurvePointType pointType, const Vector* point, void* userContext);

    /* n is order of first bezier curve segment */
    explicit BezierCurve(const Vector& startPoint, const Vector& endPoint, size_t n) :
        m_pointsCount(n + 1)
    {
        //Segment initSegment(&m_startPoint, endPoint, n);
        //m_segments.push_back(initSegment);
        auto dir = endPoint - startPoint;

        Segment initSegment_0(
            std::make_shared<StrongPoint>(startPoint.coords, BezierCurveStrongPointType::Normal),
            startPoint + dir * 0.5, n, BezierCurveStrongPointType::Normal, &m_combintaionsCache);

        Segment initSegment_1(initSegment_0.RPoints().RLastPoint(), endPoint, n,
            BezierCurveStrongPointType::Normal, &m_combintaionsCache);

        m_segments.push_back(std::move(initSegment_0));
        m_segments.push_back(std::move(initSegment_1));

        m_pointsCount = 7;
    }

    StrongPointPtr& RFirstPoint()
    {
        assert(!m_segments.empty());
        return m_segments.front().RPoints().RFirstPoint();
    }

    const StrongPointPtr& RFirstPoint() const
    {
        assert(!m_segments.empty());
        return m_segments.front().RPoints().RFirstPoint();
    }

    size_t GetPointsCount() const { return m_pointsCount; }

    // Get the reference to the first point in the curve
    StrongPointPtr& RLastPoint()
    {
        return m_segments.back().RPoints().RLastPoint();
    }

    // Get the const reference to the first point in the curve
    const StrongPointPtr& RLastPoint() const
    {
        return m_segments.back().RPoints().RLastPoint();
    }

    // Get curve segments
    const std::vector<Segment>& GetSegments() const
    {
        return m_segments;
    }
    
    // Find appropriate segment and insert point there
    template<BezierCurvePointType pointType>
    bool InsertPoint(const Vector& p)
    {
        if (p.coords[0] < RFirstPoint()->coords[0] ||
            p.coords[0] > RLastPoint()->coords[0])
        {
            return false;
        }

        for (size_t segmentIndex = 0; segmentIndex < m_segments.size(); ++segmentIndex)
        {
            if (m_segments[segmentIndex].RPoints().GetParameter(p.coords[0]) > 0)
            {
                InsertPointImpl<pointType>(p, segmentIndex);
                ++m_pointsCount;
                return true;
            }
        }

        return false;
    }

    // Instert strong point in segment
    template<BezierCurvePointType pointType>
    void InsertPointImpl(const Vector& p, size_t segmentIndex)
    {
        m_segments.insert(
            m_segments.begin() + segmentIndex,
            m_segments[segmentIndex].SplitInPoint(p,
                static_cast<BezierCurveStrongPointType>(pointType)));
    }

    // Instert weak point in segment
    template<>
    void InsertPointImpl<BezierCurvePointType::Weak>(const Vector& p, size_t segmentIndex)
    {
        m_segments[segmentIndex].RPoints().InsertWeakPoint(p);
    }

    // Use this curve as the function
    Vector operator()(const T& curveParameter) const
    {
        const T& sx = RFirstPoint()->coords[0];
        const T coord = sx + (RLastPoint()->coords[0] - sx) * curveParameter;

        for (auto& seg : m_segments)
        {
            const T segmentParameter = seg.RPoints().GetParameter(coord);

            if (segmentParameter >= 0)
            {
                return seg(segmentParameter);
            }
        }

        assert(false);
        return Vector();
    }

    /* Get point info by id.
     * You are not able to change point coordinates by the pointer received
     * using this function use MovePoint method to do that
     */
    bool GetPointInfo(size_t pointIndex,
        const Vector** point = nullptr,
        BezierCurvePointType* pointType = nullptr,
        size_t* segmentIndex = nullptr,
        size_t* indexInSegment = nullptr) const
    {
        const size_t segmentsCount = m_segments.size();
        for (size_t _segmentIndex = 0; _segmentIndex < segmentsCount; ++_segmentIndex)
        {
            const Segment& segment = m_segments[_segmentIndex];
            const Segment::Points& points = segment.RPoints();
            const size_t pointsCount = points.PointsCount();

            if (pointIndex < pointsCount)
            {
                if (segmentIndex != nullptr)
                {
                    *segmentIndex = _segmentIndex;
                }

                if (indexInSegment != nullptr)
                {
                    *indexInSegment = pointIndex;
                }

                return points.GetPointInfo(pointIndex, point, pointType);
            }

            pointIndex -= pointsCount - 1;
        }

        return false;
    }

    /* Find point by coodinates and tolerance and get it's info.
     * You are not able to change point coordinates by the pointer received
     * using this function use MovePoint method to do that
     */
    bool GetPointInfo(const Vector& p, const T& tol,
        size_t* segmentIndex = nullptr,
        BezierCurvePointType* pointType = nullptr,
        size_t* indexInSegment = nullptr,
        size_t* indexInCurve = nullptr) const
    {
        size_t prevIndices = 0;

        for(size_t segIdx = 0; segIdx < m_segments.size(); ++segIdx)
        {
            const Segment& segment = m_segments[segIdx];

            const int fnd = segment.RPoints().FindPoint(p, tol, pointType);

            if (fnd != -1)
            {
                if (segmentIndex != nullptr)
                {
                    *segmentIndex = segIdx;
                }

                if (indexInSegment != nullptr)
                {
                    *indexInSegment = fnd;
                }

                if(indexInCurve != nullptr)
                {
                    *indexInCurve = prevIndices + fnd;
                }

                return true;
            }

            prevIndices += segment.RPoints().PointsCount() - 1;
        }

        return false;
    }

    /* Remove point by it's global index */
    bool RemovePoint(size_t index)
    {
        if(index == 0 || index >= m_pointsCount - 1)
        {
            return false;
        }

        const size_t segmentsCount = m_segments.size();

        size_t segmentIndex, indexInSegment;
        if(!GetPointInfo(index, nullptr, nullptr, &segmentIndex, &indexInSegment))
        {
            return false;
        }

        auto& segmentPoints = m_segments[segmentIndex].RPoints();

        if(indexInSegment == 0)
        {
            // Merge this segment points in the pervious one
            m_segments[segmentIndex - 1].RPoints().Merge(segmentPoints);

            // Erase merged segment
            auto eraseIt = m_segments.begin();
            std::advance(eraseIt, segmentIndex);
            m_segments.erase(eraseIt);
        }
        else if(indexInSegment == segmentPoints.PointsCount() - 1)
        {
            // Merge next segment points into this segment
            segmentPoints.Merge(m_segments[segmentIndex + 1].RPoints());

            // Erase merged segment
            auto eraseIt = m_segments.begin();
            std::advance(eraseIt, segmentIndex + 1);
            m_segments.erase(eraseIt);
        }
        else
        {
            // Just remove point from segment
            segmentPoints.RemoveWeakPoint(indexInSegment);
        }

        --m_pointsCount;

        return true;
    }

private:
    // Context for WalkThroughThePoints
    struct WalkContext
    {
        size_t prevPointsCount;
        const Segment* segment;
        void* userContext;
        WalkThroughThePointsFn userFn;
    };

    // Context for MovePoint
    struct MovePointContext
    {
        size_t targetPointIndex;

        struct
        {
            Vector* point;
            BezierCurvePointType type;
        } points[3];
    };

    static bool WalkHelper(size_t pointIndex, BezierCurvePointType pointType, const Vector* point, void* context)
    {
        WalkContext* walkContext = reinterpret_cast<WalkContext*>(context);

        return walkContext->userFn(
            pointIndex,
            walkContext->prevPointsCount + pointIndex,
            const_cast<Segment*>(walkContext->segment),
            pointType,
            const_cast<Vector*>(point),
            walkContext->userContext);
    }

    static bool MovePointHelper(size_t /*indexInSegment*/, size_t indexInCurve, Segment* /*segment*/,
        BezierCurvePointType pointType, Vector* point, void* _context)
    {
        MovePointContext* context = reinterpret_cast<MovePointContext*>(_context);

        int ctxIdx = -1;

        if (context->targetPointIndex == indexInCurve)
        {
            ctxIdx = 1;
        }
        else if (context->targetPointIndex == indexInCurve + 1)
        {
            ctxIdx = 0;
        }
        else if (context->targetPointIndex + 1 == indexInCurve)
        {
            ctxIdx = 2;
        }

        if (ctxIdx != -1)
        {
            context->points[ctxIdx].point = point;
            context->points[ctxIdx].type = pointType;
        }

        return indexInCurve > context->targetPointIndex + 1;
    }

public:
    /* Call caustom function for each point in curve (const overload).
     * Returns true in case when lambda has returned 'true' at least one time
     */
    bool WalkThroughThePoints(WalkThroughThePointsConstFn f, void* context) const
    {
        return const_cast<BezierCurve*>(this)->WalkThroughThePoints(
            reinterpret_cast<WalkThroughThePointsFn>(f), context);
    }

    /* Call caustom function for each point in curve (non-const overload).
     * Returns true in case when lambda has returned 'true' at least one time
     */
    bool WalkThroughThePoints(WalkThroughThePointsFn fn, void* context)
    {
        WalkContext walkContext;
        walkContext.userContext = context;
        walkContext.userFn = fn;
        walkContext.prevPointsCount = 0;

        const size_t segmentsCount = m_segments.size();
        for (size_t segmentIndex = 0; segmentIndex < segmentsCount; ++segmentIndex)
        {
            walkContext.segment = &m_segments[segmentIndex];

            if (walkContext.segment->RPoints().WalkThroughThePoints(
                segmentIndex == 0 ? 0 : 1, WalkHelper, &walkContext))
            {
                return true;
            }

            walkContext.prevPointsCount += walkContext.segment->RPoints().GetOrder();
        }

        return false;
    }

private:

    enum Axis { x, y, z };
    enum Dir { prev, targ, next };
    
    template<Dir>
    struct SyncHelper;

    template<>
    struct SyncHelper<Dir::prev>
    {
        enum { Sign = 1 };
    };

    template<>
    struct SyncHelper<Dir::next>
    {
        enum { Sign = -1 };
    };

    template<bool sync, Dir dir>
    struct IndexHelper
    {
        static bool IndexIsValid(size_t, size_t)
        {
            return false;
        }
    };

    template<>
    struct IndexHelper<true, Dir::prev>
    {
        static bool IndexIsValid(size_t idx, size_t)
        {
            return idx > 1;
        }
    };

    template<>
    struct IndexHelper<true, Dir::next>
    {
        static bool IndexIsValid(size_t idx, size_t cnt)
        {
            return idx + 1 < cnt;
        }
    };

    /* Synchronize neighbors around strong smoothed point
     */
    template<bool sync, Dir dir>
    void SyncNeighbors(const MovePointContext& context)
    {
        constexpr const int sign = SyncHelper<dir>::Sign;

        if(IndexHelper<sync, dir>::IndexIsValid(context.targetPointIndex, m_pointsCount)
            && context.points[dir].type == BezierCurvePointType::Strong_Smooth)
        {
            BezierCurvePointType pointType;
            size_t segmentIndex;
            size_t indexInSegment;

            const size_t ptToSmooth = context.targetPointIndex - sign * 2;

            if(GetPointInfo(ptToSmooth, nullptr, &pointType, &segmentIndex, &indexInSegment) &&
                pointType == BezierCurvePointType::Weak)
            {
                // A is moved point
                const auto& a = *context.points[targ].point;

                // B is strong smooth point
                const auto& b = *context.points[dir].point;

                // C is the point for update
                auto c = m_segments[segmentIndex].RPoints()[indexInSegment];

                /* dAB and dBC are the vectors from a to b and from b to c respectively.
                 * The direction for those vetors determines my compile-time variable 'sign',
                 * but it always looks to the point A
                 */
                const auto dAB = (a - b) * sign;
                const auto dBC = (b - c) * sign;

                for(int coord = 1; coord < dimensions; ++coord)
                {
                    c.coords[coord] = b.coords[coord] - sign * (dBC.coords[x] * dAB.coords[coord]) / dAB.coords[x];
                }

                // Implicit recursive step
                MovePointImpl<true, false>(ptToSmooth, c);
            }
        }
    }

    /* Template helper for MovePoint and SyncNeighbors methods
     */
    template<bool updatePrev, bool updateNext>
    bool MovePointImpl(size_t pointIndex, const Vector& coords)
    {
        MovePointContext context;
        std::memset(&context, 0, sizeof(MovePointContext));
        context.targetPointIndex = pointIndex;

        WalkThroughThePoints(MovePointHelper, &context);

        if (context.points[targ].point == nullptr)
        {
            return false;
        }

        const T tol = 3.0;

        if ((context.points[prev].point != nullptr && coords.coords[x] > context.points[prev].point->coords[x]) &&
            (context.points[next].point != nullptr && coords.coords[x] < context.points[next].point->coords[x]))
        {
            context.points[targ].point->coords[x] = coords.coords[x];
        }

        StrongPointPtr range[] =
        {
            RFirstPoint(),
            RLastPoint()
        };

        for (int i = 1; i < dimensions; ++i)
        {
            auto& to = context.points[targ].point->coords[i];
            to = coords.coords[i];

            if(to < range[0]->coords[i])
            {
                to = range[0]->coords[i];
            }
            else if(to > range[1]->coords[i])
            {
                to = range[1]->coords[i];
            }
        }

        SyncNeighbors<updatePrev, Dir::prev>(context);

        SyncNeighbors<updateNext, Dir::next>(context);

        return true;
    }

public:
    /* Set new coordinates for point by it's index
     */
    bool MovePoint(size_t pointIndex, const Vector& coords)
    {
        return MovePointImpl<true, true>(pointIndex, coords);
    }

private:

    CombinationsCache m_combintaionsCache;
    size_t m_pointsCount;
    std::vector<Segment> m_segments;
};