#pragma once

#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include <sstream>

inline size_t Factorial(size_t val)
{
	if (val < 2) return 1;

	return val * (Factorial(val - 1));
}

/* The number of combinations from n by i */
inline size_t Combinations(size_t n, size_t i)
{
	// n! / (i!(n-i)!)
	return Factorial(n) / (Factorial(i) * Factorial(n - i));
}

template<class T>
inline T BernsteinPolynomial(size_t n, size_t i, const T& t)
{
	return
		Combinations(n, i) *
		std::pow(t, i) *
		std::pow(1 - t, n - i);
}

template<class T, size_t dimensions>
class Vector
{
public:
	explicit Vector(const T(&_coords)[dimensions])
	{
		for (int i = 0; i < dimensions; ++i)
		{
			coords[i] = _coords[i];
		}
	}

	explicit Vector()
	{
	}

	Vector operator* (const T& val) const
	{
		Vector result = *this;

		for (int i = 0; i < dimensions; ++i)
		{
			result.coords[i] *= val;
		}

		return result;
	}

	Vector& operator*= (const T& val)
	{
		for (int i = 0; i < dimensions; ++i)
		{
			coords[i] *= val;
		}

		return *this;
	}

	Vector operator+ (const Vector& p) const
	{
		Vector result = *this;

		for (int i = 0; i < dimensions; ++i)
		{
			result.coords[i] += p.coords[i];
		}

		return result;
	}

	Vector& operator+= (const Vector& p)
	{
		for (int i = 0; i < dimensions; ++i)
		{
			coords[i] += p.coords[i];
		}

		return *this;
	}

	Vector operator- (const Vector& p) const
	{
		Vector result(*this);

		for (int i = 0; i < dimensions; ++i)
		{
			result.coords[i] -= p.coords[i];
		}

		return result;
	}

	Vector operator/ (const T& val) const
	{
		Vector result(*this);

		for (int i = 0; i < dimensions; ++i)
		{
			result.coords[i] /= val;
		}

		return result;
	}

	T DistanceTo(const Vector& p) const
	{
		T result = 0;

		for (int i = 0; i < dimensions; ++i)
		{
			auto dt = coords[i] - p.coords[i];
			result += dt * dt;
		}

		return std::sqrt(result);
	}

	T coords[dimensions];
};

template<class T, size_t dimensions>
class BezierPoints
{
public:
	using Vector = Vector<T, dimensions>;

	explicit BezierPoints(Vector* sp, const Vector& ep, size_t order) :
		startPoint(sp),
		m_endPoint(std::make_unique<Vector>(ep)),
		m_order(order)
	{
		assert(sp != nullptr);

		if (order < 2)
		{
			return;
		}

		const size_t innerPoints = order - 1;

		m_points.resize(innerPoints);
		auto dir = ep - *sp;

		for (size_t i = 0; i < innerPoints; ++i)
		{
			m_points[i] = (dir * (i + 1)) / (innerPoints + 1) + *sp;
			m_points[i].coords[1] += dir.coords[1] / 5;
		}
	}
	BezierPoints(const BezierPoints&) = delete;
	BezierPoints(BezierPoints&& uref) :
		m_order(uref.m_order),
		startPoint(uref.startPoint),
		m_points(std::move(uref.m_points)),
		m_endPoint(std::move(uref.m_endPoint))
	{
		uref.startPoint = nullptr;
	}

	const size_t PointsCount() const { return m_order + 1; }
	const size_t InnerPoints() const { return m_order - 1; }
	const size_t Order() const { return m_order; }

	T GetParameter(const T& coord) const
	{
		const T range[]
		{
			startPoint->coords[0],
			m_endPoint->coords[0]
		};

		if (coord > range[1])
		{
			return -1;
		}

		return (coord - range[0]) / (range[1] - range[0]);
	}

	Vector& FirstPoint() { return *startPoint; }
	const Vector& FirstPoint() const { return *startPoint; }
	Vector& LastPoint() { return *m_endPoint; }
	const Vector& LastPoint() const { return *m_endPoint; }

	Vector& operator[](size_t idx)
	{
		const size_t pointsCount = PointsCount();
		assert(idx < pointsCount);

		if (idx == 0)
		{
			return *startPoint;
		}

		if (idx == pointsCount - 1)
		{
			return *m_endPoint;
		}

		return m_points[idx - 1];
	}

	const Vector& operator[](size_t idx) const
	{
		return const_cast<BezierPoints*>(this)->operator[](idx);
	}

	BezierPoints& operator=(const BezierPoints&) = delete;

	BezierPoints& operator=(BezierPoints&& uref)
	{

		m_order = uref.m_order;
		startPoint = uref.startPoint;
		m_points = std::move(uref.m_points);
		m_endPoint = std::move(uref.m_endPoint);

		return *this;
	}

	bool InsertWeakPoint(const Vector& p)
	{
		for(size_t pointIndex = 0; pointIndex < m_points.size(); ++pointIndex)
		{
			const Vector& point = m_points[pointIndex];

			if (p.coords[0] < point.coords[0])
			{
				m_points.insert(m_points.begin() + pointIndex, p);
				++m_order;
				return true;
			}
		}

		if (p.coords[0] < m_endPoint->coords[0])
		{
			m_points.push_back(p);
			++m_order;
			return true;
		}

		return false;
	}

	int FindPoint(const Vector& p, const T& tol, bool* isWeak) const
	{
		assert(startPoint != nullptr && m_endPoint != nullptr);

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

		if (isWeak != nullptr)
		{
			*isWeak = resIdx > 0 && resIdx < m_order;
		}

		return static_cast<int>(resIdx);
	}

	BezierPoints SplitInPoint(const Vector& p)
	{
		BezierPoints retVal(startPoint, p, 1);
		startPoint = retVal.m_endPoint.get();

		size_t pointIndex = 0;

		for (; pointIndex < m_points.size(); ++pointIndex)
		{
			const Vector& point = m_points[pointIndex];

			if (point.coords[0] > p.coords[0])
			{
				break;
			}

			retVal.m_points.push_back(point);
		}

		retVal.m_order = 1 + retVal.m_points.size();

		m_points.erase(m_points.begin(), m_points.begin() + pointIndex);
		m_order = 1 + m_points.size();

		return retVal;
	}

protected:
	BezierPoints() :
		m_endPoint(std::make_unique<Vector>())
	{}

	size_t m_order;
	Vector* startPoint;
	std::vector<Vector> m_points;

	/* This point is allocated in heap becuse is can be referenced 
	 * by next segment so its address must not be changed
	 */
	std::unique_ptr<Vector> m_endPoint;
};

/* Where n is Bezier curve order */
template<class T, size_t dimensions>
class BezierSegment
{
public:
	using Points = BezierPoints<T, dimensions>;
	using Vector = typename Points::Vector;

	explicit BezierSegment(Vector* start, const Vector& end, size_t order) :
		m_points(start, end, order)
	{}
	BezierSegment(const BezierSegment&) = delete;
	BezierSegment(BezierSegment&& uref) :
		m_points(std::move(uref.m_points))
	{}

	Points& RPoints() { return m_points; }
	const Points& RPoints() const { return m_points; }
	const Points& GetPoints() const { return m_points; }

	BezierSegment SplitInPoint(const Vector& p)
	{
		return BezierSegment(m_points.SplitInPoint(p));
	}

	Vector operator() (const T& t) const
	{
		Vector result;
		result *= 0;

		const size_t pointsCount = m_points.PointsCount();
		const size_t order = m_points.Order();

		for (size_t i = 0; i < pointsCount; ++i)
		{
			result += m_points[i] * BernsteinPolynomial<T>(order, i, t);
		}

		return result;
	}

	BezierSegment& operator=(const BezierSegment&) = delete;
	BezierSegment& operator=(BezierSegment&& uref)
	{
		m_points = std::move(uref.m_points);
		return *this;
	}

private:
	BezierSegment(Points&& points) :
		m_points(std::move(points))
	{}

	Points m_points;
};

/* Where n is Bezier curve order */
template<class T, size_t dimensions>
class BezierCurve
{
public:
	using Segment = BezierSegment<T, dimensions>;
	using Vector = typename Segment::Vector;

	/* n is order of first bezier curve segment */
	explicit BezierCurve(const Vector& startPoint, const Vector& endPoint, size_t n) :
		m_startPoint(startPoint),
		m_pointsCount(n + 1)
	{
		//Segment initSegment(&m_startPoint, endPoint, n);
		//m_segments.push_back(initSegment);
		auto dir = endPoint - startPoint;

		Segment initSegment_0(&m_startPoint, startPoint + dir * 0.5, n);
		Segment initSegment_1(&initSegment_0.RPoints().LastPoint(), endPoint, n);
		m_segments.push_back(std::move(initSegment_0));
		m_segments.push_back(std::move(initSegment_1));

		m_pointsCount = 7;
	}

	Vector& RStartPoint() { return m_startPoint; }
	const Vector& RStartPoint() const { return m_startPoint; }
	void SetStartPoint(const Vector& val) { m_startPoint = val; }
	const Vector& GetStartPoint() const { return m_startPoint; }

	Vector& REndPoint() { return m_segments.back().RPoints().LastPoint(); }
	const Vector& REndPoint() const { return m_segments.back().RPoints().LastPoint(); }
	void SetEndPoint(const Vector& val) { m_endPoint = val; }
	const Vector& GetEndPoint() const { return m_endPoint; }

	const std::vector<Segment>& GetSegments() const { return m_segments; }
	size_t GetPointsCount() const { return m_pointsCount; }

	bool GetPoint(size_t index, bool* isWeak, const Vector** p) const
	{
		assert(index < m_pointsCount);
	
		size_t minIndex = 0;
	
		const size_t segmentsCount = m_segments.size();
		for (size_t segmentIdx = 0; segmentIdx < segmentsCount; ++segmentIdx)
		{
			const Segment& segment = m_segments[segmentIdx];
			auto& points = segment.RPoints();
	
			const size_t order = points.Order();
			const size_t maxIndex = minIndex + order;
	
			if (index >= minIndex && index <= maxIndex)
			{
				const size_t segmentIndex = index - minIndex;
	
				if (isWeak != nullptr)
				{
					*isWeak = !(segmentIndex == 0 || segmentIndex == order);
				}
	
				*p = &points[segmentIndex];
				return true;
			}
	
			minIndex = maxIndex;
		}
	
		return false;
	}

	bool GetPoint(size_t index, bool* isWeak, Vector** p)
	{
		const BezierCurve* c = this;
		return c->GetPoint(index, isWeak, const_cast<const Vector**>(p));
	}

	bool InsertPoint(const Vector& p, bool isWeak)
	{
		if (p.coords[0] < RStartPoint().coords[0] ||
			p.coords[0] > REndPoint().coords[0])
		{
			return false;
		}

		for (size_t segmentIndex = 0; segmentIndex < m_segments.size(); ++segmentIndex)
		{
			auto& segment = m_segments[segmentIndex];

			if (segment.RPoints().GetParameter(p.coords[0]) > 0)
			{
				if (isWeak)
				{
					segment.RPoints().InsertWeakPoint(p);
				}
				else
				{
					m_segments.insert(
						m_segments.begin() + segmentIndex,
						segment.SplitInPoint(p));
				}

				++m_pointsCount;
				return true;
			}
		}

		return false;
	}

	Vector operator()(const T& curveParameter) const
	{
		const T& sx = RStartPoint().coords[0];
		const T coord = sx + (REndPoint().coords[0] - sx) * curveParameter;

		for (auto& seg : m_segments)
		{
			const T segmentParameter = seg.GetPoints().GetParameter(coord);

			if (segmentParameter >= 0)
			{
				return seg(segmentParameter);
			}
		}

		assert(false);
		return Vector();
	}

	bool GetPointInfo(const Vector& p, const T& tol,
		size_t* segmentIndex = nullptr,
		bool* isWeak = nullptr,
		size_t* indexInSegment = nullptr) const
	{
		for(size_t segIdx = 0; segIdx < m_segments.size(); ++segIdx)
		{
			const Segment& segment = m_segments[segIdx];

			const int fnd = segment.RPoints().FindPoint(p, tol, isWeak);

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

				return true;
			}
		}

		return false;
	}

	bool RemovePoint(size_t index) const
	{
		return true;
	}

private:

	Vector m_startPoint;
	size_t m_pointsCount;
	std::vector<Segment> m_segments;
};