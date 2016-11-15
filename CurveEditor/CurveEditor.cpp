#include <cassert>
#include <QtGui/qevent.h>
#include <QtGui/qpainter.h>
#include "CurveEditor.h"
#include "BezierCurve.h"
#include <sstream>
#include <Windows.h>
#include <QtWidgets/qmenu.h>

namespace
{

enum { X = 0, Y = 1 };

}

class CurveEditor::TransformCache
{
public:

	explicit TransformCache() :
		visited(false)
	{}

	bool visited;
	int userRange[2];
	int screenRange[2];
	QPoint screenMin;
	QPoint screenMax;
};

class CurveEditor::Impl
{
public:
	enum { DefaultCurveOrder = 3 };
	enum { Dimensions = 2 };
	using Number = double;
	using Curve = BezierCurve<double, Dimensions>;
	using Vector = Curve::Vector;

	explicit Impl(
		float minX, float maxX,
		float minY, float maxY) :
		movePtIdx(-1),
		normalBrush(QBrush(QColor(255, 255, 255))),
		highlightBrush(QBrush(QColor(0, 255, 0))),
		moveBrush(QBrush(QColor(255, 255, 0))),
		weakControlPoint(QBrush(QColor(255, 255, 0))),
		strongControlPoint(QBrush(QColor(255, 0, 0))),
		curve(
			Vector({ minX, minY }),
			Vector({ maxX, maxY }),
			DefaultCurveOrder),
		min(Vector({ minX, minY })),
		max(Vector({ maxX, maxY }))
	{}

	int movePtIdx;
	Curve curve;

	Vector min;
	Vector max;

	QBrush moveBrush;
	QBrush normalBrush;
	QBrush highlightBrush;

	QBrush weakControlPoint;
	QBrush strongControlPoint;
};

CurveEditor::CurveEditor(
	int minX, int maxX,
	int minY, int maxY,
	QWidget* parent) :
	QWidget(parent),
	m_d(std::make_unique<Impl>(
		minX, maxX, minY, maxY))
{
	setMouseTracking(true);
}

CurveEditor::~CurveEditor() = default;

/*bool CurveEditor::AddPoint(int x, int y)
{
	if (x < m_d->min.x || y < m_d->min.y ||
		x > m_d->max.x || y > m_d->max.y)
	{
		return false;
	}

	for (auto& p : m_d->points)
	{
		if (p.x == x)
		{
			return false;
		}
	}

	m_d->points.emplace_back(x, y);

	std::sort(m_d->points.begin(), m_d->points.end(),
		[](const Point& p1, const Point& p2)
	{
		return p1.x < p2.x;
	});

	return true;
}*/

void CurveEditor::paintEvent(QPaintEvent * event)
{
	assert(event != nullptr);

	QPainter painter(this);
	QPoint cursor = mapFromGlobal(QCursor::pos());

	auto r = contentsRect();
	auto start = m_d->min;
	auto end = m_d->max;

	int x1 = start.coords[X];
	int y1 = start.coords[Y];

	int x2 = start.coords[X];

	TransformCache cache;
	QPoint prevPoint;
	UserPointToScreenPoint(x1, y1, &prevPoint, &cache);

	while (++x2 < end.coords[X])
	{
		auto _p2 = m_d->curve(x2);

		std::stringstream stream;
		stream << "X : " << _p2.coords[X] << "; " <<
			"Y : " << _p2.coords[Y] << std::endl;
		OutputDebugStringA(stream.str().c_str());

		QPoint nextPoint;
		UserPointToScreenPoint(_p2.coords[X], _p2.coords[Y], &nextPoint, &cache);

		painter.drawLine(prevPoint, nextPoint);

		prevPoint = nextPoint;
	}

	auto drawPoint = [&painter](const QPoint& p)
	{
		painter.drawEllipse(p, 5, 5);
	};

	const size_t pointsCount = m_d->curve.GetPointsCount();

	auto covnertPoint = [this, &cache]
		(QPoint& res, const Impl::Curve::Vector& p)
	{
		UserPointToScreenPoint(p.coords[X], p.coords[Y], &res, &cache);
	};

	QPoint prevPt;
	covnertPoint(prevPt, m_d->curve.GetStartPoint());

	painter.setPen(QColor(127, 127, 127));
	for (size_t i = 0; i < pointsCount; ++i)
	{

		Impl::Curve::Vector* point;
		bool weak;
		m_d->curve.GetPoint(i, &weak, &point);

		painter.setBrush(weak ? m_d->weakControlPoint : m_d->strongControlPoint);

		QPoint p;
		covnertPoint(p, *point);
		
		painter.drawEllipse(p, 5, 5);

		if (i != 0)
		{
			painter.drawLine(prevPt, p);
			prevPt = p;
		}
	}
}

void CurveEditor::mouseMoveEvent(QMouseEvent * event)
{
	if (m_d->movePtIdx == -1)
	{
		return;
	}

	int x, y;
	TransformCache cache;
	ScreenPointToUserPoint(event->pos(), &x, &y, &cache);

	Impl::Vector* p = nullptr;
	m_d->curve.GetPoint(m_d->movePtIdx, nullptr, &p);

	if (m_d->movePtIdx == 0 || m_d->movePtIdx == m_d->curve.GetPointsCount() - 1)
	{
		p->coords[Y] = y;
	}
	else
	{
		Impl::Vector* prev = nullptr;
		m_d->curve.GetPoint(m_d->movePtIdx - 1, nullptr, &prev);

		Impl::Vector* next = nullptr;
		m_d->curve.GetPoint(m_d->movePtIdx + 1, nullptr, &next);

		if (x > prev->coords[X] && x < next->coords[X])
		{
			p->coords[X] = x;
			p->coords[Y] = y;
		}
	}

	update();
}

void CurveEditor::mousePressEvent(QMouseEvent * event)
{
	if (event->button() != Qt::MouseButton::LeftButton)
	{
		return;
	}

	assert(m_d->movePtIdx == -1);

	int x, y;
	TransformCache cache;
	ScreenPointToUserPoint(event->pos(), &x, &y, &cache);

	const double tol = 7.0;

	const size_t pointsCount = m_d->curve.GetPointsCount();
	for (size_t pointIndex = 0; pointIndex < pointsCount; ++pointIndex)
	{
		Impl::Vector* p;
		m_d->curve.GetPoint(pointIndex, nullptr, &p);

		if (x > p->coords[X] - tol && x < p->coords[X] + tol &&
			y > p->coords[Y] - tol && y < p->coords[Y] + tol)
		{
			m_d->movePtIdx = static_cast<int>(pointIndex);
			break;
		}
	}

	update();
}

void CurveEditor::mouseReleaseEvent(QMouseEvent * event)
{
	switch (event->button())
	{
		case Qt::MouseButton::RightButton:
			if (contentsRect().contains(event->pos()))
			{
				QMenu contextMenu(this);

				QAction insertStrong(tr("Insert strong point"), this);
				contextMenu.addAction(&insertStrong);

				TransformCache cache;
				int userCoords[CurveEditor::Impl::Dimensions];
				ScreenPointToUserPoint(event->pos(), userCoords, userCoords + 1, &cache);

				const Impl::Curve::Vector p({ double(userCoords[X]), double(userCoords[Y]) });

				QObject::connect(&insertStrong, &QAction::triggered, [this, &p]()
				{
					m_d->curve.InsertPoint(p, false);
				});

				QAction insertWeak(tr("Insert weak point"), this);
				contextMenu.addAction(&insertWeak);

				QObject::connect(&insertWeak, &QAction::triggered, [this, &p]()
				{
					m_d->curve.InsertPoint(p, true);
				});

				contextMenu.exec(mapToGlobal(event->pos()));

				update();
			}
			break;

		case Qt::MouseButton::LeftButton:
			m_d->movePtIdx = -1;
			update();
			break;

	}
}

void CurveEditor::mouseDoubleClickEvent(QMouseEvent * event)
{
	int coords[2];
	TransformCache cache;
	ScreenPointToUserPoint(event->pos(), coords, coords + 1, &cache);

	//AddPoint(coords[0], coords[1]);
	update();
}

void CurveEditor::CreateTransforCache(TransformCache* cache)
{
	assert(cache != nullptr);

	if (cache->visited)
	{
		return;
	}

	auto r = contentsRect();
	cache->screenMin = r.bottomLeft();
	cache->screenMax = r.topRight();

	cache->userRange[X] = m_d->max.coords[X] - m_d->min.coords[X];
	cache->userRange[Y] = m_d->max.coords[Y] - m_d->min.coords[Y];
	
	cache->screenRange[X] = cache->screenMax.x() - cache->screenMin.x();
	cache->screenRange[Y] = cache->screenMin.y() - cache->screenMax.y();
}

void CurveEditor::ScreenPointToUserPoint(const QPoint& p, int* x, int* y, TransformCache* cache)
{
	assert(x != nullptr && y != nullptr && cache != nullptr);

	if (!cache->visited)
	{
		CreateTransforCache(cache);
	}

	const float percents[]
	{
		float(p.x() - cache->screenMin.x()) / cache->screenRange[X],
		float(cache->screenMin.y() - p.y()) / cache->screenRange[Y]
	};

	*x = cache->userRange[X] * percents[X];
	*y = cache->userRange[Y] * percents[Y];
}

void CurveEditor::UserPointToScreenPoint(int x, int y, QPoint* p, TransformCache* cache)
{
	assert(cache != nullptr && p != nullptr);

	if (!cache->visited)
	{
		CreateTransforCache(cache);
	}

	const float percents[]
	{
		float(x - m_d->min.coords[X]) / cache->userRange[X],
		float(y - m_d->min.coords[Y]) / cache->userRange[Y]
	};

	*p = QPoint(
		cache->screenMin.x() + cache->screenRange[X] * percents[0],
		cache->screenMin.y() - cache->screenRange[Y] * percents[1]);
}

#include "moc/CurveEditor.h"