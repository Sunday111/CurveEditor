#include <cassert>
#include <QtGui/qevent.h>
#include <QtGui/qpainter.h>
#include "CurveEditor.h"
#include "BezierCurve.h"
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

void CurveEditor::paintEvent(QPaintEvent * event)
{
	assert(event != nullptr);

	QPainter painter(this);
	painter.setPen(QPen(QColor(0,0,0), 2,
		Qt::PenStyle::SolidLine,
		Qt::PenCapStyle::RoundCap,
		Qt::PenJoinStyle::RoundJoin));
	QPoint cursor = mapFromGlobal(QCursor::pos());

	auto r = contentsRect();
	auto start = m_d->min;

	const double delta = 0.001;
	double parameter = 0.0;

	TransformCache cache;
	QPoint prevPoint;
	UserPointToScreenPoint(start.coords[X], start.coords[Y], &prevPoint, &cache);

	while ((parameter += delta) < 1)
	{
		auto _p2 = m_d->curve(parameter);

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

	if (m_d->movePtIdx > 0 && m_d->movePtIdx < m_d->curve.GetPointsCount() - 1)
	{
		Impl::Vector* prev = nullptr;
		m_d->curve.GetPoint(m_d->movePtIdx - 1, nullptr, &prev);

		Impl::Vector* next = nullptr;
		m_d->curve.GetPoint(m_d->movePtIdx + 1, nullptr, &next);

		if (x > prev->coords[X] && x < next->coords[X])
		{
			p->coords[X] = x;
		}
	}

	p->coords[Y] = y;

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

				QMenu* insertPointMenu = contextMenu.addMenu(tr("Insert point"));

				QAction insertStrong(tr("Strong"), this);
				insertPointMenu->addAction(&insertStrong);

				TransformCache cache;
				int userCoords[CurveEditor::Impl::Dimensions];
				ScreenPointToUserPoint(event->pos(), userCoords, userCoords + 1, &cache);

				const Impl::Curve::Vector p({ double(userCoords[X]), double(userCoords[Y]) });

				QObject::connect(&insertStrong, &QAction::triggered, [this, &p]()
				{
					m_d->curve.InsertPoint(p, false);
				});

				QAction insertWeak(tr("Weak"), this);
				insertPointMenu->addAction(&insertWeak);

				QObject::connect(&insertWeak, &QAction::triggered, [this, &p]()
				{
					m_d->curve.InsertPoint(p, true);
				});

				QAction removePoint(tr("Remove"), this);

				bool isWeak;
				size_t pointIndex;
				size_t segmentIndex;
				if (m_d->curve.GetPointInfo(p, 3.0, &segmentIndex, &isWeak, &pointIndex))
				{
					QObject::connect(&removePoint, &QAction::triggered, [this, &p]()
					{
						assert(false);
					});

					// The cursor is over the point
					contextMenu.addAction(&removePoint);
				}

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

void CurveEditor::CreateTransformCache(TransformCache* cache)
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
		CreateTransformCache(cache);
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
		CreateTransformCache(cache);
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