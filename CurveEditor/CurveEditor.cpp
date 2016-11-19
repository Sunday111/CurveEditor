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
        strongNormalControlPoint(QBrush(QColor(255, 0, 0))),
        strongSmoothControlPoint(QBrush(QColor(0, 0, 255))),
        curve(
            Vector({ minX, minY }),
            Vector({ maxX, maxY }),
            DefaultCurveOrder),
        min(Vector({ minX, minY })),
        max(Vector({ maxX, maxY }))
    {}

    const QBrush& GetBrush(BezierCurvePointType pointType) const
    {
        switch (pointType)
        {
        case BezierCurvePointType::Strong_Normal: return strongNormalControlPoint;
        case BezierCurvePointType::Strong_Smooth: return strongSmoothControlPoint;
        case BezierCurvePointType::Weak:          return weakControlPoint;
        }

        assert(false);
        return normalBrush;
    }

    int movePtIdx;
    Curve curve;

    Vector min;
    Vector max;

    QBrush moveBrush;
    QBrush normalBrush;
    QBrush highlightBrush;

    QBrush weakControlPoint;
    QBrush strongNormalControlPoint;
    QBrush strongSmoothControlPoint;

    struct DrawPointsContext
    {
        CurveEditor* _this;
        TransformCache cache;
        QPainter* painter;
        QPoint prevPoint;
    };
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

    auto start = m_d->min;

    const double delta = 0.001;
    double parameter = 0.0;

    //Prepare drawing context
    Impl::DrawPointsContext drawContext;
    drawContext.painter = &painter;
    drawContext._this = this;

    UserPointToScreenPoint(start.coords[X], start.coords[Y], &drawContext.prevPoint, &drawContext.cache);

    /* Draw approximated curve */
    while ((parameter += delta) < 1)
    {
        auto _p2 = m_d->curve(parameter);

        QPoint nextPoint;
        UserPointToScreenPoint(_p2.coords[X], _p2.coords[Y], &nextPoint, &drawContext.cache);

        painter.drawLine(drawContext.prevPoint, nextPoint);

        drawContext.prevPoint = nextPoint;
    }

    /* Draw lines between control points */
    painter.setPen(QColor(127, 127, 127));
    m_d->curve.WalkThroughThePoints(
        [](size_t /*indexInSegment*/, size_t indexInCurve, const Impl::Curve::Segment* /*segment*/,
            BezierCurvePointType /*pointType*/, const Impl::Vector* point, void* userContext)
    {
        Impl::DrawPointsContext* context = reinterpret_cast<Impl::DrawPointsContext*>(userContext);
    
        if (indexInCurve == 0)
        {
            context->_this->UserPointToScreenPoint(point->coords[X], point->coords[Y], &context->prevPoint, &context->cache);
        }
        else
        {
            QPoint p;
            context->_this->UserPointToScreenPoint(point->coords[X], point->coords[Y], &p, &context->cache);
            context->painter->drawLine(context->prevPoint, p);
            context->prevPoint = p;
        }
    
        return false;
    }, &drawContext);

    /* Draw control points */
    m_d->curve.WalkThroughThePoints(
        [](size_t /*indexInSegment*/, size_t /*indexInCurve*/, const Impl::Curve::Segment* /*segment*/,
            BezierCurvePointType pointType, const Impl::Vector* point, void* userContext)
    {
        Impl::DrawPointsContext* context = reinterpret_cast<Impl::DrawPointsContext*>(userContext);

        context->_this->UserPointToScreenPoint(point->coords[X], point->coords[Y], &context->prevPoint, &context->cache);
        context->painter->setBrush(context->_this->m_d->GetBrush(pointType));
        context->painter->drawEllipse(context->prevPoint, 5, 5);

        return false;
    }, &drawContext);
    
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
                QMenu* insertStrongPointMenu = insertPointMenu->addMenu(tr("Strong"));

                QAction insertStrongNormalPoint(tr("Normal"), this);
                QAction insertStrongSmoothPoint(tr("Smooth"), this);
                QAction insertWeak(tr("Weak"), this);
                QAction removePoint(tr("Remove"), this);

                insertStrongPointMenu->addAction(&insertStrongNormalPoint);
                insertStrongPointMenu->addAction(&insertStrongSmoothPoint);
                insertPointMenu->addAction(&insertWeak);

                TransformCache cache;
                int userCoords[CurveEditor::Impl::Dimensions];
                ScreenPointToUserPoint(event->pos(), userCoords, userCoords + 1, &cache);
                const Impl::Curve::Vector p({ double(userCoords[X]), double(userCoords[Y]) });

                QObject::connect(&insertStrongNormalPoint, &QAction::triggered, [this, &p]()
                {
                    m_d->curve.InsertPoint<BezierCurvePointType::Strong_Normal>(p);
                });

                QObject::connect(&insertStrongSmoothPoint, &QAction::triggered, [this, &p]()
                {
                    m_d->curve.InsertPoint<BezierCurvePointType::Strong_Smooth>(p);
                });

                QObject::connect(&insertWeak, &QAction::triggered, [this, &p]()
                {
                    m_d->curve.InsertPoint<BezierCurvePointType::Weak>(p);
                });

                size_t indexInCurve;
                BezierCurvePointType pointType;
                /* The cursor is over the point and it is not first or last strong point */
                if (m_d->curve.GetPointInfo(p, 3.0, nullptr, &pointType, nullptr, &indexInCurve) &&
                    indexInCurve > 0 || indexInCurve < m_d->curve.GetPointsCount() - 1)
                {
                    QMenu* editPointMenu = contextMenu.addMenu(tr("Edit point"));
                    editPointMenu->addAction(&removePoint);

                    QObject::connect(&removePoint, &QAction::triggered, [this, indexInCurve]()
                    {
                        m_d->curve.RemovePoint(indexInCurve);
                    });
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