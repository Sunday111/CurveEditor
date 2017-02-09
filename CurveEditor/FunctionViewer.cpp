#include <cassert>
#include <QtGui/qevent.h>
#include <QtGui/qpainter.h>

#include "FunctionViewer.h"

class FunctionViewer::Impl
{
public:
    std::vector<std::unique_ptr<Function>> fns;
    size_t approximationLevel;
};

FunctionViewer::FunctionViewer(
    size_t approximationLevel,
    QWidget* parent) :
    QWidget(parent),
    m_d(std::make_unique<Impl>())
{
    m_d->approximationLevel = approximationLevel;
}
FunctionViewer::~FunctionViewer() = default;

void FunctionViewer::paintEvent(QPaintEvent* e)
{
    assert(e != nullptr);

    if(m_d->fns.empty())
    {
        return;
    }

    QPainter painter(this);
    painter.setPen(QPen(QColor(0, 0, 0), 2,
        Qt::PenStyle::SolidLine,
        Qt::PenCapStyle::RoundCap,
        Qt::PenJoinStyle::RoundJoin));

    const QRect rect = contentsRect();

    const auto delta = 1.0 / m_d->approximationLevel;

    for(auto& fn : m_d->fns)
    {
        auto t = 0.0;

        QPoint prevPoint;
        UserPosToScreenPos(rect, fn->Evaluate(t), prevPoint);

        for(size_t i = 1; i <= m_d->approximationLevel; ++i)
        {
            t += delta;

            QPoint currPoint;
            UserPosToScreenPos(rect, fn->Evaluate(t), currPoint);

            painter.drawLine(prevPoint, currPoint);
            prevPoint = currPoint;
        }
    }
}

void FunctionViewer::SetApproximationLevel(size_t val)
{
    assert(m_d != nullptr);
    m_d->approximationLevel = val;
}

bool FunctionViewer::AddFunction(std::unique_ptr<Function> fn)
{
    assert(m_d != nullptr);

    if(fn == nullptr)
    {
        return false;
    }
    
    if(std::find(m_d->fns.begin(), m_d->fns.end(), fn) !=
        m_d->fns.end())
    {
        return false;
    }

    m_d->fns.push_back(std::move(fn));

    return true;
}

size_t FunctionViewer::GetFunctionsCount() const
{
    assert(m_d != nullptr);
    return m_d->fns.size();
}

FunctionViewer::Function* FunctionViewer::GetFunction(size_t index) const
{
    assert(m_d != nullptr);
    return m_d->fns[index].get();
}

void FunctionViewer::UserPosToScreenPos(
    const QRect& screenRect,
    const Vector& inUserPos,
    QPoint& outScreenPos)
{
    outScreenPos.setX(
        screenRect.left() +
        inUserPos.coords[0] * screenRect.width()
    );

    outScreenPos.setY(
        screenRect.bottom() -
        inUserPos.coords[1] * screenRect.height()
    );
}

void FunctionViewer::ScreenPosToUserPos(
    const QRect& screenRect,
    const QPointF& inScreenPos,
    Vector& outUserPos)
{
    outUserPos.coords[0] =
        static_cast<T>(inScreenPos.x() - screenRect.left()) /
        static_cast<T>(screenRect.width());

    outUserPos.coords[1] =
        static_cast<T>(screenRect.bottom() - inScreenPos.y()) /
        static_cast<T>(screenRect.width());
}

#include "moc/FunctionViewer.h"