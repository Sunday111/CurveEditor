#include <QtGui/qevent.h>
#include <QtGui/qpainter.h>
#include "TunableParametricFunction.h"

#include "FunctionEditor.h"

class FunctionEditor::Impl
{
public:
    int pickedFunction;
    int pickedPoint;
};

FunctionEditor::FunctionEditor(size_t approximationLevel, QWidget* parent) :
    FunctionViewer(approximationLevel, parent),
    m_d(std::make_unique<Impl>())
{
    m_d->pickedPoint = -1;
    m_d->pickedFunction = -1;
}

FunctionEditor::~FunctionEditor() = default;

void FunctionEditor::paintEvent(QPaintEvent* event)
{
    FunctionViewer::paintEvent(event);

    const QRect rect = contentsRect();
    QPainter painter(this);
    painter.setPen(QPen(QColor(0, 0, 0), 2,
        Qt::PenStyle::SolidLine,
        Qt::PenCapStyle::RoundCap,
        Qt::PenJoinStyle::RoundJoin));

    const size_t fns = GetFunctionsCount();
    for(size_t i = 0; i < fns; ++i)
    {
        if(auto casted = dynamic_cast<TunableParametricFunction<T, 2>*>(GetFunction(i)))
        {
            const size_t pts = casted->GetControlPointsCount();

            for(int j = 0; j < pts; ++j)
            {
                const Vector* pPt = nullptr;

                if(casted->GetControlPoint(j, &pPt))
                {
                    QPoint point;
                    UserPosToScreenPos(rect, *pPt, point);
                    painter.drawEllipse(point, 3, 3);
                }
            }
        }
    }
}

void FunctionEditor::mouseMoveEvent(QMouseEvent* e)
{
    if(m_d->pickedFunction == -1)
    {
        return;
    }

    const QRect rect = contentsRect();
    const QPointF pos = e->localPos();

    if(auto casted = dynamic_cast<TunableParametricFunction<T, 2>*>(GetFunction(m_d->pickedFunction)))
    {
        Vector newUserPos;
        ScreenPosToUserPos(rect, pos, newUserPos);

        casted->SetControlPoint(m_d->pickedPoint, newUserPos);
    }

    update();
}

void FunctionEditor::mousePressEvent(QMouseEvent* e)
{
    const QRect rect = contentsRect();
    const auto pos = e->localPos();

    const size_t fns = GetFunctionsCount();
    for(size_t i = 0; i < fns; ++i)
    {
        if(auto casted = dynamic_cast<TunableParametricFunction<T, 2>*>(GetFunction(i)))
        {
            const size_t pts = casted->GetControlPointsCount();

            for(int j = 0; j < pts; ++j)
            {
                const Vector* pPt = nullptr;

                if(casted->GetControlPoint(j, &pPt))
                {
                    QPoint point;
                    UserPosToScreenPos(rect, *pPt, point);

                    constexpr const size_t pickTolPx = 10;

                    if(std::abs(point.x() - pos.x()) < pickTolPx &&
                        std::abs(point.x() - pos.x()) < pickTolPx)
                    {
                        m_d->pickedFunction = i;
                        m_d->pickedPoint = j;

                        return;
                    }
                }
            }
        }
    }
}

void FunctionEditor::mouseReleaseEvent(QMouseEvent*)
{
    m_d->pickedPoint = -1;
    m_d->pickedFunction = -1;
}

#include "moc/FunctionEditor.h"