#pragma once

#include <memory>

#include <QtWidgets/qwidget.h>
#include "Vector.h"
#include "ParametricFunction.h"

class FunctionViewer :
    public QWidget
{
    Q_OBJECT
public:
    using T = double;
    using Vector = Vector<T, 2>;
    using Function = ParametricFunction<T, 2>;

    FunctionViewer(
        size_t approximationLevel = 100,
        QWidget* parent = nullptr);
    ~FunctionViewer();

    void paintEvent(QPaintEvent *) override;
    void SetApproximationLevel(size_t val);
    bool AddFunction(std::unique_ptr<Function> fn);

    size_t GetFunctionsCount() const;
    Function* GetFunction(size_t index) const;

protected:
    static void UserPosToScreenPos(
        const QRect& screenRect,
        const Vector& inUserPos,
        QPoint& outScreenPos);

    static void ScreenPosToUserPos(
        const QRect& screenRect,
        const QPointF& inScreenPos,
        Vector& outUserPos);

private:
    class Impl;
    std::unique_ptr<Impl> m_d;
};