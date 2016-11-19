#pragma once

#include <memory>

#include <QtWidgets/qwidget.h>

class CurveEditor :
    public QWidget
{
    Q_OBJECT
public:
    CurveEditor(
        double minX, double maxX,
        double minY, double maxY,
        QWidget* parent = nullptr);
    ~CurveEditor();

    void paintEvent(QPaintEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;

private:
    class TransformCache;
    void CreateTransformCache(TransformCache* cache);
    void ScreenPointToUserPoint(const QPoint& p, double (&coords)[2], TransformCache* cache);
    void UserPointToScreenPoint(double const (&coords)[2], QPoint* p, TransformCache* cache);

    class Impl;
    std::unique_ptr<Impl> m_d;
};