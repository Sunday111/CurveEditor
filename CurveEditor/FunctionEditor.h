#pragma once

#include <memory>
#include "FunctionViewer.h"

class FunctionEditor :
    public FunctionViewer
{
    Q_OBJECT
public:
    FunctionEditor(
        size_t approximationLevel = 100,
        QWidget* parent = nullptr);
    ~FunctionEditor();

    void paintEvent(QPaintEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;

private:
    class Impl;
    std::unique_ptr<Impl> m_d;
};