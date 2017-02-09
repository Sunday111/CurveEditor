#include <QtWidgets/qapplication.h>
#include <QtWidgets/qmainwindow.h>

#include "CurveEditor.h"
#include "FunctionEditor.h"

#include "Ellipse.h"
#include "Line.h"
#include "Hyperbola.h"
#include "Parabola.h"

#include <cassert>

int main(int argc, char** argv)
{
    QApplication a(argc, argv);
    QMainWindow mainWindow;

    using T = double;
    constexpr const size_t dimensions = 2;

    auto functionViewer = new FunctionEditor(100, &mainWindow);

    using EllipseFunction = Ellipse<T, dimensions>;
    functionViewer->AddFunction(
        std::make_unique<EllipseFunction>(
            0.20, 0.10,
            EllipseFunction::Vector({ 0.75, 0.25 }),
            EllipseFunction::Vector({ 1.0 , 0.0 }),
            EllipseFunction::Vector({ 0.0 , 1.0 })));

    //using HyperbolaFunction = Hyperbola<T, dimensions>;
    //functionViewer->AddFunction(
    //    std::make_unique<HyperbolaFunction>(
    //        0.10, 0.10,
    //        HyperbolaFunction::Vector({ 0.5, 0.5 }),
    //        HyperbolaFunction::Vector({ 1.0 , 0.0 }),
    //        HyperbolaFunction::Vector({ 0.0 , 1.0 })));

    using LineFunction = Line<T, dimensions>;
    functionViewer->AddFunction(
        std::make_unique<LineFunction>(
            LineFunction::Vector({ 0.05, 0.05 }),
            LineFunction::Vector({ 0.45, 0.45 })));

    using ParabolaFunction = Parabola<T, dimensions>;
    functionViewer->AddFunction(
        std::make_unique<ParabolaFunction>(
            0.10,
            ParabolaFunction::Vector({ 0.5, 0.5 }),
            ParabolaFunction::Vector({ 1.0 , 0.0 }),
            ParabolaFunction::Vector({ 0.0 , 1.0 })));

    mainWindow.setCentralWidget(functionViewer);
    mainWindow.setFixedSize(1000, 1000);
    mainWindow.show();
    return a.exec();

    //QApplication a(argc, argv);
    //QMainWindow mainWindow;
    //mainWindow.setCentralWidget(
    //    new CurveEditor(0, 1000, 0, 1000, &mainWindow));
    //mainWindow.show();
    //return a.exec();
}