#include <QtWidgets/qapplication.h>
#include <QtWidgets/qmainwindow.h>

#include "CurveEditor.h"
#include "FunctionEditor.h"

#include "Ellipse.h"
#include "Line.h"

#include <cassert>

int main(int argc, char** argv)
{
    QApplication a(argc, argv);
    QMainWindow mainWindow;

    using T = double;
    constexpr const size_t dimensions = 2;
    using LineFunction = Line<T, dimensions>;
    using EllipseFunction = Ellipse<T, dimensions>;

    auto functionViewer = new FunctionEditor(100, &mainWindow);

    functionViewer->AddFunction(
        std::make_unique<EllipseFunction>(
            0.20, 0.10,
            EllipseFunction::Vector({ 0.75, 0.25 }),
            EllipseFunction::Vector({ 1.0 , 0.0 }),
            EllipseFunction::Vector({ 0.0 , 1.0 })));

    functionViewer->AddFunction(
        std::make_unique<EllipseFunction>(
            0.20, 0.20,
            EllipseFunction::Vector({ 0.25, 0.75 }),
            EllipseFunction::Vector({ 1.0 , 0.0 }),
            EllipseFunction::Vector({ 0.0 , 1.0 })));

    functionViewer->AddFunction(
        std::make_unique<LineFunction>(
            LineFunction::Vector({ 0.05, 0.05 }),
            LineFunction::Vector({ 0.45, 0.45 })));

    mainWindow.setCentralWidget(functionViewer);
    mainWindow.setFixedSize(500, 500);
    mainWindow.show();
    return a.exec();

    //QApplication a(argc, argv);
    //QMainWindow mainWindow;
    //mainWindow.setCentralWidget(
    //    new CurveEditor(0, 1000, 0, 1000, &mainWindow));
    //mainWindow.show();
    //return a.exec();
}