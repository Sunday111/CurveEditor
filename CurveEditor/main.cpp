#include <QtWidgets/qapplication.h>
#include <QtWidgets/qmainwindow.h>

#include "CurveEditor.h"

int main(int argc, char** argv)
{
    QApplication a(argc, argv);
    QMainWindow mainWindow;
    mainWindow.setCentralWidget(
        new CurveEditor(0, 1000, 0, 1000, &mainWindow));
    mainWindow.show();
    return a.exec();
}