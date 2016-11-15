#include "CMainWindow.h"
#include <QtWidgets/qfiledialog.h>

#include "CurveEditor.h"

CMainWindow::CMainWindow()
    : QMainWindow(nullptr)
{
	this->setCentralWidget(new CurveEditor(0, 1000, 0, 1000, this));
}

#include "moc/CMainWindow.h"