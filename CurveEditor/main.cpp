#include "BezierCurve.h"
#include <QtWidgets/qapplication.h>
#include "CMainWindow.h"

int main(int argc, char** argv)
{
	//struct
	//{
	//	double coords[2];
	//}
	//points[4]
	//{
	//	{ 0.0, 0.0 },
	//	{ 1.0 / 3.0, 1.0 / 3.0 },
	//	{ 2.0 / 3.0, 2.0 / 3.0 },
	//	{ 1.0, 1.0 }
	//};
	//
	//double t = 0.01;
	//std::vector<std::pair<double, double>> results;
	//
	//do
	//{
	//	double res[]{ 0.0, 0.0 };
	//
	//	const double k1[]
	//	{
	//		std::pow(1 - t, 3.0),
	//		3 * std::pow(1 - t, 2.0) * t,
	//		3 * (1 - t) * t * t,
	//		t * t * t
	//	};
	//
	//	const double k[]
	//	{
	//		BernsteinPolynomial<double>(3, 0, t),
	//		BernsteinPolynomial<double>(3, 1, t),
	//		BernsteinPolynomial<double>(3, 2, t),
	//		BernsteinPolynomial<double>(3, 3, t)
	//	};
	//
	//	for (int i = 0; i < 4; ++i)
	//	{
	//		for (int j = 0; j < 2; ++j)
	//		{
	//			res[j] += k[i] * points[i].coords[j];
	//		}
	//	}
	//
	//	results.emplace_back(res[0], res[1]);
	//}
	//while ((t += 0.01) < 1.0);


    QApplication a(argc, argv);
    CMainWindow mainWindow;
    mainWindow.show();
    return a.exec();
}