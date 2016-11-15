#pragma once

#include <memory>

#include <QtWidgets/qwidget.h>

class CurveEditor :
	public QWidget
{
	Q_OBJECT
public:
	CurveEditor(
		int minX, int maxX,
		int minY, int maxY,
		QWidget* parent = nullptr);
	~CurveEditor();

	//bool AddPoint(int x, int y);
	void paintEvent(QPaintEvent *) override;
	void mouseMoveEvent(QMouseEvent *) override;
	void mousePressEvent(QMouseEvent *) override;
	void mouseReleaseEvent(QMouseEvent *) override;
	void mouseDoubleClickEvent(QMouseEvent *) override;

private:
	class TransformCache;
	void CreateTransforCache(TransformCache* cache);
	void ScreenPointToUserPoint(const QPoint& p, int* x, int* y, TransformCache* cache);
	void UserPointToScreenPoint(int x, int y, QPoint* p, TransformCache* cache);

	class Impl;
	std::unique_ptr<Impl> m_d;
};