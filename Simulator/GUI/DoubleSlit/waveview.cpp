#include "waveview.h"

WaveView::WaveView(QWidget *parent) : QWidget(parent)
{
    setMinimumSize(340, 400);
}

void WaveView::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    p.fillRect(rect(), Qt::black);
    if (img_.isNull()) return;

    QSize s = img_.size().scaled(size(), Qt::KeepAspectRatio);
    QRect r(QPoint(0, 0), s);
    r.moveCenter(rect().center());
    p.drawImage(r, img_);
}
