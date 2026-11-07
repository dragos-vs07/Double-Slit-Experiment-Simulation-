#ifndef WAVEVIEW_H
#define WAVEVIEW_H

#include <QWidget>
#include <QImage>
#include <QPainter>
#include "simulation.h"
#include <algorithm>

class WaveView : public QWidget
{
    Q_OBJECT

public:

    explicit WaveView(QWidget *parent = nullptr);

    void showFrame(const Simulation &sim)
    {
        int w = sim.nx();
        int h = sim.ny();
        const auto &d = sim.density(); // gets the density vector containing all the |psi|^2 calculated values

        if (img_.size() != QSize(w, h))
            img_ = QImage(w, h, QImage::Format_RGB32);

        double maxd = *std::max_element(d.begin(),d.end());

        for(int i=0;i<h;i++)
            for(int j=0;j<w;j++)
            {
                double v = d[i * w + j];
                double t = std::sqrt(v / maxd);      // 0..1; sqrt brightens the faint parts
                int g = int(255 * t);
                img_.setPixel(j, h - 1 - i, qRgb(g, g, g));
            }

        update();
    }

protected:

    void paintEvent(QPaintEvent *event) override;
private:

    QImage img_;

};

#endif // WAVEVIEW_H
