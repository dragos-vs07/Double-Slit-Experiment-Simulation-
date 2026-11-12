#ifndef WAVEVIEW_H
#define WAVEVIEW_H

#pragma once

#include <QWidget>
#include <QImage>
#include <QPainter>
#include <algorithm>
#include "simulation.h"

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
                int brightness_offest = 40;
                double v = d[i * w + j];
                double t = std::min(1.0, std::sqrt(v / maxd) * 1.5);
                int g = int(255 * t);
                img_.setPixel(j, h - 1 - i, qRgb(0, g / 2, g));
            }

        update();
    }

protected:

    void paintEvent(QPaintEvent *event) override;
private:

    QImage img_;

};

#endif // WAVEVIEW_H
