#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#pragma once

#include <QMainWindow>
#include <QTimer>
#include <QPushButton>
#include "waveview.h"


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    int frame_ = 0;
    int MaxFrames_ = 480;
    void resetSim();

private:
    Simulation sim_;
    WaveView *view_;

    QTimer *timer_;
    QPushButton *playBtn_;
    void tick();
    Params current_params_;
};
#endif // MAINWINDOW_H
