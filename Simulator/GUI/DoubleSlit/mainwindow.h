#ifndef MAINWINDOW_H
#define MAINWINDOW_H

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

private:
    Simulation sim_;
    WaveView *view_;

    QTimer *timer_;
    QPushButton *playBtn_;
    void tick();
};
#endif // MAINWINDOW_H
