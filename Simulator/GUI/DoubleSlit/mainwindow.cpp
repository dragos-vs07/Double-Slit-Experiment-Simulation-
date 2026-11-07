#include "mainwindow.h"
#include "simulation.h"
#include <QDebug>
#include <QElapsedTimer>
#include <QVBoxLayout>

void MainWindow::tick()
{
    sim_.step(3);
    view_->showFrame(sim_);
}
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    view_ = new WaveView(this);

    playBtn_ = new QPushButton("Play", this);

    auto *central = new QWidget(this);
    auto *layout = new QVBoxLayout(central);
    layout->addWidget(view_, 1);      // the 1 means "this one takes the extra space"
    layout->addWidget(playBtn_);

    setCentralWidget(central);
    sim_.reset(Params{});
    view_->showFrame(sim_);

    timer_ = new QTimer(this);
    timer_->setInterval(16);   // ms between frames
    connect(timer_, &QTimer::timeout, this, &MainWindow::tick);

    connect(playBtn_, &QPushButton::clicked, this, [this]() {
        if (timer_->isActive()) { timer_->stop();  playBtn_->setText("Play"); }
        else                    { timer_->start(); playBtn_->setText("Pause"); }
    });
}

MainWindow::~MainWindow()
{
}
