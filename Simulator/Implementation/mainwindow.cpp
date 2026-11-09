#pragma once
#include "mainwindow.h"
#include "simulation.h"
#include <QDebug>
#include <QElapsedTimer>
#include <QVBoxLayout>
#include <QDoubleSpinBox>
#include <QLabel>

void MainWindow::tick()
{
    sim_.step(3);
    view_->showFrame(sim_);

    if(++frame_ >= MaxFrames_)
        resetSim();
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    {

        auto *central = new QWidget(this);
        auto *main_layout = new QHBoxLayout(central);

        auto *animation_area = new QVBoxLayout();

        view_ = new WaveView(central);

        playBtn_ = new QPushButton("Play", central);

        animation_area->addWidget(view_, 1);      // the 1 means "this one takes the extra space"
        animation_area->addWidget(playBtn_);

        sim_.reset(Params{});
        view_->showFrame(sim_);

        timer_ = new QTimer(this);
        timer_->setInterval(16);   // ms between frames

        connect(timer_, &QTimer::timeout, this, &MainWindow::tick);

        connect(playBtn_, &QPushButton::clicked, this, [this]() {
            if (timer_->isActive()) { timer_->stop();  playBtn_->setText("Play"); }
            else                    { timer_->start(); playBtn_->setText("Pause"); }
        });

        auto *params_area = new QVBoxLayout();

        /*
            struct Params{
            double delta_x = 0.2, delta_y = 0.2;    // real distance between cells
            double delta_t = 0.005; // time per timestamp
            double kx = 3, ky = 0;  // gaussian packet momentum
            double x_zero = 6.0,y_zero = 20.0;  // origin of the packet ( position of electron gun )
            double sigma_x = 2, sigma_y = 4; // dimensions of the packet

            int slit_w = 5;        // slit width (cells)
            int slit_center = 100; // midpoint between the two slits (y, in cells)
            int slit_sep = 40;     // distance between the slit centers (your c1=80, c2=120)
            };
         */
        QDoubleSpinBox   *kx_input_box = new QDoubleSpinBox(), *ky_input_box = new QDoubleSpinBox(),
                    *delta_x_input_box = new QDoubleSpinBox(), *delta_y_input_box = new QDoubleSpinBox(),
                    *sigma_x_input_box = new QDoubleSpinBox(), *sigma_y_input_box = new QDoubleSpinBox();

        Params p_default;

        kx_input_box->setValue(p_default.kx);
        ky_input_box->setValue(p_default.ky);
        delta_x_input_box->setValue(p_default.delta_x);
        delta_y_input_box->setValue(p_default.delta_y);
        sigma_x_input_box->setValue(p_default.sigma_x);
        sigma_y_input_box->setValue(p_default.sigma_y);

        auto    *kx_label = new QLabel("kx"), *ky_label = new QLabel("ky"),
                *delta_x_label = new QLabel("delta x"), *delta_y_label = new QLabel("delta y"),
                *sigma_x_label = new QLabel("sigma x"), *sigma_y_label = new QLabel("sigma y");

        auto    *kx_input = new QVBoxLayout(), *ky_input = new QVBoxLayout(),
                *delta_x_input = new QVBoxLayout(), *delta_y_input = new QVBoxLayout(),
                *sigma_x_input = new QVBoxLayout(), *sigma_y_input = new QVBoxLayout();

        kx_input -> addWidget(kx_label);
        kx_input -> addWidget(kx_input_box);

        ky_input -> addWidget(ky_label);
        ky_input -> addWidget(ky_input_box);

        sigma_x_input -> addWidget(sigma_x_label);
        sigma_x_input -> addWidget(sigma_x_input_box);

        sigma_y_input -> addWidget(sigma_y_label);
        sigma_y_input -> addWidget(sigma_y_input_box);

        delta_x_input -> addWidget(delta_x_label);
        delta_x_input -> addWidget(delta_x_input_box);

        delta_y_input -> addWidget(delta_y_label);
        delta_y_input -> addWidget(delta_y_input_box);

        QHBoxLayout *k_inputs = new QHBoxLayout();
        QHBoxLayout *delta_inputs = new QHBoxLayout();
        QHBoxLayout *sigma_inputs = new QHBoxLayout();

        k_inputs -> addLayout(kx_input);
        k_inputs -> addLayout(ky_input);
        delta_inputs -> addLayout(delta_x_input);
        delta_inputs -> addLayout(delta_y_input);
        sigma_inputs -> addLayout(sigma_x_input);
        sigma_inputs -> addLayout(sigma_y_input);

        params_area -> addLayout(k_inputs);
        params_area -> addLayout(delta_inputs);
        params_area -> addLayout(sigma_inputs);

        auto set_params_btn = new QPushButton("Submit", central);
        connect(set_params_btn, &QPushButton::clicked, this, [this, kx_input_box, ky_input_box,
                                                              delta_x_input_box, delta_y_input_box,
                                                              sigma_x_input_box, sigma_y_input_box]() {

            timer_->stop();
            playBtn_->setText("Play");

            current_params_.kx = kx_input_box->value();
            current_params_.ky = ky_input_box->value();
            current_params_.delta_x = delta_x_input_box->value();
            current_params_.delta_y = delta_y_input_box->value();
            current_params_.sigma_x = sigma_x_input_box->value();
            current_params_.sigma_y = sigma_y_input_box->value();

            resetSim();
        });

        params_area -> addWidget(set_params_btn);

        main_layout -> addLayout(params_area, 1);
        main_layout -> addLayout(animation_area, 3);


        setCentralWidget(central);
    }

void MainWindow::resetSim()
{
    frame_ = 0;
    sim_.reset(current_params_);
    view_ -> showFrame(sim_);
}

MainWindow::~MainWindow()
{}
