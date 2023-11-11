#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QtCharts>
#include <QChartView>
#include <QLineSeries>

#include <complex.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>

#define M_PI 3.14159265358979323846
#define N 128
#define F0 1

typedef std::vector<std::complex<double>> compexVector ;

void FFT(compexVector &analysValue)
{
    int n = analysValue.size();
    if (n == 1) return;

    compexVector W(n);
    double alpha;

    for (int i = 0; i < n; i ++)
    {
        alpha = 2 * M_PI * i / n;
        W[i] = std::complex<double>(cos(alpha), sin(alpha));
    }

    compexVector A(n / 2), B(n / 2);
    for (int i = 0; i < n / 2; i++)
    {
        A[i] = analysValue[i * 2];
        B[i] = analysValue[i * 2 + 1];
    }
    FFT(A);
    FFT(B);
    for ( int i = 0; i < n; i++)
    {
        analysValue[i] = A[i % (n / 2)] + W[i] * B[i % (n / 2)];
    }
}


//void createChart(QLineSeries &series,compexVector &value)
//{

//    for (int i = 0; i < value.size(); i++)
//    {
//        series.append(std::sqrt(value[i].real() * value[i].real() + value[i].imag() * value[i].imag()),
//                       std::atan2(value[i].imag(), value[i].real()));
//    }

//}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    compexVector harmonic;
    double tempVar;

    for (int i = 0; i < N; i++)
    {
        tempVar = sin(F0 * i);
//        std::cout << tempVar << " ";
        harmonic.push_back(tempVar);
    }

//    std::cout << "\n\n";

    FFT(harmonic);
//    createChart(&series, harmonic);

    FFTChart = new QChart;

    ui->m_fftChart->setChart(FFTChart);
    FFTChart->setTitle("FFT Chart");

    QValueAxis *axisX = new QValueAxis;
    axisX->setRange(-10, 10);
    axisX->setTickCount(1);
    axisX->setLabelFormat("hui");

    QValueAxis *axisY = new QValueAxis;
    axisY->setRange(0, 60);
    axisY->setTickCount(1);
    axisY->setLabelFormat("Hui");

    QLineSeries *series = new QLineSeries();

    for(int i = 0; i < harmonic.size(); i++)
    {
        series->append(std::atan2(harmonic[i].imag(), harmonic[i].real()),
                       std::sqrt(harmonic[i].real() * harmonic[i].real() + harmonic[i].imag() * harmonic[i].imag()));
    }

    FFTChart->addSeries(series);
    FFTChart->setAxisX(axisX, series);
    FFTChart->setAxisY(axisY, series);

//    for (int i = 0; i < N; i++)
//    {
//        std::cout <<" " << harmonic[i] << " ";
//    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

