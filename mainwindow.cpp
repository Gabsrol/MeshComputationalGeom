#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent),
                                          ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_checkBox_Vertices_stateChanged(int arg)
{
    // When the checkbox is checked, this function is appealed with arg = 2
    // When it is uncheked, this function is appealed with arg = 0
    // The boolean show_vertecis is updated considering arg.
    ui->rendering->displayVertices = arg == 2;
}

void MainWindow::on_checkBox_Edges_stateChanged(int arg)
{
    // When the checkbox is checked, this function is appealed with arg = 2
    // When it is uncheked, this function is appealed with arg = 0
    // The boolean displayEdges is updated considering arg.
    ui->rendering->displayEdges = arg == 2;
}

void MainWindow::on_checkBox_Faces_stateChanged(int arg)
{
    // When the checkbox is checked, this function is appealed with arg = 2
    // When it is uncheked, this function is appealed with arg = 0
    // The boolean displayFaces is updated considering arg.
    ui->rendering->displayFaces = arg == 2;
}

void MainWindow::on_checkBox_Crust_stateChanged(int arg)
{
    // When the checkbox is checked, this function is appealed with arg = 2
    // When it is uncheked, this function is appealed with arg = 0
    // The boolean displayCrust is updated considering arg.
    ui->rendering->displayCrust = arg == 2;
    ui->rendering->displayEdges = false;
}


void MainWindow::on_pushButton_released()
{
    ui->rendering->add_random_vertex();
}

void MainWindow::on_voronoiButton_released()
{
    ui->rendering->add_voronoi_centers();
}
