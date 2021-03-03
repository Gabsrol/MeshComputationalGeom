#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui{class MainWindow;}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_checkBox_Vertices_stateChanged(int arg);
    void on_checkBox_Edges_stateChanged(int arg);
    void on_checkBox_Faces_stateChanged(int arg);
    void on_checkBox_Crust_stateChanged(int arg);
    void on_pushButton_released();
    void on_voronoiButton_released();

private:
    Ui::MainWindow *ui;
};

#endif
