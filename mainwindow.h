#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui
{
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_checkBox_Vertices_stateChanged(int arg);  // Called when the "show vertices" checkbox is checked
    void on_checkBox_Edges_stateChanged(int arg);     // Called when the "show edges" checkbox is checked
    void on_checkBox_Faces_stateChanged(int arg);     // Called when the "show faces" checkbox is checked

    void on_pushButton_released();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
