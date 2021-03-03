#ifndef GLDISPLAYWIDGET_H
#define GLDISPLAYWIDGET_H

#include <QGLWidget>
#include <QtWidgets>
#include <QTimer>
#include "mesh.h"

class GLDisplayWidget : public QGLWidget
{
public:
    explicit GLDisplayWidget(QWidget *parent = 0);

    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

    Mesh _mesh;

    bool displayVertices;
    bool displayEdges;
    bool displayFaces;
    bool displayCrust;

    void add_random_vertex();
    void add_voronoi_centers();

protected:
    // Mouse Management
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

    void drawVertices();  // Display vertices of _mesh when appealed
    void drawEdges(bool crust=false);     // Same for edges
    void drawFaces();     // Same for faces

private:
    // timer, camera position and mouse position
    QTimer _timer;
    float _X, _Y, _Z;
    float _angle;
    QPoint _lastPosMouse;
};

#endif
