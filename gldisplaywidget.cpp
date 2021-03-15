#include "gldisplaywidget.h"
#ifdef __APPLE__
#include <glu.h>
#else
#include <GL/glu.h>
#endif
#include <iostream>
#include "QDebug"


void glVertexDraw(const Vertex &p)
{
    glVertex3f(p.x(), p.y(), p.z());
}

GLDisplayWidget::GLDisplayWidget(QWidget *parent) : QGLWidget(parent)
{
    // Initial position of camera
    _X = 0;
    _Y = 0;
    _Z = 0;
    _angle = 0;

    // To display vertices, edges or faces
    displayVertices = false;
    displayEdges = false;
    displayFaces = false;

    // Update the scene
    connect(&_timer, SIGNAL(timeout()), this, SLOT(updateGL()));
    _timer.start(16);
}

void GLDisplayWidget::initializeGL()
{
    // Background color
    glClearColor(0.2, 0.2, 0.2, 1);

    // Shader
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);


    // Path to off file : put your path here

    //char path_folder[256] = "C:\\Users\\briss\\OneDrive\\Bureau\\mesh_computation\\MeshComputationalGeom\\off_files\\";
    char path_folder[256] = "/Users/gabin/Ordinateur/Documents/Centrale_Lyon/3A/Secteur/Calcul_Geometrique/Mesh_Computationnal_Geometry/off_files/";
    char off_file[32] = "lapin.off";
    char * path_off_file;
    path_off_file = strcat(path_folder, off_file);


    // building mesh : comment/uncomment the options

    // ---------------------------------------------------------
    // option 1 : 3D structures like queen.off :
    // ---------------------------------------------------------

    _mesh.parseFile(path_off_file);
    _mesh.sew();


    // ---------------------------------------------------------
    // option 2 : triangulation (if your file contains only vertices coordinates)
    // ---------------------------------------------------------

    //_mesh.parseTriFile(path_off_file);
    //_mesh.triangulationFromVertices();

}

void GLDisplayWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Center the camera
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);

    // Translation
    glTranslated(_X, _Y, _Z);

    // Rotation
    glRotatef(_angle, 1.0f, 1.0f, 0.0f);

    // Color for your mesh
    glColor3f(0, 1, 0);

    // drawings with respect to options chosen in the interactive window
    if (displayVertices){drawVertices();}
    if (displayEdges){drawEdges();}
    if (displayFaces){drawFaces();}
    if (displayCrust){drawEdges(true);}
}

void GLDisplayWidget::resizeGL(int width, int height)
{
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (GLfloat)width / (GLfloat)height, 0.1f, 100.0f);
    updateGL();
}

// To draw vertices
void GLDisplayWidget::drawVertices()
{
    for (int ix_vertex = 0; ix_vertex < _mesh.n_vertices; ix_vertex++)
    {

        // draw only non infinite points
        if(_mesh.vertices[ix_vertex].is_a_to_draw_point){

            glBegin(GL_POINTS);
            if(_mesh.vertices[ix_vertex].is_a_voronoi_center){
                glColor3f(0, 1, 0);
                glVertexDraw(_mesh.vertices[ix_vertex]);
            }
            else{
                glColor3f(1, 1, 1);
                glVertexDraw(_mesh.vertices[ix_vertex]);
            }
            glEnd();
        }
    }
}

//To draw edges
void GLDisplayWidget::drawEdges(bool crust)
{
    for (int ix_face = 0; ix_face < _mesh.n_faces; ix_face++)
    {
        Face face = _mesh.faces[ix_face];

        for (int ix_vertex_in_triangle = 0; ix_vertex_in_triangle < 3; ix_vertex_in_triangle++)
        {
            if(_mesh.vertices[face.ix_vertex[(ix_vertex_in_triangle + 1) % 3]].is_a_to_draw_point
                    && _mesh.vertices[face.ix_vertex[ix_vertex_in_triangle % 3]].is_a_to_draw_point){

                // only plot edges that aren't "infinite"

                if(!crust){

                    glColor3f(1, 0, 0);
                    glBegin(GL_LINE_STRIP);
                    glVertexDraw(_mesh.vertices[face.ix_vertex[ix_vertex_in_triangle % 3]]);
                    glVertexDraw(_mesh.vertices[face.ix_vertex[(ix_vertex_in_triangle + 1) % 3]]);
                    glEnd();

                }
                else if(!_mesh.vertices[face.ix_vertex[(ix_vertex_in_triangle + 1) % 3]].is_a_voronoi_center
                            && !_mesh.vertices[face.ix_vertex[ix_vertex_in_triangle % 3]].is_a_voronoi_center){

                    glColor3f(0, 0, 1);
                    glBegin(GL_LINE_STRIP);
                    glVertexDraw(_mesh.vertices[face.ix_vertex[ix_vertex_in_triangle % 3]]);
                    glVertexDraw(_mesh.vertices[face.ix_vertex[(ix_vertex_in_triangle + 1) % 3]]);
                    glEnd();
                }
            }
        }
    }
}

//To draw faces
void GLDisplayWidget::drawFaces()
{
    glColor3f(0, 1, 0);

    for (int ix_face = 0; ix_face < _mesh.n_faces; ix_face++)
    {
        Face face = _mesh.faces[ix_face];

        if(_mesh.vertices[face.ix_vertex[0]].is_a_to_draw_point
                && _mesh.vertices[face.ix_vertex[1]].is_a_to_draw_point
                && _mesh.vertices[face.ix_vertex[2]].is_a_to_draw_point){
            glBegin(GL_TRIANGLES);
            for (int ix_vertex_in_triangle = 0; ix_vertex_in_triangle < 3; ix_vertex_in_triangle++)
            {
                glVertexDraw(_mesh.vertices[face.ix_vertex[ix_vertex_in_triangle]]);
            }
            glEnd();
        }
    }
}


void GLDisplayWidget::add_random_vertex(){
    float LO = -3;
    float HI = 3;
    float rand_x = LO + static_cast <float> (std::rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
    float rand_y = LO + static_cast <float> (std::rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
    std::cout << "add random vertex : " << std::endl;
    std::cout << rand_x << std::endl;
    std::cout << rand_y << std::endl;
    Vertex v(rand_x,rand_y,0);
    _mesh.add_vertex(v);
}

void GLDisplayWidget::add_voronoi_centers(){
    _mesh.addVoronoiCentersToTriangulation();
}


// - - - - - - - - - - - - Mouse Management  - - - - - - - - - - - - - - - -
// When you click, the position of your mouse is saved
void GLDisplayWidget::mousePressEvent(QMouseEvent *event)
{
    if (event != NULL){
        _lastPosMouse = event->pos();
//        QPoint _a = event->globalPos();
//        QPointF _b = event->localPos();
//        QPoint p = QWidget::mapFromGlobal(_a);
//        QPointF _c = event->windowPos();
//        QPointF _d = event->screenPos();



//        std::cout << "last coords :" << std::endl;
//        std::cout << _lastPosMouse.x() << " xy " << _lastPosMouse.y() << std::endl;
//        std::cout << p.x() << " xy " << p.y() << std::endl;
//        std::cout << _a.x() << " xy " << _a.y() << std::endl;
//        std::cout << _b.x() << " xy " << _b.y() << std::endl;
//        std::cout << _c.x() << " xy " << _c.y() << std::endl;
//        std::cout << _d.x() << " xy " << _d.y() << std::endl;
    }
}


// Mouse movement management
void GLDisplayWidget::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - _lastPosMouse.x();
    // int dy = event->y() - lastPosMouse.y();

    if (event != NULL)
    {
        _angle += dx;
        _lastPosMouse = event->pos();

        updateGL();
    }
}

// Mouse Management for the zoom
void GLDisplayWidget::wheelEvent(QWheelEvent *event)
{
    QPoint numDegrees = event->angleDelta();
    float stepZoom = 0.1;
    if (!numDegrees.isNull())
    {
        _Z = (numDegrees.x() > 0 || numDegrees.y() > 0) ? _Z + stepZoom : _Z - stepZoom;
    }
}

