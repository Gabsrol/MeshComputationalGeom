#ifndef MESH_H
#define MESH_H

#include <QGLWidget>
#include <iostream>
#include <cmath>

class Vertex
{

private:

    // coordonates

    float _x,_y,_z;

public:

    // constructors

    Vertex();
    Vertex(float x_, float y_, float z_);

    // getters

    float x() const;
    float y() const;
    float z() const;

    // some useful attributes

    int ix_incident_face; // index of an incident face
    bool is_a_to_draw_point = true; // to prevent "infinite" points from being drawn
    bool is_a_voronoi_center = false; // used to implement CRUST algorithm

    // some useful methods or operators

    float getNorm() const { return std::sqrt(_x * _x + _y * _y + _z * _z); } // The L2 norm of the vertex, seen as a vector.
    float operator*(const Vertex &p) const { return _x * p._x + _y * p._y + _z * p._z; }       // Scalar product with another vector
    Vertex operator*(float a) const { return Vertex(_x * a, _y * a, _z * a); }                 // Product with a scalar
    Vertex operator/(float a) const { return Vertex(_x / a, _y / a, _z / a); }                 // Division with a scalar !=0
    Vertex operator+(const Vertex &p) const { return Vertex(_x + p._x, _y + p._y, _z + p._z); } // Addition with another vectors
    Vertex operator-(const Vertex &p) const { return Vertex(_x - p._x, _y - p._y, _z - p._z); } // Substraction with another vectors
    bool operator!=(const Vertex other_vertex) const;
    Vertex cross(const Vertex &p) const; // Cross product
};

class Face
{
public:

    // constructors

    Face();
    Face(int ix_vertex0_, int ix_vertex1_, int ix_vertex2_);

    // attributes for sewing

    int ix_vertex[3];       // indices of the 3 vertices forming the face
    int adjacent_faces[3]; // indices of the 3 adjacent faces

};



class Mesh
{
public:

    // constructor

    Mesh();

    // attributes

    std::vector<Vertex> vertices; // list of vertices
    std::vector<Face> faces; // list of faces
    int n_vertices, n_faces; // number of vertices and faces
    float x_min = -3;
    float x_max = 3;
    float y_min = -3;
    float y_max = 3;

    // file parsing

    void parseFile(const char file[]);
    void parseTriFile(const char file[]);

    // sewing, triangulation and voronoi

    void sew();
    void triangulationFromVertices();
    void addVoronoiCentersToTriangulation();

    // other useful methods

    void add_vertex(Vertex v);
    float vertexInCircumscribingCircle(Face face, Vertex P);
    bool isDelaunay(int ix_face1, int ix_vertex_oppose_1);
    QList<std::pair<int, int>> flipEdge(int ix_face1, int ix_vertex_oppose_1_initial);
    float orientationTest(Vertex A, Vertex B, Vertex C);
    float inTriangleTest(Face face, Vertex P);
    void insertionInFace(int i_P, int ix_face);
    bool infinitEdge(int ix_face, int ix_vertex);
    void lawsonAroundVertex(int i_P);
    std::vector<float> voronoiCenter(int ix_face);
};

#endif
