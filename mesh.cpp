#include "mesh.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>


//-----------------------------------------------------------------
// Vertex
//-----------------------------------------------------------------

// Constructors

Vertex::Vertex(){_x = 0;_y = 0;_z = 0;}
Vertex::Vertex(float x_, float y_, float z_){_x = x_;_y = y_;_z = z_;}

// Getter

float Vertex::x() const{return _x;}
float Vertex::y() const{return _y;}
float Vertex::z() const{return _z;}

// Operators

Vertex Vertex::cross(const Vertex &other_vertex) const
{
    float result_x = _y * other_vertex.z() - _z * other_vertex.y();
    float result_y = _z * other_vertex.x() - _x * other_vertex.z();
    float result_z = _x * other_vertex.y() - _y * other_vertex.x();
    return Vertex(result_x, result_y, result_z);
}

//-----------------------------------------------------------------
// Face
//-----------------------------------------------------------------

// Constructors

Face::Face(){ix_vertex[0] = 0;ix_vertex[1] = 0;ix_vertex[2] = 0;}
Face::Face(int ix_vertex_1_, int ix_vertex_2_, int ix_vertex_3_){ix_vertex[0] = ix_vertex_1_;ix_vertex[1] = ix_vertex_2_;ix_vertex[2] = ix_vertex_3_;}

//-----------------------------------------------------------------
// Mesh
//-----------------------------------------------------------------

// Constructors

Mesh::Mesh(){}

void Mesh::parseFile(const char file_name[])
{

    FILE *pFile;
    pFile = fopen(file_name, "r");

    if (pFile != NULL){std::cout << "opening " << file_name << " succeeded" << std::endl;}
    else{std::cout << "opening " << file_name << " failed" << std::endl;}

    vertices.clear();
    int n_edges = 0;
    faces.clear();
    n_faces = 0;

    // first line : n_vertices, n_faces, n_edges
    fscanf(pFile, "%d %d %d\n", &n_vertices, &n_faces, &n_edges);
    std::cout << "nb points: " << n_vertices << ", nb faces: " << n_faces << std::endl;


    vertices.reserve(n_vertices);
    faces.reserve(n_faces);

    float x, y, z;

    for (int ix_vertex = 0; ix_vertex < n_vertices; ix_vertex++)
    {
        fscanf(pFile, "%f %f %f\n", &x, &y, &z);
        vertices.push_back(Vertex(x, y, z));

    }
    int n_face, ix_vertex_1, ix_vertex_2, ix_vertex_3;
    for (int ix_face = 0; ix_face < n_faces; ix_face++)
    {
            fscanf(pFile, "%d %d %d %d\n", &n_face, &ix_vertex_1, &ix_vertex_2, &ix_vertex_3);
            faces.push_back(Face(ix_vertex_1, ix_vertex_2, ix_vertex_3));
    }

    fclose(pFile);

    std::cout << "end of reading" << std::endl;

}

// To sew the mesh

void Mesh::sew()
{
    std::cout << "sewing" << std::endl;

    bool is_sewed[n_vertices];
    for (int ix_vertex = 0; ix_vertex < n_vertices; ix_vertex++)
    {
        is_sewed[ix_vertex] = false;
    }

    std::map<std::pair<int, int>, std::pair<int, int>> pair_map;

    for (int ix_face = 0; ix_face < n_faces; ix_face++) // for each face
    {
        for (int ix_vertex_of_face = 0; ix_vertex_of_face < 3; ix_vertex_of_face++)
        {
            int ix_vertex = faces[ix_face].ix_vertex[ix_vertex_of_face];
            if (!is_sewed[ix_vertex])
            {
                vertices[ix_vertex].ix_incident_face = ix_face;
                is_sewed[ix_vertex] = true;
            };
        }

        // setting adjacent faces on each edge
        std::pair<int, int> edge;
        std::pair<int, int> pair_face_vertex;

        // Order vertex
        for (int ix_vertex_of_face = 0; ix_vertex_of_face < 3; ix_vertex_of_face++)
        {
            if (faces[ix_face].ix_vertex[(ix_vertex_of_face + 1) % 3] < faces[ix_face].ix_vertex[(ix_vertex_of_face + 2) % 3])
            {
                edge = {faces[ix_face].ix_vertex[(ix_vertex_of_face + 1) % 3], faces[ix_face].ix_vertex[(ix_vertex_of_face + 2) % 3]};
            }
            else
            {
                edge = {faces[ix_face].ix_vertex[(ix_vertex_of_face + 2) % 3], faces[ix_face].ix_vertex[(ix_vertex_of_face + 1) % 3]};
            }

            if (pair_map.find(edge) == pair_map.end()) // if pair does not exist yet
            {
                pair_map[edge] = {ix_face, ix_vertex_of_face}; // (ix_vertex_1, ix_vertex_2) : (ix_face, ix_vertex)
            }
            else
            {
                pair_face_vertex = pair_map.at(edge);
                int i_other_face = pair_face_vertex.first;
                int i_other_vertex = pair_face_vertex.second;

                faces[ix_face].adjacent_faces[ix_vertex_of_face] = i_other_face;
                faces[i_other_face].adjacent_faces[i_other_vertex] = ix_face;
            }
        }
    }
    std::cout << "end sewing" << std::endl;
}

// to add vertex in a mesh
void Mesh::add_vertex(Vertex v){
    vertices.push_back(v);
    n_vertices+=1;

    int ix_vertex = vertices.size()-1;
    for (int ix_face = 0; ix_face < n_faces; ix_face++)
    {
        if (inTriangleTest(faces[ix_face], vertices[ix_vertex]) > 0)
        {
            insertionInFace(ix_vertex, ix_face);
            lawsonAroundVertex(ix_vertex);
            break;
        }
    };
}

//is vertex P on the circumscribed circle ?
// = 0 on edge
// > 0 yes
// < 0 no
float Mesh::vertexInCircumscribingCircle(Face face, Vertex V)
{

    Vertex v1 = vertices[face.ix_vertex[0]];
    Vertex v2 = vertices[face.ix_vertex[1]];
    Vertex v3 = vertices[face.ix_vertex[2]];

    Vertex v1_sqrd = Vertex(v1.x(), v1.y(), v1.x() * v1.x() + v1.y() * v1.y());
    Vertex v2_sqrd = Vertex(v2.x(), v2.y(), v2.x() * v2.x() + v2.y() * v2.y());
    Vertex v3_sqrd = Vertex(v3.x(), v3.y(), v3.x() * v3.x() + v3.y() * v3.y());

    Vertex V_sqrd = Vertex(V.x(), V.y(), V.x() * V.x() + V.y() * V.y());

    Vertex v12_sqrd = v2_sqrd - v1_sqrd;
    Vertex v13_sqrd = v3_sqrd - v1_sqrd;
    Vertex v1V_sqrd = V_sqrd - v1_sqrd;

    Vertex crossprod = v12_sqrd.cross(v13_sqrd);
    float res = -(crossprod * v1V_sqrd);

    return res;
}

// test if an edge is locally Delaunay
bool Mesh::isDelaunay(int ix_face_1, int ix_opposit_vertex_1)
{

    Face face_1 = faces[ix_face_1];
    Face face_2;

    if (face_1.ix_vertex[0] == ix_opposit_vertex_1)
    {
        face_2 = faces[face_1.adjacent_faces[0]];
    }
    else if (face_1.ix_vertex[1] == ix_opposit_vertex_1)
    {
        face_2 = faces[face_1.adjacent_faces[1]];
    }
    else
    {
        face_2 = faces[face_1.adjacent_faces[2]];
    };

    int ix_opposit_vertex_2;

    if (face_2.adjacent_faces[0] == ix_face_1)
    {
        ix_opposit_vertex_2 = face_2.ix_vertex[0];
    }
    else if (face_2.adjacent_faces[1] == ix_face_1)
    {
        ix_opposit_vertex_2 = face_2.ix_vertex[1];
    }
    else
    {
        ix_opposit_vertex_2 = face_2.ix_vertex[2];
    };

    float test1 = vertexInCircumscribingCircle(face_1, vertices[ix_opposit_vertex_2]);
    float test2 = vertexInCircumscribingCircle(face_2, vertices[ix_opposit_vertex_1]);

    return ((test1 <= 0) && (test2 <= 0));
}

// to flip edge
QList<std::pair<int, int>> Mesh::flipEdge(int ix_face_1, int ix_initial_opposit_vertex_1)
{
    Face &face_1 = faces[ix_face_1];
    int ix_vertex_1 = ix_initial_opposit_vertex_1;

    int ix_face_2 = -1;
    int ix_face_3 = -1;
    int ix_face_4 = -1;
    int ix_face_5 = -1;
    int ix_face_6 = -1;

    int ix_vertex_2 = -1;
    int ix_vertex_3 = -1;
    int ix_vertex_4 = -1;

    if (ix_initial_opposit_vertex_1 == face_1.ix_vertex[0])
    {
        ix_face_2 = face_1.adjacent_faces[0];
        ix_face_4 = face_1.adjacent_faces[1];
        ix_face_3 = face_1.adjacent_faces[2];

        ix_vertex_4 = face_1.ix_vertex[1];
        ix_vertex_2 = face_1.ix_vertex[2];

    }
    else if (ix_initial_opposit_vertex_1 == face_1.ix_vertex[1])
    {
        ix_face_2 = face_1.adjacent_faces[1];
        ix_face_3 = face_1.adjacent_faces[0];
        ix_face_4 = face_1.adjacent_faces[2];

        ix_vertex_4 = face_1.ix_vertex[2];
        ix_vertex_2 = face_1.ix_vertex[0];
    }
    else
    {
        ix_face_2 = face_1.adjacent_faces[2];
        ix_face_3 = face_1.adjacent_faces[1];
        ix_face_4 = face_1.adjacent_faces[0];


        ix_vertex_4 = face_1.ix_vertex[0];
        ix_vertex_2 = face_1.ix_vertex[1];
    };

    Face &face_2 = faces[ix_face_2];

    if (ix_face_1 == face_2.adjacent_faces[0])
    {
        ix_face_5 = face_2.adjacent_faces[2];
        ix_face_6 = face_2.adjacent_faces[1];

        ix_vertex_3 = face_2.ix_vertex[0];

    }
    else if (ix_face_1 == face_2.adjacent_faces[1])
    {
        ix_face_5 = face_2.adjacent_faces[0];
        ix_face_6 = face_2.adjacent_faces[2];

        ix_vertex_3 = face_2.ix_vertex[1];
    }
    else
    {
        ix_face_5 = face_2.adjacent_faces[1];
        ix_face_6 = face_2.adjacent_faces[0];

        ix_vertex_3 = face_2.ix_vertex[2];

    };

    // update of sewing after flipping the edge

    face_1.ix_vertex[0] = ix_vertex_1;
    face_1.ix_vertex[1] = ix_vertex_4;
    face_1.ix_vertex[2] = ix_vertex_3;

    face_1.adjacent_faces[0] = ix_face_6;
    face_1.adjacent_faces[1] = ix_face_2;
    face_1.adjacent_faces[2] = ix_face_3;

    face_2.ix_vertex[0] = ix_vertex_3;
    face_2.ix_vertex[1] = ix_vertex_2;
    face_2.ix_vertex[2] = ix_vertex_1;

    face_2.adjacent_faces[0] = ix_face_4;
    face_2.adjacent_faces[1] = ix_face_1;
    face_2.adjacent_faces[2] = ix_face_5;

    for (int i = 0; i < 3; i++)
    {
        if (ix_face_4 >= 0 && faces[ix_face_4].adjacent_faces[i] == ix_face_1)
        {
            faces[ix_face_4].adjacent_faces[i] = ix_face_2;
        };
        if (ix_face_6 >= 0 && faces[ix_face_6].adjacent_faces[i] == ix_face_2)
        {
            faces[ix_face_6].adjacent_faces[i] = ix_face_1;
        };
    };

    // Vertices

    vertices[ix_vertex_1].ix_incident_face = ix_face_1;
    vertices[ix_vertex_4].ix_incident_face = ix_face_1;
    vertices[ix_vertex_3].ix_incident_face = ix_face_2;
    vertices[ix_vertex_2].ix_incident_face = ix_face_2;

    // get the new edges to test in lawson algorithm

    QList<std::pair<int, int>> quad_edges;
    quad_edges.push_back({ix_face_1, ix_vertex_1});
    quad_edges.push_back({ix_face_2, ix_vertex_1});

    return quad_edges;
}

// test if 3 points are oriented positively
float Mesh::orientationTest(Vertex v1, Vertex v2, Vertex v3)
{

    Vertex v21 = v2 - v1;
    Vertex v32 = v3 - v2;

    Vertex p1 = Vertex(v21.x(), v21.y(), 0);
    Vertex p2 = Vertex(v32.x(), v32.y(), 0);

    Vertex crossprod = p1.cross(p2);

    return crossprod.z();
}

// test if a vertex is inside a face
// 1 : inside
// 2 : on the edge
// 3 : outside the edge
float Mesh::inTriangleTest(Face face, Vertex vertex)
{

    Vertex v1 = vertices[face.ix_vertex[0]];
    Vertex v2 = vertices[face.ix_vertex[1]];
    Vertex v3 = vertices[face.ix_vertex[2]];

    float test1 = orientationTest(v1, v2, vertex);
    float test2 = orientationTest(v2, v3, vertex);
    float test3 = orientationTest(v3, v1, vertex);

    if ((test1 > 0) && (test2 > 0) && (test3 > 0))
    {
        return 1;
    }
    else if ((test1 >= 0) && (test2 >= 0) && (test3 >= 0))
    {
        return 0;
    }
    else
    {
        return -1;
    };
}


// insertion point i_P in face ix_face
void Mesh::insertionInFace(int i_P, int ix_face)
{
    Face old_face = faces[ix_face];

    int i_A = old_face.ix_vertex[0];
    int i_B = old_face.ix_vertex[1];
    int i_C = old_face.ix_vertex[2];

    int ix_face_A = old_face.adjacent_faces[0];
    int ix_face_B = old_face.adjacent_faces[1];
    int ix_face_C = old_face.adjacent_faces[2];

    Vertex &A = vertices[i_A];
    Vertex &B = vertices[i_B];
    Vertex &C = vertices[i_C];
    Vertex &P = vertices[i_P];

    // ABP
    faces[ix_face] = Face(i_A, i_B, i_P);
    int i_ABP = ix_face;

    // BCP
    faces.push_back(Face(i_B, i_C, i_P));
    int i_BCP = n_faces;

    // CAP
    faces.push_back(Face(i_C, i_A, i_P));
    int i_CAP = n_faces + 1;

    Face &ABP = faces[i_ABP];
    Face &BCP = faces[i_BCP];
    Face &CAP = faces[i_CAP];

    // Add faces and update adjacent faces around

    n_faces += 2;

    ABP.adjacent_faces[0] = i_BCP;
    ABP.adjacent_faces[1] = i_CAP;
    ABP.adjacent_faces[2] = ix_face_C;

    BCP.adjacent_faces[0] = i_CAP;
    BCP.adjacent_faces[1] = i_ABP;
    BCP.adjacent_faces[2] = ix_face_A;

    CAP.adjacent_faces[0] = i_ABP;
    CAP.adjacent_faces[1] = i_BCP;
    CAP.adjacent_faces[2] = ix_face_B;


    for (int i = 0; i < 3; i++)
    {
        if (ix_face_A >= 0 && faces[ix_face_A].adjacent_faces[i] == ix_face)
        {
            faces[ix_face_A].adjacent_faces[i] = i_BCP;
        };
        if (ix_face_B >= 0 && faces[ix_face_B].adjacent_faces[i] == ix_face)
        {
            faces[ix_face_B].adjacent_faces[i] = i_CAP;
        };
        if (ix_face_C >= 0 && faces[ix_face_C].adjacent_faces[i] == ix_face)
        {
            faces[ix_face_C].adjacent_faces[i] = i_ABP;
        };
    }

    // modify incident faces for the vertices

    A.ix_incident_face = i_ABP;
    B.ix_incident_face = i_BCP;
    C.ix_incident_face = i_CAP;
    P.ix_incident_face = i_ABP;
}



// test if edge opposite to ix_vertex is an "infinite" edge i.e, if it's an edge between two infinite points.
bool Mesh::infinitEdge(int ix_face, int ix_vertex)
{
    if (faces[ix_face].ix_vertex[0] == ix_vertex){return faces[ix_face].adjacent_faces[0] < 0;}
    else if (faces[ix_face].ix_vertex[1] == ix_vertex){return faces[ix_face].adjacent_faces[1] < 0;}
    else if (faces[ix_face].ix_vertex[2] == ix_vertex){ return faces[ix_face].adjacent_faces[2] < 0;}
    else
    {return true;};
}

// After inserting a vertex into a face, flip around the vertex.
void Mesh::lawsonAroundVertex(int i_P)
{

    QList<std::pair<int, int>> toFlipEdges;
    for (int ix_face = 0; ix_face < n_faces; ix_face++)
    {
        for (int i = 0; i < 3; i++)
        {
            if (faces[ix_face].ix_vertex[i] == i_P)
            {
                toFlipEdges.push_back({ix_face, i_P});
            }
        }
    }

    while (!toFlipEdges.isEmpty())
    {
        std::pair<int, int> edge = toFlipEdges.takeFirst();
        int ix_face = edge.first;
        int ix_vertex = edge.second;

        // we don't flip infinite edges
        if (!infinitEdge(ix_face, ix_vertex))
        {
            if (!isDelaunay(ix_face, ix_vertex))
            {
                QList<std::pair<int, int>> newEdgesToFlip = flipEdge(ix_face, ix_vertex);
                while (!newEdgesToFlip.isEmpty())
                {

                    std::pair<int, int> newEdge = newEdgesToFlip.takeFirst();
                    int ix_face_newEdge = newEdge.first;
                    int ix_vertex_newEdge = newEdge.second;
                    if (!infinitEdge(ix_face_newEdge, ix_vertex_newEdge))
                    {
                        if (!isDelaunay(ix_face_newEdge, ix_vertex_newEdge))
                        {
                            toFlipEdges.push_back(newEdge);
                        };
                    };
                };
            };
        };
    };
}


// ----------------------------------------------
// --------Crust Algorithm--------
// ----------------------------------------------


//Parse a file with vertices and store the data into vertices.
void Mesh::parseTriFile(const char file_name[])
{
    FILE *pFile;
    pFile = fopen(file_name, "r");

    if (pFile != NULL)
    {
        std::cout << "file " << file_name << " opened" << std::endl;
    }
    else
    {
        std::cout << "file " << file_name << " not opened" << std::endl;
    }

    vertices.clear();
    faces.clear();

    fscanf(pFile, "%d\n", &n_vertices);
    std::cout << "nb points: " << n_vertices << std::endl;

    n_faces = 0;
    vertices.reserve(n_vertices);

    float x, y, z;

    fscanf(pFile, "%f %f %f\n", &x, &y, &z);
    vertices.push_back(Vertex(x, y, z));

    for (int ix_vertex = 1; ix_vertex < n_vertices; ix_vertex++)
    {
        fscanf(pFile, "%f %f %f\n", &x, &y, &z);
        vertices.push_back(Vertex(x, y, z));
    }
    fclose(pFile);
    std::cout << "end of reading" << std::endl;

}

// Draw delaunay triangulation from vertices without triangles.
void Mesh::triangulationFromVertices()
{
    x_min = vertices[0].x();
    x_max = vertices[0].x();
    y_min = vertices[0].y();
    y_max = vertices[0].y();

    for (int ix_vertex = 1; ix_vertex < n_vertices; ix_vertex++)
    {
        if (vertices[ix_vertex].x() < x_min)
        {x_min = vertices[ix_vertex].x();};
        if (vertices[ix_vertex].x() > x_max)
        {x_max = vertices[ix_vertex].x();};
        if (vertices[ix_vertex].y() < y_min)
        {y_min = vertices[ix_vertex].y();};
        if (vertices[ix_vertex].y() > y_max)
        {y_max = vertices[ix_vertex].y();};
    };

    // add the square containing all the other vertices
    float x1 = x_min - (x_max-x_min);
    float x2 = x_max + (x_max-x_min);
    float y1 = y_min - (y_max-y_min);
    float y2 = y_max + (y_max-y_min);

    vertices.push_back(Vertex(x1, y1, 0));
    vertices.push_back(Vertex(x1, y2, 0));
    vertices.push_back(Vertex(x2, y1, 0));
    vertices.push_back(Vertex(x2, y2, 0));

    n_vertices += 4;

    Vertex &inf_vrtx_1 = vertices[n_vertices - 4];
    Vertex &inf_vrtx_2 = vertices[n_vertices - 3];
    Vertex &inf_vrtx_3 = vertices[n_vertices - 2];
    Vertex &inf_vrtx_4 = vertices[n_vertices - 1];

    faces.push_back(Face(n_vertices - 4, n_vertices - 2, n_vertices - 3));
    faces.push_back(Face(n_vertices - 1, n_vertices - 3, n_vertices - 2));

    n_faces += 2;

    Face &triangle_inf_vrtx_1 = faces[n_faces - 2];
    Face &triangle_inf_vrtx_4 = faces[n_faces - 1];

    inf_vrtx_1.ix_incident_face = n_faces - 2;
    inf_vrtx_2.ix_incident_face = n_faces - 2;
    inf_vrtx_3.ix_incident_face = n_faces - 1;
    inf_vrtx_4.ix_incident_face = n_faces - 1;

    // do not draw these points
    inf_vrtx_1.is_a_to_draw_point = false;
    inf_vrtx_2.is_a_to_draw_point = false;
    inf_vrtx_3.is_a_to_draw_point = false;
    inf_vrtx_4.is_a_to_draw_point = false;

    triangle_inf_vrtx_1.adjacent_faces[0] = n_faces - 1;
    triangle_inf_vrtx_1.adjacent_faces[1] = -1;
    triangle_inf_vrtx_1.adjacent_faces[2] = -1;

    triangle_inf_vrtx_4.adjacent_faces[0] = n_faces - 2;
    triangle_inf_vrtx_4.adjacent_faces[1] = -1;
    triangle_inf_vrtx_4.adjacent_faces[2] = -1;

    for (int ix_vertex = 0; ix_vertex < (n_vertices - 4); ix_vertex++)
    {
        for (int ix_face = 0; ix_face < n_faces; ix_face++)
        {
            if (inTriangleTest(faces[ix_face], vertices[ix_vertex]) > 0)
            {
                insertionInFace(ix_vertex, ix_face);
                lawsonAroundVertex(ix_vertex);
                break;
            }
        };
    };
}


//Return the Voronoi point coordinates associated to a Delaunay Triangle
std::vector<float> Mesh::voronoiCenter(int ix_face)
{
    std::vector<float> coordinates;

    Face &triangle = faces[ix_face];

    Vertex &A = vertices[triangle.ix_vertex[0]];
    Vertex &B = vertices[triangle.ix_vertex[1]];
    Vertex &C = vertices[triangle.ix_vertex[2]];

    std::cout << " orientation: " << orientationTest(A,B,C) << std::endl;

    // vectors
    Vertex AB = B.operator-(A);
    Vertex BA = A.operator-(B);
    Vertex BC = C.operator-(B);
    Vertex CB = B.operator-(C);
    Vertex AC = C.operator-(A);
    Vertex CA = A.operator-(C);

    // tangentes
    float tan_A = AB.cross(AC).z()/(AB.operator*(AC));
    float tan_B = BC.cross(BA).z()/(BC.operator*(BA));
    float tan_C = CA.cross(CB).z()/(CA.operator*(CB));

    // barycentric coords
    float alpha = tan_B+tan_C;
    float beta = tan_A + tan_C;
    float gamma = tan_B + tan_A;
    float sum_ainf_vrtx_1 = alpha + beta + gamma;

    float x_voronoi = (alpha*A.x() + beta*B.x() + gamma*C.x())/sum_ainf_vrtx_1;
    coordinates.push_back(x_voronoi);
    float y_voronoi = (alpha*A.y() + beta*B.y() + gamma*C.y())/sum_ainf_vrtx_1;
    coordinates.push_back(y_voronoi);
    float z_voronoi = (alpha*A.z() + beta*B.z() + gamma*C.z())/sum_ainf_vrtx_1;
    coordinates.push_back(z_voronoi);

    return coordinates;
}


void Mesh::addVoronoiCentersToTriangulation(){

    int nb_voronoi_vertices = 0;

    for (int ix_face = 0; ix_face < n_faces; ix_face++)
    // we calculate voronoi dual point for each faces that is not an infinite face
    {
        if(vertices[faces[ix_face].ix_vertex[0]].is_a_to_draw_point
                && vertices[faces[ix_face].ix_vertex[1]].is_a_to_draw_point
                && vertices[faces[ix_face].ix_vertex[2]].is_a_to_draw_point){

            std::vector<float> voronoiCoordinates = voronoiCenter(ix_face);
            std::cout << " vorCoord0: " << voronoiCoordinates[0] << std::endl;
            std::cout << " vorCoord1: " << voronoiCoordinates[1] << std::endl;
            std::cout << " vorCoord2: " << voronoiCoordinates[2] << std::endl;

            vertices.push_back(Vertex(voronoiCoordinates[0],voronoiCoordinates[1],voronoiCoordinates[2]));
            n_vertices+=1;
            nb_voronoi_vertices+=1;

            // define this point as a voronoi point
            vertices[n_vertices-1].is_a_voronoi_center = true;
        }
    }

    // now we insert precedent points into the mesh structure
    for (int i_voronoix_vertex = 0; i_voronoix_vertex<nb_voronoi_vertices; i_voronoix_vertex++){


        for (int ix_face = 0; ix_face < n_faces; ix_face++){

            if (inTriangleTest(faces[ix_face], vertices[n_vertices-nb_voronoi_vertices+i_voronoix_vertex]) > 0)
            {
                insertionInFace(n_vertices-nb_voronoi_vertices+i_voronoix_vertex, ix_face);
                lawsonAroundVertex(n_vertices-nb_voronoi_vertices+i_voronoix_vertex);
                break;
            }
        };
    };
}


// To Do next : CRUST --> une fonction qui calcule les centres de voronoi de chaque triangle,
// puis qui refait la triangulation et puis qui dessine que les arÃªtes composÃ©es seulement des sommets initiaux
