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

bool Vertex::operator!=(const Vertex other_vertex) const
{
    return ((x() != other_vertex.x()) | (y() != other_vertex.y()) | (z() != other_vertex.z()));
}

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
Face::Face(int ix_vertex0_, int ix_vertex1_, int ix_vertex2_){ix_vertex[0] = ix_vertex0_;ix_vertex[1] = ix_vertex1_;ix_vertex[2] = ix_vertex2_;}

//-----------------------------------------------------------------
// Mesh
//-----------------------------------------------------------------

// Constructors

Mesh::Mesh(){}

void Mesh::parseFile(const char file_name[])
{
    // reading file

    FILE *pFile;
    pFile = fopen(file_name, "r");

    // print success or failure

    if (pFile != NULL){std::cout << "opening " << file_name << " succeeded" << std::endl;}
    else{std::cout << "opening " << file_name << " failed" << std::endl;}

    // Initialize vertices and number of edges

    vertices.clear();
    int n_edges = 0;

    // first line : n_vertices, n_faces, n_edges
    fscanf(pFile, "%d %d %d\n", &n_vertices, &n_faces, &n_edges);
    std::cout << "nb points: " << n_vertices << ", nb faces: " << n_faces << std::endl;

    // Reserve the right size for vertices and faces.
    vertices.reserve(n_vertices);
    faces.reserve(n_faces);

    // Initialize working variables
    float x, y, z;
    float x_min, x_max, y_min, y_max, z_min, z_max;

    // 1st point :
    fscanf(pFile, "%f %f %f\n", &x, &y, &z); // Stockage de la premiere ligne
    vertices.push_back(Vertex(x, y, z));     // Ajout du point dans le vecteur vertices

    x_min = x;
    x_max = x;
    y_min = y;
    y_max = y;
    z_min = z;
    z_max = z;

    // Tous les autres points :
    for (int ix_vertex = 1; ix_vertex < n_vertices; ix_vertex++)
    {
        fscanf(pFile, "%f %f %f\n", &x, &y, &z); // Stockage de la ligne lue
        vertices.push_back(Vertex(x, y, z));     // Ajout du point dans le vecteur vertices

        if (x < x_min)
        {
            x_min = x;
        };
        if (x > x_max)
        {
            x_max = x;
        };
        if (y < y_min)
        {
            y_min = y;
        };
        if (y > y_max)
        {
            y_max = y;
        };
        if (z < z_min)
        {
            z_min = z;
        };
        if (z > z_max)
        {
            z_max = z;
        };
    }

    // On remet Ã  jour la liste en centrant cette fois les coordonnees
    float x_middle = (x_max + x_min) / 2;
    float y_middle = (y_max + y_min) / 2;
    float z_middle = (z_max + z_min) / 2;

    // Mise Ã  jour de tous les vertex en les remplaÃ§ant par leur Ã©quivalent centrÃ© en 0,0,0.
    for (int ix_vertex = 0; ix_vertex < n_vertices; ix_vertex++)
    {
        Vertex vertex = vertices[ix_vertex];
        vertices[ix_vertex] = Vertex(vertex.x() - x_middle, vertex.y() - y_middle, vertex.z() - z_middle);
    }

    // Stocke toutes les faces dans faces
    int n_face, ix_vertex0, ix_vertex1, ix_vertex2;
    for (int ix_face = 0; ix_face < n_faces; ix_face++)
    {
        fscanf(pFile, "%d %d %d %d\n", &n_face, &ix_vertex0, &ix_vertex1, &ix_vertex2);
        faces.push_back(Face(ix_vertex0, ix_vertex1, ix_vertex2));
    }

    // Logs the ending of the process
    fclose(pFile);
    std::cout << "end of reading" << std::endl;
}

/**
 * Connects the index of one adjacent face for every vertices of the mesh.
 * Connects the 3 adjacent faces for every faces of the mesh.
 */
void Mesh::sew()
{
    std::cout << "begin sewing" << std::endl;

    // vert_done est un liste de boolean. Si vert_done[5] == true, Ã§a signifie que le vertex 5 Ã  Ã©tÃ© traitÃ©
    // Initilisation de vert_done
    bool vert_done[n_vertices];
    for (int ix_vertex = 0; ix_vertex < n_vertices; ix_vertex++)
    {
        vert_done[ix_vertex] = false;
    }

    // Map est un dict. ici : {(int, int) : (int, int), ...}
    std::map<std::pair<int, int>, std::pair<int, int>> myMap;

    for (int ix_face = 0; ix_face < n_faces; ix_face++)
    {
        // On se balade dans toutes les faces dÃ©finie, et on rÃ©cupÃ¨re les indices des trois points(=vertex)
        Face *face = &faces[ix_face]; // Un pointeur qui pointe vers la i_eme face

        // On regarde si on a dÃ©ja traitÃ© les vertex(=points). Si non, on dit qu'on les a traitÃ©
        // puis on ajoute l'indice de la face comme face incidente au vertex(=point)
        for (int ix_vertex_in_triangle = 0; ix_vertex_in_triangle < 3; ix_vertex_in_triangle++)
        {
            int ix_vertex = face->ix_vertex[ix_vertex_in_triangle];
            if (!vert_done[ix_vertex])
            {
                vertices[ix_vertex].ix_incident_face = ix_face; // On ajoute Ã  l'objet vertex une face inscidente
                vert_done[ix_vertex] = true;                     // On le considÃ¨re traitÃ©
            };
        }

        // setting adjacent faces on each edge:
        std::pair<int, int> edge;                            //  ArrÃªtes
        std::pair<int, int> i_other_face_and_i_other_vertex; //  {face id whose edge was'nt yet register in the map,  its opposed vertex id (0,1 or 2) }

        // Order vertex
        for (int ix_vertex_in_triangle = 0; ix_vertex_in_triangle < 3; ix_vertex_in_triangle++)
        {
            if (face->ix_vertex[(ix_vertex_in_triangle + 1) % 3] < face->ix_vertex[(ix_vertex_in_triangle + 2) % 3])
            {
                edge = {face->ix_vertex[(ix_vertex_in_triangle + 1) % 3], face->ix_vertex[(ix_vertex_in_triangle + 2) % 3]};
            }
            else
            {
                edge = {face->ix_vertex[(ix_vertex_in_triangle + 2) % 3], face->ix_vertex[(ix_vertex_in_triangle + 1) % 3]};
            }
            // edge est la pair d'indice dans l'ordre croissant de l'arrÃ¨te

            if (myMap.find(edge) == myMap.end())
            { // Si l'arrÃ¨te n'est pas dans la map
                // On l'ajoute
                myMap[edge] = {ix_face, ix_vertex_in_triangle}; // (ix_vertex0, ix_vertex1) : (ix_face, ix_vertex) (c'est le vertex opposÃ©)
            }
            else
            {                                                     // L'arrÃªte appartient dÃ©ja Ã  un autre triangle de la map
                i_other_face_and_i_other_vertex = myMap.at(edge); //(ix_face de l'autre triangle, ix_face dans l'autre triangle (0,1,2))
                int i_other_face = i_other_face_and_i_other_vertex.first;
                int i_other_vertex = i_other_face_and_i_other_vertex.second;

                face->adjacent_faces[ix_vertex_in_triangle] = i_other_face; // On stocke l'indice de la face adjacente opposÃ© au vertex(=point) nÂ°2
                // On fait la mÃªme chose Ã  l'autre triangle
                faces[i_other_face].adjacent_faces[i_other_vertex] = ix_face;
            }
        }
    }
    std::cout << "end sewing" << std::endl;
}

void Mesh::add_vertex(Vertex v){
    vertices.push_back(v);
    n_vertices+=1;

    int ix_vertex = vertices.size()-1;
    for (int ix_face = 0; ix_face < n_faces; ix_face++)
    {
        if (inTriangleTest(faces[ix_face], vertices[ix_vertex]) > 0)
        {
            insertionTriangle(ix_vertex, ix_face);
            lawsonAroundVertex(ix_vertex);
            break;
        }
        else if (inTriangleTest(faces[ix_face], vertices[ix_vertex]) == 0)
        {
            insertionInEdge(ix_face, ix_vertex); // pb ici probablement
            lawsonAroundVertex(ix_vertex);
            break;
        };
        // Il manque le cas inTriangleTest == 0 (ie sur une arete; cf insertion arete qui marche pas (?) pr le moment)
    };
}

/**
 * Tests if a vertex is within the circumscribed circle of the face.
 *
 * @param face for the test.
 * @param P the vertex for which we want to know if it's in the circumscribed circle.
 * @return float = 0 if P is on the circumscribed circle.
 *                > 0 if P is inside the circumscribed circle.
 *                < 0 if P is outside the circumscribed circle.
 */
float Mesh::vertexInCircumscribingCircle(Face face, Vertex P)
{
    // On analyse seulement les dimensions x, y
    // On regarde en projetant sur le paraboloÃ¯de z = x*x + y*y

    Vertex A = vertices[face.ix_vertex[0]];
    Vertex B = vertices[face.ix_vertex[1]];
    Vertex C = vertices[face.ix_vertex[2]];

    Vertex A_hyper = Vertex(A.x(), A.y(), A.x() * A.x() + A.y() * A.y());
    Vertex B_hyper = Vertex(B.x(), B.y(), B.x() * B.x() + B.y() * B.y());
    Vertex C_hyper = Vertex(C.x(), C.y(), C.x() * C.x() + C.y() * C.y());
    Vertex P_hyper = Vertex(P.x(), P.y(), P.x() * P.x() + P.y() * P.y());

    Vertex AB_hyper = B_hyper - A_hyper;
    Vertex AC_hyper = C_hyper - A_hyper;
    Vertex AP_hyper = P_hyper - A_hyper;

    Vertex pdtVect = AB_hyper.cross(AC_hyper);
    float pdtScal = -(pdtVect * AP_hyper);

    return pdtScal;
}

/**
 * Tests whether an edge is Delaunay.
 *
 * @param ix_face1 a face that contains the edge.
 * @param ix_vertex_oppose_1 the index of the vertex opposite the edge
 * @returns true if Delaunay, false otherwise.
 */
bool Mesh::isDelaunay(int ix_face1, int ix_vertex_oppose_1)
{

    Face face1 = faces[ix_face1];
    Face face2;
    if (face1.ix_vertex[0] == ix_vertex_oppose_1)
    {
        face2 = faces[face1.adjacent_faces[0]];
    }
    else if (face1.ix_vertex[1] == ix_vertex_oppose_1)
    {
        face2 = faces[face1.adjacent_faces[1]];
    }
    else
    {
        face2 = faces[face1.adjacent_faces[2]];
    };
    int ix_vertex_oppose_2;
    if (face2.adjacent_faces[0] == ix_face1)
    {
        ix_vertex_oppose_2 = face2.ix_vertex[0];
    }
    else if (face2.adjacent_faces[1] == ix_face1)
    {
        ix_vertex_oppose_2 = face2.ix_vertex[1];
    }
    else
    {
        ix_vertex_oppose_2 = face2.ix_vertex[2];
    };

    float test1 = vertexInCircumscribingCircle(face1, vertices[ix_vertex_oppose_2]);
    float test2 = vertexInCircumscribingCircle(face2, vertices[ix_vertex_oppose_1]);

    // On demande Ã  ce que les points soient sur la bordure ou strictement Ã  l'extÃ©rieur du cercle
    return ((test1 <= 0) && (test2 <= 0));
}

/**
 * Flips the edge.
 * It is therefore assumed that the edge is not on the contour.
 *
 * @param ix_face1 the index of the first triangle to which the edge belongs.
 * @param vertex_oppose_1_initial the opposite vertex, belonging to the
 *     second triangle whose edge is common.
 * @returns the 2 edges that frame the flipped edge, and not being in
 *     the same initial triangle as vertex_oppose_1_initial.
 */
QList<std::pair<int, int>> Mesh::flipEdge(int ix_face1, int ix_vertex_oppose_1_initial)
{
    Face &face1 = faces[ix_face1];

    int ix_vertex_a = ix_vertex_oppose_1_initial;
    int ix_vertex_b = -1;
    int ix_vertex_c = -1;
    int ix_vertex_d = -1;
    int ix_face2 = -1;
    int ix_face3 = -1;
    int ix_face4 = -1;
    int ix_face5 = -1;
    int ix_face6 = -1;

    if (ix_vertex_oppose_1_initial == face1.ix_vertex[0])
    {
        ix_face2 = face1.adjacent_faces[0];
        ix_vertex_d = face1.ix_vertex[1];
        ix_vertex_b = face1.ix_vertex[2];
        ix_face4 = face1.adjacent_faces[1];
        ix_face3 = face1.adjacent_faces[2];
    }
    else if (ix_vertex_oppose_1_initial == face1.ix_vertex[1])
    {
        ix_face2 = face1.adjacent_faces[1];
        ix_vertex_d = face1.ix_vertex[2];
        ix_vertex_b = face1.ix_vertex[0];
        ix_face4 = face1.adjacent_faces[2];
        ix_face3 = face1.adjacent_faces[0];
    }
    else
    {
        ix_face2 = face1.adjacent_faces[2];
        ix_vertex_d = face1.ix_vertex[0];
        ix_vertex_b = face1.ix_vertex[1];
        ix_face4 = face1.adjacent_faces[0];
        ix_face3 = face1.adjacent_faces[1];
    };
    Face &face2 = faces[ix_face2];

    if (ix_face1 == face2.adjacent_faces[0])
    {

        ix_vertex_c = face2.ix_vertex[0];
        ix_face6 = face2.adjacent_faces[1];
        ix_face5 = face2.adjacent_faces[2];
    }
    else if (ix_face1 == face2.adjacent_faces[1])
    {

        ix_vertex_c = face2.ix_vertex[1];
        ix_face6 = face2.adjacent_faces[2];
        ix_face5 = face2.adjacent_faces[0];
    }
    else
    {

        ix_vertex_c = face2.ix_vertex[2];
        ix_face6 = face2.adjacent_faces[0];
        ix_face5 = face2.adjacent_faces[1];
    };

    // Actualisation des donnees apres le flip

    // Faces 1 et 2

    face1.ix_vertex[0] = ix_vertex_a;
    face1.ix_vertex[1] = ix_vertex_d;
    face1.ix_vertex[2] = ix_vertex_c;

    face1.adjacent_faces[0] = ix_face6;
    face1.adjacent_faces[1] = ix_face2;
    face1.adjacent_faces[2] = ix_face3;

    face2.ix_vertex[0] = ix_vertex_c;
    face2.ix_vertex[1] = ix_vertex_b;
    face2.ix_vertex[2] = ix_vertex_a;

    face2.adjacent_faces[0] = ix_face4;
    face2.adjacent_faces[1] = ix_face1;
    face2.adjacent_faces[2] = ix_face5;

    // Faces 3 4 5 6
    // On actualise la face adj partant du vertex Ã  l'oppose de notre quadrilatere

    for (int i = 0; i < 3; i++)
    {
        if (ix_face4 >= 0 && faces[ix_face4].adjacent_faces[i] == ix_face1)
        {
            faces[ix_face4].adjacent_faces[i] = ix_face2;
        };
        if (ix_face6 >= 0 && faces[ix_face6].adjacent_faces[i] == ix_face2)
        {
            faces[ix_face6].adjacent_faces[i] = ix_face1;
        };
    };

    // Vertices

    vertices[ix_vertex_a].ix_incident_face = ix_face1;
    vertices[ix_vertex_d].ix_incident_face = ix_face1;
    vertices[ix_vertex_c].ix_incident_face = ix_face2;
    vertices[ix_vertex_b].ix_incident_face = ix_face2;

    QList<std::pair<int, int>> aretes_quadrilatere;
    aretes_quadrilatere.push_back({ix_face1, ix_vertex_a}); // arete f6
    aretes_quadrilatere.push_back({ix_face2, ix_vertex_a}); // arete f5

    return aretes_quadrilatere;
}

/**
 * Tests if A, B and C are oriented in the trigonometric direction.
 *
 * @param A, B, C vertices for the orientation test.
 * @returns float which is positive if A B C are oriented in the
 *     trigonometric direction, zero if aligned, negative otherwise.
*/
float Mesh::orientationTest(Vertex A, Vertex B, Vertex C)
{
    // On travaille dans le plan z = 0
    Vertex AB = B - A;
    Vertex BC = C - B;

    Vertex AB_proj = Vertex(AB.x(), AB.y(), 0);
    Vertex BC_proj = Vertex(BC.x(), BC.y(), 0);

    Vertex vector = AB_proj.cross(BC_proj);
    float value = vector.z(); //*n_vertices;
    return value;
}

/**
 * Tests if a vertex is inside a face.
 *
 * @param face used for the test.
 * @param vertex we want to know if it's inside the face.
 * @returns float which is positive if the vertex is inside the
 *     triangle, 0 if it is on an edge, negative if it is outside.
*/
float Mesh::inTriangleTest(Face face, Vertex vertex)
{
    // On travaille dans le plan z = 0
    Vertex A = vertices[face.ix_vertex[0]];
    Vertex B = vertices[face.ix_vertex[1]];
    Vertex C = vertices[face.ix_vertex[2]];

    float test1 = orientationTest(A, B, vertex);
    float test2 = orientationTest(B, C, vertex);
    float test3 = orientationTest(C, A, vertex);

    if ((test1 > 0) && (test2 > 0) && (test3 > 0))
    {
        // Renvoi val strictement positive si toutes les orientation sont strictement positive
        // Ou 0 si au moins une des orientations est nulle (ie vertex sur un bord)
        return 1;
    }
    else if ((test1 >= 0) && (test2 >= 0) && (test3 >= 0))
    {
        return 0;
    }
    else
    {
        // On est forcÃ©ment en dehors du triangle
        // On renvoi une valeur nÃ©gative
        return -1;
    };
}

/**
 * Insert a point in a face.
 *
 * @param i_P index of the point to be inserted.
 * @param ix_face index of the face in which we want to insert the vertex.
 */
void Mesh::insertionTriangle(int i_P, int ix_face)
{
    Face old_face = faces[ix_face];

    // Creer nouveaux triangles

    int i_A = old_face.ix_vertex[0];
    int i_B = old_face.ix_vertex[1];
    int i_C = old_face.ix_vertex[2];

    // On prend les faces autour du triangle si les voisins existent
    // Si l'indice  = -1, c'est qu'il n'y a pas de face voisine (cf plus tard)

    int ix_face_A = old_face.adjacent_faces[0];
    int ix_face_B = old_face.adjacent_faces[1];
    int ix_face_C = old_face.adjacent_faces[2];

    Vertex &A = vertices[i_A];
    Vertex &B = vertices[i_B];
    Vertex &C = vertices[i_C];
    Vertex &P = vertices[i_P];

    // InsÃ©rer ces triangles dans le tableau faces et supprimer l'ancien

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

    // Actualiser paramÃ¨tres du Mesh

    n_faces += 2;

    // Ajouter les adjacences des nouveaux triangles

    ABP.adjacent_faces[0] = i_BCP;
    ABP.adjacent_faces[1] = i_CAP;
    ABP.adjacent_faces[2] = ix_face_C;

    BCP.adjacent_faces[0] = i_CAP;
    BCP.adjacent_faces[1] = i_ABP;
    BCP.adjacent_faces[2] = ix_face_A;

    CAP.adjacent_faces[0] = i_ABP;
    CAP.adjacent_faces[1] = i_BCP;
    CAP.adjacent_faces[2] = ix_face_B;

    // Adapter les adjacences des triangles autour

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

    // Modifier les ix_incident_face des Vertices

    A.ix_incident_face = i_ABP;
    B.ix_incident_face = i_BCP;
    C.ix_incident_face = i_CAP;
    P.ix_incident_face = i_ABP;
}

/**
 * Inserts a vertex on an edge
 *
 * @param ix_face1 index of the face containing the edge.
 * @param i_P index of the vertex to be inserted.
 */
void Mesh::insertionInEdge(int ix_face1, int i_P)
{
    Face face1 = faces[ix_face1];

    int i_A = -1;
    int i_B = -1;
    int i_C = -1;
    int i_D = -1;
    int ix_face2 = -1;
    int ix_face3 = -1;
    int ix_face4 = -1;
    int ix_face5 = -1;
    int ix_face6 = -1;

    if (orientationTest(vertices[face1.ix_vertex[0]], vertices[face1.ix_vertex[1]], vertices[i_P]) == 0)
    {
        i_A = face1.ix_vertex[2];
        i_B = face1.ix_vertex[0];
        i_D = face1.ix_vertex[1];
        ix_face2 = face1.adjacent_faces[2];
        ix_face6 = face1.adjacent_faces[0];
        ix_face3 = face1.adjacent_faces[1];
    }
    else if (orientationTest(vertices[face1.ix_vertex[1]], vertices[face1.ix_vertex[2]], vertices[i_P]) == 0)
    {
        i_A = face1.ix_vertex[0];
        i_B = face1.ix_vertex[1];
        i_D = face1.ix_vertex[2];
        ix_face2 = face1.adjacent_faces[0];
        ix_face6 = face1.adjacent_faces[1];
        ix_face3 = face1.adjacent_faces[2];
    }
    else if (orientationTest(vertices[face1.ix_vertex[2]], vertices[face1.ix_vertex[0]], vertices[i_P]) == 0)
    {
        i_A = face1.ix_vertex[1];
        i_B = face1.ix_vertex[2];
        i_D = face1.ix_vertex[0];
        ix_face2 = face1.adjacent_faces[1];
        ix_face6 = face1.adjacent_faces[2];
        ix_face3 = face1.adjacent_faces[0];
    };

    for (int i = 0; i < 3; i++)
    {
        if (faces[ix_face2].adjacent_faces[i] == ix_face1)
        {
            i_C = faces[ix_face2].ix_vertex[i];
            ix_face4 = faces[ix_face2].adjacent_faces[(i + 1) % 3];
            ix_face5 = faces[ix_face2].adjacent_faces[(i + 2) % 3];
        };
    };

    Vertex &A = vertices[i_A];
    Vertex &B = vertices[i_B];
    Vertex &C = vertices[i_C];
    Vertex &D = vertices[i_D];
    Vertex &P = vertices[i_P];

    // InsÃ©rer ces triangles dans le tableau faces et supprimer l'ancien

    faces[ix_face1] = Face(i_D, i_A, i_P); // DAP
    faces[ix_face2] = Face(i_C, i_D, i_P); // CDP

    faces.push_back(Face(i_A, i_B, i_P)); // ABP
    faces.push_back(Face(i_B, i_C, i_P)); // BCP

    int i_DAP = ix_face1;
    int i_CDP = ix_face2;
    int i_ABP = n_faces;
    int i_BCP = n_faces + 1;

    Face &DAP = faces[i_DAP];
    Face &CDP = faces[i_CDP];
    Face &ABP = faces[i_ABP];
    Face &BCP = faces[i_BCP];

    // Actualiser paramÃ¨tres du Mesh

    n_faces += 2;

    // Ajouter les adjacences des nouveaux triangles

    DAP.adjacent_faces[0] = i_ABP;
    DAP.adjacent_faces[1] = i_CDP;
    DAP.adjacent_faces[2] = ix_face6;

    CDP.adjacent_faces[0] = i_DAP;
    CDP.adjacent_faces[1] = i_BCP;
    CDP.adjacent_faces[2] = ix_face5;

    ABP.adjacent_faces[0] = i_BCP;
    ABP.adjacent_faces[1] = i_DAP;
    ABP.adjacent_faces[2] = ix_face3;

    BCP.adjacent_faces[0] = i_CDP;
    BCP.adjacent_faces[1] = i_ABP;
    BCP.adjacent_faces[2] = ix_face4;

    // Adapter les adjacences des triangles autour

    for (int i = 0; i < 3; i++)
    {
        if (ix_face3 >= 0 && faces[ix_face3].adjacent_faces[i] == ix_face1)
        {
            faces[ix_face3].adjacent_faces[i] = i_ABP;
        };
        if (ix_face4 >= 0 && faces[ix_face4].adjacent_faces[i] == ix_face2)
        {
            faces[ix_face4].adjacent_faces[i] = i_BCP;
        };
        if (ix_face5 >= 0 && faces[ix_face5].adjacent_faces[i] == ix_face2)
        {
            faces[ix_face5].adjacent_faces[i] = i_CDP;
        };
        if (ix_face6 >= 0 && faces[ix_face6].adjacent_faces[i] == ix_face1)
        {
            faces[ix_face6].adjacent_faces[i] = i_DAP;
        };
    }

    // Modifier les ix_incident_face des Vertices

    A.ix_incident_face = i_ABP;
    B.ix_incident_face = i_BCP;
    C.ix_incident_face = i_CDP;
    D.ix_incident_face = i_DAP;
    P.ix_incident_face = i_ABP;
}

/**
 * Tests whether the edge opposite the vertex is on the contour.
 *
 * @param ix_face index of the face containing the edge.
 * @param ix_vertex index of the vertex opposite the edge.
 * @returns true if the edge is on the contour, false otherwise.
 */
bool Mesh::infinitEdge(int ix_face, int ix_vertex)
{
    if (faces[ix_face].ix_vertex[0] == ix_vertex)
    {
        return faces[ix_face].adjacent_faces[0] < 0;
    }
    else if (faces[ix_face].ix_vertex[1] == ix_vertex)
    {
        return faces[ix_face].adjacent_faces[1] < 0;
    }
    else if (faces[ix_face].ix_vertex[2] == ix_vertex)
    {
        return faces[ix_face].adjacent_faces[2] < 0;
    }
    else
    {
        return true;
    };
}

/**
 * After inserting a vertex into a face, flip around the vertex until the triangle is Delaunay.
 *
 * @param i_P index of the vertex
 */
void Mesh::lawsonAroundVertex(int i_P)
{
    // On part d'un vertex i_P (dans ix_face) qui vient d'Ãªtre insere, et on fait des flips recursifs pour qu'a la fin il soit insere et tout soit Delaunay
    // On recupere les trois ou quatre aretes autour du P dans une file
    QList<std::pair<int, int>> atraiter;
    for (int ix_face = 0; ix_face < n_faces; ix_face++)
    {
        for (int i = 0; i < 3; i++)
        {
            if (faces[ix_face].ix_vertex[i] == i_P)
            {
                atraiter.push_back({ix_face, i_P});
            }
        }
    }
    // On lance la boucle while et on remplie et traite la file
    while (!atraiter.isEmpty())
    {
        std::pair<int, int> face_et_vertex = atraiter.takeFirst();
        int ix_face = face_et_vertex.first;
        int ix_vertex = face_et_vertex.second;

        if (!infinitEdge(ix_face, ix_vertex))
        {
            if (!isDelaunay(ix_face, ix_vertex))
            {
                QList<std::pair<int, int>> nouvelle_queue = flipEdge(ix_face, ix_vertex); // On fait le flip et rÃ©cupÃ¨re les arete Ã  retester
                while (!nouvelle_queue.isEmpty())
                {

                    std::pair<int, int> face_et_vertex_quadrilatere = nouvelle_queue.takeFirst();
                    int ix_face_quadrilatere = face_et_vertex_quadrilatere.first;
                    int ix_vertex_quadrilatere = face_et_vertex_quadrilatere.second;
                    if (!infinitEdge(ix_face_quadrilatere, ix_vertex_quadrilatere))
                    {
                        if (!isDelaunay(ix_face_quadrilatere, ix_vertex_quadrilatere))
                        {
                            atraiter.push_back(face_et_vertex_quadrilatere);
                        };
                    };
                };
            };
        };
    };
}


// ----------------------------------------------
// --------Crust Algorithm-------- MON CODE COMMENCE A PARTIR DE LA
// ----------------------------------------------


/*
 * Parse a file with vertices and store the data into vertices.

 * params:
  - file_name the path of the file.
*/
void Mesh::parseTriFile(const char file_name[])
{
    // Lecture du fichier, et stockage dans vertices
    FILE *pFile;
    pFile = fopen(file_name, "r");

    // Logs the success of the opening operation
    if (pFile != NULL)
    {
        std::cout << "file " << file_name << " opened" << std::endl;
    }
    else
    {
        std::cout << "file " << file_name << " not opened" << std::endl;
    }

    // Initialize vertices and n_edgess
    vertices.clear();
    faces.clear();

    // Read the first line, that gives the number of vertices,
    // the number of faces and the number of edges
    fscanf(pFile, "%d\n", &n_vertices);
    // Logs the result
    std::cout << "nb points: " << n_vertices << std::endl;
    n_faces = 0;

    // Reserve the right size for vertices and faces.
    vertices.reserve(n_vertices);

    // Initialize working variables
    float x, y, z;
    float x_min, x_max, y_min, y_max, z_min, z_max;

    // 1st point :
    fscanf(pFile, "%f %f %f\n", &x, &y, &z); // Stockage de la premiere ligne
    vertices.push_back(Vertex(x, y, z));     // Ajout du point dans le vecteur vertices

    x_min = x;
    x_max = x;
    y_min = y;
    y_max = y;
    z_min = z;
    z_max = z;

    // Tous les autres points :
    for (int ix_vertex = 1; ix_vertex < n_vertices; ix_vertex++)
    {
        fscanf(pFile, "%f %f %f\n", &x, &y, &z); // Stockage de la ligne lue
        vertices.push_back(Vertex(x, y, z));     // Ajout du point dans le vecteur vertices

        if (x < x_min)
        {
            x_min = x;
        };
        if (x > x_max)
        {
            x_max = x;
        };
        if (y < y_min)
        {
            y_min = y;
        };
        if (y > y_max)
        {
            y_max = y;
        };
        if (z < z_min)
        {
            z_min = z;
        };
        if (z > z_max)
        {
            z_max = z;
        };
    }

    // On remet Ã  jour la liste en centrant cette fois les coordonnees
    float x_middle = (x_max + x_min) / 2;
    float y_middle = (y_max + y_min) / 2;
    float z_middle = (z_max + z_min) / 2;

    // Mise Ã  jour de tous les vertex en les remplaÃ§ant par leur Ã©quivalent centrÃ© en 0,0,0.
    for (int ix_vertex = 0; ix_vertex < n_vertices; ix_vertex++)
    {
        Vertex vertex = vertices[ix_vertex];
        vertices[ix_vertex] = Vertex(vertex.x() - x_middle, vertex.y() - y_middle, vertex.z() - z_middle);
    }

    fclose(pFile);
    std::cout << "end of reading" << std::endl;
}

/*
 * Draw delaunay triangulation from vertices without triangles.

 */
void Mesh::triangulationFromVertices()
{
    // On ne prend pas en compte la dimension z
    // On part d'un maillage sans triangle, seulement des points

    // Tout d'abord, on crée le rectangle qui va contenir tous les points

    float x_min = vertices[0].x();
    float x_max = vertices[0].x();
    float y_min = vertices[0].y();
    float y_max = vertices[0].y();

    for (int ix_vertex = 1; ix_vertex < n_vertices; ix_vertex++)
    {
        if (vertices[ix_vertex].x() < x_min)
        {
            x_min = vertices[ix_vertex].x();
        };
        if (vertices[ix_vertex].x() > x_max)
        {
            x_max = vertices[ix_vertex].x();
        };
        if (vertices[ix_vertex].y() < y_min)
        {
            y_min = vertices[ix_vertex].y();
        };
        if (vertices[ix_vertex].y() > y_max)
        {
            y_max = vertices[ix_vertex].y();
        };
    };

    float x_min2 = x_min - (x_max - x_min) / 2;
    float x_max2 = x_max + (x_max - x_min) / 2;
    float y_min2 = y_min - (y_max - y_min) / 2;
    float y_max2 = y_max + (y_max - y_min) / 2;

    vertices.push_back(Vertex(x_min2, y_min2, 0)); // bas gauche
    vertices.push_back(Vertex(x_min2, y_max2, 0)); // haut gauche
    vertices.push_back(Vertex(x_max2, y_min2, 0)); // bas droit
    vertices.push_back(Vertex(x_max2, y_max2, 0)); // haut droit

    n_vertices += 4;

    Vertex &bg = vertices[n_vertices - 4];
    Vertex &hg = vertices[n_vertices - 3];
    Vertex &bd = vertices[n_vertices - 2];
    Vertex &hd = vertices[n_vertices - 1];

    faces.push_back(Face(n_vertices - 4, n_vertices - 2, n_vertices - 3)); // Triangle bas gauche
    faces.push_back(Face(n_vertices - 1, n_vertices - 3, n_vertices - 2)); // Triangle haut droit

    n_faces += 2;

    Face &triangle_bg = faces[n_faces - 2];
    Face &triangle_hd = faces[n_faces - 1];

    bg.ix_incident_face = n_faces - 2;
    hg.ix_incident_face = n_faces - 2;
    bd.ix_incident_face = n_faces - 1;
    hd.ix_incident_face = n_faces - 1;

    // do not draw these points
    bg.is_a_to_draw_point = false;
    hg.is_a_to_draw_point = false;
    bd.is_a_to_draw_point = false;
    hd.is_a_to_draw_point = false;

    triangle_bg.adjacent_faces[0] = n_faces - 1;
    triangle_bg.adjacent_faces[1] = -1;
    triangle_bg.adjacent_faces[2] = -1;

    triangle_hd.adjacent_faces[0] = n_faces - 2;
    triangle_hd.adjacent_faces[1] = -1;
    triangle_hd.adjacent_faces[2] = -1;

    for (int ix_vertex = 0; ix_vertex < (n_vertices - 4); ix_vertex++)
    {
        for (int ix_face = 0; ix_face < n_faces; ix_face++)
        {
            if (inTriangleTest(faces[ix_face], vertices[ix_vertex]) > 0)
            {
                insertionTriangle(ix_vertex, ix_face);
                lawsonAroundVertex(ix_vertex);
                break;
            }
            else if (inTriangleTest(faces[ix_face], vertices[ix_vertex]) == 0)
            {
                insertionInEdge(ix_face, ix_vertex); // pb ici probablement
                lawsonAroundVertex(ix_vertex);
                break;
            };
        };
    };
}


/*
 * Return the Voronoi point coordinates associated to a Delaunay Triangle
 *
*/
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
    float sum_abg = alpha + beta + gamma;

    float x_voronoi = (alpha*A.x() + beta*B.x() + gamma*C.x())/sum_abg;
    coordinates.push_back(x_voronoi);
    float y_voronoi = (alpha*A.y() + beta*B.y() + gamma*C.y())/sum_abg;
    coordinates.push_back(y_voronoi);
    float z_voronoi = (alpha*A.z() + beta*B.z() + gamma*C.z())/sum_abg;
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
            std::cout << "nb points: " << n_faces << " vorCoord0: " << voronoiCoordinates[0] << std::endl;
            std::cout << "nb points: " << n_vertices << " vorCoord1: " << voronoiCoordinates[1] << std::endl;
            std::cout << "nb points: " << n_vertices << " vorCoord2: " << voronoiCoordinates[2] << std::endl;

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
                insertionTriangle(n_vertices-nb_voronoi_vertices+i_voronoix_vertex, ix_face);
                lawsonAroundVertex(n_vertices-nb_voronoi_vertices+i_voronoix_vertex);
                break;
            }
            else if (inTriangleTest(faces[ix_face], vertices[n_vertices-nb_voronoi_vertices+i_voronoix_vertex]) == 0)
            {
                insertionInEdge(n_vertices-nb_voronoi_vertices+i_voronoix_vertex, ix_face); // pb ici probablement
                lawsonAroundVertex(n_vertices-nb_voronoi_vertices+i_voronoix_vertex);
                break;
            };
        };
    };
}


// To Do next : CRUST --> une fonction qui calcule les centres de voronoi de chaque triangle,
// puis qui refait la triangulation et puis qui dessine que les arÃªtes composÃ©es seulement des sommets initiaux
