#include "mesh.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

// ----------------------------------------------
// -------------------Functions------------------
// ----------------------------------------------

/**
 * Compute cotan between two vectors (described with vertices).
 *
 * @param vertex1, vertex2 the vertex describing the vectors.
 * @returns result of cotan operation.
 */
double cotan(const Vertex vertex1, const Vertex vertex2)
{
    double cosalpha = (vertex1 * vertex2) / (vertex1.getNorm() * vertex2.getNorm());
    double result = cosalpha / sqrt(1 - cosalpha * cosalpha);
    return result;
}

//-----------------------------------------------
// --------------------Vertices------------------
// ----------------------------------------------

/**
 * Default constructor. Set the coordinates to 0.
 */
Vertex::Vertex()
{
    _x = 0;
    _y = 0;
    _z = 0;
}

/**
 * Constructor with explicit data.
 *
 * @param x_, y_, z_ the coordinates of the vertex.
 */
Vertex::Vertex(double x_, double y_, double z_)
{
    _x = x_;
    _y = y_;
    _z = z_;
}

/**
 * Getter of x.
 *
 * @returns x coordinate of the vertex.
 */
double Vertex::x() const
{
    return _x;
}

/**
 * Getter of y.
 *
 * @returns y coordinate of the vertex.
 */
double Vertex::y() const
{
    return _y;
}

/**
 * Getter of z.
 *
 * @returns z coordinate of the vertex.
 */
double Vertex::z() const
{
    return _z;
}

/**
 * Normalize the vector.
 *
 * @returns the normalized vector.
 */
Vertex Vertex::getNormalized() const
{
    double norm = getNorm();
    return Vertex(_x / norm, _y / norm, _z / norm);
}

/**
 * Operator !=. true if the other vertex has the same coordinates, false otherwise.
 *
 * @param other_vertex the vertex we want to compare it to.
 * @return result the boolean result.
 */
bool Vertex::operator!=(const Vertex other_vertex) const
{
    return ((x() != other_vertex.x()) | (y() != other_vertex.y()) | (z() != other_vertex.z()));
}

/**
 * Cross product (or vector product) with another vector.
 *
 * @param other_vertex the vector with which we want to do the cross product.
 * @return result the resulting vertex.
 */
Vertex Vertex::cross(const Vertex &other_vertex) const
{
    double result_x = _y * other_vertex.z() - _z * other_vertex.y();
    double result_y = _z * other_vertex.x() - _x * other_vertex.z();
    double result_z = _x * other_vertex.y() - _y * other_vertex.x();
    return Vertex(result_x, result_y, result_z);
}

// ----------------------------------------------
// ------------------------Face------------------
// ----------------------------------------------

/**
 * Default constructor. 
 * 
 * The face is described with three vertices. These vertices are
 * store into face object with their index. By default the three indexes are 0.
 */
Face::Face()
{
    i_vertex[0] = 0;
    i_vertex[1] = 0;
    i_vertex[2] = 0;
}

/**
 * Constructor with explicit data.
 * 
 * The face is described with three vertices. These vertices are
 * store into face object with their index
 * 
 * @param i_vertex0_, i_vertex1_, i_vertex2_ the indexes of the three vertex
 *     describing the triangle.
 */
Face::Face(int i_vertex0_, int i_vertex1_, int i_vertex2_)
{
    i_vertex[0] = i_vertex0_;
    i_vertex[1] = i_vertex1_;
    i_vertex[2] = i_vertex2_;
}

// ----------------------------------------------
// ----------------------Mesh--------------------
// ----------------------------------------------

/**
 * Default constructor.
 */
Mesh::Mesh()
{
}

/**
 * Parse an off file and store the data into verticesTab and facesTab.
 * 
 * @param file_name the path of the file.
 */
void Mesh::parseFile(const char file_name[])
{
    // Lecture du fichier, et stockage dans verticesTab et facesTab
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

    // Initilize verticesTab and nb_edges
    verticesTab.clear();
    int nb_edge = 0;

    // Read the first line, that gives the number of vertices,
    // the number of faces and the number of edges
    fscanf(pFile, "%d %d %d\n", &nb_vertex, &nb_faces, &nb_edge);
    // Logs the result
    std::cout << "nb points: " << nb_vertex << ", nb faces: " << nb_faces << std::endl;

    // Reserve the right size for verticesTab and facesTab.
    verticesTab.reserve(nb_vertex);
    facesTab.reserve(nb_faces);

    // Initialize working variables
    double x, y, z;
    double x_min, x_max, y_min, y_max, z_min, z_max;

    // 1st point :
    fscanf(pFile, "%lf %lf %lf\n", &x, &y, &z); // Stockage de la premiere ligne
    verticesTab.push_back(Vertex(x, y, z));     // Ajout du point dans le vecteur verticesTab

    x_min = x;
    x_max = x;
    y_min = y;
    y_max = y;
    z_min = z;
    z_max = z;

    // Tous les autres points :
    for (int i_vertex = 1; i_vertex < nb_vertex; i_vertex++)
    {
        fscanf(pFile, "%lf %lf %lf\n", &x, &y, &z); // Stockage de la ligne lue
        verticesTab.push_back(Vertex(x, y, z));     // Ajout du point dans le vecteur verticesTab

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

    // On remet à jour la liste en centrant cette fois les coordonnees
    double x_middle = (x_max + x_min) / 2;
    double y_middle = (y_max + y_min) / 2;
    double z_middle = (z_max + z_min) / 2;

    // Mise à jour de tous les vertex en les remplaçant par leur équivalent centré en 0,0,0.
    for (int i_vertex = 0; i_vertex < nb_vertex; i_vertex++)
    {
        Vertex vertex = verticesTab[i_vertex];
        verticesTab[i_vertex] = Vertex(vertex.x() - x_middle, vertex.y() - y_middle, vertex.z() - z_middle);
    }

    // Stocke toutes les faces dans facesTab
    int n_face, i_vertex0, i_vertex1, i_vertex2;
    for (int i_face = 0; i_face < nb_faces; i_face++)
    {
        fscanf(pFile, "%d %d %d %d\n", &n_face, &i_vertex0, &i_vertex1, &i_vertex2);
        facesTab.push_back(Face(i_vertex0, i_vertex1, i_vertex2));
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

    // vert_done est un liste de boolean. Si vert_done[5] == true, ça signifie que le vertex 5 à été traité
    // Initilisation de vert_done
    bool vert_done[nb_vertex];
    for (int i_vertex = 0; i_vertex < nb_vertex; i_vertex++)
    {
        vert_done[i_vertex] = false;
    }

    // Map est un dict. ici : {(int, int) : (int, int), ...}
    std::map<std::pair<int, int>, std::pair<int, int>> myMap;

    for (int i_face = 0; i_face < nb_faces; i_face++)
    {
        // On se balade dans toutes les faces définie, et on récupère les indices des trois points(=vertex)
        Face *face = &facesTab[i_face]; // Un pointeur qui pointe vers la i_eme face

        // On regarde si on a déja traité les vertex(=points). Si non, on dit qu'on les a traité
        // puis on ajoute l'indice de la face comme face incidente au vertex(=point)
        for (int i_vertex_in_triangle = 0; i_vertex_in_triangle < 3; i_vertex_in_triangle++)
        {
            int i_vertex = face->i_vertex[i_vertex_in_triangle];
            if (!vert_done[i_vertex])
            {
                verticesTab[i_vertex].i_incident_face = i_face; // On ajoute à l'objet vertex une face inscidente
                vert_done[i_vertex] = true;                     // On le considère traité
            };
        }

        // setting adjacent faces on each edge:
        std::pair<int, int> edge;                            //  Arrêtes
        std::pair<int, int> i_other_face_and_i_other_vertex; //  {face id whose edge was'nt yet register in the map,  its opposed vertex id (0,1 or 2) }

        // Order vertex
        for (int i_vertex_in_triangle = 0; i_vertex_in_triangle < 3; i_vertex_in_triangle++)
        {
            if (face->i_vertex[(i_vertex_in_triangle + 1) % 3] < face->i_vertex[(i_vertex_in_triangle + 2) % 3])
            {
                edge = {face->i_vertex[(i_vertex_in_triangle + 1) % 3], face->i_vertex[(i_vertex_in_triangle + 2) % 3]};
            }
            else
            {
                edge = {face->i_vertex[(i_vertex_in_triangle + 2) % 3], face->i_vertex[(i_vertex_in_triangle + 1) % 3]};
            }
            // edge est la pair d'indice dans l'ordre croissant de l'arrète

            if (myMap.find(edge) == myMap.end())
            { // Si l'arrète n'est pas dans la map
                // On l'ajoute
                myMap[edge] = {i_face, i_vertex_in_triangle}; // (i_vertex0, i_vertex1) : (i_face, i_vertex) (c'est le vertex opposé)
            }
            else
            {                                                     // L'arrête appartient déja à un autre triangle de la map
                i_other_face_and_i_other_vertex = myMap.at(edge); //(i_face de l'autre triangle, i_face dans l'autre triangle (0,1,2))
                int i_other_face = i_other_face_and_i_other_vertex.first;
                int i_other_vertex = i_other_face_and_i_other_vertex.second;

                face->adjacent_faces[i_vertex_in_triangle] = i_other_face; // On stocke l'indice de la face adjacente opposé au vertex(=point) n°2
                // On fait la même chose à l'autre triangle
                facesTab[i_other_face].adjacent_faces[i_other_vertex] = i_face;
            }
        }
    }
    std::cout << "end sewing" << std::endl;
}
/*
Vertex *Mesh::vertices_begin()
{
    return verticesTab.begin();
}

Vertex *Mesh::vertices_past_the_end()
{
    return verticesTab.end();
}

Face *Mesh::faces_begin()
{
    return facesTab.begin();
}

Face *Mesh::faces_past_the_end()
{
    return facesTab.end();
}

Circulator_on_faces Mesh::incident_faces(int i_vertex)
{
    return Circulator_on_faces(this, i_vertex);
}

Circulator_on_vertices Mesh::adjacent_vertices(int i_vertex)
{
    return Circulator_on_vertices(this, i_vertex);
}
*/
/*
void Mesh::computeLaplacian()
{
    std::cout << "begin compute laplacian" << std::endl;
    for (int i_vertex = 0; i_vertex < nb_vertex; i_vertex++)
    {
        double area = 0;
        Vertex laplacian(0, 0, 0);

        Circulator_on_vertices circulator_on_vertices(this, i_vertex);

        Vertex vertexO = verticesTab[i_vertex];
        Vertex vertexA;
        Vertex vertexB;
        Vertex vertexC;

        Vertex vectorAtoO;
        Vertex vectorAtoB;
        Vertex vectorCtoO;
        Vertex vectorCtoB;
        Vertex vectorOtoB;

        vertexA = *circulator_on_vertices;
        ++circulator_on_vertices;
        vertexB = *circulator_on_vertices;
        ++circulator_on_vertices;
        vertexC = *circulator_on_vertices;

        do
        {
            vectorAtoO = vertexO - vertexA;
            vectorAtoB = vertexB - vertexA;
            vectorCtoO = vertexO - vertexC;
            vectorCtoB = vertexB - vertexC;
            vectorOtoB = vertexB - vertexO;

            laplacian = laplacian + vectorOtoB * (cotan(vectorAtoO, vectorAtoB) + cotan(vectorCtoO, vectorCtoB));
            area += vectorAtoO.cross(vectorAtoB).getNorm() / 2 / 3;

            vertexA = vertexB;
            vertexB = vertexC;
            ++circulator_on_vertices;
            vertexC = *circulator_on_vertices;

        } while (vertexA != *circulator_on_vertices.begin());

        laplacianTab.push_back(laplacian / (2 * area));
    }
    std::cout << "laplacian computed" << std::endl;
}
*/

/**
 * Tests if a vertex is within the circumscribed circle of the face.
 * 
 * @param face for the test.
 * @param P the vertex for which we want to know if it's in the circumscribed circle.
 * @return double = 0 if P is on the circumscribed circle.
 *                > 0 if P is inside the circumscribed circle.
 *                < 0 if P is outside the circumscribed circle.
 */
double Mesh::vertexInCircumscribingCircle(Face face, Vertex P)
{
    // On analyse seulement les dimensions x, y
    // On regarde en projetant sur le paraboloïde z = x*x + y*y

    Vertex A = verticesTab[face.i_vertex[0]];
    Vertex B = verticesTab[face.i_vertex[1]];
    Vertex C = verticesTab[face.i_vertex[2]];

    Vertex A_hyper = Vertex(A.x(), A.y(), A.x() * A.x() + A.y() * A.y());
    Vertex B_hyper = Vertex(B.x(), B.y(), B.x() * B.x() + B.y() * B.y());
    Vertex C_hyper = Vertex(C.x(), C.y(), C.x() * C.x() + C.y() * C.y());
    Vertex P_hyper = Vertex(P.x(), P.y(), P.x() * P.x() + P.y() * P.y());

    Vertex AB_hyper = B_hyper - A_hyper;
    Vertex AC_hyper = C_hyper - A_hyper;
    Vertex AP_hyper = P_hyper - A_hyper;

    Vertex pdtVect = AB_hyper.cross(AC_hyper);
    double pdtScal = -(pdtVect * AP_hyper);

    return pdtScal;
}

/**
 * Tests whether an edge is Delaunay.
 * 
 * @param i_face1 a face that contains the edge.
 * @param i_vertex_oppose_1 the index of the vertex opposite the edge
 * @returns true if Delaunay, false otherwise.
 */
bool Mesh::isDelaunay(int i_face1, int i_vertex_oppose_1)
{

    Face face1 = facesTab[i_face1];
    Face face2;
    if (face1.i_vertex[0] == i_vertex_oppose_1)
    {
        face2 = facesTab[face1.adjacent_faces[0]];
    }
    else if (face1.i_vertex[1] == i_vertex_oppose_1)
    {
        face2 = facesTab[face1.adjacent_faces[1]];
    }
    else
    {
        face2 = facesTab[face1.adjacent_faces[2]];
    };
    int i_vertex_oppose_2;
    if (face2.adjacent_faces[0] == i_face1)
    {
        i_vertex_oppose_2 = face2.i_vertex[0];
    }
    else if (face2.adjacent_faces[1] == i_face1)
    {
        i_vertex_oppose_2 = face2.i_vertex[1];
    }
    else
    {
        i_vertex_oppose_2 = face2.i_vertex[2];
    };

    double test1 = vertexInCircumscribingCircle(face1, verticesTab[i_vertex_oppose_2]);
    double test2 = vertexInCircumscribingCircle(face2, verticesTab[i_vertex_oppose_1]);

    // On demande à ce que les points soient sur la bordure ou strictement à l'extérieur du cercle
    return ((test1 <= 0) && (test2 <= 0));
}

/**
 * Flips the edge.
 * It is therefore assumed that the edge is not on the contour.
 * 
 * @param i_face1 the index of the first triangle to which the edge belongs.
 * @param vertex_oppose_1_initial the opposite vertex, belonging to the
 *     second triangle whose edge is common.
 * @returns the 2 edges that frame the flipped edge, and not being in
 *     the same initial triangle as vertex_oppose_1_initial.
 */
QList<std::pair<int, int>> Mesh::flipEdge(int i_face1, int i_vertex_oppose_1_initial)
{
    Face &face1 = facesTab[i_face1];

    int i_vertex_a = i_vertex_oppose_1_initial;
    int i_vertex_b = -1;
    int i_vertex_c = -1;
    int i_vertex_d = -1;
    int i_face2 = -1;
    int i_face3 = -1;
    int i_face4 = -1;
    int i_face5 = -1;
    int i_face6 = -1;

    if (i_vertex_oppose_1_initial == face1.i_vertex[0])
    {
        i_face2 = face1.adjacent_faces[0];
        i_vertex_d = face1.i_vertex[1];
        i_vertex_b = face1.i_vertex[2];
        i_face4 = face1.adjacent_faces[1];
        i_face3 = face1.adjacent_faces[2];
    }
    else if (i_vertex_oppose_1_initial == face1.i_vertex[1])
    {
        i_face2 = face1.adjacent_faces[1];
        i_vertex_d = face1.i_vertex[2];
        i_vertex_b = face1.i_vertex[0];
        i_face4 = face1.adjacent_faces[2];
        i_face3 = face1.adjacent_faces[0];
    }
    else
    {
        i_face2 = face1.adjacent_faces[2];
        i_vertex_d = face1.i_vertex[0];
        i_vertex_b = face1.i_vertex[1];
        i_face4 = face1.adjacent_faces[0];
        i_face3 = face1.adjacent_faces[1];
    };
    Face &face2 = facesTab[i_face2];

    if (i_face1 == face2.adjacent_faces[0])
    {

        i_vertex_c = face2.i_vertex[0];
        i_face6 = face2.adjacent_faces[1];
        i_face5 = face2.adjacent_faces[2];
    }
    else if (i_face1 == face2.adjacent_faces[1])
    {

        i_vertex_c = face2.i_vertex[1];
        i_face6 = face2.adjacent_faces[2];
        i_face5 = face2.adjacent_faces[0];
    }
    else
    {

        i_vertex_c = face2.i_vertex[2];
        i_face6 = face2.adjacent_faces[0];
        i_face5 = face2.adjacent_faces[1];
    };

    // Actualisation des donnees apres le flip

    // Faces 1 et 2

    face1.i_vertex[0] = i_vertex_a;
    face1.i_vertex[1] = i_vertex_d;
    face1.i_vertex[2] = i_vertex_c;

    face1.adjacent_faces[0] = i_face6;
    face1.adjacent_faces[1] = i_face2;
    face1.adjacent_faces[2] = i_face3;

    face2.i_vertex[0] = i_vertex_c;
    face2.i_vertex[1] = i_vertex_b;
    face2.i_vertex[2] = i_vertex_a;

    face2.adjacent_faces[0] = i_face4;
    face2.adjacent_faces[1] = i_face1;
    face2.adjacent_faces[2] = i_face5;

    // Faces 3 4 5 6
    // On actualise la face adj partant du vertex à l'oppose de notre quadrilatere

    for (int i = 0; i < 3; i++)
    {
        if (i_face4 >= 0 && facesTab[i_face4].adjacent_faces[i] == i_face1)
        {
            facesTab[i_face4].adjacent_faces[i] = i_face2;
        };
        if (i_face6 >= 0 && facesTab[i_face6].adjacent_faces[i] == i_face2)
        {
            facesTab[i_face6].adjacent_faces[i] = i_face1;
        };
    };

    // Vertices

    verticesTab[i_vertex_a].i_incident_face = i_face1;
    verticesTab[i_vertex_d].i_incident_face = i_face1;
    verticesTab[i_vertex_c].i_incident_face = i_face2;
    verticesTab[i_vertex_b].i_incident_face = i_face2;

    QList<std::pair<int, int>> aretes_quadrilatere;
    aretes_quadrilatere.append({i_face1, i_vertex_a}); // arete f6
    aretes_quadrilatere.append({i_face2, i_vertex_a}); // arete f5

    return aretes_quadrilatere;
}

/**
 * Tests if A, B and C are oriented in the trigonometric direction.
 * 
 * @param A, B, C vertices for the orientation test.
 * @returns double which is positive if A B C are oriented in the
 *     trigonometric direction, zero if aligned, negative otherwise.
*/
double Mesh::orientationTest(Vertex A, Vertex B, Vertex C)
{
    // On travaille dans le plan z = 0
    Vertex AB = B - A;
    Vertex BC = C - B;

    Vertex AB_proj = Vertex(AB.x(), AB.y(), 0);
    Vertex BC_proj = Vertex(BC.x(), BC.y(), 0);

    Vertex vector = AB_proj.cross(BC_proj);
    double value = vector.z(); //*nb_vertex;
    return value;
}

/**
 * Tests if a vertex is inside a face.
 * 
 * @param face used for the test.
 * @param vertex we want to know if it's inside the face. 
 * @returns double which is positive if the vertex is inside the
 *     triangle, 0 if it is on an edge, negative if it is outside.
*/
double Mesh::inTriangleTest(Face face, Vertex vertex)
{
    // On travaille dans le plan z = 0
    Vertex A = verticesTab[face.i_vertex[0]];
    Vertex B = verticesTab[face.i_vertex[1]];
    Vertex C = verticesTab[face.i_vertex[2]];

    double test1 = orientationTest(A, B, vertex);
    double test2 = orientationTest(B, C, vertex);
    double test3 = orientationTest(C, A, vertex);

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
        // On est forcément en dehors du triangle
        // On renvoi une valeur négative
        return -1;
    };
}

/**
 * Insert a point in a face.
 * 
 * @param i_P index of the point to be inserted.
 * @param i_face index of the face in which we want to insert the vertex.
 */
void Mesh::insertionTriangle(int i_P, int i_face)
{
    Face old_face = facesTab[i_face];

    // Creer nouveaux triangles

    int i_A = old_face.i_vertex[0];
    int i_B = old_face.i_vertex[1];
    int i_C = old_face.i_vertex[2];

    // On prend les faces autour du triangle si les voisins existent
    // Si l'indice  = -1, c'est qu'il n'y a pas de face voisine (cf plus tard)

    int i_face_A = old_face.adjacent_faces[0];
    int i_face_B = old_face.adjacent_faces[1];
    int i_face_C = old_face.adjacent_faces[2];

    Vertex &A = verticesTab[i_A];
    Vertex &B = verticesTab[i_B];
    Vertex &C = verticesTab[i_C];
    Vertex &P = verticesTab[i_P];

    // Insérer ces triangles dans le tableau faces et supprimer l'ancien

    facesTab[i_face] = Face(i_A, i_B, i_P); // ABP
    facesTab.append(Face(i_B, i_C, i_P));   // BCP
    facesTab.append(Face(i_C, i_A, i_P));   // CAP

    int i_ABP = i_face;
    int i_BCP = nb_faces;
    int i_CAP = nb_faces + 1;

    Face &ABP = facesTab[i_ABP];
    Face &BCP = facesTab[i_BCP];
    Face &CAP = facesTab[i_CAP];

    // Actualiser paramètres du Mesh

    nb_faces += 2;

    // Ajouter les adjacences des nouveaux triangles

    ABP.adjacent_faces[0] = i_BCP;
    ABP.adjacent_faces[1] = i_CAP;
    ABP.adjacent_faces[2] = i_face_C;

    BCP.adjacent_faces[0] = i_CAP;
    BCP.adjacent_faces[1] = i_ABP;
    BCP.adjacent_faces[2] = i_face_A;

    CAP.adjacent_faces[0] = i_ABP;
    CAP.adjacent_faces[1] = i_BCP;
    CAP.adjacent_faces[2] = i_face_B;

    // Adapter les adjacences des triangles autour

    for (int i = 0; i < 3; i++)
    {
        if (i_face_A >= 0 && facesTab[i_face_A].adjacent_faces[i] == i_face)
        {
            facesTab[i_face_A].adjacent_faces[i] = i_BCP;
        };
        if (i_face_B >= 0 && facesTab[i_face_B].adjacent_faces[i] == i_face)
        {
            facesTab[i_face_B].adjacent_faces[i] = i_CAP;
        };
        if (i_face_C >= 0 && facesTab[i_face_C].adjacent_faces[i] == i_face)
        {
            facesTab[i_face_C].adjacent_faces[i] = i_ABP;
        };
    }

    // Modifier les i_incident_face des Vertices

    A.i_incident_face = i_ABP;
    B.i_incident_face = i_BCP;
    C.i_incident_face = i_CAP;
    P.i_incident_face = i_ABP;
}

/**
 * Inserts a vertex on an edge
 * 
 * @param i_face1 index of the face containing the edge.
 * @param i_P index of the vertex to be inserted.
 */
void Mesh::insertionInArete(int i_face1, int i_P)
{
    Face face1 = facesTab[i_face1];

    int i_A = -1;
    int i_B = -1;
    int i_C = -1;
    int i_D = -1;
    int i_face2 = -1;
    int i_face3 = -1;
    int i_face4 = -1;
    int i_face5 = -1;
    int i_face6 = -1;

    if (orientationTest(verticesTab[face1.i_vertex[0]], verticesTab[face1.i_vertex[1]], verticesTab[i_P]) == 0)
    {
        i_A = face1.i_vertex[2];
        i_B = face1.i_vertex[0];
        i_D = face1.i_vertex[1];
        i_face2 = face1.adjacent_faces[2];
        i_face6 = face1.adjacent_faces[0];
        i_face3 = face1.adjacent_faces[1];
    }
    else if (orientationTest(verticesTab[face1.i_vertex[1]], verticesTab[face1.i_vertex[2]], verticesTab[i_P]) == 0)
    {
        i_A = face1.i_vertex[0];
        i_B = face1.i_vertex[1];
        i_D = face1.i_vertex[2];
        i_face2 = face1.adjacent_faces[0];
        i_face6 = face1.adjacent_faces[1];
        i_face3 = face1.adjacent_faces[2];
    }
    else if (orientationTest(verticesTab[face1.i_vertex[2]], verticesTab[face1.i_vertex[0]], verticesTab[i_P]) == 0)
    {
        i_A = face1.i_vertex[1];
        i_B = face1.i_vertex[2];
        i_D = face1.i_vertex[0];
        i_face2 = face1.adjacent_faces[1];
        i_face6 = face1.adjacent_faces[2];
        i_face3 = face1.adjacent_faces[0];
    };

    for (int i = 0; i < 3; i++)
    {
        if (facesTab[i_face2].adjacent_faces[i] == i_face1)
        {
            i_C = facesTab[i_face2].i_vertex[i];
            i_face4 = facesTab[i_face2].adjacent_faces[(i + 1) % 3];
            i_face5 = facesTab[i_face2].adjacent_faces[(i + 2) % 3];
        };
    };

    Vertex &A = verticesTab[i_A];
    Vertex &B = verticesTab[i_B];
    Vertex &C = verticesTab[i_C];
    Vertex &D = verticesTab[i_D];
    Vertex &P = verticesTab[i_P];

    // Insérer ces triangles dans le tableau faces et supprimer l'ancien

    facesTab[i_face1] = Face(i_D, i_A, i_P); // DAP
    facesTab[i_face2] = Face(i_C, i_D, i_P); // CDP

    facesTab.append(Face(i_A, i_B, i_P)); // ABP
    facesTab.append(Face(i_B, i_C, i_P)); // BCP

    int i_DAP = i_face1;
    int i_CDP = i_face2;
    int i_ABP = nb_faces;
    int i_BCP = nb_faces + 1;

    Face &DAP = facesTab[i_DAP];
    Face &CDP = facesTab[i_CDP];
    Face &ABP = facesTab[i_ABP];
    Face &BCP = facesTab[i_BCP];

    // Actualiser paramètres du Mesh

    nb_faces += 2;

    // Ajouter les adjacences des nouveaux triangles

    DAP.adjacent_faces[0] = i_ABP;
    DAP.adjacent_faces[1] = i_CDP;
    DAP.adjacent_faces[2] = i_face6;

    CDP.adjacent_faces[0] = i_DAP;
    CDP.adjacent_faces[1] = i_BCP;
    CDP.adjacent_faces[2] = i_face5;

    ABP.adjacent_faces[0] = i_BCP;
    ABP.adjacent_faces[1] = i_DAP;
    ABP.adjacent_faces[2] = i_face3;

    BCP.adjacent_faces[0] = i_CDP;
    BCP.adjacent_faces[1] = i_ABP;
    BCP.adjacent_faces[2] = i_face4;

    // Adapter les adjacences des triangles autour

    for (int i = 0; i < 3; i++)
    {
        if (i_face3 >= 0 && facesTab[i_face3].adjacent_faces[i] == i_face1)
        {
            facesTab[i_face3].adjacent_faces[i] = i_ABP;
        };
        if (i_face4 >= 0 && facesTab[i_face4].adjacent_faces[i] == i_face2)
        {
            facesTab[i_face4].adjacent_faces[i] = i_BCP;
        };
        if (i_face5 >= 0 && facesTab[i_face5].adjacent_faces[i] == i_face2)
        {
            facesTab[i_face5].adjacent_faces[i] = i_CDP;
        };
        if (i_face6 >= 0 && facesTab[i_face6].adjacent_faces[i] == i_face1)
        {
            facesTab[i_face6].adjacent_faces[i] = i_DAP;
        };
    }

    // Modifier les i_incident_face des Vertices

    A.i_incident_face = i_ABP;
    B.i_incident_face = i_BCP;
    C.i_incident_face = i_CDP;
    D.i_incident_face = i_DAP;
    P.i_incident_face = i_ABP;
}

/**
 * A naive mesh is created from points by inserting them successively.
 */
void Mesh::naiveInsertion()
{
    // On ne prend pas en compte la dimension z
    // On part d'un maillage sans triangle, seulement des points

    // Tout d'abord, on crée le rectangle qui va contenir tous les points

    double x_min = verticesTab[0].x();
    double x_max = verticesTab[0].x();
    double y_min = verticesTab[0].y();
    double y_max = verticesTab[0].y();

    for (int i_vertex = 1; i_vertex < nb_vertex; i_vertex++)
    {
        if (verticesTab[i_vertex].x() < x_min)
        {
            x_min = verticesTab[i_vertex].x();
        };
        if (verticesTab[i_vertex].x() > x_max)
        {
            x_max = verticesTab[i_vertex].x();
        };
        if (verticesTab[i_vertex].y() < y_min)
        {
            y_min = verticesTab[i_vertex].y();
        };
        if (verticesTab[i_vertex].y() > y_max)
        {
            y_max = verticesTab[i_vertex].y();
        };
    };

    x_min -= (x_max - x_min) / 20;
    x_max += (x_max - x_min) / 20;
    y_min -= (y_max - y_min) / 20;
    y_max += (y_max - y_min) / 20;

    verticesTab.append(Vertex(x_min, y_min, 0)); // bas gauche
    verticesTab.append(Vertex(x_min, y_max, 0)); // haut gauche
    verticesTab.append(Vertex(x_max, y_min, 0)); // bas droit
    verticesTab.append(Vertex(x_max, y_max, 0)); // haut droit

    nb_vertex += 4;

    Vertex &bg = verticesTab[nb_vertex - 4];
    Vertex &hg = verticesTab[nb_vertex - 3];
    Vertex &bd = verticesTab[nb_vertex - 2];
    Vertex &hd = verticesTab[nb_vertex - 1];

    facesTab.append(Face(nb_vertex - 4, nb_vertex - 2, nb_vertex - 3)); // Triangle bas gauche
    facesTab.append(Face(nb_vertex - 1, nb_vertex - 3, nb_vertex - 2)); // Triangle haut droit

    nb_faces += 2;

    Face &triangle_bg = facesTab[nb_faces - 2];
    Face &triangle_hd = facesTab[nb_faces - 1];

    bg.i_incident_face = nb_faces - 2;
    hg.i_incident_face = nb_faces - 2;
    bd.i_incident_face = nb_faces - 1;
    hd.i_incident_face = nb_faces - 1;

    triangle_bg.adjacent_faces[0] = nb_faces - 1;
    triangle_bg.adjacent_faces[1] = -1;
    triangle_bg.adjacent_faces[2] = -1;

    triangle_hd.adjacent_faces[0] = nb_faces - 2;
    triangle_hd.adjacent_faces[1] = -1;
    triangle_hd.adjacent_faces[2] = -1;

    for (int i_vertex = 0; i_vertex < (nb_vertex - 4); i_vertex++)
    {
        for (int i_face = 0; i_face < nb_faces; i_face++)
        {

            if (inTriangleTest(facesTab[i_face], verticesTab[i_vertex]) > 0)
            {
                insertionTriangle(i_vertex, i_face);
                break;
            };
            // Il manque le cas inTriangleTest == 0 (ie sur une arete; cf insertion arete qui marche pas (?) pr le moment)
        };
    };
}

/**
 * Tests whether the edge opposite the vertex is on the contour.
 * 
 * @param i_face index of the face containing the edge.
 * @param i_vertex index of the vertex opposite the edge.
 * @returns true if the edge is on the contour, false otherwise.
 */
bool Mesh::areteEnBordure(int i_face, int i_vertex)
{
    if (facesTab[i_face].i_vertex[0] == i_vertex)
    {
        return facesTab[i_face].adjacent_faces[0] < 0;
    }
    else if (facesTab[i_face].i_vertex[1] == i_vertex)
    {
        return facesTab[i_face].adjacent_faces[1] < 0;
    }
    else if (facesTab[i_face].i_vertex[2] == i_vertex)
    {
        return facesTab[i_face].adjacent_faces[2] < 0;
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
    // On part d'un vertex i_P (dans i_face) qui vient d'être insere, et on fait des flips recursifs pour qu'a la fin il soit insere et tout soit Delaunay
    // On recupere les trois ou quatre aretes autour du P dans une file
    QList<std::pair<int, int>> atraiter;
    for (int i_face = 0; i_face < nb_faces; i_face++)
    {
        for (int i = 0; i < 3; i++)
        {
            if (facesTab[i_face].i_vertex[i] == i_P)
            {
                atraiter.append({i_face, i_P});
            }
        }
    }
    // On lance la boucle while et on remplie et traite la file
    while (!atraiter.isEmpty())
    {
        std::pair<int, int> face_et_vertex = atraiter.takeFirst();
        int i_face = face_et_vertex.first;
        int i_vertex = face_et_vertex.second;

        if (!areteEnBordure(i_face, i_vertex))
        {
            if (!isDelaunay(i_face, i_vertex))
            {
                QList<std::pair<int, int>> nouvelle_queue = flipEdge(i_face, i_vertex); // On fait le flip et récupère les arete à retester
                while (!nouvelle_queue.isEmpty())
                {

                    std::pair<int, int> face_et_vertex_quadrilatere = nouvelle_queue.takeFirst();
                    int i_face_quadrilatere = face_et_vertex_quadrilatere.first;
                    int i_vertex_quadrilatere = face_et_vertex_quadrilatere.second;
                    if (!areteEnBordure(i_face_quadrilatere, i_vertex_quadrilatere))
                    {
                        if (!isDelaunay(i_face_quadrilatere, i_vertex_quadrilatere))
                        {
                            atraiter.append(face_et_vertex_quadrilatere);
                        };
                    };
                };
            };
        };
    };
}

/**
 * Similar to naiveInsertion, ensuring Delaunay after each insertion.
 */
void Mesh::naiveInsertionAndLawson()
{
    // On ne prend pas en compte la dimension z
    // On part d'un maillage sans triangle, seulement des points

    // Tout d'abord, on crée le rectangle qui va contenir tous les points

    double x_min = verticesTab[0].x();
    double x_max = verticesTab[0].x();
    double y_min = verticesTab[0].y();
    double y_max = verticesTab[0].y();

    for (int i_vertex = 1; i_vertex < nb_vertex; i_vertex++)
    {
        if (verticesTab[i_vertex].x() < x_min)
        {
            x_min = verticesTab[i_vertex].x();
        };
        if (verticesTab[i_vertex].x() > x_max)
        {
            x_max = verticesTab[i_vertex].x();
        };
        if (verticesTab[i_vertex].y() < y_min)
        {
            y_min = verticesTab[i_vertex].y();
        };
        if (verticesTab[i_vertex].y() > y_max)
        {
            y_max = verticesTab[i_vertex].y();
        };
    };

    double x_min2 = x_min - (x_max - x_min) / 20;
    double x_max2 = x_max + (x_max - x_min) / 20;
    double y_min2 = y_min - (y_max - y_min) / 20;
    double y_max2 = y_max + (y_max - y_min) / 20;

    verticesTab.append(Vertex(x_min2, y_min2, 0)); // bas gauche
    verticesTab.append(Vertex(x_min2, y_max2, 0)); // haut gauche
    verticesTab.append(Vertex(x_max2, y_min2, 0)); // bas droit
    verticesTab.append(Vertex(x_max2, y_max2, 0)); // haut droit

    nb_vertex += 4;

    Vertex &bg = verticesTab[nb_vertex - 4];
    Vertex &hg = verticesTab[nb_vertex - 3];
    Vertex &bd = verticesTab[nb_vertex - 2];
    Vertex &hd = verticesTab[nb_vertex - 1];

    facesTab.append(Face(nb_vertex - 4, nb_vertex - 2, nb_vertex - 3)); // Triangle bas gauche
    facesTab.append(Face(nb_vertex - 1, nb_vertex - 3, nb_vertex - 2)); // Triangle haut droit

    nb_faces += 2;

    Face &triangle_bg = facesTab[nb_faces - 2];
    Face &triangle_hd = facesTab[nb_faces - 1];

    bg.i_incident_face = nb_faces - 2;
    hg.i_incident_face = nb_faces - 2;
    bd.i_incident_face = nb_faces - 1;
    hd.i_incident_face = nb_faces - 1;

    triangle_bg.adjacent_faces[0] = nb_faces - 1;
    triangle_bg.adjacent_faces[1] = -1;
    triangle_bg.adjacent_faces[2] = -1;

    triangle_hd.adjacent_faces[0] = nb_faces - 2;
    triangle_hd.adjacent_faces[1] = -1;
    triangle_hd.adjacent_faces[2] = -1;

    for (int i_vertex = 0; i_vertex < (nb_vertex - 4); i_vertex++)
    {
        for (int i_face = 0; i_face < nb_faces; i_face++)
        {
            if (inTriangleTest(facesTab[i_face], verticesTab[i_vertex]) > 0)
            {
                insertionTriangle(i_vertex, i_face);
                lawsonAroundVertex(i_vertex);
                break;
            }
            else if (inTriangleTest(facesTab[i_face], verticesTab[i_vertex]) == 0)
            {
                insertionInArete(i_face, i_vertex); // pb ici probablement
                lawsonAroundVertex(i_vertex);
                break;
            };
        };
    };
}

// ----------------------------------------------
// --------Useful function for circulator--------
// ----------------------------------------------

/**
 * Finds the index (0, 1 or 2) of a vertex belonging to a face
 * 
 * @param face we're working on
 * @param i_vertex index of the vertex in the mesh
 * @return 0, 1 or 2, the index of the vertex within the face
 */
int find_i_vertex_in_face(Face *face, int i_vertex)
{
    int i_vertex_in_face = 0;
    while (face->i_vertex[i_vertex_in_face] != i_vertex)
    {
        i_vertex_in_face++;
    }
    return i_vertex_in_face;
}

// ----------------------------------------------
// --------Crust Algorithm--------
// ----------------------------------------------


/**
 * Parse a file with vertices and store the data into verticesTab.
 *
 * @param file_name the path of the file.
 */
void Mesh::parseTriFile(const char file_name[])
{
    // Lecture du fichier, et stockage dans verticesTab
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

    // Initilize verticesTab and nb_edges
    verticesTab.clear();

    // Read the first line, that gives the number of vertices,
    // the number of faces and the number of edges
    fscanf(pFile, "%d\n", &nb_vertex);
    // Logs the result
    std::cout << "nb points: " << nb_vertex << std::endl;

    // Reserve the right size for verticesTab and facesTab.
    verticesTab.reserve(nb_vertex);
    int nb_edges = 0;

    // Initialize working variables
    double x, y, z;
    double x_min, x_max, y_min, y_max, z_min, z_max;

    // 1st point :
    fscanf(pFile, "%lf %lf %lf\n", &x, &y, &z); // Stockage de la premiere ligne
    verticesTab.push_back(Vertex(x, y, z));     // Ajout du point dans le vecteur verticesTab

    x_min = x;
    x_max = x;
    y_min = y;
    y_max = y;
    z_min = z;
    z_max = z;

    // Tous les autres points :
    for (int i_vertex = 1; i_vertex < nb_vertex; i_vertex++)
    {
        fscanf(pFile, "%lf %lf %lf\n", &x, &y, &z); // Stockage de la ligne lue
        verticesTab.push_back(Vertex(x, y, z));     // Ajout du point dans le vecteur verticesTab

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

    // On remet à jour la liste en centrant cette fois les coordonnees
    double x_middle = (x_max + x_min) / 2;
    double y_middle = (y_max + y_min) / 2;
    double z_middle = (z_max + z_min) / 2;

    // Mise à jour de tous les vertex en les remplaçant par leur équivalent centré en 0,0,0.
    for (int i_vertex = 0; i_vertex < nb_vertex; i_vertex++)
    {
        Vertex vertex = verticesTab[i_vertex];
        verticesTab[i_vertex] = Vertex(vertex.x() - x_middle, vertex.y() - y_middle, vertex.z() - z_middle);
    }

    fclose(pFile);
    std::cout << "end of reading" << std::endl;
}

/**
 * Draw delaunay triangulation from vertices.
 *
 */
void Mesh::triangulationFromVertices()
{
    // On ne prend pas en compte la dimension z
    // On part d'un maillage sans triangle, seulement des points

    // Tout d'abord, on crée le rectangle qui va contenir tous les points

    double x_min = verticesTab[0].x();
    double x_max = verticesTab[0].x();
    double y_min = verticesTab[0].y();
    double y_max = verticesTab[0].y();

    for (int i_vertex = 1; i_vertex < nb_vertex; i_vertex++)
    {
        if (verticesTab[i_vertex].x() < x_min)
        {
            x_min = verticesTab[i_vertex].x();
        };
        if (verticesTab[i_vertex].x() > x_max)
        {
            x_max = verticesTab[i_vertex].x();
        };
        if (verticesTab[i_vertex].y() < y_min)
        {
            y_min = verticesTab[i_vertex].y();
        };
        if (verticesTab[i_vertex].y() > y_max)
        {
            y_max = verticesTab[i_vertex].y();
        };
    };

    double x_min2 = x_min - (x_max - x_min) / 2;
    double x_max2 = x_max + (x_max - x_min) / 2;
    double y_min2 = y_min - (y_max - y_min) / 2;
    double y_max2 = y_max + (y_max - y_min) / 2;

    verticesTab.append(Vertex(x_min2, y_min2, 0)); // bas gauche
    verticesTab.append(Vertex(x_min2, y_max2, 0)); // haut gauche
    verticesTab.append(Vertex(x_max2, y_min2, 0)); // bas droit
    verticesTab.append(Vertex(x_max2, y_max2, 0)); // haut droit

    nb_vertex += 4;

    Vertex &bg = verticesTab[nb_vertex - 4];
    Vertex &hg = verticesTab[nb_vertex - 3];
    Vertex &bd = verticesTab[nb_vertex - 2];
    Vertex &hd = verticesTab[nb_vertex - 1];

    facesTab.append(Face(nb_vertex - 4, nb_vertex - 2, nb_vertex - 3)); // Triangle bas gauche
    facesTab.append(Face(nb_vertex - 1, nb_vertex - 3, nb_vertex - 2)); // Triangle haut droit

    nb_faces += 2;

    Face &triangle_bg = facesTab[nb_faces - 2];
    Face &triangle_hd = facesTab[nb_faces - 1];

    bg.i_incident_face = nb_faces - 2;
    hg.i_incident_face = nb_faces - 2;
    bd.i_incident_face = nb_faces - 1;
    hd.i_incident_face = nb_faces - 1;

    triangle_bg.adjacent_faces[0] = nb_faces - 1;
    triangle_bg.adjacent_faces[1] = -1;
    triangle_bg.adjacent_faces[2] = -1;

    triangle_hd.adjacent_faces[0] = nb_faces - 2;
    triangle_hd.adjacent_faces[1] = -1;
    triangle_hd.adjacent_faces[2] = -1;

    for (int i_vertex = 0; i_vertex < (nb_vertex - 4); i_vertex++)
    {
        for (int i_face = 0; i_face < nb_faces; i_face++)
        {
            if (inTriangleTest(facesTab[i_face], verticesTab[i_vertex]) > 0)
            {
                insertionTriangle(i_vertex, i_face);
                lawsonAroundVertex(i_vertex);
                break;
            }
            else if (inTriangleTest(facesTab[i_face], verticesTab[i_vertex]) == 0)
            {
                insertionInArete(i_face, i_vertex); // pb ici probablement
                lawsonAroundVertex(i_vertex);
                break;
            };
        };
    };
}
