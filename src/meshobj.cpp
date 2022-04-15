#include "meshobj.h"

#include <math.h>
#include <cassert>

MeshObj::MeshObj(const CubComplex &K)
    : m_K(K), m_r(0.05), m_s(0)
{
}


/**
 * @brief Compute the mesh of a cubical complex
 */
void MeshObj::set_parameters(double r, double s)
{
    m_r = r;
    m_s = s;
    if (0 > m_s || s > m_r || m_r >= 1)
    {
        std::cerr << "Error in MeshObj::setFromCubComplex(): wrong parameters. Choose 0 <= s <= r < 1" << std::endl;
        exit(EXIT_FAILURE);
    }
}


/**
 * @brief Write the CubComplex mesh in an obj file
 */
void MeshObj::write_complex(const std::string &filename)
{
    for (int q = 0; q <= 3; q++)
    {
        const std::string label = std::to_string(q) + "cubes";

        std::list<std::vector<int>> qfaces;
        for (int c = 0; c < m_K.size_space(); ++c)
        {
            if (m_K.at(c) && m_K.dim(c) == q)
            {
                add_cube(c, qfaces);
            }
        }
        m_faces[label] = qfaces;
    }
    write(filename);
}


void MeshObj::add_chain(const Chain &x, const std::string &p_label)
{
    std::list<std::vector<int>> faces;
    for (int c : x)
    {
        add_cube(c, faces);
    }
    m_faces[p_label] = faces;
}


void MeshObj::add_cube(int cube, std::list<std::vector<int> > &faces)
{
    int n = m_vertices.size();
    std::vector<int> coor = m_K.coordinates(cube);
    std::vector<double> u(3), v(3); // corners of the cubical cell
    for (int j = 0; j < 3; j++)
    {
        if (coor.at(j) % 2 == 0)
        {
            u.at(j) = coor.at(j)/2 - (m_r - m_s);
            v.at(j) = coor.at(j)/2 + (m_r - m_s);
        }
        else
        {
            u.at(j) = coor.at(j)/2 + 0 + (m_r + m_s);
            v.at(j) = coor.at(j)/2 + 1 - (m_r + m_s);
        }
    }
    std::vector<double> vertex(3);
    vertex[0] = u[0]; vertex[1] = u[1]; vertex[2] = u[2]; m_vertices.push_back(vertex);   // (-,-,-)
    vertex[0] = v[0]; vertex[1] = u[1]; vertex[2] = u[2]; m_vertices.push_back(vertex);   // (+,-,-)
    vertex[0] = u[0]; vertex[1] = v[1]; vertex[2] = u[2]; m_vertices.push_back(vertex);   // (-,+,-)
    vertex[0] = v[0]; vertex[1] = v[1]; vertex[2] = u[2]; m_vertices.push_back(vertex);   // (+,+,-)
    vertex[0] = u[0]; vertex[1] = u[1]; vertex[2] = v[2]; m_vertices.push_back(vertex);   // (-,-,+)
    vertex[0] = v[0]; vertex[1] = u[1]; vertex[2] = v[2]; m_vertices.push_back(vertex);   // (+,-,+)
    vertex[0] = u[0]; vertex[1] = v[1]; vertex[2] = v[2]; m_vertices.push_back(vertex);   // (-,+,+)
    vertex[0] = v[0]; vertex[1] = v[1]; vertex[2] = v[2]; m_vertices.push_back(vertex);   // (+,+,+)
    std::vector<int> face(4);
    face[0] = n+2; face[1] = n+4; face[2] = n+8; face[3] = n+6; faces.push_back(face);
    face[0] = n+1; face[1] = n+5; face[2] = n+7; face[3] = n+3; faces.push_back(face);
    face[0] = n+3; face[1] = n+7; face[2] = n+8; face[3] = n+4; faces.push_back(face);
    face[0] = n+1; face[1] = n+2; face[2] = n+6; face[3] = n+5; faces.push_back(face);
    face[0] = n+5; face[1] = n+6; face[2] = n+8; face[3] = n+7; faces.push_back(face);
    face[0] = n+1; face[1] = n+3; face[2] = n+4; face[3] = n+2; faces.push_back(face);
}


void MeshObj::reset()
{
    m_index.clear();
    m_vertices.clear();
    m_faces.clear();
}


int MeshObj::addVertex(const std::vector<double> &vertex)
{
    std::map<std::vector<double>,int>::iterator it = m_index.find(vertex);
    if (it == m_index.end())	// if the vertex was not created before
    {
        m_vertices.push_back(vertex);
        m_index[vertex] = m_vertices.size();
        return m_vertices.size();
    }
    return it->second;
}


void MeshObj::addFace(const std::string &label, const std::vector<int> &face)
{
    const auto it = m_faces.find(label);
    if (it == m_faces.end())   // if that label does not exist
    {
        std::list<std::vector<int>> L;
        L.push_back(face);
        m_faces[label] = L;
    }
    else
    {
        it->second.push_back(face);
    }
}


/**
 * @brief Write the mesh in an obj file
 */
void MeshObj::write(const std::string &filename) const
{
    /* Create the file */
    std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!file)
    {
        std::cerr << "Error in MeshObj::write(" << filename << "): impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    file << "# Obj file created by aldom (Aldo Gonzalez-Lorenzo)" << std::endl;
    file << "# " << m_vertices.size() << " vertices" << std::endl;
    file << "# " << m_faces.size() << " square faces" << std::endl;
    for (auto it = m_vertices.cbegin(); it != m_vertices.cend(); ++it)
    {
        file << "v " << it->at(0) << " " << it->at(1) << " " << it->at(2) << std::endl;
    }
    for (auto it_m = m_faces.cbegin(); it_m != m_faces.cend(); ++it_m)
    {
        file << "o " << it_m->first << std::endl;    // print the label, then the faces
        for (auto it_l = it_m->second.cbegin(); it_l != it_m->second.cend(); ++it_l)
        {
            file << "f";
            for (auto it_v = it_l->cbegin(); it_v != it_l->cend(); ++it_v)
                file << " " << *it_v;
            file << std::endl;
        }
    }
    file.close();
}
