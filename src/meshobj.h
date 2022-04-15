#ifndef MESHOBJ_H
#define MESHOBJ_H

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <string>
#include "cubcomplex.h"

/**
 * @class MeshObj
 * @brief A class for making obj files. There is a list of vertices in R^3 and
 * a list of faces (triangles or rectangles).
 *
 * The parameter m_r measure the width of the cubical cells. It must be between
 * 0 and 1. A good value is 0.1
 * The parameter m_s measures the gap between the cells. It must be smaller
 * than m_r, and it can be zero.
 *
 */
class MeshObj
{
public:
    MeshObj(const CubComplex &K);

    void set_parameters(double r, double s);
    void write_complex(const std::string &filename);
    void add_chain(const Chain &x, const std::string &p_label);
    void add_cube(int cube, std::list<std::vector<int>> &faces);
    void reset();




    int addVertex(const std::vector<double> &vertex);
    void addFace(const std::string &label, const std::vector<int> &face);

    void write(const std::string &filename) const;

    inline int get_nVertices() const { return m_vertices.size(); }

private:
    const CubComplex &m_K;
    double m_r;   /// radius of each q-cube
    double m_s;   /// separation between the cubes

    std::map<std::vector<double>,int> m_index;                    /// index[(z,y,z)] gives the index of the vertex (x,y,x). Inverse of "vertices"
    /* TODO vector or lists? */
    std::vector<std::vector<double>> m_vertices;                  /// vertices in the order to print
    std::map<std::string, std::list<std::vector<int>>> m_faces;   /// faces ordered by label
};

#endif // MESHOBJ_H
