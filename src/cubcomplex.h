#ifndef CUBCOMPLEX_H
#define CUBCOMPLEX_H

#include "digobject.h"
#include <set>
#include <string>
#include <vector>
#include <list>
#include <set>

/**
 * @class CubComplex
 * @brief Class for a cubical complex. It can be the cubical complex associated
 * to a digital object (with 26-connectivity) or a full cubical complex.
 *
 * The cubical complex can also be full.
 */

typedef std::set<int> Chain;

class CubComplex
{
public:
    CubComplex(const DigObject & obj);
    void set_from_object();
    void set_from_sdt();
    void only_q_cells(int q);

    Chain d(int cube) const;
    Chain cod(int cube) const;

    std::vector<int> coordinates(int cube) const;
    int index_of_cube(int x, int y, int z) const;
    int dim(int cube) const;
    int translate(int cube, int x, int y, int z) const;
    int cube_of_voxel(int voxel) const;
    std::list<int> all_faces(int cube) const;
    std::vector<int> boundary(int cube) const;
    std::vector<int> coboundary(int cube) const;
    std::list<int> bottom_faces(int cube) const;
    std::list<int> incident_voxels(int cube) const;
    std::list<int> incident(int cube) const;
    std::list<int> coincident(int cube) const;
    bool in_space(int x, int y, int z) const;
    double distance(int c1, int c2) const;


//    std::list<Cube> allFaces(const std::list<Cube> &cubes) const;
//    std::list<Cube> closure_filtered(const Cube &c) const;
//    std::list<Cube> allCofaces(const Cube &c) const;
//    std::list<Cube> allCofaces(const std::list<Cube> &cubes) const;
//    std::list<Cube> cofaces_n(const Cube &c) const;
//    bool isInComplex(const Cube &c, int &face) const;
//    bool isInComplex(const Cube &c) const;

//    int addCell(const Cube &c);
//    int addCellWithFaces(const Cube &c);

//    /* Coordinates */
//    int num_of_cube(const Cube &c) const;
//    Cube cube_with_num(int n) const;
//    std::list<Cube> cube_with_num(const std::list<int> &L) const;
//    Cube voxel2cube(int voxel) const;
//    int cube2voxel(const Cube &c) const;
//    bool isInBoundingBox(const Cube &c) const;
    
    /* Input/Output */
    void print() const;
    void printCoord(const std::vector<int> &coord) const;

    /* Getters */
//    inline Cube get_complex(int face)   const { return m_complex.at(face);   }
    inline int size_space() const { return m_size.at(0)*m_size.at(1)*m_size.at(2); }
    inline int size(int i) const { return m_size.at(i); }
    inline int size() const { return m_nb_cubes; }
    inline bool at(int cube) const { return cube >= 0 && m_cube.at(cube); }
    std::vector<int> cubes_of_dim(int q) const;
    inline std::string filename() const { return m_obj.filename(); }

private:
    void set_full();
    void reset();
    
private:
    const DigObject &m_obj;             /// Object from which we build the cubical complex (with 26-connectivity)
    std::vector<int> m_size;            /// Size in every direction
    std::vector<bool> m_cube;           /// m_cube[c] == true if the cube belongs to the object
    int m_nb_cubes;             /// number of cubes
//    std::map<Cube, int> m_complex_map;  /// Map between the cells and its indices in complex: cube -> face */
//    std::vector<Cube> m_complex;        /// Vector of the cubes in the complex: face -> cube

};

#endif // CUBCOMPLEX_H
