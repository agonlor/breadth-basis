#include "cubcomplex.h"

#include <sstream>	// std::istringstream
#include <queue>
#include <math.h>	// pow()
#include <cassert>

/**
 * @brief CubComplex::CubComplex Generates the cubical complex associated to a
 * digital object. It does not insert the cubes.
 * @param obj
 */
CubComplex::CubComplex(const DigObject &obj)
    : m_obj(obj), m_nb_cubes(0)
{
    m_size.resize(3);
    for (int i = 0; i < 3; i++)
    {
        m_size.at(i) = 2 * m_obj.size(i) + 1;
    }
}


/**
 * @brief CubComplex::set_from_image Set the cubical complex from its digital
 * object
 */
void CubComplex::set_from_object()
{
    reset();
    m_cube.resize(m_size.at(0)*m_size.at(1)*m_size.at(2), false);
    for (int i = 0; i < m_obj.size(); i++)
    {
        if (m_obj.at(i))
        {
            const std::list<int> faces = all_faces(cube_of_voxel(i));
            for (auto c : faces)
            {
                if (m_cube.at(c) == false)
                {
                    m_nb_cubes++;
                    m_cube.at(c) = true;
                }
            }
        }
    }
    std::clog << "Cubical complex build with " << m_nb_cubes << " cubes" << std::endl;
}

/**
 * @brief CubComplex::set_from_sdt
 * Set the full cubical complex
 */
void CubComplex::set_from_sdt()
{
    set_full();
}

/**
 * @brief CubComplex::only_q_cells
 * Remove all but the (q-1)-cells, q-cells and (q+1)-cells.
 * I do this because I am usually interested in only one homology group
 */
void CubComplex::only_q_cells(int q)
{
    for (std::size_t c = 0; c < m_cube.size(); ++c)
    {
        if (m_cube.at(c) && abs(q - dim(c)) > 1)
        {
            m_cube.at(c) = false;
            m_nb_cubes--;
        }
    }
    std::clog << "Cubical complex pruned with " << m_nb_cubes << " cubes" << std::endl;

}

/**
 * @brief CubComplex::set_full All the cubes belong to the cubical complex
 * @todo Do this
 */
void CubComplex::set_full()
{
    m_cube.clear();
    m_cube.resize(m_size.at(0) * m_size.at(1) * m_size.at(2), true);
    m_nb_cubes = m_cube.size();
}


void CubComplex::reset()
{
    m_cube.clear();
    m_nb_cubes = 0;

}


///////////////////////////////////////////////////////////////////////////////




/**
 * @brief CubComplex::coordinates
 * @param cube The cube defined as an integer
 * @return The coordinates of a cube in a cubical complex given its index
 */
std::vector<int> CubComplex::coordinates(int cube) const
{
    std::vector<int> coord(3);
    coord[0] = cube % m_size[0];
    cube /= m_size[0];
    coord[1] = cube % m_size[1];
    cube /= m_size[1];
    coord[2] = cube;
    return coord;
}

/**
 * @brief CubComplex::index_of_cube
 * @return The index of a cube in a cubical complex given its coordinates
 */
int CubComplex::index_of_cube(int x, int y, int z) const
{
    assert(0 <= x && x < m_size[0]);
    assert(0 <= y && y < m_size[1]);
    assert(0 <= z && z < m_size[2]);

    int index = 0;
    index += x;
    index += y * m_size[0];
    index += z * m_size[0]*m_size[1];
    return index;
}

/**
 * @brief CubComplex::dim_of_cube
 * @param cube
 * @return The dimension of the cube given its index
 */
int CubComplex::dim(int cube) const
{
    const std::vector<int> coord = coordinates(cube);
    return coord[0]%2 + coord[1]%2 + coord[2]%2;
}


/**
 * @brief CubComplex::translate Obtain the cube after doing a translation.
 * I use this to obtain the vertices around a vertex
 * @return The index of the translated cube or -1 if it is outside of the bounding box
 */
int CubComplex::translate(int cube, int x, int y, int z) const
{
    std::vector<int> coord = coordinates(cube);
    coord.at(0) += x;
    coord.at(1) += y;
    coord.at(2) += z;
    if (!(0 <= coord.at(0) && coord.at(0) < m_size.at(0) &&
          0 <= coord.at(1) && coord.at(1) < m_size.at(1) &&
          0 <= coord.at(2) && coord.at(2) < m_size.at(2)))
        return -1;
    return index_of_cube(coord.at(0), coord.at(1), coord.at(2));
}


/**
 * @brief CubComplex::cube_of_voxel
 * @param index
 * @return The index of the 3-cube associated to a voxel, defined by its index
 */
int CubComplex::cube_of_voxel(int voxel) const
{
    std::vector<int> coord = m_obj.coordinates(voxel);
    return index_of_cube(2*coord.at(0) + 1, 2*coord.at(1) + 1, 2*coord.at(2) + 1);
}

/**
 * @brief CubComplex::all_faces
 * @return The indices of the faces of the 3-cube associated to the voxels with
 * index @param voxel. The cubes are sorted by dimension
 *
 * @todo Do not use the function index_of_cube and compute the indices locally
 * @note The function does not make much sense because I took it from the code for
 * TB pairs and I did not want to rewrite everything
 */
std::list<int> CubComplex::all_faces(int cube) const
{
    const std::vector<int> coord = coordinates(cube); // voxel coordinates
    assert(dim(cube) == 3);
    const int x = coord[0] / 2;
    const int y = coord[1] / 2;
    const int z = coord[2] / 2;

    std::list<int> f;
//    f.push_back(translate(cube, -1, -1, -1)); // dim 0
    f.push_back(index_of_cube(2*x  , 2*y  , 2*z  )); // dim 0
    f.push_back(index_of_cube(2*x+2, 2*y  , 2*z  )); // dim 0
    f.push_back(index_of_cube(2*x+1, 2*y  , 2*z  )); // dim 1
    f.push_back(index_of_cube(2*x  , 2*y+2, 2*z  )); // dim 0
    f.push_back(index_of_cube(2*x  , 2*y+1, 2*z  )); // dim 1
    f.push_back(index_of_cube(2*x+2, 2*y+2, 2*z  )); // dim 0
    f.push_back(index_of_cube(2*x+2, 2*y+1, 2*z  )); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y+2, 2*z  )); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y+1, 2*z  )); // dim 2
    f.push_back(index_of_cube(2*x  , 2*y  , 2*z+2)); // dim 0
    f.push_back(index_of_cube(2*x  , 2*y  , 2*z+1)); // dim 1
    f.push_back(index_of_cube(2*x+2, 2*y  , 2*z+2)); // dim 0
    f.push_back(index_of_cube(2*x+2, 2*y  , 2*z+1)); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y  , 2*z+2)); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y  , 2*z+1)); // dim 2
    f.push_back(index_of_cube(2*x  , 2*y+2, 2*z+2)); // dim 0
    f.push_back(index_of_cube(2*x  , 2*y+2, 2*z+1)); // dim 1
    f.push_back(index_of_cube(2*x  , 2*y+1, 2*z+2)); // dim 1
    f.push_back(index_of_cube(2*x  , 2*y+1, 2*z+1)); // dim 2
    f.push_back(index_of_cube(2*x+2, 2*y+2, 2*z+2)); // dim 0
    f.push_back(index_of_cube(2*x+2, 2*y+2, 2*z+1)); // dim 1
    f.push_back(index_of_cube(2*x+2, 2*y+1, 2*z+2)); // dim 1
    f.push_back(index_of_cube(2*x+2, 2*y+1, 2*z+1)); // dim 2
    f.push_back(index_of_cube(2*x+1, 2*y+2, 2*z+2)); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y+2, 2*z+1)); // dim 2
    f.push_back(index_of_cube(2*x+1, 2*y+1, 2*z+2)); // dim 2
    f.push_back(index_of_cube(2*x+1, 2*y+1, 2*z+1)); // dim 3
    return f;
}


/**
 * @brief CubComplex::boundary
 * @return The list of the faces of @param cube_index of codimension 1
 */
std::vector<int> CubComplex::boundary(int cube) const
{
    assert(at(cube)); // Call it only on a cube of the cubical complex

    std::vector<int> faces;
    const std::vector<int> coord = coordinates(cube);
    int offset = 1;
    for (int i = 0; i < 3; i++)
    {
        if (coord.at(i) % 2 == 1)
        {
            faces.push_back(cube - offset);
            faces.push_back(cube + offset);
        }
        offset *= m_size[i];
    }
    return faces;
}


/**
 * @brief CubComplex::coboundary
 * @param cube
 * @return The cofaces of a cube that belong to the cubical complex
 */
std::vector<int> CubComplex::coboundary(int cube) const
{
    assert(at(cube)); // Call it only on a cube of the cubical complex

    std::vector<int> cofaces;
    const std::vector<int> coord = coordinates(cube);
    int offset = 1;
    int c2;
    for (int i = 0; i < 3; i++)
    {
        if (coord.at(i) % 2 == 0)
        {
            c2 = cube - offset;
            if (coord.at(i)-1 >= 0      && at(c2)) cofaces.push_back(c2);
            c2 = cube + offset;
            if (coord.at(i)+1 < size(i) && at(c2)) cofaces.push_back(c2);
        }
        offset *= m_size[i];
    }
    return cofaces;
}


/**
 * @brief CubComplex::bottom_faces
 * @param cube
 * @return The 0-dimensional faces of the cube
 */
std::list<int> CubComplex::bottom_faces(int cube) const
{
    std::list<int> pointels;
    const std::vector<int> coord = coordinates(cube);
    const int x = coord.at(0);
    const int y = coord.at(1);
    const int z = coord.at(2);
    if (x%2 == 1 && y%2 == 1 && z%2 == 1)
    {
        pointels.push_back(index_of_cube(x-1, y-1, z-1));
        pointels.push_back(index_of_cube(x+1, y-1, z-1));
        pointels.push_back(index_of_cube(x-1, y+1, z-1));
        pointels.push_back(index_of_cube(x+1, y+1, z-1));
        pointels.push_back(index_of_cube(x-1, y-1, z+1));
        pointels.push_back(index_of_cube(x+1, y-1, z+1));
        pointels.push_back(index_of_cube(x-1, y+1, z+1));
        pointels.push_back(index_of_cube(x+1, y+1, z+1));
    }
    else if (x%2 == 1 && y%2 == 1 && z%2 == 0)
    {
        pointels.push_back(index_of_cube(x-1, y-1, z  ));
        pointels.push_back(index_of_cube(x+1, y-1, z  ));
        pointels.push_back(index_of_cube(x-1, y+1, z  ));
        pointels.push_back(index_of_cube(x+1, y+1, z  ));
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 1)
    {
        pointels.push_back(index_of_cube(x-1, y  , z-1));
        pointels.push_back(index_of_cube(x+1, y  , z-1));
        pointels.push_back(index_of_cube(x-1, y  , z+1));
        pointels.push_back(index_of_cube(x+1, y  , z+1));
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 1)
    {
        pointels.push_back(index_of_cube(x , y-1, z-1));
        pointels.push_back(index_of_cube(x , y+1, z-1));
        pointels.push_back(index_of_cube(x , y-1, z+1));
        pointels.push_back(index_of_cube(x , y+1, z+1));
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 0)
    {
        pointels.push_back(index_of_cube(x-1, y  , z  ));
        pointels.push_back(index_of_cube(x+1, y  , z  ));
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 0)
    {
        pointels.push_back(index_of_cube(x  , y-1, z  ));
        pointels.push_back(index_of_cube(x  , y+1, z  ));
    }
    else if (x%2 == 0 && y%2 == 0 && z%2 == 1)
    {
        pointels.push_back(index_of_cube(x  , y  , z-1));
        pointels.push_back(index_of_cube(x  , y  , z+1));
    }
    else
    {
        pointels.push_back(cube);
    }
    return pointels;
}


/**
 * @brief CubComplex::incident_voxels
 * @param cube
 * @return The list of the indices of the voxels that are adjacent to a given cube.
 */
std::list<int> CubComplex::incident_voxels(int cube) const
{
    std::list<int> voxels;
    const std::vector<int> coord = coordinates(cube);
    const int x = coord[0];
    const int y = coord[1];
    const int z = coord[2];
    if (x%2 == 0 && y%2 == 0 && z%2 == 0)
    {
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y+1)/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y+1)/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y-1)/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y-1)/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y+1)/2, (z-1)/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y+1)/2, (z-1)/2));
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y-1)/2, (z-1)/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y-1)/2, (z-1)/2));
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 0)
    {
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y+1)/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y-1)/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y+1)/2, (z-1)/2));
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y-1)/2, (z-1)/2));
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 0)
    {
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y  )/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y  )/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y  )/2, (z-1)/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y  )/2, (z-1)/2));
    }
    else if (x%2 == 0 && y%2 == 0 && z%2 == 1)
    {
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y+1)/2, (z  )/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y+1)/2, (z  )/2));
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y-1)/2, (z  )/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y-1)/2, (z  )/2));
    }
    else if (x%2 == 1 && y%2 == 1 && z%2 == 0)
    {
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y  )/2, (z+1)/2));
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y  )/2, (z-1)/2));
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 1)
    {
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y+1)/2, (z  )/2));
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y-1)/2, (z  )/2));
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 1)
    {
        voxels.push_back(m_obj.index_of_voxel((x+1)/2, (y  )/2, (z  )/2));
        voxels.push_back(m_obj.index_of_voxel((x-1)/2, (y  )/2, (z  )/2));
    }
    else
    {
        voxels.push_back(m_obj.index_of_voxel((x  )/2, (y  )/2, (z  )/2));
    }
    return voxels;
}


/**
 * @brief CubComplex::incident
 * @param cube
 * @return List of incident cubes, that is, the cofaces of the faces
 */
std::list<int> CubComplex::incident(int cube) const
{
    std::list<int> cubes;
    const std::vector<int> coord = coordinates(cube);
    const int x = coord.at(0);
    const int y = coord.at(1);
    const int z = coord.at(2);
    int c;
    if (x%2 == 1 && y%2 == 1 && z%2 == 1)
    {
        c = translate(cube, -1,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0,  0, +1); if (at(c)) cubes.push_back(c);
    }
    else if (x%2 == 1 && y%2 == 1 && z%2 == 0)
    {
        c = translate(cube, -1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, -2,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +2,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -2,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +2,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, +1); if (at(c)) cubes.push_back(c);
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 1)
    {
//        pointels.push_back(index_of_cube(x-1, y  , z-1));
//        pointels.push_back(index_of_cube(x+1, y  , z-1));
//        pointels.push_back(index_of_cube(x-1, y  , z+1));
//        pointels.push_back(index_of_cube(x+1, y  , z+1));
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 1)
    {
//        pointels.push_back(index_of_cube(x , y-1, z-1));
//        pointels.push_back(index_of_cube(x , y+1, z-1));
//        pointels.push_back(index_of_cube(x , y-1, z+1));
//        pointels.push_back(index_of_cube(x , y+1, z+1));
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 0)
    {
        c = translate(cube, -2,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +2,  0,  0); if (at(c)) cubes.push_back(c);
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 0)
    {
        c = translate(cube,  0, -2,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +2,  0); if (at(c)) cubes.push_back(c);
    }
    else if (x%2 == 0 && y%2 == 0 && z%2 == 1)
    {
        c = translate(cube,  0,  0, -2); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0,  0, +2); if (at(c)) cubes.push_back(c);
    }
    return cubes;
}


/**
 * @brief CubComplex::coincident
 * @param cube
 * @return List of coincident cubes, that is the faces of the cofaces
 */
//std::list<int> CubComplex::coincident(int cube) const
//{
//    std::list<int> cubes;
//    int c;
//    const std::vector<int> coord = coordinates(cube);
//    const int x = coord.at(0);
//    const int y = coord.at(1);
//    const int z = coord.at(2);
//    if (x%2 == 1 && y%2 == 1 && z%2 == 1)
//    {
////        c = translate(cube, -1,  0,  0); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube, +1,  0,  0); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0, -1,  0); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0, +1,  0); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0,  0, -1); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0,  0, +1); if (c >= 0 && at(c)) cubes.push_back(c);
//    }
//    else if (x%2 == 1 && y%2 == 1 && z%2 == 0)
//    {
////        c = translate(cube, -1,  0, -1); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube, +1,  0, -1); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0, -1, -1); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0, +1, -1); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube, -2,  0,  0); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube, +2,  0,  0); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0, -2,  0); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0, +2,  0); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube, -1,  0, +1); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube, +1,  0, +1); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0, -1, +1); if (c >= 0 && at(c)) cubes.push_back(c);
////        c = translate(cube,  0, +1, +1); if (c >= 0 && at(c)) cubes.push_back(c);
//    }
//    else if (x%2 == 1 && y%2 == 0 && z%2 == 1)
//    {
////        pointels.push_back(index_of_cube(x-1, y  , z-1));
////        pointels.push_back(index_of_cube(x+1, y  , z-1));
////        pointels.push_back(index_of_cube(x-1, y  , z+1));
////        pointels.push_back(index_of_cube(x+1, y  , z+1));
//    }
//    else if (x%2 == 0 && y%2 == 1 && z%2 == 1)
//    {
////        pointels.push_back(index_of_cube(x , y-1, z-1));
////        pointels.push_back(index_of_cube(x , y+1, z-1));
////        pointels.push_back(index_of_cube(x , y-1, z+1));
////        pointels.push_back(index_of_cube(x , y+1, z+1));
//    }
//    else if (x%2 == 1 && y%2 == 0 && z%2 == 0)
//    {
//        if (at(translate(cube,  0, -1,  0))) // z == 0
//        {
//            c = translate(cube, -1, -1,  0); if (at(c)) cubes.push_back(c);
//            c = translate(cube, +1, -1,  0); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0, -2,  0); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube,  0, +1,  0))) // z == 0
//        {
//            c = translate(cube, -1, +1,  0); if (at(c)) cubes.push_back(c);
//            c = translate(cube, +1, +1,  0); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0, +2,  0); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube,  0,  0, -1))) // y == 0
//        {
//            c = translate(cube, -1,  0, -1); if (at(c)) cubes.push_back(c);
//            c = translate(cube, +1,  0, -1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0,  0, -2); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube,  0,  0, +1))) // y == 0
//        {
//            c = translate(cube, -1,  0, +1); if (at(c)) cubes.push_back(c);
//            c = translate(cube, +1,  0, +1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0,  0, +2); if (at(c)) cubes.push_back(c);
//        }
//    }
//    else if (x%2 == 0 && y%2 == 1 && z%2 == 0)
//    {
//        if (at(translate(cube, -1,  0,  0))) // z == 0
//        {
//            c = translate(cube, -1, -1,  0); if (at(c)) cubes.push_back(c);
//            c = translate(cube, -1, +1,  0); if (at(c)) cubes.push_back(c);
//            c = translate(cube, -2,  0,  0); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube, +1,  0,  0))) // z == 0
//        {
//            c = translate(cube, +1, -1,  0); if (at(c)) cubes.push_back(c);
//            c = translate(cube, +1, +1,  0); if (at(c)) cubes.push_back(c);
//            c = translate(cube, +2,  0,  0); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube,  0,  0, -1))) // x == 0
//        {
//            c = translate(cube,  0, -1, -1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0, +1, -1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0,  0, -2); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube,  0,  0, +1))) // x == 0
//        {
//            c = translate(cube,  0, -1, +1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0, +1, +1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0,  0, +2); if (at(c)) cubes.push_back(c);
//        }
//    }
//    else if (x%2 == 0 && y%2 == 0 && z%2 == 1)
//    {
//        if (at(translate(cube, -1,  0,  0))) // y == 0
//        {
//            c = translate(cube, -1,  0, -1); if (at(c)) cubes.push_back(c);
//            c = translate(cube, -1,  0, +1); if (at(c)) cubes.push_back(c);
//            c = translate(cube, -2,  0,  0); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube, +1,  0,  0))) // y == 0
//        {
//            c = translate(cube, +1,  0, -1); if (at(c)) cubes.push_back(c);
//            c = translate(cube, +1,  0, +1); if (at(c)) cubes.push_back(c);
//            c = translate(cube, +2,  0,  0); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube,  0, -1,  0))) // x == 0
//        {
//            c = translate(cube,  0, -1, -1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0, -1, +1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0, -2,  0); if (at(c)) cubes.push_back(c);
//        }
//        if (at(translate(cube,  0, +1,  0))) // x == 0
//        {
//            c = translate(cube,  0, +1, -1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0, +1, +1); if (at(c)) cubes.push_back(c);
//            c = translate(cube,  0, +2,  0); if (at(c)) cubes.push_back(c);
//        }
//    }
//    return cubes;
//}


/**
 * @brief CubComplex::incident
 * @param cube
 * @return List of coincident cubes, that is the faces of the cofaces
 */
std::list<int> CubComplex::coincident(int cube) const
{
    std::list<int> cubes;
    int c;
    const std::vector<int> coord = coordinates(cube);
    const int x = coord.at(0);
    const int y = coord.at(1);
    const int z = coord.at(2);
    if (x%2 == 1 && y%2 == 1 && z%2 == 1)
    {
//        c = translate(cube, -1,  0,  0); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube, +1,  0,  0); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0, -1,  0); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0, +1,  0); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0,  0, -1); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0,  0, +1); if (c >= 0 && at(c)) cubes.push_back(c);
    }
    else if (x%2 == 1 && y%2 == 1 && z%2 == 0)
    {
//        c = translate(cube, -1,  0, -1); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube, +1,  0, -1); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0, -1, -1); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0, +1, -1); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube, -2,  0,  0); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube, +2,  0,  0); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0, -2,  0); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0, +2,  0); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube, -1,  0, +1); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube, +1,  0, +1); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0, -1, +1); if (c >= 0 && at(c)) cubes.push_back(c);
//        c = translate(cube,  0, +1, +1); if (c >= 0 && at(c)) cubes.push_back(c);
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 1)
    {
//        pointels.push_back(index_of_cube(x-1, y  , z-1));
//        pointels.push_back(index_of_cube(x+1, y  , z-1));
//        pointels.push_back(index_of_cube(x-1, y  , z+1));
//        pointels.push_back(index_of_cube(x+1, y  , z+1));
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 1)
    {
//        pointels.push_back(index_of_cube(x , y-1, z-1));
//        pointels.push_back(index_of_cube(x , y+1, z-1));
//        pointels.push_back(index_of_cube(x , y-1, z+1));
//        pointels.push_back(index_of_cube(x , y+1, z+1));
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 0)
    {
        c = translate(cube, -1, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -2,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +2,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0,  0, -2); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0,  0, +2); if (at(c)) cubes.push_back(c);
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 0)
    {
        c = translate(cube, -1, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, -2,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1, -1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1, +1,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +2,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0,  0, -2); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0,  0, +2); if (at(c)) cubes.push_back(c);
    }
    else if (x%2 == 0 && y%2 == 0 && z%2 == 1)
    {
        c = translate(cube, -1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, -1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube, -2,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +1,  0, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube, +2,  0,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, -2,  0); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, -1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +1, +1); if (at(c)) cubes.push_back(c);
        c = translate(cube,  0, +2,  0); if (at(c)) cubes.push_back(c);
    }
    return cubes;
}



/**
 * @brief Boundary of a cube. We assume that it is a cube of the cubical complex
 * @param face
 * @return
 */
Chain CubComplex::d(int cube) const
{
    Chain x;
    const std::vector<int> faces = boundary(cube);
    for (int face : faces)
        x.insert(face);
    return x;
}

/**
 * @brief Coboundary of a face
 */
Chain CubComplex::cod(int cube) const
{
    Chain x;
    const std::vector<int> cofaces = coboundary(cube);
    for (int coface : cofaces)
        x.insert(coface);
    return x;
}


///////////////////////////////////////////////////////////////////////////////


/**
 * @brief Faces of any dimension of the cubes in L, even themselves
 * @param cubes is a list of cubes
 * @return
 */
//std::list<Cube> CubComplex::allFaces(const std::list<Cube> &cubes) const
//{
////    std::set<Cube> S;
////    for (auto it = cubes.cbegin(); it != cubes.cend(); ++it)
////    {
////        const std::list<Cube> L2 = it->allFaces();
////        S.insert(L2.cbegin(), L2.cend());
////    }
////    std::list<Cube> L2(S.cbegin(), S.cend());
////    return L2;
//}

///**
// * @brief Closure of c (all its faces, event itself) ordered by dimension.
// * The resulting list is a filtration
// * @param c
// * @return
// */
//std::list<Cube> CubComplex::closure_filtered(const Cube &c) const
//{
////    const std::list<Cube> L = c.allFaces();
////    /* Now we sort the cubes by dimension */
////    std::list<Cube> L2;
////    for (int q = 0; q <= m_dimension; q++)
////    {
////        for (auto it = L.cbegin(); it != L.cend(); ++it)
////        {
////            if (it->dim() == q)
////                L2.push_back(*it);
////        }
////    }
////    return L2;
////    /* TODO: I think it is enough to revert L */
////    assert(false); // because this isn't finished
//}

///**
// * @brief Cofaces of any dimension of the cube c in the complex
// * @param c
// * @return
// */
//std::list<Cube> CubComplex::allCofaces(const Cube &c) const
//{
////    const std::list<Cube> L(1, c);
////    return allCofaces(L);
//}


///**
// * @brief Cofaces of any dimension of the list of cubes L in the complex
// * @param L
// * @return
// */
//std::list<Cube> CubComplex::allCofaces(const std::list<Cube> &cubes) const
//{
////    /* Put all the cofaces of L in a set S */
////    std::set<Cube> S;
////    for (auto it = cubes.cbegin(); it != cubes.cend(); ++it)
////    {
////        const std::list<Cube> L2 = it->allCofaces();
////        S.insert(L2.cbegin(), L2.cend());     // S <- L2
////    }
////    /* Take only those that belong to the complex */
////    std::list<Cube> L2;
////    for (auto it = S.cbegin(); it != S.cend(); ++it)
////    {
////        if (isInComplex(*it))
////        {
////            L2.push_back(*it);
////        }
////    }
////    return L2;
//}

///**
// * @brief Maximal cofaces of c in K.
// * I just use it in Filtration::computeFromMap()
// * @param c
// * @return
// */
//std::list<Cube> CubComplex::cofaces_n(const Cube &c) const
//{
////    std::list<Cube> L;
////    const std::list<Cube> cofaces = c.cofaces_n();
////    for (auto it = cofaces.cbegin(); it != cofaces.cend(); ++it)
////    {
////        if (isInComplex(*it))
////            L.push_back(*it);
////    }
////    return L;
//}

///**
// * @brief Check if the cell belong to the complex and get its index
// * @param c
// * @param face
// * @return
// */
//bool CubComplex::isInComplex(const Cube &c, int &face) const
//{
////    std::map<Cube,int>::const_iterator it = m_complex_map.find(c);
////    if (it == m_complex_map.end())
////        return false;
////    face = it->second;
////    return true;
//}

///**
// * @brief Check if the cell belong to the complex
// * @param c
// * @return
// */
//bool CubComplex::isInComplex(const Cube &c) const
//{
////    std::map<Cube,int>::const_iterator it = m_complex_map.find(c);
////    return !(it == m_complex_map.end());
//}

///**
// * @brief Add a single cube to the complex (not its faces)
// * @param c
// * @return The index of this cube
// */
//int CubComplex::addCell(const Cube &c)
//{
////    int i;

////    if (!isInComplex(c, i)) // add it only if it is not there
////    {
////        i = m_complex.size();
////        m_complex.push_back(c);
////        m_complex_map[c] = i;
////        m_numberCells++;
////        m_dimension = std::max(m_dimension, c.dim());
////    }
////    return i;
//}

///**
// * @brief Add a cube and all its faces
// * @param c is the cube to add
// * @return @note I don't know
// */
//int CubComplex::addCellWithFaces(const Cube &c)
//{
////    const std::list<Cube> L = c.allFaces(); // The list of all the cells that we will add

////    for (auto it = L.cbegin(); it != L.cend(); ++it)
////    {
////        addCell(*it);
////    }
////    return m_complex_map[c];
//}


///////////////////////////////////////////////////////////////////


///**
// * @brief Get the numerical identifier of the cube c
// * @param c
// * @return
// */
//int CubComplex::num_of_cube(const Cube &c) const
//{
////    int num = 0;
////    const std::vector<int> coord = c.get_coord();
////    for (std::size_t i = 0; i < coord.size(); i++)
////    {
////        num += coord.at(i) * m_accSize.at(i);
////    }
////    return num;
//}


///**
// * @brief Get the cube with numerical identifier n
// * @param n
// * @return
// */
//Cube CubComplex::cube_with_num(int n) const
//{
////    std::vector<int> coord(m_dimensionEmbedding);

////    for (int i = 0; i < m_dimensionEmbedding; i++)
////    {
////        coord.at(i) = n % m_size.at(i);
////        n /= m_size.at(i);
////    }
////    const Cube c(coord);
////    return c;
//}

///**
// * @brief Get the list of numbers of a list of cubes
// * @param L
// * @return
// */
//std::list<Cube> CubComplex::cube_with_num(const std::list<int> &L) const
//{
////    std::list<Cube> cubes;
////    for (auto it = L.cbegin(); it != L.cend(); ++it)
////    {
////        cubes.push_back(cube_with_num(*it));
////    }
////    return cubes;
//}

///**
// * @brief Get the cube corresponding to the given voxel
// * It can be a n-cube if we built the primal complex, or a 0-cubes otherwise
// *
// * @warning I have some issues with this, because it forces the class to "remember"
// * the image.
// *
// * @param voxel
// * @return
// */
//Cube CubComplex::voxel2cube(int voxel) const
//{
////    std::vector<int> coord = Im->coord_of_voxel(voxel);
////    for (std::size_t i = 0; i < coord.size(); i++)
////    {
////        coord[i] *= 2;
////        if (m_primal > 0)
////            coord[i]++;
////    }
////    const Cube c(coord);
////    return c;
//}

///**
// * @brief  Get the voxel associated to one of the 3-cofaces in the full cubical
// * complex.
// * @warning This is not done for the dual associated cubical complex to an image
// * @param c
// * @return
// */
//int CubComplex::cube2voxel(const Cube &c) const
//{
////    std::queue<Cube> Q;
////    Q.push(c);
////    while (!Q.empty())
////    {
////        const Cube c1 = Q.front();
////        Q.pop();
////        const std::vector<int> coord = c1.get_coord();
////        bool b = true;  // this avoids repeated cubes in the queue
////        for (std::size_t i = 0; i < coord.size() && b; i++)
////        {
////            Cube c2;
////            if (coord.at(i) % 2 == 0)
////            {
////                c2 = c1.nextCube(+1, i);    if (isInBoundingBox(c2)) {Q.push(c2);}
////                c2 = c1.nextCube(-1, i);    if (isInBoundingBox(c2)) {Q.push(c2);}
////                b = false;
////            }
////        }
////        if (b)
////        {
////            std::vector<int> v(coord);
////            for (auto it = v.begin(); it != v.end(); ++it)
////                *it /= 2;
////            return Im->coord2voxel(v);
////        }
////    }
////    assert(false); // I should never arrive here
//}


//bool CubComplex::isInBoundingBox(const Cube &c) const
//{
////    const std::vector<int> coord = c.get_coord();
////    for (std::size_t i = 0; i < coord.size(); i++)
////    {
////        if (coord[i] < 0 || coord[i] >= m_size.at(i))
////        {
////            return false;
////        }
////    }
////    return true;
//}


////////////////////////////////////////////////////////////////

/**
 * @brief Print all the cells with their coordinates
 */
void CubComplex::print() const
{
//	int i = 0;
//    for (auto it = m_complex.cbegin(); it != m_complex.cend(); ++it)
//    {
//        std::cout << i << " : " << *it << std::endl;
//        i++;
//    }
}

/**
 * @brief Print the coordinates
 * @param coord
 */
void CubComplex::printCoord(const std::vector<int> &coord) const
{
    std::cout << "(";
    for (std::size_t i = 0; i < coord.size()-1; i++)
	{
        std::cout << coord.at(i) << ", ";
	}
    std::cout << coord.at(coord.size()-1) << ")";
}

std::vector<int> CubComplex::cubes_of_dim(int q) const
{
    std::vector<int> cubes;
    for (std::size_t i = 0; i < m_cube.size(); i++)
        if (m_cube.at(i) && dim(i) == q)
            cubes.push_back(i);
    return cubes;
}

/**
 * @brief CubComplex::distance
 * @param c1
 * @param c2
 * @return The Euclidean distance between two cubes of dimension 0
 */
double CubComplex::distance(int c1, int c2) const
{
//    assert(dim(c1) == 0);
//    assert(dim(c2) == 0);
    const std::vector<int> coor1 = coordinates(c1);
    const std::vector<int> coor2 = coordinates(c2);
    double d2 = 0;
    // Euclidean distance
    for (int k = 0; k < 3; k++)
        d2 += (coor1.at(k) - coor2.at(k))*(coor1.at(k) - coor2.at(k));
    return sqrt(d2);
    // L1 distance
//    for (int k = 0; k < 3; k++)
//        d2 += abs(coor1.at(k) - coor2.at(k));
    return d2;
}
