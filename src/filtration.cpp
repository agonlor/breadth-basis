#include "filtration.h"
//#include "discgeodesic.h"
#include "meshobj.h"

#include <cassert>
#include <numeric>      // std::iota

Filtration::Filtration(const CubComplex &complex)
    : m_K(complex)
{
    m_pos_of_cube.resize(m_K.size_space(), -1);
    m_cube_at_pos.resize(m_K.size(), -1);
    m_value_of_pos.resize(m_K.size());
}


/**
 * @brief Filtration::geodesic
 * Make the filtration defined by the geodesic distance from a 0-cube of the complex.
 * @param c A 0-dimensional cube
 */
void Filtration::geodesic(int c)
{
    assert(m_K.dim(c) == 0);

    // compute the geodesic distance of the 0-cubes
    std::map<int,int> dist; // map from 0-cubes to real (actually only positive integers)
    std::queue<std::pair<int, int>> fifo;
    fifo.push(std::make_pair(c, 0));
    while (!fifo.empty())
    {
        const std::pair<int, int> cur_pair = fifo.front(); fifo.pop();
        if (dist.count(cur_pair.first) == 0)
        {
            const int cube = cur_pair.first;
            const int d = cur_pair.second;
            dist[cube] = d;
            int c2;
            c2 = m_K.translate(cube, -2, 0, 0); if (m_K.at(c2)) fifo.push(std::make_pair(c2, d+1));
            c2 = m_K.translate(cube, +2, 0, 0); if (m_K.at(c2)) fifo.push(std::make_pair(c2, d+1));
            c2 = m_K.translate(cube, 0, -2, 0); if (m_K.at(c2)) fifo.push(std::make_pair(c2, d+1));
            c2 = m_K.translate(cube, 0, +2, 0); if (m_K.at(c2)) fifo.push(std::make_pair(c2, d+1));
            c2 = m_K.translate(cube, 0, 0, -2); if (m_K.at(c2)) fifo.push(std::make_pair(c2, d+1));
            c2 = m_K.translate(cube, 0, 0, +2); if (m_K.at(c2)) fifo.push(std::make_pair(c2, d+1));
        }
    }

    // assign to each cube the maximal value of its 0-cubes and sort them
    std::vector<std::tuple<int,int,int>> items; // tuples (value, dimension, cube)
    for (int c = 0; c < m_K.size_space(); c++)
    {
        if (m_K.at(c))
        {
            const std::list<int> pointels = m_K.bottom_faces(c);
            int max_value = 0;
            for (const int c2 : pointels)
                max_value = std::max(max_value, dist.at(c2));
            items.push_back(std::make_tuple(max_value, m_K.dim(c), c));
        }
    }
    std::sort(items.begin(), items.end());
    int pos = 0;
    for (const auto tup : items)
    {
        m_pos_of_cube.at(std::get<2>(tup)) = pos;
        m_cube_at_pos.at(pos) = std::get<2>(tup);
        m_value_of_pos.at(pos) = std::get<0>(tup);
        pos++;
    }
}


/**
 * @brief Filtration::sdt
 * Compute the filtration based on the signed distance transform (sdt) of the
 * input digital object
 */
void Filtration::sdt(const DigObject &DO)
{
    // sort voxels by signed distance transform value
    std::vector<std::pair<double, int>> voxels(DO.size());
    for (std::size_t i = 0; i < voxels.size(); i++)
        voxels.at(i) = std::make_pair(DO.sdt(i), i);
    std::sort(voxels.begin(), voxels.end());
//    std::clog << "voxels sorted  [" << DO.elapsed_sec() << "s]" << std::endl;

    // build filtered cubical complex
    assert(m_K.size_space() == m_K.size());
    int pos = 0;
    for (const auto &voxel : voxels)
    {
        // add faces of this voxels, sorted by dimension, unless they are already added
        const std::list<int> cubes = m_K.all_faces(m_K.cube_of_voxel(voxel.second));
        for (int cube : cubes)
        {
            if (m_K.at(cube) && m_pos_of_cube.at(cube) < 0) // if that cube is not already added
            {
                m_pos_of_cube.at(cube) = pos;
                m_cube_at_pos.at(pos) = cube;
                m_value_of_pos.at(pos) = voxel.first;
                pos++;
            }
        }
    }
    std::clog << "Filtration sdt done [" << DO.elapsed_sec() << "s]" << std::endl;
}

/**
 * @brief Filtration::distance
 * Make a filtration of the cubical complex associated to the digital object DO
 * based on the distance to the voxel associated to the cube tau
 * @param DO
 * @param tau
 */
void Filtration::distance_voxel(const DigObject &DO, int tau)
{
    // find the voxel that made the cube tau
    std::list<int> inc_voxels = m_K.incident_voxels(tau);
    double min_sdt = std::numeric_limits<double>::max();
    int center_voxel = 0;
    for (int voxel : inc_voxels)
        if (voxel >= 0 && DO.sdt(voxel) < min_sdt)
        {
            min_sdt = DO.sdt(voxel);
            center_voxel = voxel;
        }
    assert(min_sdt != std::numeric_limits<double>::max());

    // sort voxels of the object by distance to center_voxel
    std::vector<std::pair<double, double>> voxels;
    for (int i = 0; i < DO.size(); i++)
        if (DO.at(i))
            voxels.push_back(std::make_pair(DO.distance(center_voxel, i), i));
    std::sort(voxels.begin(), voxels.end());

    // build filtered cubical complex
    int pos = 0;
    for (auto it = voxels.cbegin(); it != voxels.cend(); ++it)
    {
        // add faces of this voxels, sorted by dimension, unless they are already added
        const std::list<int> cubes = m_K.all_faces(m_K.cube_of_voxel(it->second));
        for (int cube : cubes)
        {
            if (m_K.at(cube) && m_pos_of_cube.at(cube) < 0) // if that cube is not already added
            {
                m_pos_of_cube.at(cube) = pos;
                m_cube_at_pos.at(pos) = cube;
                m_value_of_pos.at(pos) = it->first;
            }
        }
    }
}

/**
 * @brief Filtration::distance_vertex
 * Make a filtration of the cubical complex associated to the digital object DO
 * based on the distance from vertices to the cubical cell tau
 * @param DO
 * @param tau
 */
void Filtration::distance_vertex(int tau)
{
    std::list<int> tau_v = m_K.bottom_faces(tau); // vertices incident to tau
    std::map<int,int> dist; // map from 0-cubes to real (actually only positive integers)
    for (int c = 0; c < m_K.size_space(); c++)
    {
        if (m_K.at(c) && m_K.dim(c) == 0)
        {
            double min_dist = std::numeric_limits<double>::max();
            for (int c2 : tau_v)
                min_dist = std::min(min_dist, m_K.distance(c, c2));
            dist[c] = min_dist;
        }
    }

    // assign to each cube the maximal value of its 0-cubes and sort them
    std::vector<std::tuple<int,int,int>> items; // tuples (value, dimension, cube)
    for (int c = 0; c < m_K.size_space(); c++)
    {
        if (m_K.at(c))
        {
            const std::list<int> pointels = m_K.bottom_faces(c);
            int max_value = 0;
            for (const int c2 : pointels)
                max_value = std::max(max_value, dist.at(c2));
            items.push_back(std::make_tuple(max_value, m_K.dim(c), c));
        }
    }
    std::sort(items.begin(), items.end());

    int pos = 0;
    for (const auto tup : items)
    {
        m_pos_of_cube.at(std::get<2>(tup)) = pos;
        m_cube_at_pos.at(pos) = std::get<2>(tup);
        m_value_of_pos.at(pos) = std::get<0>(tup);
        pos++;
    }
}


/**
 * @brief Filtration::geodesic_vertex
 * Find the closest 0-cube of K to any of the 0-faces of tau.
 * Then, do a geodesic filtration
 * @param tau
 * @todo I think that this search can be optimized
 */
void Filtration::geodesic_vertex(int tau)
{
    std::list<int> tau_v = m_K.bottom_faces(tau); // vertices incident to tau
    double min_dist = std::numeric_limits<double>::max();
    int cube = -1;
    for (int c = 0; c < m_K.size_space(); c++)
        if (m_K.at(c) && m_K.dim(c) == 0)
            for (int c2 : tau_v)
                if (m_K.distance(c, c2) < min_dist)
                {
                    min_dist = m_K.distance(c, c2);
                    cube = c;
                }
    geodesic(cube);
}


void Filtration::trivial()
{
    for (int c = 0; c < m_K.size_space(); c++)
        if (m_K.at(c) && m_K.dim(c) == 0)
        {
            geodesic(c);
            return;
        }



//    int pos = 0;
//    for (int q = 0; q <= 3; q++)
//        for (int c = 0; c < m_K.size_space(); c++)
//            if (m_K.at(c) && m_K.dim(c) == q)
//            {
//                m_pos_of_cube.at(c) = pos;
//                m_cube_at_pos.at(pos) = c;
//                pos++;
//            }
}
