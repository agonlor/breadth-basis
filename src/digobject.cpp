#include "digobject.h"

#include "DGtal/base/BasicFunctors.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/images/IntervalForegroundPredicate.h"

#include <string>

typedef DGtal::ImageSelector<DGtal::Z3i::Domain, int>::Type Image3D;

DigObject::DigObject(const std::string &filename)
    : m_filename(filename),
      m_benchmark(false),
      m_start_time(std::chrono::steady_clock::now())
{
}


void DigObject::set_benchmark()
{
    m_benchmark = true;
    std::cout << "Filename & voxels & sdt & matrix & pairs" << std::endl;
}


/**
 * @brief Read a PGM (3D) file
 * Compute the distance transform of the foreground and the background
 * Write the shifted (so the values start by 0) SDT of the picture
 */
void DigObject::set_from_pgm()
{
    const Image3D volume = DGtal::GenericReader<Image3D>::import(m_filename + ".pgm");
    const DGtal::Z3i::Domain dom = volume.domain();

    if (!m_benchmark)
    {
        std::cout << "Name: " << m_filename + ".pgm" << std::endl;
        std::cout << "Size of the image: "
                  << dom.myUpperBound[0] + 1 << " x "
                  << dom.myUpperBound[1] + 1 << " x "
                  << dom.myUpperBound[2] + 1
                  << " [" << dom.size() << "]" << std::endl;
    }

    m_size.resize(3);
    for (int i = 0; i < 3; i++)
        m_size[i] = dom.myUpperBound[i] + 1;

    m_voxel.resize(m_size.at(0)*m_size.at(1)*m_size.at(2), false);
    std::size_t i = 0;
    for (DGtal::Z3i::Domain::Iterator it = dom.begin(); it != dom.end(); it++)
    {
       m_voxel.at(i) = (volume(*it) > 0);
       i++;
    }

    prune_small_components();

    //    for (Z3i::Domain::Iterator it = dom.begin(); it != dom.end(); it++) // print the values of the voxels for debugging
    //        std::cout << *it << " value = " << volume(*it) << std::endl;
    //    std::cout << "Volume loaded from " << fileName << "." << std::endl;

    typedef DGtal::functors::IntervalForegroundPredicate<Image3D> Binarizer;
    const Binarizer b_in(volume, 0, 255);    // the interval ]0,255] is considered true
    typedef DGtal::DistanceTransformation<DGtal::Z3i::Space, Binarizer, DGtal::Z3i::L2Metric> DTL2;   // the Distance Transform for L2 (Euclidean)
    const DTL2 dt_in(&volume.domain(), &b_in, &DGtal::Z3i::l2Metric);   // compute the distance transform in dt_in

    const Binarizer b_out(volume, -1, 0);
    const DTL2 dt_out(&volume.domain(), &b_out, &DGtal::Z3i::l2Metric);

    m_sdt.resize(dom.size());
    i = 0;
    for (DGtal::Z3i::Domain::ConstIterator it = dom.begin(); it != dom.end(); it++)
    {
        m_sdt.at(i) = dt_out(*it) - dt_in(*it);
        i++;
    }

    if (m_benchmark)
    {
        std::cout << m_filename << " & "
                  << m_size[0]*m_size[1]*m_size[2] << " & "
                  << elapsed_sec() << " & ";
    }
//    else
//        std::clog << "Signed distance transform computed [" << elapsed_sec() << "s]" << std::endl;
}

/**
 * @brief DigObject::index_of_voxel
 * @return The index of voxel (x,y,z) or -1 if it is outside the box
 */
int DigObject::index_of_voxel(int x, int y, int z) const
{
    if (x < 0 || x >= (int)m_size.at(0)) return -1;
    if (z < 0 || y >= (int)m_size.at(1)) return -1;
    if (x < 0 || z >= (int)m_size.at(2)) return -1;

    int index = 0;
    index += x;
    index += y * m_size[0];
    index += z * m_size[0]*m_size[1];
    return index;
}

/**
 * @brief DigObject::coordinates
 * @param index
 * @return The coordinates of a voxel in a discrete object given its index
 */
std::vector<int> DigObject::coordinates(int index) const
{
    std::vector<int> coord(3);
    coord[0] = index % m_size[0];
    index /= m_size[0];
    coord[1] = index % m_size[1];
    index /= m_size[1];
    coord[2] = index;
    return coord;
}


/**
 * @brief DigObject::value_of_cube
 * @param cube
 * @return Given the index of a cube, get its incident voxels (3-cubes) and get the minimum value
 * @note This is not used anymore
 */
//double DigObject::value_of_cube(int cube) const
//{
//    const std::list<int> voxels = incident_voxels(cube);
//    double min_sdt = std::numeric_limits<double>::max();
//    for (int voxel : voxels)
//        if (voxel >= 0)
//            min_sdt = std::min(min_sdt, m_sdt.at(voxel));
//    assert(min_sdt != std::numeric_limits<double>::max());
//    return min_sdt;
//}


/**
 * @brief DigObject::voxel_of_cube
 * @param cube
 * @return Given a cube in the filter, return the index of the voxels that adds
 * it to the filter. That is, the adjacent voxel with least signed distance
 * transform value
 */
//int DigObject::voxel_of_cube(int cube) const
//{
//    const std::list<int> voxels = incident_voxels(cube);
//    double min_sdt = std::numeric_limits<double>::max();
//    int first_voxel = 0;
//    for (int voxel : voxels)
//        if (voxel >= 0 && m_sdt.at(voxel) < min_sdt)
//        {
//            min_sdt = m_sdt.at(voxel);
//            first_voxel = voxel;
//        }
//    assert(min_sdt != std::numeric_limits<double>::max());
//    return first_voxel;
//}


double DigObject::elapsed_sec() const
{
    auto cur_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = cur_time - m_start_time;
    return elapsed_seconds.count();
}

/**
 * @brief DigObject::distance
 * @param i
 * @param j
 * @return Euclidean distance between two voxels;
 */
double DigObject::distance(int i, int j) const
{
    const std::vector<int> coor1 = coordinates(i);
    const std::vector<int> coor2 = coordinates(j);
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


/**
 * @brief DigObject::prune_small_components
 * Remove all the connected components except dor the largest one.
 * warning
 */
void DigObject::prune_small_components()
{
    std::vector<std::vector<int>> components;
    std::set<int> visited;
    for (std::size_t v = 0; v < m_voxel.size(); v++)
        if (m_voxel.at(v) && visited.count(v) == 0)
        {
            std::vector<int> component;
            std::queue<int> fifo; // queue of voxels to visit
            fifo.push(v);
            while (!fifo.empty())
            {
                const int cur_v = fifo.front(); fifo.pop();
                if (visited.count(cur_v) == 0)
                {
                    component.push_back(cur_v);
                    visited.insert(cur_v);
                    int next_v;
                    next_v = translate(cur_v, -1, -1, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v,  0, -1, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1, -1, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, -1,  0, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v,  0,  0, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1,  0, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, -1, +1, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v,  0, +1, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1, +1, -1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, -1, -1,  0); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v,  0, -1,  0); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1, -1,  0); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, -1,  0,  0); if (at(next_v)) fifo.push(next_v);
//                    next_v = translate(cur_v,  0,  0,  0); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1,  0,  0); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, -1, +1,  0); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v,  0, +1,  0); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1, +1,  0); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, -1, -1, +1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v,  0, -1, +1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1, -1, +1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, -1,  0, +1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v,  0,  0, +1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1,  0, +1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, -1, +1, +1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v,  0, +1, +1); if (at(next_v)) fifo.push(next_v);
                    next_v = translate(cur_v, +1, +1, +1); if (at(next_v)) fifo.push(next_v);
                }
            }
            components.push_back(component);
        }
    if (components.size() == 1)
        return;

    // sort the components by size
    std::sort(components.begin(), components.end(), [](const std::vector<int> &c1, const std::vector<int> &c2) {
        return c1.size() < c2.size();
    });
    int count = 0;
    for (std::size_t i = 0; i < components.size() - 1; ++i)
    {
        for (int v : components.at(i))
            m_voxel.at(v) = false;
        count += components.at(i).size();
    }
    std::clog << count << " voxels removed to clear small components" << std::endl;
}


int DigObject::translate(int voxel, int x, int y, int z) const
{
    std::vector<int> coord = coordinates(voxel);
    coord.at(0) += x;
    coord.at(1) += y;
    coord.at(2) += z;
    if (!(0 <= coord.at(0) && coord.at(0) < size(0) &&
          0 <= coord.at(1) && coord.at(1) < size(1) &&
          0 <= coord.at(2) && coord.at(2) < size(2)))
        return -1;
    return index_of_voxel(coord.at(0), coord.at(1), coord.at(2));
}
