#ifndef DIGOBJECT_H
#define DIGOBJECT_H

#include <string>
#include <fstream>
#include <vector>
#include <chrono>

// DGtal
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/images/ImageSelector.h>

/**
 * @brief The DigObject class
 * Class for the binary object (2D image or 3D volume) and the computation of
 * its thickness-breadth pairs
 */
class DigObject
{
public:
    DigObject(const std::string &filename);
    void set_benchmark();
    void set_from_pgm();

    std::vector<int> coordinates(int index) const;
    int index_of_voxel(int x, int y, int z) const;
//    double value_of_cube(int cube) const;
//    int voxel_of_cube(int cube) const;
    double distance(int i, int j) const;

    inline int size(int i) const { return m_size.at(i); }
    inline int size() const { return m_size.at(0)*m_size.at(1)*m_size.at(2); }
    inline bool at(int voxel) const { return voxel >= 0 && m_voxel.at(voxel); }
    inline std::string filename() const { return m_filename; }
    inline double sdt(int i) const { return m_sdt.at(i); }

    double elapsed_sec() const;

private:
    void prune_small_components();
    int translate(int voxel, int x, int y, int z) const;

private:
    std::string m_filename;         /// name of the file containing the object
    std::vector<unsigned> m_size;   /// number of voxels in the object along each axis
    std::vector<double> m_sdt;      /// signed distance transform
    std::vector<bool> m_voxel;      /// m_voxel[v] == true if the voxel belongs to the object
    bool m_benchmark;               /// output running times in one line
    std::chrono::steady_clock::time_point m_start_time;
};

#endif // DIGOBJECT_H
