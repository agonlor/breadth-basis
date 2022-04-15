#ifndef FILTRATION_H
#define FILTRATION_H

#include "cubcomplex.h"

#include <string>
#include <fstream>  // Pour lire des fichier de texte.
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <list>
#include <climits>

/**
* @class Filtration
* @brief Class for computing a filtration on a cubical complex
*
* @note Is it any better to make m_pos_of_cube a std::vector<int>?
*/
class Filtration
{
public:
    Filtration(const CubComplex & complex);

    void geodesic(int c);
    void sdt(const DigObject &DO);
    void distance_voxel(const DigObject &DO, int tau);
    void distance_vertex(int tau);
    void geodesic_vertex(int tau);
    void trivial();

    inline int pos_of_cube(int cube) const { return m_pos_of_cube.at(cube); }
    inline int cube_at_pos(int pos)  const { return m_cube_at_pos.at(pos);  }
    inline int value(int cube)       const { return m_value_of_pos.at(m_pos_of_cube.at(cube)); }
    inline std::size_t size()        const { return m_cube_at_pos.size();   }

private:
    const CubComplex &m_K;
//    std::map<int,int> m_pos_of_cube;    /// position of a cube (given its index) in the filter
    std::vector<int> m_pos_of_cube;    /// position of a cube (given its index) in the filter
    std::vector<int> m_cube_at_pos;     /// the inverse: the (index of the) cube at a given position in the filter
    std::vector<int> m_value_of_pos;    /// value by the filtration of the cube at a given position in the filter
};

#endif // FILTRATION_H
