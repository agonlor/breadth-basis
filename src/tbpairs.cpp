#include "tbpairs.h"

#include "filtration.h"

#include <numeric>      // std::iota

// PHAT
#include "../include/phat/compute_persistence_pairs.h" // wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include "../include/phat/representations/default_representations.h" // main data structure (choice affects performance)
#include "../include/phat/algorithms/standard_reduction.h" // algorithm (choice affects performance)
#include "../include/phat/algorithms/chunk_reduction.h"
#include "../include/phat/algorithms/row_reduction.h"
#include "../include/phat/algorithms/twist_reduction.h"

TBPairs::TBPairs(const DigObject &DO)
    : m_DO(DO), m_K(DO)
{
    m_K.set_from_sdt();
}

/**
 * @brief TBPairs::compute
 * @param q
 * @return The TB pairs of cells
 */
std::vector<std::pair<int,int>> TBPairs::compute(int q)
{
    // make the filter using the signed distance transform
    Filtration F(m_K);
    F.sdt(m_DO);

    // make the boundary matrix
    // first define a boundary matrix with the chosen internal representation
    phat::boundary_matrix<phat::bit_tree_pivot_column> boundary_matrix;
    boundary_matrix.set_num_cols(m_K.size_space());
    std::vector<phat::index> temp_col; // current column of the matrix
    for (int i = 0; i < m_K.size_space(); i++)
    {
        boundary_matrix.set_dim(i, m_K.dim(F.cube_at_pos(i)));
        const Chain f = m_K.d(F.cube_at_pos(i)); // list of faces of the cube
        temp_col.resize(f.size());
        int j = 0;
        for (int face : f)
        {
            temp_col.at(j) = F.pos_of_cube(face);
            j++;
        }
        std::sort(temp_col.begin(), temp_col.end()); // @note This is very important!
        boundary_matrix.set_col(i, temp_col);
    }

//    std::clog << "The boundary matrix has " << boundary_matrix.get_num_cols() << " columns and " << boundary_matrix.get_num_entries() << " entries." << std::endl;
//    std::clog << "Boundary matrix built [" << m_DO.elapsed_sec() << "s]" << std::endl;

    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs<phat::chunk_reduction>(pairs, boundary_matrix);
//    std::clog << "Persistent homology computed [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // get the thickness-breadth pairs
    std::vector<std::pair<int,int>> tb_pairs;
    for (phat::index idx = 0; idx < pairs.get_num_pairs(); idx++)
    {
//        std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
        const int t_cell = F.cube_at_pos(pairs.get_pair(idx).first);
        const int b_cell = F.cube_at_pos(pairs.get_pair(idx).second);
        const double t_value = F.value(t_cell);
        const double b_value = F.value(b_cell);
        if (m_K.dim(t_cell) == q && t_value < 0 && b_value > 0)
        {
            tb_pairs.push_back(std::make_pair(t_cell, b_cell));
        }
    }
    return tb_pairs;
}
