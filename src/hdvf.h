#ifndef HDVF_H
#define HDVF_H

#include "cubcomplex.h"
#include "filtration.h"

#include <vector>
#include <string>
#include <map>
#include <set>

typedef std::map<int, Chain> mapset;

enum class CellType { Primary, Secondary, Critical };

struct Annotation
{
    int dimension = 0;
    std::vector<std::pair<int,int>> pairs;
    std::map<int,Chain> f;

    bool coeff(const Chain &x, int gamma)
    {
        bool nonzero = false;
        for (int c : x)
            if (f[gamma].count(c) > 0 || c == gamma)
                nonzero = !nonzero;
        return nonzero;
    }

    Chain image(const Chain &x)
    {
        Chain y;
        for (const auto &it : f)
            if (coeff(x, it.first))
                y.insert(it.first);
        return y;
    }

    void update(const Chain &x, int gamma)
    {
        const Chain F2 = f.at(gamma);
        const Chain D2 = image(x);
        for (int gamma2 : D2)
        {
            Chain f_gamma2 = f.at(gamma2);
            for (int sigma : F2)
            {
                std::pair<Chain::iterator,bool> ret = f_gamma2.insert(sigma);
                if (ret.second == false) // the chqin already contains the cell
                    f_gamma2.erase(ret.first);
            }
            f[gamma2] = f_gamma2;
        }
        f.erase(gamma);
    }
};

struct HDVF_Options
{
    bool all_dimensions = true;
    int one_dimension = -1;
    bool h = true;
    bool f = true;
    bool g = true;
    bool d = true;
};

/**
 * @class HDVF
 * @brief This is the class that deals with the HDVF.
 *
 * It is similar to HDVF, but it uses maps instead of matrices.
 * The code is uglier, but the class HDVF is too slow
 *
 * Note that:
 * 1. If we read a V file, we only need to compute H at each step. The rest can be done at the end.
 * 2. If we are building a HDVF, D must be computedat each step. But, if we want the reduction, we also need H.
 *    We can do F and G at the end.
 * 3. There are three types: [1] with (d), [2] with (d,h) and [3] with (d,h,f,g)
 *
 *
 */
class HDVF
{
public:
    HDVF(const CubComplex &K);
    void set_options(const HDVF_Options &options);
    HDVF(const CubComplex &K, const HDVF_Options &opt);

    /* Methods for establishing a HDVF */
    void set_from_pairing(const std::string &filename);

    void homology();
    void homology_rand();
    void compute_filtration(const Filtration &F);
    Annotation compute_tb_filtration(const Filtration &F, int q);

    void add(int s, int t);

    Chain critical_cells(int q) const;

    Chain get(const mapset &m, int k) const;
    int coeff(const Chain &x, int face) const;
    void add_to_chain(Chain &x, int face) const;
    Chain get_Fr(int c) const;
    Chain get_Gc(int c) const;

    Chain get_d(int c) const;
    Chain get_dstar(int c) const;
    Chain get_h(int c) const;
    Chain get_fstar(int c) const;
    Chain f(const Chain &x) const;
    Chain g(int c) const;
    void write_fstar(const std::string &filename) const;
    void write_g(const std::string &filename) const;

    std::vector<Chain> homology_cycles(int q) const;


    /* Algebraic operators computed using the HDVF */

    bool isCritical(int i) const;
    void numberCriticalCells() const;

    /* I/O */
    void print_Matrix(char c, int q) const;
    void save(const std::string &filename) const;

    inline CellType type(int s) const { return m_partition.at(s); }
    inline int nb_criticals(int q) const { return m_Dc.at(q).size(); }
    void reset();

private:
    void update_reduction(int s, int t);
    bool find_pair(int &s, int &t) const;




private:
    const CubComplex & m_K;    /// Pointer to the CubComplex where we define the HDVF
    std::map<int, CellType> m_partition;    /// A map sending each cube to its cell type (P, S or C)
    std::vector<mapset> m_Hc;  /// Hc[q](c) = h(c), where c is a face of K
    std::vector<mapset> m_Fr;  /// Fr[q](c) = f*(c)
    std::vector<mapset> m_Gc;  /// Gc[q](c) = g(c)
    std::vector<mapset> m_Dc;  /// Dc[q](c) = d(c)
    std::vector<mapset> m_Dr;  /// Dr[q](c) = d*(c)
    HDVF_Options m_options;
};


#endif // HDVF_H
