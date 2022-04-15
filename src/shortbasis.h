#ifndef SHORTBASIS_H
#define SHORTBASIS_H

#include "digobject.h"
#include "cubcomplex.h"
#include "filtration.h"
#include "hdvf.h"



enum class Algorithm { VoxelDistance, VertexDistance, Geodesic, Chen, Dey, Deyq1 };


struct Options
{
    std::string filename = "none";
    Algorithm algorithm = Algorithm::Geodesic;
    int seed = 0;
    bool benchmark = false;

    std::string comment() const
    {
        std::string str;
        if (algorithm == Algorithm::VoxelDistance)
            str = "tb_voxels_distance";
        else if (algorithm == Algorithm::VertexDistance)
            str = "tb_vertex_distance";
        else if (algorithm == Algorithm::Geodesic)
            str = "tb";
        else if(algorithm == Algorithm::Chen)
            str = "chen_freedman_2010";
        else if(algorithm == Algorithm::Dey)
            str = "deylw_2018";
        else if(algorithm == Algorithm::Deyq1)
            str = "deylw_2018_q1";
        else
            str = "undefined";
        return str;
    }
};

class ShortBasis
{
public:
    ShortBasis(const DigObject &DO, const Options &options);

private:
    void chen(int q);
    void dey(int q);
    void dey_sampling(int q, int k = 1);
    void dey_q1();
    void tb_old(int q);
    void tb(int q);
    void write(const CubComplex &K) const;
    void write_obj(const CubComplex &K) const;
    void statistics() const;
    void remove_sealed_holes(Chain &criticals, const HDVF &X, const Filtration &F) const;
    int min_radius_generator(const Filtration &F, const Chain &criticals) const;
    std::vector<int> sorted_criticals(const CubComplex &K, const Filtration &F, const HDVF &X, int q) const;
    Chain sum_chains(const Chain &x1, const Chain &x2) const;
    void add_cycles_from_point(const CubComplex &K, int q, int cube, std::vector<Chain> &cycles) const;
    std::vector<Chain> shortest_cycles(const CubComplex &K, int c) const;


private:
    const DigObject &m_DO;
    const Options &m_opt;
    std::vector<Chain> m_basis;
};

#endif // SHORTBASIS_H
