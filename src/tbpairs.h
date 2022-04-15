#ifndef TBPAIRS_H
#define TBPAIRS_H

#include "digobject.h"
#include "cubcomplex.h"

class TBPairs
{
public:
    TBPairs(const DigObject &DO);
    std::vector<std::pair<int,int>> compute(int q);

private:
    const DigObject &m_DO;
    CubComplex m_K;
};

#endif // TBPAIRS_H
