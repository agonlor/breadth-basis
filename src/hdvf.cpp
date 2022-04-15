#include "hdvf.h"

#include <sstream>      // std::ostringstream
#include <cassert>


/**
 * @brief Initialize a HDVF from a cubical complex
 */
HDVF::HDVF(const CubComplex &K)
    : m_K(K)
{
    assert(m_partition.empty() && m_Dc.empty() && m_Dc.empty());
    m_Dc.resize(3 + 1);
    m_Dr.resize(3 + 1);

    /* Initialize the reduced boundary */
    for (int c = 0; c < m_K.size_space(); c++)
    {
        if (m_K.at(c))
        {
            m_partition[c] = CellType::Critical;
            const int q = m_K.dim(c);
            if (q > 0)
            {
                const Chain d_c = m_K.d(c);
                m_Dc[q][c] = d_c;
            }
            if (q < 3)
            {
                const Chain cod_c = m_K.cod(c);
                m_Dr[q+1][c] = cod_c;
            }
        }
    }

    /* Initialize the reduction (h, f, g) */
    assert(m_Hc.empty() && m_Fr.empty() && m_Gc.empty());
    m_Hc.resize(3 + 1);
    m_Fr.resize(3 + 1);
    m_Gc.resize(3 + 1);
}


void HDVF::set_options(const HDVF_Options &options)
{
    m_options = options;
}


/**
 * @brief ...
 * It takes a list of pairs of cubes, using the numeration of cubes instead of the faces of the CubComplex.
 * Note that here K is considered to be a cubical complex.
 * This file comes from somewhere, notably the persistent homology software
 */
void HDVF::set_from_pairing(const std::string &filename)
{
//    /* Open the file */
//    std::ifstream file(filename.c_str(), std::ios::in);
//    if (!file.good())
//    {
//        std::cerr << "Error in HDVF::setFromV(" << filename << "): impossible to open." << std::endl;
//        exit(EXIT_FAILURE);
//    }

//    std::string line;
//    while (getline(file, line))
//    {
//        std::istringstream iss(line);
//        int s, t;
//        double fs, ft;
//        iss >> s >> t >> fs >> ft;        // Read the line
//        if (file.eof()) break;   // This must be done like this. It means that if we are at the end, we stop
//        //std::cout << "s=" << s << " "; m_K.printCoord(s); std::cout << " t=" << t << " "; m_K.printCoord(t);  std::cout << std::endl; std::cout << m_K.dimCube(s) << " -- " << m_K.dimCube(t) << std::endl;
//        if (m_Kc == 0)
//        {
//            std::cerr << "Error in HDVF::setFromV(): this is not a cubical complex." << std::endl;
//            exit(EXIT_FAILURE);
//        }
//        const Cube c_s = m_Kc->cube_with_num(s);
//        const Cube c_t = m_Kc->cube_with_num(t);
//        if (m_Kc->isInComplex(c_s) && m_Kc->isInComplex(c_t))    // Sometimes V contains arrows not in the complex. There are also unpaired cells
//        {
//            if (c_t.dim() == c_s.dim() + 1)
//            {
//                const int f_s = m_Kc->get_face(c_s);
//                const int f_t = m_Kc->get_face(c_t);
//                m_V.push_back( std::pair<int,int>(f_s, f_t) ); // We put the two faces of K in V
//            }
//            else
//            {
//                std::cerr << "Error in HDVF::setFromV(" << filename << "): dim(" << t << ") != dim(" << s << ")+1" << std::endl;
//                exit(EXIT_FAILURE);
//            }
//        }
//    }
//    file.close();  // We close the file
}


/**
 * @brief Compute a HDVF until there's no 1 in the reduced boundary
 * @todo For some reason this is very slow. Why?
 */
void HDVF::homology()
{
    int s, t;
    while (find_pair(s, t))
    {
        add(s, t);
    }
}


/**
 * @brief Compute a HDVF until there's no 1 in the reduced boundary
 * Here we traverse the cells in a random way, shuffling the indices
 */
void HDVF::homology_rand()
{
//    /* Shuffle the indices */
//    const int n = m_K.get_numberCells();
//    std::vector<int> v(n);
//    for (int i = 0; i < n; i++) v[i] = i;
//    for (int i = 0; i < n-1; i++)
//    {
//        int j = rand() % (n - i);
//        std::swap(v[i], v[i+j]);
//    }
//    /* Find a perfect (?) HDVF */
//    bool idem;
//    do
//    {
//        idem = true;
//        int s, t;
//        for (int i = 0; i < n && idem; i++)  // for each cell
//        {
//            const int q = m_K.dim(v[i]);
//            if (q > 0 && m_partition[v[i]] == CellType::Critical)
//            {
//                const Chain x = get(m_Dc[q], v[i]);
//                for (auto it = x.cbegin(); it != x.cend() && idem; ++it)
//                {
//                    s = *it;
//                    t = v[i];
//                    idem = false;
//                }
//            }
//        }
//        if (!idem)
//        {
//            /* Add it */
//            idem = false;
//            add(s, t);
//        }
//    } while (!idem);
}


/**
 * @brief Search a pair of critical cells (s,t) such that we can add them to the
 * HDVF.
 * This means that <d(t), s> = 1.
 * @return true if a pair is found, false otherwise
 */
bool HDVF::find_pair(int &s, int &t) const
{
    for (int q = 1; q <= 3; q++) // for each dimension q
        for (const auto &kv : m_Dc[q])   // for each q-cell t
            if (!kv.second.empty())
            {
                t = kv.first;
                s = *(kv.second.begin());
                return true;
            }
    return false;
}



void HDVF::compute_filtration(const Filtration &F)
{
    for (std::size_t i = 0; i < F.size(); i++)
    {
        const int t = F.cube_at_pos(i);
        const int q = m_K.dim(t);
        if (q > 0)
        {
            int max_pos = -1;
            int s = -1;
            for (int c : m_Dc[q].at(t))   // for each (q-1)-cell in d(t)
            {
                if (F.pos_of_cube(c) > max_pos)
                {
                    max_pos = F.pos_of_cube(c);
                    s = c;
                }
            }
            if (s > -1)
                add(s, t);
        }
    }
}

/**
 * @brief HDVF::compute_tb_filtration
 * Compute the map f of half of this filtration plus the thickness-breadth pairs.
 * We only take dimension q
 * @param f
 *
 * @note This could go in a derived class for the problem
 */
Annotation HDVF::compute_tb_filtration(const Filtration &F, int q)
{
    std::size_t i = 0;
    // first, we compute the filtration for the negative values
    while (F.value(F.cube_at_pos(i)) < 0)
    {
        const int t = F.cube_at_pos(i);
        const int q = m_K.dim(t);
        if (q > 0)
        {
            int max_pos = -1;
            int s = -1;
            for (int c : m_Dc[q].at(t))   // for each (q-1)-cell in d(t)
            {
                if (F.pos_of_cube(c) > max_pos)
                {
                    max_pos = F.pos_of_cube(c);
                    s = c;
                }
            }
            if (s > -1)
                add(s, t);
        }
        i++;
    }
    // save the current values of the map f
    const mapset map_f = m_Fr.at(q);
    std::clog << "beta_" << q << " = " << map_f.size() << std::endl;
    Annotation annot;
    annot.dimension = q;
    // continue computing the filtration an get the thickness-breadth pairs
    // note that we stop when all the original holes are filled
    while (i < F.size() && annot.pairs.size() < map_f.size())
    {
        const int t = F.cube_at_pos(i);
        const int qt = m_K.dim(t);
        if (qt > 0)
        {
            int max_pos = -1;
            int s = -1;
            for (int c : m_Dc[qt].at(t))   // for each (q-1)-cell in d(t)
            {
                if (F.pos_of_cube(c) > max_pos)
                {
                    max_pos = F.pos_of_cube(c);
                    s = c;
                }
            }
            if (s > -1)
            {
                add(s, t);
                // if the paired critical cell is one of the critical cells of the object
                // we keep this pair of cells
                if (m_K.dim(s) == q && map_f.count(s) > 0)
                {
                    annot.pairs.push_back(std::make_pair(s, t));
                    annot.f[s] = map_f.at(s); // without the critical cell
                }
            }
        }
        i++;
    }
    return annot;
}

/**
 * @brief Add the pair (s,t) of faces of K to the HDVF
 * @param s
 * @param t
 * @param trace
 */
void HDVF::add(int s, int t)
{
//    std::clog << "Adding cells " << s << ", " << t << " to the HDVF (q=" << m_K.dim(s) << ")" << std::endl;
    assert(m_partition[s] == CellType::Critical);
    assert(m_partition[t] == CellType::Critical);
    assert(m_K.dim(s) + 1 == m_K.dim(t));
    assert(get_d(t).count(s) > 0);

    update_reduction(s, t);  // Update the reduction: H_q, F_q, G_{q+1}, D_{q+1}

    m_partition[s] = CellType::Primary; // Update the partition
    m_partition[t] = CellType::Secondary;

    // if I want to see the matrices
    if (false)
    {
        const int q = m_K.dim(s);
        print_Matrix('D', q+1);
        print_Matrix('H', q);
        print_Matrix('F', q);
        print_Matrix('G', q+1);
    }
}


/**
 * @brief Update the reduction after adding the pair of faces (s,t).
 * Sometimes a boundary is 0 but we do not remove the column from m_Dc
 * @todo Why?
 *
 * @param q is the dimension of s
 * @param s
 * @param t
 */
void HDVF::update_reduction(int s, int t)
{
    /* The submatrices that I need */
    const int q = m_K.dim(s);
    const Chain F21 = get_Fr(s);
    const Chain G12 = get_Gc(t);
    Chain D12 = get(m_Dc[q+1], t); D12.erase(s);
    Chain D21 = get(m_Dr[q+1], s); D21.erase(t);

    /* Compute D'_q */
    if (q > 0)
    {
        for (int ci : m_Dc[q][s])   // remove the col d(s)
            m_Dr[q][ci].erase(s);
        m_Dc[q].erase(s);
    }

    /* Compute D'_{q+1} */
    for (int ci : m_Dc[q+1][t])   // remove the col d(t)
        m_Dr[q+1][ci].erase(t);
    m_Dc[q+1].erase(t);
    for (int cj : m_Dr[q+1][s])   // remove the row d*(s)
        m_Dc[q+1][cj].erase(s);
    m_Dr[q+1].erase(s);
    for (int ci : D12)   // D'11
        for (int cj : D21)
        {
            add_to_chain(m_Dc[q+1][cj], ci);
            add_to_chain(m_Dr[q+1][ci], cj);
        }

    /* Compute D'_{q+2} */
    if (q+2 <= 3)
    {
        for (int cj : m_Dr[q+2][t])   // remove the row d*(t)
            m_Dc[q+2][cj].erase(t);
        m_Dr[q+2].erase(t);
    }

    /* Compute H'_q */
    if (m_options.h)
    {
        for (int ci : G12)   // H_11
            for (int cj : F21)
                add_to_chain(m_Hc[q][cj], ci);
        m_Hc[q][s] = G12;  // H_12
        for (int c : F21)   // H_21
            add_to_chain(m_Hc[q][c], t);
        add_to_chain(m_Hc[q][s], t);    // H_22
    }

    /* Compute F'_q */
    if (m_options.f)
    {
        m_Fr[q].erase(s);   // remove f*(s)
        for (int ci : D12)  // F_11
            for (int cj : F21)
                add_to_chain(m_Fr[q][ci], cj);
        for (int ci : D12)   // F_12
            add_to_chain(m_Fr[q][ci], s);
    }

    /* Compute F'_{q+1} (easy) */
    if (m_options.f)
        m_Fr[q+1].erase(t);

    /* Compute G'_q (easy) */
    if (m_options.g)
        m_Gc[q].erase(s);

    /* Compute G'_{q+1} */
    if (m_options.g)
    {
        m_Gc[q+1].erase(t);   // Remove g(tau)
        for (int ci : G12)   // G_11
            for (int cj : D21)
                add_to_chain(m_Gc[q+1][cj], ci);
        for (int cj : D21)   // G_21
            add_to_chain(m_Gc[q+1][cj], t);
    }
}


/**
 * @brief HDVF::critical_cells
 * @param q
 * @return A chain with all the critical cells of dimension q
 */
Chain HDVF::critical_cells(int q) const
{
    Chain x;
    for (const auto &pair : m_partition)
        if (pair.second == CellType::Critical && m_K.dim(pair.first) == q)
            add_to_chain(x, pair.first);
    return x;
}

/**
 * @brief Get a chain in a mapset.
 * If it isn't there, it returns an empty chain
 */
Chain HDVF::get(const mapset &m, int k) const
{
    Chain x;
    mapset::const_iterator it = m.find(k);
    if (it != m.end())
        x = it->second;
    return x;
}


/**
 * @brief Coefficient of a face in a chain
 */
int HDVF::coeff(const Chain &x, int face) const
{
    Chain::const_iterator it = x.find(face);
    if (it == x.end())
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


/**
 * @brief Add face to a chain
 */
void HDVF::add_to_chain(Chain &x, int face) const
{
    std::pair<Chain::iterator,bool> ret = x.insert(face);
    if (ret.second == false) // if x already contains that face
        x.erase(ret.first);
}


/**
 * @brief F* of the face c
 * That is, the row of the matrix F corresponding to the face c
 */
Chain HDVF::get_Fr(int c) const
{
    assert(m_partition.at(c) == CellType::Critical);

    const int q = m_K.dim(c);
    if (m_Fr[q].count(c) > 0)
        return m_Fr[q].at(c);
    return Chain();
}


/**
 * @brief G of the face c
 */
Chain HDVF::get_Gc(int c) const
{
    assert(m_partition.at(c) == CellType::Critical);

    const int q = m_K.dim(c);
    if (m_Gc[q].count(c) > 0)
        return m_Gc[q].at(c);
    return Chain();
}




/**
 * @brief Reduced boundary of the face c of K
 */
Chain HDVF::get_d(int c) const
{
    assert(m_partition.at(c) == CellType::Critical);

    const int q = m_K.dim(c);
    const Chain &x = m_Dc[q].at(c);
    return x;
}

/**
 * @brief Reduced coboundary of the face c of K
 * @param c
 * @return
 */
Chain HDVF::get_dstar(int c) const
{
    assert(m_partition.at(c) == CellType::Critical);

    const int q = m_K.dim(c);
    const Chain &x = m_Dr[q+1].at(c);
    return x;
}


Chain HDVF::get_h(int c) const
{
    if (m_partition.at(c) == CellType::Primary)
    {
        const int q = m_K.dim(c);
        return m_Hc[q].at(c);
    }
    else
    {
        return Chain();
    }
}


/**
 * @brief f* of the face "c"
 */
Chain HDVF::get_fstar(int c) const
{
    assert(m_partition.at(c) == CellType::Critical);
    Chain x = get_Fr(c);
    x.insert(c);
    return x;
}

/**
 * @brief HDVF::f
 * @param x
 * @return The chain f(x)
 */
Chain HDVF::f(const Chain &x) const
{
    Chain y; // y = f(x)
    const int q = m_K.dim(*x.cbegin());

    // first, we compute f for the primary cells. Note that, since f is defined
    // on rows, this is a bit convoluted
    for (auto it = m_Fr.at(q).cbegin(); it != m_Fr.at(q).cend(); ++it)
    {
        const Chain f_row = get_fstar(it->first);
        for (int cx : x)
            if (f_row.count(cx))
                add_to_chain(y, it->first);
    }
    return y;
}


Chain HDVF::g(int c) const
{
    Chain x = get_Gc(c);
    x.insert(c);
    return x;
}


void HDVF::write_fstar(const std::string &filename) const
{
    /* Create the file */
    std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!file)
    {
        std::cerr << "Error in HDVF::write_fstar(" << filename << "): impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    /* Add the values of f* */
    for (const auto & pair : m_partition)
    {
        if (pair.second == CellType::Critical)
        {
            const Chain x = get_fstar(pair.first);
            file << pair.first << " ";
            for (auto it = x.cbegin(); it != x.cend(); ++it)
                file << *it << " ";
            file << std::endl;
        }
    }
    file.close();
}


void HDVF::write_g(const std::string &filename) const
{
    /* Create the file */
    std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!file)
    {
        std::cerr << "Error in HDVF::write_g(" << filename << "): impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    /* Add the values of g */
    for (const auto & pair : m_partition)
    {
        if (pair.second == CellType::Critical)
        {
            const Chain x = g(pair.first);
            file << pair.first << " ";
            for (auto it = x.cbegin(); it != x.cend(); ++it)
                file << *it << " ";
            file << std::endl;
        }
    }
    file.close();
}




/**
 * @brief Tells if the cell is critical
 */
bool HDVF::isCritical(int i) const
{
    return (m_partition.at(i) == CellType::Critical);
}


/**
 * @brief Print the number of critical cells in each dimension
 */
void HDVF::numberCriticalCells() const
{
    /* get the number of critical cells */
    std::vector<int> betti(3 + 1);
    for (const auto & pair : m_partition)
    {
        if (pair.second == CellType::Critical)
            betti.at(m_K.dim(pair.first))++;
    }

    /* display it */
    std::cout << "beta = ";
    for (auto it = betti.cbegin(); it != betti.cend(); ++it)
    {
        std::cout << *it << " ";
    }
    std::cout << "__" << std::endl;
}


/**
 * @brief Print the reduction.
 * I actually don't print the coefficients, but who cares?
 */
void HDVF::print_Matrix(char c, int q) const
{
    if (c == 'H')
    {
        std::cout << "\n H_" << q << std::endl;
        for (auto it = m_Hc[q].cbegin(); it != m_Hc[q].cend(); ++it)
        {
            std::cout << "h(" << it->first << ") =";
            //                m_K.printChain(it->second);
            std::cout << std::endl;
        }
    }
    else if (c == 'F')
    {
        std::cout << "\n F_" << q << std::endl;
        for (auto it = m_Fr[q].cbegin(); it != m_Fr[q].cend(); ++it)
        {
            std::cout << "f*(" << it->first << ") =";
            //                m_K.printChain(it->second);
            std::cout << std::endl;
        }
    }
    else if (c == 'G')
    {
        std::cout << "\n G_" << q << std::endl;
        for (auto it = m_Gc[q].cbegin(); it != m_Gc[q].cend(); ++it)
        {
            std::cout << "g(" << it->first << ") =";
            //                m_K.printChain(it->second);
            std::cout << std::endl;
        }
    }
    else if (c == 'D')
    {
        std::cout << "\n D_" << q << std::endl;
        for (auto it = m_Dc[q].begin(); it != m_Dc[q].end(); ++it)
        {
            std::cout << "\td(" << it->first << ") =";
//            m_K.printChain(it->second);
            std::cout << std::endl;
        }
    }
}

/**
 * @brief Save the actual HDVF in a file.
 * It can be recovered by void HDVF_fromFile() (not written)
 *
 * @todo It isn't finished, and I do not need it right now
 */
void HDVF::save(const std::string &filename) const
{
//    /* Create the file */
//    std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
//    if (!file)
//    {
//        std::cerr << "Error in HDVF::write_fstar(" << filename << "): impossible to create the output text file." << std::endl;
//        exit(EXIT_FAILURE);
//    }

//    for (int q = 0; q < listP.size(); q++)
//    {
//        for (int i = 0; i < listP[q].size(); i++)
//        {
//            file << listP[q].at(i) << " " << listS[q+1].at(i) << std::endl;
//        }
//    }
//    file.close();
}

/**
 * @brief HDVF::homology_cycles
 * @param q
 * @return The cycles g_q(c) for each q-critical cell
 */
std::vector<Chain> HDVF::homology_cycles(int q) const
{
    if (!m_options.g)
    {
        std::cerr << "The map g has not been computed for this HDVF" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::vector<Chain> cycles;
    for (const auto &pair : m_Gc.at(q))
        cycles.push_back(g(pair.first));
    return cycles;
}


void HDVF::reset()
{
    m_partition.clear();
    m_Hc.clear();
    m_Fr.clear();
    m_Gc.clear();
    m_Dc.clear();
    m_Dr.clear();
}
