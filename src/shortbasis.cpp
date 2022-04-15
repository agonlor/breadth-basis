#include "shortbasis.h"

#include "cubcomplex.h"
#include "filtration.h"
#include "hdvf.h"
#include "tbpairs.h"
#include "meshobj.h"

#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

#include <random>       // std::default_random_engine
#include <iomanip>      // std::put_time

ShortBasis::ShortBasis(const DigObject &DO, const Options &options)
    : m_DO(DO), m_opt(options)
{
    std::cout << "Computing cycles for " << m_opt.filename
              << " with algorithm " << m_opt.comment() << ":" << m_opt.seed << std::endl;
    if (m_opt.algorithm == Algorithm::Chen)
        chen(1);
    else if (m_opt.algorithm == Algorithm::Dey)
        dey_sampling(1); // dey(1);
    else if (m_opt.algorithm == Algorithm::Deyq1)
        dey_q1();
    else
    {
//        tb_old(1);
        tb(1);
    }
}

/**
 * @brief ShortBasis::chen
 * Compute the minimum homology basis (w.r.t. the radius) of dimension q of
 * a digital object DO
 * @param q
 */
void ShortBasis::chen(int q)
{
    assert(m_opt.algorithm == Algorithm::Chen);
    CubComplex K(m_DO);
    K.set_from_object();
    K.only_q_cells(q);
    // We iteratively find the minimum non-trivial cycle, add it to the basis
    // and close that hole (symbolically)
    while (true)
    {
        int min_radius = K.size_space(); // this is an upper bound
        Chain min_cycle;
        const std::vector<int> pointels = K.cubes_of_dim(0);
        // for each vertex of the complex, we set the geodesic distance as a filtration
        // and look for the first homology cycle we found
        for (std::size_t i = 0; i < pointels.size(); ++i)
        {
            std::clog << "-- Try pointel " << pointels.at(i) << " (" << 100*i/pointels.size() << "%)\r";
            Filtration F(K);
            F.geodesic(pointels.at(i));
            HDVF X(K);
            HDVF_Options opt; opt.h = false;
            X.set_options(opt);
            X.compute_filtration(F);
            Chain criticals = X.critical_cells(q);
            remove_sealed_holes(criticals, X, F);
            if (criticals.empty())
            {
                write(K);
                write_obj(K);
                statistics();
                return;
            }
            const int cube = min_radius_generator(F, criticals);
            if (F.value(cube) < min_radius)
            {
                min_radius = F.value(cube);
                min_cycle = X.g(cube);
            }
        }
        m_basis.push_back(min_cycle);
        std::clog << "cycle #" << m_basis.size() << " with radius " << min_radius << " and " << min_cycle.size() << " cubes [" << m_DO.elapsed_sec() << "s]" << std::endl;
    }
}

/**
 * @brief ShortBasis::dey
 * Compute a minimum q-dimensional homology basis (w.r.t radius) using the
 * algorithm in Dey, Li and Wang (2018)
 * @param q
 */
void ShortBasis::dey(int q)
{
    assert(m_opt.algorithm == Algorithm::Dey);
    // compute the cycles in each geodesic filtration
    CubComplex K(m_DO);
    K.set_from_object();
    K.only_q_cells(q);
    const std::vector<int> pointels = K.cubes_of_dim(0);
    std::vector<Chain> cycles;
    for (std::size_t i = 0; i < pointels.size(); ++i)
    {
        std::clog << "-- Try pointel " << pointels.at(i) << " (" << 100*i/pointels.size() << "%)\r";
        add_cycles_from_point(K, q, pointels.at(i), cycles);
    }

    // sort the cycles by length (not by radius)
    std::sort(cycles.begin(), cycles.end(), [](const Chain &x1, const Chain &x2) {
        return x1.size() < x2.size();
    });

    // compute an annotation and apply it to the cycles
    HDVF X(K);
    X.homology();
    std::vector<Chain> f_cycles(cycles.size());
    for (std::size_t i = 0; i < f_cycles.size(); ++i)
        f_cycles.at(i) = X.f(cycles.at(i));

    // extract the earliest basis
    for (std::size_t i = 0; i < f_cycles.size(); ++i)
    {
        if (!f_cycles.at(i).empty())
        {
            m_basis.push_back(cycles.at(i));
            std::clog << "cycle #" << m_basis.size() << " of length " << cycles.at(i).size() << " found [" << m_DO.elapsed_sec() << "s]" << std::endl;
            const int pivot = *(f_cycles.at(i).begin());
            for (std::size_t j = i+1; j < f_cycles.size(); ++j)
                if (f_cycles.at(j).count(pivot) > 0)
                    f_cycles.at(j) = sum_chains(f_cycles.at(j), f_cycles.at(i));
        }
    }
    std::clog << "Basis found [" << m_DO.elapsed_sec() << "s]" << std::endl;
    write(K);
    write_obj(K);
    statistics();
}

/**
 * @brief ShortBasis::dey_sampling
 * Like ShortBasis::dey, but we randomly choose k vertices
 * @param q
 */
void ShortBasis::dey_sampling(int q, int k)
{
    assert(m_opt.algorithm == Algorithm::Dey);
    assert(k > 0);

    // compute an annotation
    CubComplex K(m_DO);
    K.set_from_object();
    K.only_q_cells(q);
    HDVF X(K);
    HDVF_Options opt; opt.h = false; opt.g = false;
    X.set_options(opt);
    Filtration F(K); F.trivial();
    X.compute_filtration(F);
    std::clog << "beta_" << q << " = " << X.nb_criticals(q) << std::endl;

    // compute the cycles in each geodesic filtration
    std::vector<int> pointels = K.cubes_of_dim(0);
    std::default_random_engine engine(m_opt.seed);
    std::shuffle(pointels.begin(), pointels.end(), engine);
    std::vector<Chain> cycles;
    for (int i = 0; i < k*X.nb_criticals(q); ++i)
    {
//        std::clog << "-- Try pointel " << pointels.at(i) << " (" << 100*i/(k*X.nb_criticals(q)) << "%)\r";
        std::clog << "-- Do filtration from c=" << pointels.at(i) << " (" << 100*i/(k*X.nb_criticals(q)) << "%) [" << m_DO.elapsed_sec() << "s]" << std::endl;
        add_cycles_from_point(K, q, pointels.at(i), cycles);
    }
    std::clog << cycles.size() << " cycles found [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // sort the cycles by length (not by radius)
    std::sort(cycles.begin(), cycles.end(), [](const Chain &x1, const Chain &x2) {
        return x1.size() < x2.size();
    });

    // apply the annotation to the cycles
    std::vector<Chain> f_cycles(cycles.size());
    for (std::size_t i = 0; i < f_cycles.size(); ++i)
        f_cycles.at(i) = X.f(cycles.at(i));
    std::clog << "cycles annotated [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // extract the earliest basis
    for (std::size_t i = 0; i < f_cycles.size(); ++i)
    {
        if (!f_cycles.at(i).empty())
        {
            m_basis.push_back(cycles.at(i));
//            std::clog << "cycle #" << m_basis.size() << " of length " << cycles.at(i).size() << " found [" << m_DO.elapsed_sec() << "s]" << std::endl;
            const int pivot = *(f_cycles.at(i).begin());
            for (std::size_t j = i+1; j < f_cycles.size(); ++j)
                if (f_cycles.at(j).count(pivot) > 0)
                    f_cycles.at(j) = sum_chains(f_cycles.at(j), f_cycles.at(i));
        }
    }
    std::clog << "Basis found [" << m_DO.elapsed_sec() << "s]" << std::endl;
    write(K);
    write_obj(K);
    statistics();
}

/**
 * @brief ShortBasis::dey_q1
 * Simplified version of the algorithm by Dey, Li and Wang (2018) to compute a
 * minimal homology basis of dimension 1.
 * We compute all the shortest cycles, sort them, annotate them, and extract the
 * earliest basis.
 * @param q
 */
void ShortBasis::dey_q1()
{
    assert(m_opt.algorithm == Algorithm::Deyq1);

    // compute an annotation
    CubComplex K(m_DO);
    K.set_from_object();
    K.only_q_cells(1);
    HDVF X(K);
    HDVF_Options opt; opt.h = false; opt.g = false;
    X.set_options(opt);
    Filtration F(K); F.trivial();
    X.compute_filtration(F);
    std::clog << "annotation computed [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // compute the shortest cycles from each vertex
    const std::vector<int> pointels = K.cubes_of_dim(0);
    std::vector<Chain> cycles;
    for (std::size_t i = 0; i < pointels.size(); ++i)
    {
//        std::clog << "-- Try pointel " << pointels.at(i) << " (" << 100*i/pointels.size() << "%)[" << m_DO.elapsed_sec() << "s]\r";
        std::clog << "-- Try pointel " << pointels.at(i) << " (" << 100*i/pointels.size() << "%)\r";
        const std::vector<Chain> cur_cycles = shortest_cycles(K, pointels.at(i));
        for (const Chain &x : cur_cycles)
            if (!X.f(x).empty())
                cycles.push_back(x);
    }
    std::clog << cycles.size() << " non-trivial cycles found [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // sort the cycles by length
    std::sort(cycles.begin(), cycles.end(), [](const Chain &x1, const Chain &x2) {
        return x1.size() < x2.size();
    });

    // annotate the cycles
    std::vector<Chain> f_cycles(cycles.size());
    for (std::size_t i = 0; i < f_cycles.size(); ++i)
        f_cycles.at(i) = X.f(cycles.at(i));
    std::clog << "cycles annotated [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // extract the earliest basis
    for (std::size_t i = 0; i < f_cycles.size(); ++i)
    {
        if (!f_cycles.at(i).empty())
        {
            m_basis.push_back(cycles.at(i));
            std::clog << "cycle #" << m_basis.size() << " of length " << cycles.at(i).size() << " found [" << m_DO.elapsed_sec() << "s]" << std::endl;
            const int pivot = *(f_cycles.at(i).begin());
            for (std::size_t j = i+1; j < f_cycles.size(); ++j)
                if (f_cycles.at(j).count(pivot) > 0)
                    f_cycles.at(j) = sum_chains(f_cycles.at(j), f_cycles.at(i));
        }
    }
    std::clog << "Basis found [" << m_DO.elapsed_sec() << "s]" << std::endl;
    write(K);
    write_obj(K);
    statistics();
}

/**
 * @brief ShortBasis::tb_old
 * Compute a homology basis using the heuristic based on the TB balls
 * @param q
 */
void ShortBasis::tb_old(int q)
{
    CubComplex K(m_DO);
    K.set_from_sdt();
    K.only_q_cells(q);
    Filtration F(K);
    F.sdt(m_DO);
    std::clog << "First filtration done [" << m_DO.elapsed_sec() << "s]" << std::endl;
    HDVF X(K);
    HDVF_Options opt; opt.h = false; opt.g = false;
    X.set_options(opt);
    Annotation annotation = X.compute_tb_filtration(F, q);
    X.reset();
    std::clog << "First HDVF done [" << m_DO.elapsed_sec() << "s]" << std::endl;

    K.set_from_object();
    K.only_q_cells(q);
    for (const auto &tb_pair : annotation.pairs) // for each hole
    {
        Filtration Fd(K);
        if (m_opt.algorithm == Algorithm::VoxelDistance)
            Fd.distance_voxel(m_DO, tb_pair.second);
        else if (m_opt.algorithm == Algorithm::VertexDistance)
            Fd.distance_vertex(tb_pair.second);
        else
            Fd.geodesic_vertex(tb_pair.second);
        HDVF X(K);
        opt.g = true; opt.f = false;
        X.set_options(opt);
        X.compute_filtration(Fd);
        std::vector<int> crits = sorted_criticals(K, Fd, X, q);
        for (int s : crits)
        {
            const Chain g_s = X.g(s);
            if (annotation.coeff(g_s, tb_pair.first))
            {
                m_basis.push_back(g_s);
                std::clog << "cycle #" << m_basis.size() << " of length " << g_s.size() << " found [" << m_DO.elapsed_sec() << "s]" << std::endl;
                annotation.update(g_s, tb_pair.first);
                break;
            }
        }
    }
    std::clog << "Basis found [" << m_DO.elapsed_sec() << "s]" << std::endl;
    write(K);
    write_obj(K);
    statistics();
}

/**
 * @brief ShortBasis::tb
 * Compute a homology basis using the heuristic based on the TB balls
 * @param q
 */
void ShortBasis::tb(int q)
{
    assert(m_opt.algorithm == Algorithm::Geodesic);

    // compute the thickness-breadth pairs
    CubComplex K(m_DO);
    K.set_from_sdt();
    TBPairs Y(m_DO);
    const std::vector<std::pair<int,int>> tb_pairs = Y.compute(q);
    std::clog << tb_pairs.size() << " TB pairs computed [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // compute the cycles in each geodesic filtration
    K.set_from_object();
    K.only_q_cells(q);
    std::vector<Chain> cycles;
    int i = 0;
    for (const auto &tb_pair : tb_pairs)
    {
        std::clog << "-- Do filtration from c=" << tb_pair.second << " (" << 100*i/tb_pairs.size() << "%) [" << m_DO.elapsed_sec() << "s]" << std::endl;
        Filtration F(K); F.geodesic_vertex(tb_pair.second);
        HDVF X(K);
        HDVF_Options opt; opt.h = false; opt.f = false; X.set_options(opt);
//        LuckyReduction X(K);
        X.compute_filtration(F);
        const std::vector<Chain> cur_cycles = X.homology_cycles(q);
        for (const Chain &x : cur_cycles)
            cycles.push_back(x);
        i++;
    }
    std::clog << cycles.size() << " cycles found [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // sort the cycles by length (not by radius)
    std::sort(cycles.begin(), cycles.end(), [](const Chain &x1, const Chain &x2) {
        return x1.size() < x2.size();
    });

    // compute an annotation and apply it to the cycles
    Filtration F(K); F.trivial();
    HDVF X(K);
    HDVF_Options opt; opt.h = false; opt.g = false; X.set_options(opt);
    X.compute_filtration(F);
    std::clog << "annotation computed [" << m_DO.elapsed_sec() << "s]" << std::endl;
    std::vector<Chain> f_cycles(cycles.size());
    for (std::size_t i = 0; i < f_cycles.size(); ++i)
        f_cycles.at(i) = X.f(cycles.at(i));
    std::clog << "cycles annotated [" << m_DO.elapsed_sec() << "s]" << std::endl;

    // extract the earliest basis
    for (std::size_t i = 0; i < f_cycles.size(); ++i)
    {
        if (!f_cycles.at(i).empty())
        {
            m_basis.push_back(cycles.at(i));
//            std::clog << "cycle #" << m_basis.size() << " of length " << cycles.at(i).size() << " found [" << m_DO.elapsed_sec() << "s]" << std::endl;
            const int pivot = *(f_cycles.at(i).begin());
            for (std::size_t j = i+1; j < f_cycles.size(); ++j)
                if (f_cycles.at(j).count(pivot) > 0)
                    f_cycles.at(j) = sum_chains(f_cycles.at(j), f_cycles.at(i));
        }
    }
    std::clog << "Basis found [" << m_DO.elapsed_sec() << "s]" << std::endl;

    write(K);
    write_obj(K);
    statistics();
}

/**
 * @brief ShortBasis::write Write the short homology basis in a JSON file
 * @param filename
 * @todo Put metadata: time, method, statistics
 */
void ShortBasis::write(const CubComplex &K) const
{
    rapidjson::StringBuffer s;
    rapidjson::Writer<rapidjson::StringBuffer> writer(s);
    writer.StartObject();

    writer.Key("meta");
    writer.StartObject();
    writer.Key("object");
    writer.String((K.filename() + ".pgm").c_str());
    writer.Key("size");
    writer.StartArray();
    for (unsigned i = 0; i < 3; i++)
        writer.Uint(K.size(i));
    writer.EndArray();
    writer.Key("algorithm");
    writer.String(m_opt.comment().c_str());
    writer.Key("arguments");
    writer.StartArray();
    writer.Uint(m_opt.seed);
    writer.EndArray();
    writer.Key("elapsed_time");
    writer.Double(m_DO.elapsed_sec());
    writer.EndObject();

    writer.Key("cycles");
    writer.StartArray();
    for (const Chain &x : m_basis)
    {
        writer.StartArray();
        for (int cube : x)
        {
            const std::vector<int> coord = K.coordinates(cube);
            writer.StartObject();
            writer.Key("x");
            writer.Uint(coord.at(0));
            writer.Key("y");
            writer.Uint(coord.at(1));
            writer.Key("z");
            writer.Uint(coord.at(2));
            writer.EndObject();
        }
        writer.EndArray();
    }
    writer.EndArray();
    writer.EndObject();

    std::ostringstream oss;
    oss << m_opt.filename;
    std::time_t tt = std::chrono::system_clock::to_time_t (std::chrono::system_clock::now());
    struct std::tm * ptm = std::localtime(&tt);
    oss << "." << std::put_time(ptm, "%Y%m%d-%H%M%S");
    oss << ".json";
    std::string filename = oss.str();

    std::ofstream file(filename , std::ios::out | std::ios::trunc);
    if (!(file))
    {
        std::cerr << "Error in CubComplex::write_chains(): impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }
    file << s.GetString();
    std::clog << "Basis written to " << filename << std::endl;
}


/**
 * @brief ShortBasis::write_obj
 * Write the basis into a OBJ mesh file. This also make an OBJ file for the object.
 * I would prefer to do this with a Python script.
 * @param K
 * @param filename
 */
void ShortBasis::write_obj(const CubComplex &K) const
{
    MeshObj M(K);
    int i = 0;
    for (const Chain &x : m_basis)
    {
        M.add_chain(x, "cycle_" + std::to_string(++i));
    }
    M.write(K.filename() + "_cycles.obj");
}


void ShortBasis::statistics() const
{
    int count = 0;
    std::clog << m_basis.size() << " cycles of size: ";
    std::vector<int> lengths;
    for (const Chain &x : m_basis)
    {
        lengths.push_back(x.size());
        count += x.size();
    }
    std::sort(lengths.begin(), lengths.end());
    for (int n : lengths)
        std::clog << n << " ";
    std::clog << "[" << count << "]" << std::endl;
}



/**
 * @brief ShortBasis::remove_sealed_holes
 * @details Instead of sealing the holes in the complex, we remove the
 * critical cells of the cycles that we have added to the basis.
 * For this, we simply simulate the operator A applying f(x)
 * @param criticals
 * @param X
 * @param F
 */
void ShortBasis::remove_sealed_holes(Chain &criticals, const HDVF &X, const Filtration &F) const
{
    for (const Chain &y : m_basis)
    {
        const Chain fy = X.f(y);
        int max_value = -1;
        int max_cube;
        for (int c : fy)
        {
            if (max_value < F.value(c))
            {
                max_value = F.value(c);
                max_cube = c;
            }
        }
        criticals.erase(max_cube);
    }
}


/**
 * @brief ShortBasis::min_radius_generator
 * @param F
 * @param criticals
 * @return The cube in criticals that has the minimum value in the filtration F
 */
int ShortBasis::min_radius_generator(const Filtration &F, const Chain &criticals) const
{
    int cur_min_radius = std::numeric_limits<int>::max();
    int cube = -1;
    for (int c : criticals)
    {
        if (F.value(c) < cur_min_radius)
        {
            cur_min_radius = F.value(c);
            cube = c;
        }
    }
    return cube;
}


std::vector<int> ShortBasis::sorted_criticals(const CubComplex &K, const Filtration &F, const HDVF &X, int q) const
{
    std::vector<int> criticals;
    for (std::size_t i = 0; i < F.size(); i++)
    {
        const int c = F.cube_at_pos(i);
        if (X.type(c) == CellType::Critical && K.dim(c) == q)
            criticals.push_back(c);
    }
    return criticals;
}


Chain ShortBasis::sum_chains(const Chain &x1, const Chain &x2) const
{
    Chain y = x1;
    for (int c : x2)
    {
        std::pair<Chain::iterator,bool> ret = y.insert(c);
        if (ret.second == false) // if y already contains c
            y.erase(ret.first);
    }
    return y;
}


void ShortBasis::add_cycles_from_point(const CubComplex &K, int q, int cube, std::vector<Chain> &cycles) const
{
    Filtration F(K);
    F.geodesic(cube);
    HDVF X(K);
    HDVF_Options opt; opt.h = false; opt.f = false; X.set_options(opt);
    X.compute_filtration(F);
    const std::vector<Chain> cur_cycles = X.homology_cycles(q);
    for (const Chain &x : cur_cycles)
        cycles.push_back(x);
}


std::vector<Chain> ShortBasis::shortest_cycles(const CubComplex &K, int c) const
{
    assert(K.dim(c) == 0);
    // make a shortest path three
    std::map<int,int> edge;
    std::queue<std::pair<int,int>> fifo; // queue of pairs (vertex, previous edge)
    fifo.push(std::make_pair(c, -1));
    std::vector<bool> tree_edges(K.size_space(), false); //  tree_edges[c1] == true iff the edge is in the tree
    while (!fifo.empty())
    {
        const std::pair<int, int> cur_pair = fifo.front(); fifo.pop();
        if (edge.count(cur_pair.first) == 0)
        {
            const int c0 = cur_pair.first;
            const int c1 = cur_pair.second;
            edge[c0] = c1;
            if (c1 >= 0) // this is for the first element in the queue
                tree_edges.at(c1) = true;
            int next_c0, next_c1;
            next_c1 = K.translate(c0, -1, 0, 0);
            next_c0 = K.translate(c0, -2, 0, 0);
            if (K.at(next_c1)) fifo.push(std::make_pair(next_c0, next_c1));
            next_c1 = K.translate(c0, +1, 0, 0);
            next_c0 = K.translate(c0, +2, 0, 0);
            if (K.at(next_c1)) fifo.push(std::make_pair(next_c0, next_c1));
            next_c1 = K.translate(c0, 0, -1, 0);
            next_c0 = K.translate(c0, 0, -2, 0);
            if (K.at(next_c1)) fifo.push(std::make_pair(next_c0, next_c1));
            next_c1 = K.translate(c0, 0, +1, 0);
            next_c0 = K.translate(c0, 0, +2, 0);
            if (K.at(next_c1)) fifo.push(std::make_pair(next_c0, next_c1));
            next_c1 = K.translate(c0, 0, 0, -1);
            next_c0 = K.translate(c0, 0, 0, -2);
            if (K.at(next_c1)) fifo.push(std::make_pair(next_c0, next_c1));
            next_c1 = K.translate(c0, 0, 0, +1);
            next_c0 = K.translate(c0, 0, 0, +2);
            if (K.at(next_c1)) fifo.push(std::make_pair(next_c0, next_c1));
        }
    }

    // Get the cycles
    std::vector<Chain> cycles;
    for (std::size_t c1 = 0; c1 < tree_edges.size(); ++c1)
        if (K.at(c1) && !tree_edges.at(c1) && K.dim(c1) == 1)
        {
            Chain cycle;
            cycle.insert(c1);
            const std::vector<int> faces = K.boundary(c1);
            for (int c0 : faces)
                while (edge[c0] >= 0)
                {
                    auto ret = cycle.insert(edge[c0]);
                    if (ret.second == false) // if x already contains that face
                        cycle.erase(edge[c0]);
                    const std::vector<int> faces2 = K.boundary(edge[c0]);
                    c0 = (faces2.at(0) != c0) ? faces2.at(0) : faces2.at(1);
                }
            cycles.push_back(cycle);
        }
    return cycles;
}
