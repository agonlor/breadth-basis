#include "digobject.h"
#include "shortbasis.h"

#include <string>

void help(char **argv)
{
    std::cout << " == Short Homology Basis ==" << std::endl;
    std::cout << " Compute a short homology basis for a 3D digital object" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: " << argv[0] << " <input file name> [OPTIONS]" << std::endl;
    std::cout << "  <input file name>: a pgm file of type P2 (ASCII)" << std::endl;
    std::cout << "OPTIONS" << std::endl;
//    std::cout << "  --chen: use Chen and Freedman's algorithm (too long)" << std::endl;
    std::cout << "  --dey: use sampled version of Dey et al's algorithm (too long)" << std::endl;
    std::cout << "  --seed s: random seed" << std::endl;
//    std::cout << "  --geodist: use geodesic distance filtration (default)" << std::endl;
//    std::cout << "  --verdist: use vertex distance filtration" << std::endl;
//    std::cout << "  --voxdist: use voxel distance filtration" << std::endl;
//    std::cout << "  --benchmark: output running times" << std::endl;
}

void usage(int argc, char **argv)
{
    std::cout << "Usage: " << argv[0] << " <input file name> [OPTIONS]" << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "Example: " << argv[0] << " volume.pgm" << std::endl;
    std::cout << "More details: " << argv[0] << " -h" << std::endl;
}


Options read_options(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg.compare("--benchmark") == 0)
            options.benchmark = true;
        else if (arg.compare("--chen") == 0)
            options.algorithm = Algorithm::Chen;
        else if (arg.compare("--dey") == 0)
            options.algorithm = Algorithm::Dey;
        else if (arg.compare("--deyq1") == 0)
            options.algorithm = Algorithm::Deyq1;
        else if (arg.compare("--verdist") == 0)
            options.algorithm = Algorithm::VertexDistance;
        else if (arg.compare("--voxdist") == 0)
            options.algorithm = Algorithm::VoxelDistance;
        else if (arg.compare("--geodist") == 0)
            options.algorithm = Algorithm::Geodesic;
        else if (arg.compare("--seed") == 0 || arg.compare("-s") == 0)
        {
            if (i+1 < argc)
                options.seed = std::stoi(argv[++i]);
            else { std::cerr << "You must give an integer." << std::endl; exit(EXIT_FAILURE); }
        }
        else if (arg.compare("--help") == 0 || arg.compare("-h") == 0)
        {
            help(argv);
            exit(EXIT_SUCCESS);
        }
        else
        {
            options.filename = arg;
            options.filename.erase(options.filename.end()-4, options.filename.end()); // filename is argvc[1] without ".pgm"
        }
    }
    return options;
}


int main(int argc, char **argv)
{
    if (argc <  2)
    {
        usage(argc, argv);
        return 0;
    }
    const Options options = read_options(argc, argv);

    DigObject DO(options.filename);
    if (options.benchmark)
        DO.set_benchmark();
    DO.set_from_pgm();

    ShortBasis basis(DO, options);

    return EXIT_SUCCESS;
}
