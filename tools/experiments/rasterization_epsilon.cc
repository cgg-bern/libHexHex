#include <HexHex/HexHex.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Utils/FileAccessor.hh>
#include <HexHex/Utils/csvfile.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/TestMeshes.hh>

namespace HexHex::Experiments
{
struct Experiment
{
    std::string in_file;
    Config config;
};
}

int main()
{
    using namespace HexHex::Experiments;
    using namespace HexHex;

    std::string ovmb_dir = "/Users/tobiaskohler/Downloads/ovmb/";

    csvfile csv("/Users/tobiaskohler/Downloads/hexhex-experiments-rasterization-epsilon-testrun.csv");
    csv << "Model" << "Scale" << "#T" << "#H" << "Epsilon" << "#Threads" << "Vertex Extraction" << "Total" << "Memory" << "Error" << endrow;

    std::vector<std::string> models = {
        "s01b", "s04b", "s09u", "s10u", "s17c",         "s01c", "s01u", "s02b", "s03u", "s05c",
        "n03u", "n04c", "n07c", "n10u", "n12b",         "n02b", "n05u", "n08c", "n09u", "n13c",
        "i01c", "i02c", "i09u", "i18c", "i25u",          "i06c", "i10c", "i11u", "i14c", "i23c",
    };
    models = {"s17c"};

    //std::vector<Experiment> experiments = getScalingAndHashmapAndOldVertexExtractionExperiments({"s17c"});
    std::vector<Experiment> experiments;
    for (const auto& model : models)
        for (int nthreads : {8}) for (int scale : {12})
                for (int exp = -16; exp >= -16; --exp)
    {
                    double epsilon = std::pow(10, exp);
        experiments.push_back(Experiment{.in_file=model, .config=Config{.num_threads=nthreads, .igm_scaling_factor=scale, .rasterization_epsilon=epsilon}});
    }

    for (const auto& e : experiments)
    {
        double mb;

        std::cout << e.in_file << ": eps = " << e.config.rasterization_epsilon << std::endl;

        auto maybe_inputmesh = loadInputFromFile(ovmb_dir + e.in_file);
        if (!maybe_inputmesh.has_value()) {
            std::cerr << "Failed to load input mesh " << e.in_file << std::endl;
            return 1;
        }
        HexExtractor he(maybe_inputmesh->mesh, maybe_inputmesh->igm);

        Config config = e.config;
        config.verbose = true;

        {
            ScopedStopWatch _{sw::root};
            sw::root.reset();

            // omp_set_num_threads(config.num_threads);
            // he.CONFIG = config;
            // he.preprocessParametrization();
            // he.extractHexVertices();
            he.extract(config);
        }
        std::cout << Timekeeper::HierarchicalStopWatchResult(HexHex::sw::extract) << std::endl;

        // Write Result to csv file
        csv << e.in_file
            << config.igm_scaling_factor
            << he.getInputMesh().n_cells()
            << he.getOutputMesh().n_cells()
            << config.rasterization_epsilon
            << config.num_threads
            << toSeconds(HexHex::sw::extractHexVertices.elapsed())
            << toSeconds(HexHex::sw::extract.elapsed())
            << he.get_peak_memory_usage_bytes()
            << he.hasDetectedError()
            << endrow;
    }

    return 0;

}
