#include <HexHex/HexHex.hh>
#include <HexHex/Utils/FileAccessor.hh>
#include <HexHex/Utils/csvfile.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/TestMeshes.hh>


namespace HexHex::Experiments
{

static const std::vector<std::string> ALL_MODELS =
{
    "i01c", "i01u", "i02c", "i02u", "i06c", "i06u", "i09c", "i09u", "i10c", "i10u",
    "i11c", "i11u", "i12c", "i12u", "i13c", "i14b", "i14c", "i14u", "i15c", "i15u",
    "i18b", "i18c", "i18u", "i20c", "i22c", "i22u", "i23c", "i23u", "i25c", "i25u",
    "n02b", "n02c", "n02u", "n03b", "n03c", "n03u", "n04b", "n04c", "n04u",
    "n05c", "n05u", "n06c", "n06u", "n07c", "n07u", "n08c", "n08u", "n09c", "n09u",
    "n10c", "n10u", "n12b", "n12c", "n12u", "n13b", "n13c", "n13u",
    "s01b", "s01c", "s01u", "s02b", "s02c", "s02u", "s03b", "s03c", "s03u",
    "s04b", "s04c", "s04u", "s05c", "s05u", "s06b", "s06c", "s07b", "s07c", "s07u",
    "s08b", "s08c", "s08u", "s09b", "s09c", "s09u", "s10c", "s10u",
    "s11b", "s11c", "s11u", "s12b", "s12c", "s12u", "s13b", "s13c", "s13u",
    "s14b", "s14c", "s14u", "s15b", "s15c", "s15u", "s16c", "s17b", "s17c", "s17u"
};


struct Experiment
{
    std::string in_file;
    std::string features_file = "";
    Config config;
};

// Run each model single and multithreaded for different scales
std::vector<Experiment> getBaseExperiments()
{
    std::vector<Experiment> exps;

    for (const auto& model : ALL_MODELS)
    {
        for (int nthreads : {1,8})
        {
            for (int scale : {1,2,4})
            {
                Config config;
                config.num_threads = nthreads;
                config.igm_scaling_factor = scale;

                exps.push_back(Experiment{
                    .in_file=model,
                    .features_file = model,
                    .config = config});
            }
        }
    }

    return exps;
}

std::vector<Experiment> getAlgoHexExperiments()
{
    std::vector<Experiment> exps;

    for (const auto& model : {
             "Rebuttal-AlgoHex-IGM-i02c.hexex",
             "Rebuttal-AlgoHex-IGM-i25u.hexex",
             "Rebuttal-AlgoHex-IGM-s09u.hexex",
             "Rebuttal-AlgoHex-IGM-s17c.hexex"
    })
    {
        Config config;
        config.num_threads = 1;
        config.igm_scaling_factor = 1;

        exps.push_back(Experiment{
                                  .in_file=model,
                                  .features_file = model,
                                  .config = config});
    }

    return exps;
}

std::vector<Experiment> getSphereScaleExperiments()
{
    std::vector<Experiment> exps;
    for (int nthreads : {1,8})
    {
        for (int scale = 1; scale <= 16; ++scale)
        {
            Config config;
            config.num_threads = nthreads;
            config.igm_scaling_factor = scale;

            exps.push_back(Experiment{
              .in_file="s17c",
              .features_file = "s17c",
              .config = config});
        }
    }
    return exps;
}

// Run each model single and multithreaded for different scales, with ple and ple + plf
std::vector<Experiment> getPiecewiseLinearExperiments()
{
    std::vector<Experiment> exps;

    for (const auto& model : ALL_MODELS)
    {
        exps.push_back(Experiment{model, model, Config{
            .extract_piecewise_linear_edges=true, .num_threads=8, .igm_scaling_factor=1}});
        exps.push_back(Experiment{model, model, Config{
            .extract_piecewise_linear_edges=true, .extract_piecewise_linear_faces=true, .num_threads=8, .igm_scaling_factor=1}});
    }

    return exps;
}

std::vector<Experiment> getRasterizationEpsilonExperiments()
{
    std::vector<Experiment> exps;
    std::vector<double> epsilons = {
        10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,
        1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,
        1e-16,1e-17,1e-18,1e-19,1e-20,1e-21,1e-22,0
    };

    for (const auto& model : ALL_MODELS)
    {
        for (int nthreads : {8})
        {
            for (int scale : {1,2,4})
            {
                for (double epsilon : epsilons)
                {
                    Config config;
                    config.num_threads = nthreads;
                    config.igm_scaling_factor = scale;
                    config.rasterization_epsilon = epsilon;

                    exps.push_back(Experiment{
                      .in_file=model,
                      .features_file = model,
                      .config = config});
                }
            }
        }
    }

    return exps;
}

std::vector<Experiment> getRasterizationEpsilonExperiments2()
{
    std::vector<Experiment> exps;
    std::vector<double> epsilons = {
        10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,
        1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,
        1e-16,1e-17,1e-18,1e-19,1e-20,1e-21,1e-22,0
    };

    for (const auto& model : {"n07u", "n08c"})
    {
        for (int nthreads : {8})
        {
            for (int scale : {1,2,4})
            {
                for (double epsilon : epsilons)
                {
                    Config config;
                    config.num_threads = nthreads;
                    config.igm_scaling_factor = scale;
                    config.rasterization_epsilon = epsilon;

                    exps.push_back(Experiment{
                                              .in_file=model,
                                              .features_file = model,
                                              .config = config});
                }
            }
        }
    }

    return exps;
}

std::vector<Experiment> getBoundingBoxExperiments()
{
    std::vector<Experiment> exps;
    for (int nthreads : {1,8})
    {
        for (int scale = 1; scale <= 16; ++scale)
        {
            Config config;
            config.num_threads = nthreads;
            config.igm_scaling_factor = scale;

            exps.push_back(Experiment{
                                      .in_file="s17c",
                                      .features_file = "s17c",
                                      .config = config});
        }
    }
    return exps;
}

// Run each model single and multithreaded, recomputing valences and transitions
std::vector<Experiment> getPreprocessingRecomputationExperiments()
{
    std::vector<Experiment> exps;

    for (const auto& model : ALL_MODELS)
    {
        Config config;
        config.num_threads = 1;
        config.igm_scaling_factor = 1;
        config.hf_transition_prop_name = "";
        config.e_valence_prop_name = "";
        exps.push_back(Experiment{model, model, config});
    }

    return exps;
}



}

int main(int argc, char**argv)
{
    using namespace HexHex::Experiments;
    using namespace HexHex;

    // TetrahedralMesh cube;
    // IntegerGridMap igm;
    // createCube(cube, igm, 1, 1);
    // HexExtractor he(cube, igm);
    // he.saveInputTetMesh("/Users/tobiaskohler/Downloads/Cube.ovmb");
    // return 0;

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <ovmb_dir> <out.json>" << std::endl;
        return 1;
    }
    std::string ovmb_dir = argv[1]; // "/Users/tobiaskohler/Uni/HexHex/dataset/tet-ovmb/";
    std::string json_file = argv[2]; // "/Users/tobiaskohler/Uni/HexHex/evaluation/reports/hex2/TEMP.json";
    std::vector<nlohmann::json> jsons;

    std::vector<Experiment> experiments = getRasterizationEpsilonExperiments();

    // Experiment e;
    // e.in_file = "i09u";
    // e.config.igm_scaling_factor = 1;
    // e.config.rasterization_epsilon = 1e-17;
    // e.config.num_threads = 8;
    // e.config.recompute_valences_if_found = true;
    // e.config.recompute_transitions_if_found = true;
    // e.config.extract_piecewise_linear_faces = true;
    // experiments.push_back(e);

    size_t n_runs = experiments.size();
    size_t n_failures = 0;
    double max_failure_epsilon = -1;

    for (const auto& e : experiments)
    {
        double mb;

        std::cout << e.in_file << ", t = " << e.config.num_threads
                  << ", s = " << e.config.igm_scaling_factor
                  << ", epsilon = " << e.config.rasterization_epsilon
                  << std::endl;

        auto maybe_inputmesh = loadInputFromFile(ovmb_dir + e.in_file);
        if (!maybe_inputmesh.has_value()) {
            std::cerr << "Failed to load input mesh " << e.in_file << std::endl;
            return 1;
        }

        Config config = e.config;
        config.verbose = false;

        HexHexOutput res;
        {
            ScopedStopWatch _{sw::root};
            sw::root.reset();

            res = extractHexMesh(maybe_inputmesh->mesh, maybe_inputmesh->igm, config);
            if (!res.success) {
                n_failures += 1;
                max_failure_epsilon = std::max(max_failure_epsilon, config.rasterization_epsilon);
            }
        }

        // Store the json
        nlohmann::json j = res.report;
        //std::cout << std::setw(4) << j << std::endl;
        j["in-file"] = e.in_file;
        j["out-file"] = "";
        jsons.emplace_back(j);

    }

    // Write the final json
    std::ofstream o(json_file);
    nlohmann::json j;
    j["experiments"] = jsons;
    o << std::setw(4) << j << std::endl;
    o.close();

    std::cerr << n_failures << "/" << n_runs << " Failures with max epsilon = " << max_failure_epsilon << std::endl;


    return 0;

}
