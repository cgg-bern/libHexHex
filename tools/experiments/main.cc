#include <HexHex/HexHex.hh>
#include <HexHex/Utils/FileAccessor.hh>
#include <HexHex/Utils/csvfile.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/TestMeshes.hh>


namespace HexHex::Experiments
{

static const std::vector<std::string> ALL_MODELS =
{
    "i01c_m1_new_hex_igm",
    "i01u_m1_new_hex_igm",
    "i02c_m2_new_hex_igm",
    "i02u_m2_new_hex_igm",
    "i06c_m6_new_hex_igm",
    "i06u_m6_new_hex_igm",
    "i09c_m9_new_hex_igm",
    "i09u_m9_new_hex_igm",
    "i10c_simp_new_hex_igm",
    "i10u_simp_new_hex_igm",
    "i11c_s1_new_hex_igm",
    "i11u_s1_new_hex_igm",
    "i12c_s5_new_hex_igm",
    "i12u_s5_new_hex_igm",
    "i13c_s6_new_hex_igm",
    "i14b_s7_new_hex_igm",
    "i14c_s7_new_hex_igm",
    "i14u_s7_new_hex_igm",
    "i15c_s8_new_hex_igm",
    "i15u_s8_new_hex_igm",
    "i18b_s22_new_hex_igm",
    "i18c_s22_new_hex_igm",
    "i18u_s22_new_hex_igm",
    "i20c_s25_new_hex_igm",
    "i22c_s27_new_hex_igm",
    "i22u_s27_new_hex_igm",
    "i23c_s31_new_hex_igm",
    "i23u_s31_new_hex_igm",
    "i25c_s40_new_hex_igm",
    "i25u_s40_new_hex_igm",
    "n02b_skijump_anti_box_cyl_new_hex_igm",
    "n02c_skijump_anti_box_cyl_new_hex_igm",
    "n02u_skijump_anti_box_cyl_new_hex_igm",
    "n03b_skijump_box_cyl_new_hex_igm",
    "n03c_skijump_box_cyl_new_hex_igm",
    "n03u_skijump_box_cyl_new_hex_igm",
    "n04b_transition_prism_new_hex_igm",
    "n04c_transition_prism_new_hex_igm",
    "n04u_transition_prism_new_hex_igm",
    "n05c_box_min_pcyls_new_hex_igm",
    "n05u_box_min_pcyls_new_hex_igm",
    "n06c_anti_pentapyr_new_hex_igm",
    "n06u_anti_pentapyr_new_hex_igm",
    "n07c_anti_pyramid_new_hex_igm",
    "n07u_anti_pyramid_new_hex_igm",
    "n08c_pentapyr_new_hex_igm",
    "n08u_pentapyr_new_hex_igm",
    "n09c_pyramid_new_hex_igm",
    "n09u_pyramid_new_hex_igm",
    "n10c_qtorus_cyl_new_hex_igm",
    "n10u_qtorus_cyl_new_hex_igm",
    "n12b_limit_cycle_genus0_new_hex_igm",
    "n12c_limit_cycle_genus0_new_hex_igm",
    "n12u_limit_cycle_genus0_new_hex_igm",
    "n13b_acute_line_new_hex_igm",
    "n13c_acute_line_new_hex_igm",
    "n13u_acute_line_new_hex_igm",
    "s01b_cube_new_hex_igm",
    "s01c_cube_new_hex_igm",
    "s01u_cube_new_hex_igm",
    "s02b_prism_new_hex_igm",
    "s02c_prism_new_hex_igm",
    "s02u_prism_new_hex_igm",
    "s03b_pentagon_new_hex_igm",
    "s03c_pentagon_new_hex_igm",
    "s03u_pentagon_new_hex_igm",
    "s04b_tetrahedron_new_hex_igm",
    "s04c_tetrahedron_new_hex_igm",
    "s04u_tetrahedron_new_hex_igm",
    "s05c_cube_sphere_new_hex_igm",
    "s05u_cube_sphere_new_hex_igm",
    "s06b_dodecahedron_new_hex_igm",
    "s06c_dodecahedron_new_hex_igm",
    "s07b_notch_new_hex_igm",
    "s07c_notch_new_hex_igm",
    "s07u_notch_new_hex_igm",
    "s08b_cross_cyls_dr_new_hex_igm",
    "s08c_cross_cyls_dr_new_hex_igm",
    "s08u_cross_cyls_dr_new_hex_igm",
    "s09b_bridge_new_hex_igm",
    "s09c_bridge_new_hex_igm",
    "s09u_bridge_new_hex_igm",
    "s10c_cyl_cutsphere_new_hex_igm",
    "s10u_cyl_cutsphere_new_hex_igm",
    "s11b_cube_cyl_new_hex_igm",
    "s11c_cube_cyl_new_hex_igm",
    "s11u_cube_cyl_new_hex_igm",
    "s12b_cube_rounded_1_new_hex_igm",
    "s12c_cube_rounded_1_new_hex_igm",
    "s12u_cube_rounded_1_new_hex_igm",
    "s13b_cube_rounded_2_new_hex_igm",
    "s13c_cube_rounded_2_new_hex_igm",
    "s13u_cube_rounded_2_new_hex_igm",
    "s14b_cube_corner_sub_sphere_new_hex_igm",
    "s14c_cube_corner_sub_sphere_new_hex_igm",
    "s14u_cube_corner_sub_sphere_new_hex_igm",
    "s15b_cylinder_new_hex_igm",
    "s15c_cylinder_new_hex_igm",
    "s15u_cylinder_new_hex_igm",
    "s16c_torus_new_hex_igm",
    "s17b_sphere_new_hex_igm",
    "s17c_sphere_new_hex_igm",
    "s17u_sphere_new_hex_igm"
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

        auto maybe_inputmesh = loadInputFromFile(ovmb_dir + e.in_file + ".ovmb");
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
