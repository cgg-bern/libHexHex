#include <HexHex/HexHex.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/TestMeshes.hh>
#include <HexHex/Commons/FlatMap.hh>
#include <fstream>
#include <filesystem>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/IO/ovmb_write.hh>


#include "test_cube.hh"


struct RefinementBenchmarkConfig {
    int scale = 1;
    int n_boxes = 1;
    std::filesystem::path filename_out_json;
    bool pwe = false;
    bool pwl = false;
    int nthreads = 1;
};
int main(int argc, char **argv)
{
    using namespace HexHex;

    RefinementBenchmarkConfig config;
    if (argc == 4) {
        config.filename_out_json = argv[1];
        config.scale = std::atoi(argv[2]);
        config.n_boxes = std::atoi(argv[3]);
    } else {
        std::cerr << "Usage: " << argv[0] << " <out.json> <scale> <n_boxes>" << std::endl;
        exit(1);
    }

    TetrahedralMesh cube;
    //auto igm = cube.request_cell_property<FlatMap<VertexHandle, Parameter>>("Parametrization");
    auto igm = cube.request_halfface_property<Parameter>("Parametrization");
    cube.set_persistent(igm);
    //createCube(cube, igm, 1, config.scale);
    mycube(cube, igm, config.scale, config.n_boxes);
#if 0
    OpenVolumeMesh::IO::ovmb_write("cube.ovmb", cube);
    std::cout << "saved" << std::endl;
#endif

    nlohmann::json report;
    {
        ScopedStopWatch _{sw::root};
        sw::root.reset();
        auto hexhex_config = HexHex::Config{
            .extract_piecewise_linear_edges = config.pwe,
                .extract_piecewise_linear_faces = config.pwl,
                .verbose = false,
                .hf_transition_prop_name = "HalffaceTransiton",
                .e_valence_prop_name = "edge_valance",
                .sanitize_igm = true,
                .num_threads = config.nthreads,
                .igm_scaling_factor = 1,
        };
        auto out = extractHexMesh(cube, igm, hexhex_config);
        report = out.report;
    }
#if 0
    he.saveOutputHexMesh("out.ovmb");
#endif

    report["scaling"] = config.scale;

    if (!config.filename_out_json.empty()) {
        std::ofstream json{config.filename_out_json};
        json << std::setw(4) << report << std::endl;
    }

    return 0;

}
