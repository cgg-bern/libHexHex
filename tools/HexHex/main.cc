#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <HexHex/Utils/FileAccessor.hh>
#include <iostream>
#include <fstream>
#include <filesystem>

#include <HexHex/HexHex.hh>
#include <nlohmann/json.hpp>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/Utils.hh>
#include <libTimekeeper/StopWatchPrinting.hh>
#include <CLI/CLI.hpp>

namespace OVM = OpenVolumeMesh;

struct Options {
    std::filesystem::path inTetFile;
    std::filesystem::path outHexFile;
    std::optional<std::filesystem::path> outPWLFile;
    std::optional<std::filesystem::path> outReportFile;
    std::optional<std::filesystem::path> inConfigFile;
    std::optional<int> igm_scaling_factor;
    std::optional<int> num_threads;
};

int main(int argc, char* argv[])
{
    CLI::App app{"HexHex: Highspeed Extraction of Hexahedral Meshes"};
    argv = app.ensure_utf8(argv);

    Options options;
    app.add_option("-i, --in", options.inTetFile, "Input file (.ovmb, .ovm, .hexex)")->required();
    app.add_option("-o, --out-hex", options.outHexFile, "Output file (.ovmb, .ovm, .mesh)")->required();
    app.add_option("--out-pwl", options.outPWLFile, "Output file for piecewise linear mesh");
    app.add_option("--report", options.outReportFile, "Output file with details about the extraction process (.json)");
    app.add_option("--config", options.inConfigFile, "Config file (.json). Used when parameters are not explicitly set.");
    app.add_option("--scale", options.igm_scaling_factor, "Parametrization scaling factor (positive integer)")->check(CLI::PositiveNumber);
    app.add_option("--nthreads", options.num_threads, "Number of threads or nonpositive to use number of available cores");

    CLI11_PARSE(app, argc, argv);

    HexHex::Config config;
    if (options.inConfigFile)
    {
        HexHex::loadConfig(options.inConfigFile.value(), config);
    }
    if (options.outPWLFile) {
        config.extract_piecewise_linear_faces = true;
        config.extract_piecewise_linear_edges = true;
    }
    if (options.igm_scaling_factor.has_value()) {
        config.igm_scaling_factor = options.igm_scaling_factor.value();
    }
    if (options.num_threads.has_value()) {
        config.num_threads = options.num_threads.value();
    }

    HexHex::TetrahedralMesh tetmesh;
    OVM::HalfFacePropertyT<HexHex::Vec3d> igm = tetmesh.request_halfface_property<HexHex::Vec3d>("HexHex::Parametrization");
    std::cout << "Load Input Tet Mesh from " << options.inTetFile << std::endl;
    if (!HexHex::loadInputFromFile(options.inTetFile, tetmesh, igm)) {
        return 1;
    }
    auto res = HexHex::extractHexMesh(tetmesh, igm, config);

    // Hex Mesh
    if (res.hex_mesh != nullptr) {
        std::cout << "Save HexHex Hex Mesh to " << options.outHexFile << std::endl;
        HexHex::saveOutputToFile(options.outHexFile, *res.hex_mesh);
    } else {
        std::cerr << "Hex extraction failed!" << std::endl;
    }

    // Piecewise Linear Mesh
    if (res.piecewise_linear_mesh != nullptr) {
        std::cout << "Save HexHex Piecewise Linear Mesh to " << *options.outPWLFile << std::endl;
        OpenVolumeMesh::IO::ovmb_write(*options.outPWLFile, *res.piecewise_linear_mesh);
    } else if (options.outPWLFile) {
        std::cerr << "Piecewise-linear extraction failed!" << std::endl;
    }

    // Report
    if (options.outReportFile)
    {
        std::cout << "Save HexHex Report to " << *options.outReportFile << std::endl;
        nlohmann::json j = res.report;
        j["tet_mesh_filename"] = options.inTetFile.string();
        j["hex_mesh_filename"] = options.outHexFile.string();
        std::ofstream o(*options.outReportFile);
        o << std::setw(4) << j << std::endl;
    }

    return 0;
}
