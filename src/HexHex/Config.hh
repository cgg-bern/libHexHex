#pragma once

#include <HexHex/Config/Export.hh>
#include <cstddef>
#include <string>
#include <filesystem>

namespace HexHex {

struct Config
{
public:
    bool extract_piecewise_linear_edges = false; // Whether to extract picewise linear edges
    bool extract_piecewise_linear_faces = false; // Whether to extract picewise linear faces (also enables piecewise linear edges)
    bool verbose = true; // Whether to print some progress to the console
    std::string hf_transition_prop_name = "HalffaceTransiton"; // Whether to always recompute halfface transitions or use the preexisting AlgoHex Property
    std::string e_valence_prop_name = "edge_valance"; // Whether to always recompute edge valences or use the preexisting AlgoHex Property
    bool sanitize_igm = true; // Whether to enforce the Integer-Grid Map constraints. This is highly recommended as there are no guarantees if the IGM is not valid!
    int num_threads = 1; // Number of Threads
    int igm_scaling_factor = 1; // Resulting in scale^3 as many hex-cells
    double sanitization_epsilon = 1e-2; // Report error when trying to move vertex parameter by this much or more
    bool assert_igm_validity = true; // Test all IGM properties after the sanitization
    double rasterization_epsilon = 1e-6; // Tet Inflation
};

bool HEXHEX_EXPORT loadConfig(const std::filesystem::path& filename, Config& config);

} // namespace HexHex

