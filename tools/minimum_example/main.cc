#include <HexHex/HexHex.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/TestMeshes.hh>

int main()
{
    // Create a Tet Mesh Cube consisting of T*T*T*6 Tets with parameters resulting in a Hex Mesh Cube with H*H*H Hexes
    int T = 2, H = 3;
    auto test_cube = HexHex::createCube(T, H);

    // Initialize the Hex Extractor
    HexHex::HexExtractor he(test_cube.mesh, test_cube.igm);
    {
    HexHex::ScopedStopWatch _{HexHex::sw::root};

    // Extract the Hex Mesh
    he.extract(HexHex::Config{
        .extract_piecewise_linear_edges = false,
        .extract_piecewise_linear_faces = false,
        .verbose = true,
        .num_threads = 1,
        .igm_scaling_factor = 1
    });
    }
    std::cout << HexHex::sw::root << std::endl;
    return 0;

}
