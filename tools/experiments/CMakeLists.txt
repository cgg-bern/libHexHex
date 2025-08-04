
cmake_minimum_required(VERSION 3.14)

add_executable (experiments main.cc)
target_link_libraries(experiments OpenVolumeMesh::OpenVolumeMesh HexHex::HexHex)

# disable OpenVolumeMesh debug output:
target_compile_definitions(experiments PUBLIC NDEBUG)


add_executable (benchmark_refinement benchmark_refinement.cc)
target_link_libraries(benchmark_refinement OpenVolumeMesh::OpenVolumeMesh HexHex::HexHex)

add_executable (rasterization_epsilon rasterization_epsilon.cc)
target_link_libraries(rasterization_epsilon OpenVolumeMesh::OpenVolumeMesh HexHex::HexHex)

