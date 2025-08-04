cmake_minimum_required(VERSION 3.14)


add_executable (minimum_example main.cc)
target_link_libraries(minimum_example OpenVolumeMesh::OpenVolumeMesh HexHex::HexHex)

# disable OpenVolumeMesh debug output:
target_compile_definitions(minimum_example PUBLIC NDEBUG)
