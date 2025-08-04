cmake_minimum_required(VERSION 3.14)

add_executable (HexHex_cli main.cc)
target_link_libraries(HexHex_cli OpenVolumeMesh::OpenVolumeMesh HexHex::HexHex CLI11::CLI11)

target_compile_options(HexHex_cli PUBLIC)

set_property(TARGET HexHex_cli PROPERTY OUTPUT_NAME HexHex)
