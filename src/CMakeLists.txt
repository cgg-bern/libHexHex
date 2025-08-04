cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
cmake_policy(VERSION 3.25.0)

project (HexHex
    VERSION 0.5.0
    LANGUAGES CXX C
    HOMEPAGE_URL http://github.com/cgg-bern/libHexHex
    )

add_library (HexHex
    HexHex/HexHex.hh
    HexHex/HexHex.cc
    HexHex/c_api.h
    HexHex/c_api.cc
    HexHex/HexExtractor.hh
    HexHex/HexExtractor.cc
    HexHex/Config.cc

    HexHex/Utils/FileAccessor.hh
    HexHex/Utils/FileAccessor.cc

    HexHex/Utils/Stopwatches.hh
    HexHex/Utils/Typedefs.hh
    HexHex/Utils/Utils.hh
    HexHex/Utils/Utils.cc
    HexHex/Utils/TestMeshes.hh
    HexHex/Utils/csvfile.hh


    HexHex/Predicates/ExactPredicates.hh
    HexHex/Predicates/ExactPredicates.cc
    HexHex/Predicates/predicates.c

    HexHex/Commons/Direction.hh
    HexHex/Commons/Matrix4x4T.hh
    HexHex/Commons/MeshElement.hh
    HexHex/Commons/Permutation.hh
    HexHex/Commons/Transition.hh
    HexHex/Commons/Transition.cc
    HexHex/Commons/FlatMap.hh
    HexHex/Preprocessing/Parametrization.hh
    HexHex/Preprocessing/Parametrization.cc

    HexHex/VertexExtraction/HexVertexExtractor.hh
    HexHex/VertexExtraction/HexVertexExtractor.cc

    HexHex/Preprocessing/HexHexPreprocessor.hh
    HexHex/Preprocessing/HexHexPreprocessor.cc


    HexHex/PiecewiseLinear/PLMesh.cc
    HexHex/PiecewiseLinear/PLHexFaceExtractor.hh
    HexHex/PiecewiseLinear/PLHexFaceExtractor.cc
    HexHex/PiecewiseLinear/halfface_traversal_check.cc

    HexHex/LocalTopology/Propellers.hh
    HexHex/LocalTopology/ConstantIndices.hh
    HexHex/LocalTopology/HexVertexGenerator.hh
    HexHex/LocalTopology/HexVertexGenerator.cc
    HexHex/LocalTopology/Propellers.cc
    HexHex/GlobalTopology/HexEdgeExtractor.hh
    HexHex/GlobalTopology/HexEdgeExtractor.cc
    HexHex/GlobalTopology/HexCellExtractor.hh
    HexHex/GlobalTopology/HexCellExtractor.cc
    HexHex/Utils/HexHexError.hh
    HexHex/VertexExtraction/epsilon_geometry.hh
    HexHex/Config.hh
    HexHex/PiecewiseLinear/halfface_traversal_check.hh

    HexHex/Commons/MeshElement.cc
    HexHex/PiecewiseLinear/PLMesh.hh
    HexHex/Utils/TetMeshCache.hh
    HexHex/Utils/TetMeshCache.cc
    HexHex/GlobalTopology/cell_connectivity.hh
    HexHex/GlobalTopology/cell_connectivity.cc
    HexHex/Utils/Threading.hh
    HexHex/Predicates/GeneratorPredicates.hh
    HexHex/Predicates/GeneratorPredicates.cc
    HexHex/Predicates/TracingPredicates.hh
    HexHex/Predicates/TracingPredicates.cc
    HexHex/Commons/CellIGM.hh
    HexHex/Utils/TestMeshes.cc
)
add_library (HexHex::HexHex ALIAS HexHex)

target_link_libraries (HexHex PUBLIC
    OpenVolumeMesh::OpenVolumeMesh
    Timekeeper::libTimekeeper
    nlohmann_json::nlohmann_json
    absl::flat_hash_map
    absl::flat_hash_set
    absl::hash
)


if (TARGET OpenMP::OpenMP_CXX)
    target_link_libraries (HexHex PUBLIC OpenMP::OpenMP_CXX)
endif()

if (TARGET PkgConfig::jemalloc)
    message("Building with jemalloc.")
    target_link_libraries (HexHex PUBLIC PkgConfig::jemalloc)
    target_sources(HexHex PRIVATE HexHex/Utils/MemoryUsage/MemoryUsage_jemalloc.cc)
else()
    message("Building WITHOUT jemalloc.")
    target_sources(HexHex PRIVATE HexHex/Utils/MemoryUsage/MemoryUsage_dummy.cc)
endif()

include(GenerateExportHeader)
generate_export_header(HexHex
    BASE_NAME HEXHEX
    EXPORT_FILE_NAME "HexHex/Config/Export.hh"
    DEPRECATED_MACRO_NAME "CMAKE_HEXHEX_DEPRECATED"
    NO_DEPRECATED_MACRO_NAME "CMAKE_HEXHEX_NO_DEPRECATED"
)

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/HexHex/Config/Version.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/HexHex/Config/Version.hh"
)



target_include_directories(HexHex
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties (HexHex PROPERTIES
    VERSION ${HexHex_VERSION_MAJOR}.${HexHex_VERSION_MINOR}
    SOVERSION ${HexHex_VERSION_MAJOR}.${HexHex_VERSION_MINOR}
    POSITION_INDEPENDENT_CODE ON
    AUTOMOC OFF
    AUTOUIC OFF
    CXX_VISIBILITY_PRESET hidden
    )
target_compile_features(HexHex PUBLIC cxx_std_20)


target_compile_definitions(HexHex
    PUBLIC
        "$<$<CXX_COMPILER_ID:MSVC>:NOMINMAX>"
        "$<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>"
    )

target_compile_options(HexHex
    PUBLIC
    "$<$<CXX_COMPILER_ID:MSVC>:/bigobj>"
    "$<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>" # Enable updated __cplusplus macro

    PRIVATE
    "$<$<CXX_COMPILER_ID:GNU>:-Wall>"
    "$<$<CXX_COMPILER_ID:Clang>:-Wall>"
    "$<$<CXX_COMPILER_ID:AppleClang>:-Wall>"
    )

set_source_files_properties(HexHex/Predicates/predicates.c
    PROPERTIES
    # make "if MSVC" work
    COMPILE_DEFINITIONS "$<$<CXX_COMPILER_ID:MSVC>:MSVC=1>"
# Shewchuk's predicates.c use some old-style C, we don't need warnings about it:
    COMPILE_FLAGS "$<$<CXX_COMPILER_ID:MSVC>:/wd4131>"
    #COMPILE_FLAGS "$<$<CXX_COMPILER_ID:GNU>:-Wno-deprecated-non-prototype>"
    COMPILE_FLAGS "$<$<CXX_COMPILER_ID:Clang>:-Wno-deprecated-non-prototype -Wno-incompatible-pointer-types-discards-qualifiers>"
    COMPILE_FLAGS "$<$<CXX_COMPILER_ID:AppleClang>:-Wno-deprecated-non-prototype -Wno-incompatible-pointer-types-discards-qualifiers>"
    LANGUAGE C
    )

#
# In order for the exact predicates to work the compiler
# must not generate x87 FPU code as this leads to the use
# of extended precision registers which lead to
# wrong results.
#
# As SSE does not have extended precision registers,
# forcing the generation of SSE code ensures that the
# exact predicates produce correct results.
#

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES ".*Clang|GNU")
    if (("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "amd64|.*86") )
       target_compile_options(HexHex PUBLIC -msse -mfpmath=sse)
    elseif (("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "arm64") )
    else()
    message (WARNING
        "You are compiling for a ${CMAKE_SYSTEM_PROCESSOR} system."
        "If this is *not* x86, everything is fine."
        "If it is x86 based, CMakeLists.txt must be updated "
        "to specify the correct compiler flags.")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # On 64-bit platforms, MSVC is supposed to only generate x87 code for
    #   'long double' types, not for doubles.
    # For 32-bit platforms, there is code in exactinit() that should
    #   make sure that no x87 extended internal precision is used.
else ()
    message (WARNING
        "You are using an unknown compiler which might generate x87 FPU code "
        "that breaks the exact predicates. If you know how to detect this compiler "
        "and which flags "
        "ensure that this compiler produces SSE code, please submit a pull request :)")
endif ()
