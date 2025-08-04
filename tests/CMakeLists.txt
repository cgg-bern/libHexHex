if (NOT TARGET GTest::gtest_main)
    FetchContent_Declare( googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.17.x
        )
    if(WIN32)
        # avoid linking errors, cf https://stackoverflow.com/questions/12540970/how-to-make-gtest-build-mdd-instead-of-mtd-by-default-using-cmake
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif()
    FetchContent_MakeAvailable(googletest)
endif()

include(GoogleTest)

add_executable(hexhex_unittests
    cellextraction_test.cc
    common.cc
    c_api_test.cc
    gridisomorphism_test.cc
    predicate_test.cc
    navigation_test.cc
    unittest.cc
    )
target_link_libraries(hexhex_unittests PUBLIC GTest::gtest_main)
target_link_libraries(hexhex_unittests PUBLIC HexHex::HexHex)


# For the unittest we don't want the install rpath as set by acg_add_executable
set_target_properties (hexhex_unittests PROPERTIES  BUILD_WITH_INSTALL_RPATH 0 )

# Set output directory to ${BINARY_DIR}/Build
set (OUTPUT_DIR "${CMAKE_BINARY_DIR}/Build/${ACG_PROJECT_BINDIR}")
set_target_properties(hexhex_unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})

foreach(CONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG} UPCONFIG)
    set_target_properties(hexhex_unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${UPCONFIG} ${OUTPUT_DIR})
endforeach()

if ( CMAKE_GENERATOR MATCHES "^Visual Studio 11.*" )
    target_compile_definitions( _VARIADIC_MAX=10 )
endif()

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/testdata/" DESTINATION "${OUTPUT_DIR}/testdata")
file(MAKE_DIRECTORY ${OUTPUT_DIR}/Results)

gtest_add_tests(TARGET hexhex_unittests
    EXTRA_ARGS "--gtest_color=yes;--gtest_output=xml"
    WORKING_DIRECTORY "${OUTPUT_DIR}"
    )

