if (GTEST_FOUND)
    include(GoogleTest)

    FILE(GLOB UNITTEST_SRC
        cellextraction_test.cc
        common.cc
        fileaccess_test.cc
        gridisomorphism_test.cc
        navigation_test.cc
        predicate_test.cc
        unittest.cc)

    # Create unittest executable
    add_executable(unittests ${UNITTEST_SRC})
    target_link_libraries(unittests HexEx::HexEx GTest::GTest GTest::Main)
    if ( NOT WIN32)
      target_link_libraries(unittests pthread)
    endif()

    # For the unittest we don't want the install rpath as set by acg_add_executable
    set_target_properties (unittests PROPERTIES  BUILD_WITH_INSTALL_RPATH 0 )

    # Set output directory to ${BINARY_DIR}/Build
    set (OUTPUT_DIR "${CMAKE_BINARY_DIR}/Build/${ACG_PROJECT_BINDIR}")
    set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})

    foreach(CONFIG ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${CONFIG} UPCONFIG)
      set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${UPCONFIG} ${OUTPUT_DIR})
    endforeach()

    if ( CMAKE_GENERATOR MATCHES "^Visual Studio 11.*" )
        target_compile_definitions( _VARIADIC_MAX=10 )
    endif()


    if ( NOT WIN32 )
      # Set compiler flags
      if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set_target_properties(unittests PROPERTIES COMPILE_FLAGS "-g -Wno-long-long")
      else()
        set_target_properties(unittests PROPERTIES COMPILE_FLAGS "-g -pedantic -Wno-long-long")
      endif()
    else()
      # Set compiler flags
      set_target_properties(unittests PROPERTIES COMPILE_FLAGS "" )
    endif()

    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/testdata/" DESTINATION "${OUTPUT_DIR}/testdata")
    file(MAKE_DIRECTORY ${OUTPUT_DIR}/Results)

    gtest_add_tests(TARGET unittests
        EXTRA_ARGS "--gtest_color=yes;--gtest_output=xml"
        WORKING_DIRECTORY "${OUTPUT_DIR}"
        )

else(GTEST_FOUND)
    message(STATUS "Google testing framework was not found, unittests disabled.")
endif(GTEST_FOUND)
