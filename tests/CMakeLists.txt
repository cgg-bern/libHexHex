if (NOT DEFINED HEXEX_INCLUDE_DIRS)
  set(HEXEX_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/src/")
endif()

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${HEXEX_INCLUDE_DIRS}
  )

message("HEXEX_INCLUDE_DIRS: ${HEXEX_INCLUDE_DIRS}")

if (GTEST_FOUND)

    enable_testing()

    # Set correct include paths so that the compiler can find the headers
    include_directories(${GTEST_INCLUDE_DIRS})
    # set additional link directories
    link_directories(${GTEST_LIBRARY_DIR} )

    if ( CMAKE_GENERATOR MATCHES "^Visual Studio 11.*" )
      add_definitions( /D _VARIADIC_MAX=10 )
    endif()

    FILE(GLOB UNITTEST_SRC *.cc)

    # Create new target named unittests
    # Create unittest executable
    add_executable(unittests ${UNITTEST_SRC})

    target_include_directories(unittests PRIVATE ${HEXEX_INCLUDE_DIRS})

    # For the unittest we don't want the install rpath as set by acg_add_executable
    set_target_properties ( unittests PROPERTIES  BUILD_WITH_INSTALL_RPATH 0 )

    # Set output directory to ${BINARY_DIR}/Build
    set (OUTPUT_DIR "${CMAKE_BINARY_DIR}/Build/${ACG_PROJECT_BINDIR}")
    set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    foreach(CONFIG ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${CONFIG} UPCONFIG)
      set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${UPCONFIG} ${OUTPUT_DIR})
    endforeach()

    if ( NOT WIN32)
      # Link against all necessary libraries
      target_link_libraries(unittests HexEx ${LIBRARIES} ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} pthread)
    else()
      target_link_libraries(unittests HexEx ${LIBRARIES} ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES})
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

	if ( NOT WIN32 )
      #copy test files to unittest output dir
      #acg_copy_after_build(unittests ${CMAKE_CURRENT_LIST_DIR}/testdata ${OUTPUT_DIR}/testdata)
      file(MAKE_DIRECTORY ${OUTPUT_DIR}/Results)
	else()
      #copy test files to Build dir because in unittest output dir dlls are missing
      #acg_copy_after_build(unittests ${CMAKE_CURRENT_LIST_DIR}/testdata ${CMAKE_BINARY_DIR}/Build/testdata)
      file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Build/Results)
	endif()

  add_test(NAME AllTestsIn_HexEx_tests WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Unittests" COMMAND "${CMAKE_BINARY_DIR}/Unittests/unittests")

else(GTEST_FOUND)
    message(STATUS "Google testing framework was not found, unittests disabled.")
endif(GTEST_FOUND)
