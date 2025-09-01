/*
 * Copyright 2025 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */

#include <gtest/gtest.h>
#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <HexHex/Utils/FileAccessor.hh>
#include "common.hh"
#include <HexHex/HexExtractor.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

using namespace HexHex;

TEST(CellExtraction, CellExtractionTest)
{

  for (auto transition : {false, true})
    for (auto size : {1,2,3,4,5})
    {

      TetrahedralMesh mesh;
      if (transition)
        createCubeWithTransition(mesh, size);
      else
        createCube(mesh, size);

      EXPECT_EQ(8u, mesh.n_vertices());

      IGM parametrization = mesh.request_cell_property<CellIGM>("Parametrization");
      HexExtractor hexExtractor(mesh, parametrization);

      Config config;
      config.verbose = false;
      hexExtractor.extract(config);

      HexahedralMesh& hexMesh = hexExtractor.getOutputMesh();

      unsigned int expectedNumberOfCells = size*size*size;

      EXPECT_EQ(expectedNumberOfCells, hexMesh.n_cells()) << "Size: " << size;

      OpenVolumeMesh::IO::FileManager fileManager;
      fileManager.writeFile(std::string("Results/Cube_") + std::to_string(size) + ".ovm", hexMesh);

    }

}


TEST(CellExtraction, CellExtractionMasterVertexTest)
{

  TetrahedralMesh mesh;
  createMasterVertexMesh(mesh);

  IGM parametrization = mesh.request_cell_property<CellIGM>("Parametrization");
  HexExtractor hexExtractor(mesh, parametrization);

  Config config;
  config.verbose = false;
  hexExtractor.extract(config);

  HexahedralMesh& hexMesh = hexExtractor.getOutputMesh();

  unsigned int expectedNumberOfCells = 8;

  EXPECT_EQ(expectedNumberOfCells, hexMesh.n_cells());

  OpenVolumeMesh::IO::FileManager fileManager;
  fileManager.writeFile("Results/MasterVertex.ovm", hexMesh);

}

TEST(CellExtraction, CellExtractionPrismTest)
{
  TetrahedralMesh mesh;
  createPrism(mesh);
  IGM parametrization = mesh.request_cell_property<CellIGM>("Parametrization");
  HexExtractor hexExtractor(mesh, parametrization);

  Config config;
  config.verbose = false;
  hexExtractor.extract(config);

  HexahedralMesh& hexMesh = hexExtractor.getOutputMesh();

  unsigned int expectedNumberOfCells = 6;

  EXPECT_EQ(expectedNumberOfCells, hexMesh.n_cells());

  OpenVolumeMesh::IO::FileManager fileManager;
  fileManager.writeFile("Results/prism.ovm", hexMesh);

}

void test_file(const std::string& _filename, size_t _n_expected_cells)
{
  // Test with 1 and 8 threads
  for (uint nthreads : {1,8})
  {
      std::string filename = "testdata/" + _filename;
      auto maybe_pmesh = loadInputFromFile(filename);
      ASSERT_TRUE(maybe_pmesh.has_value());
      HexExtractor hexExtractor(maybe_pmesh->mesh, maybe_pmesh->igm);

      ASSERT_GT(hexExtractor.getInputMesh().n_cells(), 0) << "could not load mesh";

      Config config;
      config.verbose = false;
      config.num_threads = nthreads;
      hexExtractor.extract(config);

      HexahedralMesh& hexMesh = hexExtractor.getOutputMesh();

      EXPECT_EQ(hexMesh.n_cells(), _n_expected_cells) << "Not correctly extracted";

      OpenVolumeMesh::IO::ovmb_write(("Results/" + _filename + ".ovmb").c_str(), hexMesh);
  }
}

TEST(CellExtraction, CellExtractionS06U)
{
  test_file("s06u.hexex", 5167);
}

TEST(CellExtraction, CellExtractionN08C)
{
  test_file("n08c.hexex", 6368);
}

TEST(CellExtraction, CellExtractionI11c)
{
    test_file("i11c.hexex", 4982);
}
