#include <gtest/gtest.h>
#include "common.hh"
#include <HexHex/HexExtractor.hh>
#include <HexHex/Predicates/GeneratorPredicates.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

using namespace HexHex;

void build_tet(TetrahedralMesh& mesh,
              std::array<Vec3d,4> params
              )
{
    mesh.clear();

    mesh.add_vertex(params[0]);
    mesh.add_vertex(params[1]);
    mesh.add_vertex(params[2]);
    mesh.add_vertex(params[3]);

    mesh.add_cell(VertexHandle(0),VertexHandle(1),VertexHandle(2),VertexHandle(3));

}

TEST(Predicates, generatorIsVertexTest)
{
    TetrahedralMesh mesh;
    std::array<Parameter, 4> params = {Vec3d(0, 0, 0),Vec3d(0, 0, 1),Vec3d(1, 0, 0),Vec3d(0, 1, 0)};
    std::array<VertexHandle, 4> vhs = {VertexHandle(0),VertexHandle(1),VertexHandle(2),VertexHandle(3)};
    build_tet(mesh, params);

    TetMeshCache meshCache(mesh);
    meshCache.build();

    EXPECT_EQ(
        computeVertexGeneratorElement(meshCache, CellHandle(0), vhs, params, Parameter(0,0,0)),
        VertexHandle(0)
        );

    EXPECT_EQ(
        computeVertexGeneratorElement(meshCache, CellHandle(0), vhs, params, Parameter(0,0,1)),
        VertexHandle(1)
        );

    EXPECT_EQ(
        computeVertexGeneratorElement(meshCache, CellHandle(0), vhs, params, Parameter(1,0,0)),
        VertexHandle(2)
        );

    EXPECT_EQ(
        computeVertexGeneratorElement(meshCache, CellHandle(0), vhs, params, Parameter(0,1,0)),
        VertexHandle(3)
        );
}


TEST(Predicates, generatorIsEdgeTest)
{
    TetrahedralMesh mesh;
    std::array<Parameter, 4> params = {Vec3d(0, 0, 0),Vec3d(0, 0, 1),Vec3d(1, 0, 0),Vec3d(0, 1, 0)};
    std::array<VertexHandle, 4> vhs = {VertexHandle(0),VertexHandle(1),VertexHandle(2),VertexHandle(3)};
    build_tet(mesh, params);

    TetMeshCache meshCache(mesh);
    meshCache.build();

    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0,0,0.2)).is_edge());
    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0,0.9,0)).is_edge());
    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0.01,0,0)).is_edge());
    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0.5,0,0.5)).is_edge());
    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0.25,0.75,0)).is_edge());
    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0,0.75,0.25)).is_edge());
}

TEST(Predicates, generatorIsFaceTest)
{
    TetrahedralMesh mesh;
    std::array<Parameter, 4> params = {Vec3d(0, 0, 0),Vec3d(0, 0, 1),Vec3d(1, 0, 0),Vec3d(0, 1, 0)};
    std::array<VertexHandle, 4> vhs = {VertexHandle(0),VertexHandle(1),VertexHandle(2),VertexHandle(3)};
    build_tet(mesh, params);

    TetMeshCache meshCache(mesh);
    meshCache.build();

    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0.2,0,0.1)).is_face());
    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0,0.4,0.1)).is_face());
    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0.3,0,1e-9)).is_face());
    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0.25,0.5,0.25)).is_face());
}

TEST(Predicates, generatorIsCellTest)
{
    TetrahedralMesh mesh;
    std::array<Parameter, 4> params = {Vec3d(0, 0, 0),Vec3d(0, 0, 1),Vec3d(1, 0, 0),Vec3d(0, 1, 0)};
    std::array<VertexHandle, 4> vhs = {VertexHandle(0),VertexHandle(1),VertexHandle(2),VertexHandle(3)};
    build_tet(mesh, params);

    TetMeshCache meshCache(mesh);
    meshCache.build();

    EXPECT_TRUE(computeVertexGeneratorElement(meshCache, CellHandle(0), vhs,
                                              params, Parameter(0.3,0.2,0.1)).is_cell());
}
