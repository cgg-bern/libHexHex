/*
 * Copyright 2025 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */

#include <gtest/gtest.h>
#include "common.hh"
#include <HexHex/HexExtractor.hh>
#include <HexHex/Predicates/TracingPredicates.hh>

using namespace HexHex;

class NavigationTest : public ::testing::Test {
protected:
    virtual void SetUp() {

        cell_vhs[0] = mesh.add_vertex(Vec3d(0,0,0));
        cell_vhs[1] = mesh.add_vertex(Vec3d(1,0,0));
        cell_vhs[2] = mesh.add_vertex(Vec3d(1,0,-1));
        cell_vhs[3] = mesh.add_vertex(Vec3d(0,1,0));

        cell_hfhs[0] = mesh.halfface_handle(mesh.add_face({cell_vhs[0],cell_vhs[1],cell_vhs[2]}), 0);
        cell_hfhs[1] = mesh.halfface_handle(mesh.add_face({cell_vhs[0],cell_vhs[2],cell_vhs[3]}), 0);
        cell_hfhs[2] = mesh.halfface_handle(mesh.add_face({cell_vhs[0],cell_vhs[3],cell_vhs[1]}), 0);
        cell_hfhs[3] = mesh.halfface_handle(mesh.add_face({cell_vhs[1],cell_vhs[3],cell_vhs[2]}), 0);

        ch  = mesh.add_cell({cell_hfhs[0],cell_hfhs[1],cell_hfhs[2],cell_hfhs[3]},true);

        for (int i = 0; i < 4; ++i)
            cell_params[i] = mesh.vertex(cell_vhs[i]);

        cachedMesh = new TetMeshCache(mesh);
    }

    virtual void TearDown() { delete cachedMesh; }

    TetrahedralMesh mesh;
    TetMeshCache* cachedMesh;

    std::array<VertexHandle, 4> cell_vhs;
    std::array<Parameter, 4> cell_params;
    std::array<HalfFaceHandle, 4> cell_hfhs;
    OpenVolumeMesh::CellHandle ch;


};


TEST_F(NavigationTest, SimpleFacePiercing) {

    cell_params[0] = Parameter(  0, 0, 0);
    cell_params[1]  = Parameter(0.5,-1, 1);
    cell_params[2]  = Parameter(0.5,-1,-1);
    cell_params[3]  = Parameter(0.5, 1, 0);

    MeshElement traceThrough;
    HalfFaceHandle traceHfh = findOppositeNextHalfFace(
        *cachedMesh,
        ch, cell_vhs, cell_params,
        HalfFaceHandle(-1),
        Parameter(0,0,0), Direction(1,0,0),
        Direction(0,0,1),
        traceThrough, true
        );

    EXPECT_EQ(cell_hfhs[3], traceHfh);
    EXPECT_TRUE(traceThrough.is_face());
}

TEST_F(NavigationTest, EdgePiercing) {

    cell_params[0] = Parameter(  0, 0,-1);
    cell_params[1] = Parameter(  0, 0, 1);
    cell_params[2] = Parameter(0.5,-1, 0);
    cell_params[3] = Parameter(0.5, 1, 0);

    MeshElement traceThrough;
    HalfFaceHandle traceHfh = findOppositeNextHalfFace(
        *cachedMesh,
        ch, cell_vhs, cell_params,
        HalfFaceHandle(-1),
        Parameter(0,0,0), Direction(1,0,0),
        Direction(0,0,1),
        traceThrough, true
        );

    EXPECT_EQ(cell_hfhs[3], traceHfh);
    EXPECT_TRUE(traceThrough.is_edge());
}

TEST_F(NavigationTest, DoubleEdgePiercing) {

    cell_params[0] = Parameter(  0, 0,-1);
    cell_params[1] =  Parameter(  0, 0, 1);
    cell_params[2] =  Parameter(0.5,-1, 0);
    cell_params[3] =  Parameter(0.5, 1, 0);

    MeshElement traceThrough;
    HalfFaceHandle traceHfh = findOppositeNextHalfFace(
        *cachedMesh,
        ch, cell_vhs, cell_params,
        HalfFaceHandle(-1),
        Parameter(0,0,0), Direction(1,0,0),
        Direction(0,1,0),
        traceThrough, true
        );

    EXPECT_TRUE(traceHfh == cell_hfhs[1] || traceHfh == cell_hfhs[3]);
    EXPECT_TRUE(traceThrough.is_edge());

}


TEST_F(NavigationTest, VertexPiercingSimple) {

    cell_params[0] = Parameter(  0,-0.5, 0.5);
    cell_params[1] = Parameter(  0,-0.5,-0.5);
    cell_params[2] = Parameter(  0, 0.5,   0);
    cell_params[3] = Parameter(0.5,   0,   0);

    MeshElement traceThrough;
    HalfFaceHandle traceHfh = findOppositeNextHalfFace(
        *cachedMesh,
        ch, cell_vhs, cell_params,
        HalfFaceHandle(-1),
        Parameter(0,0,0), Direction(1,0,0),
        Direction(0,1,0),
        traceThrough, true
        );

    EXPECT_NE(cell_hfhs[0], traceHfh);
    EXPECT_TRUE(traceThrough.is_vertex());
}

TEST_F(NavigationTest, VertexPiercingWithEdge) {

    cell_params[0] = Parameter(  0,-0.5, 0.5);
    cell_params[1] = Parameter(  0,-0.5,-0.5);
    cell_params[2] = Parameter(  0, 0.5,   0);
    cell_params[3] = Parameter(0.5,   0,   0);

    MeshElement traceThrough;
    HalfFaceHandle traceHfh = findOppositeNextHalfFace(
        *cachedMesh,
        ch, cell_vhs, cell_params,
        HalfFaceHandle(-1),
        Parameter(0,0,0), Direction(1,0,0),
        Direction(0,1,0),
        traceThrough, true
        );

    EXPECT_TRUE(traceHfh == cell_hfhs[1] || traceHfh == cell_hfhs[3]);
    EXPECT_TRUE(traceThrough.is_vertex());

}

