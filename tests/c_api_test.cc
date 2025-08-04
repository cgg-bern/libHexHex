/*
 * Copyright 2025 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 *
 * This file is part of HexHex.
 */

#include <gtest/gtest.h>
#include <HexHex/c_api.h>
#include "common.hh"

using namespace HexHex;

TEST(c_api, basic)
{
    TetrahedralMesh mesh;
    createCube(mesh, 1.);

    std::vector<double> vpos(mesh.n_vertices()*3);
    std::vector<uint32_t> tet_vertices(mesh.n_cells()*4);
    std::vector<double> tet_igm(mesh.n_cells()*4*3);
    for (const auto vh: mesh.vertices()) {
        const auto &p = mesh.vertex(vh);
        vpos[vh.idx()*3] = p[0];
        vpos[vh.idx()*3+1] = p[1];
        vpos[vh.idx()*3+2] = p[2];
    }
    for (const auto ch: mesh.cells()) {
        int i = 0;
        for (const auto vh: mesh.tet_vertices(ch)) {
            tet_vertices[4*ch.idx() + i] = vh.idx();
            Vec3d uv = mesh.vertex(vh) * 3;
            tet_igm[(4*ch.idx() + i)*3 + 0 ] = uv[0];
            tet_igm[(4*ch.idx() + i)*3 + 1 ] = uv[1];
            tet_igm[(4*ch.idx() + i)*3 + 2 ] = uv[2];
            ++i;
        }
    }
    HexMesh hm = ::hexhex_extract_simple(
            mesh.n_vertices(),
            vpos.data(),
            mesh.n_cells(),
            tet_vertices.data(),
            tet_igm.data());

    EXPECT_EQ(hm.n_vertices, 4*4*4);
    EXPECT_EQ(hm.n_hexes, 3*3*3);
}

