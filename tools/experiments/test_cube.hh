#pragma once

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <HexHex/Utils/Typedefs.hh>
//#include <OpenVolumeMesh/IO/ovmb_write.hh>

const unsigned char TET_CUBE_6[6][4][3] = {
    {{1,1,1}, {1,0,1}, {0,0,1}, {1,1,0}},
    {{0,1,0}, {1,1,0}, {0,0,1}, {0,0,0}},
    {{0,1,1}, {1,1,1}, {0,0,1}, {1,1,0}},
    {{0,0,1}, {1,0,1}, {1,0,0}, {1,1,0}},
    {{0,0,0}, {0,0,1}, {1,0,0}, {1,1,0}},
    {{0,1,0}, {1,1,0}, {0,1,1}, {0,0,1}}
};

using Vec3d = OpenVolumeMesh::Vec3d;

inline void createCube(HexHex::TetrahedralMesh &mesh,
        //HexEx::PerCellOVM::VertexPropertyT<Vec3d>
        auto &parameters,
        std::array<int, 3> n_blocks,
        std::array<int, 3> hex_dims)
{
    mesh.clear();
    using namespace OpenVolumeMesh;

    /// number of vertices in each dimension
    std::array<int, 3> vdim = {n_blocks[0] + 1, n_blocks[1] + 1, n_blocks[2] + 1};

    std::vector<VertexHandle> all_vh(vdim[0] * vdim[1] * vdim[2]);
    auto vhs =[&](int x, int y, int z) -> VertexHandle& {
        return all_vh.at((vdim[1] * x + y) * vdim[2] + z);
    };

    auto vigm = mesh.request_vertex_property<Vec3d>("uvw");
    // Add vertices
    for (int x = 0; x <= n_blocks[0]; ++x) {
        for (int y = 0; y <= n_blocks[1]; ++y) {
            for (int z = 0; z <= n_blocks[2]; ++z) {
                auto param =  Vec3d(
                        double(x)/(n_blocks[0]),
                        double(y)/(n_blocks[1]),
                        double(z)/(n_blocks[2]));
                auto p = param;
#if 1
                // TODO: twirl this a bit :)
                auto alpha = p[2]/10;
                auto beta = p[0]/10;
                auto gamma = p[1]/10;
                p = Vec3d(
                        cos(alpha) * p[0] - sin(alpha) * p[1],
                        sin(alpha) * p[0] + cos(alpha) * p[1],
                        p[2]);
                p = Vec3d(
                        p[0],
                        cos(beta) * p[1] - sin(beta) * p[2],
                        sin(beta) * p[1] + cos(beta) * p[2]);
                p = Vec3d(
                        cos(-alpha) * p[0] - sin(-alpha) * p[1],
                        sin(-alpha) * p[0] + cos(-alpha) * p[1],
                        p[2]);
                p = Vec3d(
                        cos(gamma) * p[0] - sin(gamma) * p[2],
                        p[1],
                        sin(gamma) * p[0] + cos(gamma) * p[2]);
#endif
                auto vh = mesh.add_vertex(p);
                vhs(x,y,z) = vh;
                vigm[vh] = {
                    hex_dims[0] * param[0],
                    hex_dims[1] * param[1],
                    hex_dims[2] * param[2]};
            }
        }
    }

    auto add_tets = [&](int x, int y, int z) -> void
    {
        for (int i = 0; i < 6; ++i)
        {
            std::vector<VertexHandle> res; res.reserve(4);
            for (int j = 0; j < 4; ++j) res.push_back(vhs(
                        x+TET_CUBE_6[i][j][0],
                        y+TET_CUBE_6[i][j][1],
                        z+TET_CUBE_6[i][j][2]));
            mesh.add_cell(res);
        }
    };

    // Add tets
    for (int x = 0; x < n_blocks[0]; ++x) {
        for (int y = 0; y < n_blocks[1]; ++y) {
            for (int z = 0; z < n_blocks[2]; ++z) {
                add_tets(x, y, z);
            }
        }
    }

    // Set parameters
    for (auto c_it = mesh.c_iter(); c_it.is_valid(); ++c_it) {
        for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it) {
            const auto &p = mesh.vertex(*cv_it);
            parameters[mesh.vertex_opposite_halfface(*c_it, *cv_it)] = vigm[*cv_it];
        }
    }
}

inline void mycube(HexHex::TetrahedralMesh &cube,
        auto &igm,
        int xxx,
        int n = 1)
{
    auto hex_dim = std::array<int,3>{xxx,xxx/2+1,xxx/3+1};
    auto n_hexes = size_t(hex_dim[0]) * hex_dim[1] * hex_dim[2];
    std::cout << "will produce " << n_hexes << " hexes." << std::endl;
    if (n == 0) {
        double fac = 1./M_E;
        n = std::lround(std::ceil(std::pow(n_hexes * fac, 1./3)));
    }
    auto n_blocks = std::array<int,3>{n,n,n};
    std::cout << "n blocks: " << n << std::endl;
    createCube(cube, igm, n_blocks, hex_dim);
}
