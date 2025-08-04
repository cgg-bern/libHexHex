#include <HexHex/Utils/TestMeshes.hh>
#include <HexHex/Commons/Transition.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <random>

namespace HexHex
{

namespace {
    // 6 Tetrahedra build a tileable unit cube
    static const unsigned char TET_CUBE_6[6][4][3] = {
        {{1,1,1}, {1,0,1}, {0,0,1}, {1,1,0}},
        {{0,1,0}, {1,1,0}, {0,0,1}, {0,0,0}},
        {{0,1,1}, {1,1,1}, {0,0,1}, {1,1,0}},
        {{0,0,1}, {1,0,1}, {1,0,0}, {1,1,0}},
        {{0,0,0}, {0,0,1}, {1,0,0}, {1,1,0}},
        {{0,1,0}, {1,1,0}, {0,1,1}, {0,0,1}}
    };
}

void createHexCube(HexahedralMesh& mesh)
{
    mesh.clear();

    std::vector<VertexHandle> vhs;
    vhs.reserve(8);
    vhs.push_back(mesh.add_vertex(Vec3d( 0.0, 0.0, 0.0)));
    vhs.push_back(mesh.add_vertex(Vec3d( 1.0, 0.0, 0.0)));
    vhs.push_back(mesh.add_vertex(Vec3d( 1.0, 0.0, 1.0)));
    vhs.push_back(mesh.add_vertex(Vec3d( 0.0, 0.0, 1.0)));
    vhs.push_back(mesh.add_vertex(Vec3d( 0.0, 1.0, 0.0)));
    vhs.push_back(mesh.add_vertex(Vec3d( 1.0, 1.0, 0.0)));
    vhs.push_back(mesh.add_vertex(Vec3d( 1.0, 1.0, 1.0)));
    vhs.push_back(mesh.add_vertex(Vec3d( 0.0, 1.0, 1.0)));

    std::vector<HalfFaceHandle> hfhs;
    hfhs.reserve(6);
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[0],vhs[1],vhs[2],vhs[3]}), 1));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[7],vhs[6],vhs[5],vhs[4]}), 1));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[1],vhs[0],vhs[4],vhs[5]}), 1));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[2],vhs[1],vhs[5],vhs[6]}), 1));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[3],vhs[2],vhs[6],vhs[7]}), 1));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[0],vhs[3],vhs[7],vhs[4]}), 1));
    mesh.add_cell(hfhs);
}

TestMesh createCube(uint tetLength, uint hexLength)
{
    TestMesh result;
    auto& mesh = result.mesh;
    auto& parameters = result.igm;

    uint n = tetLength+1;
    std::vector<VertexHandle> all_vh(n*n*n);
    auto vhs =[&](uint x, uint y, uint z) -> VertexHandle& {
        return all_vh.at(z + n*y + n*n*x);
    };

    // Add vertices
    for (uint x = 0; x < n; ++x)
        for (uint y = 0; y < n; ++y)
            for (uint z = 0; z < n; ++z)
                vhs(x,y,z) = (mesh.add_vertex(Vec3d(x,y,z)/tetLength));

    auto add_tets = [&](uint x, uint y, uint z) -> void
    {
        for (uint i = 0; i < 6; ++i)
        {
            std::vector<VertexHandle> res; res.reserve(4);
            for (uint j = 0; j < 4; ++j) res.push_back(vhs(
                    x+TET_CUBE_6[i][j][0],
                    y+TET_CUBE_6[i][j][1],
                    z+TET_CUBE_6[i][j][2]));
            mesh.add_cell(res);
        }
    };

    // Add tets
    for (uint x = 0; x < tetLength; ++x)
        for (uint y = 0; y < tetLength; ++y)
            for (uint z = 0; z < tetLength; ++z)
                add_tets(x, y, z);

    // Set parameters
    for (auto hf_it = mesh.hf_iter(); hf_it.is_valid(); ++hf_it) {
        if (!mesh.incident_cell(*hf_it).is_valid()) {continue;}
        VertexHandle vh = mesh.halfface_opposite_vertex(*hf_it);
        parameters[*hf_it] = hexLength * mesh.vertex(vh);
    }
    return result;
}

TestMesh createCylinder(uint valence, uint hexScale)
{
    TestMesh result;
    auto& mesh = result.mesh;
    auto& parameters = result.igm;

    if (valence < 3) valence = 3;

    double h = 0.5;
    double r1 = 1.0;
    double r2 = 2.0;

    auto angle = [&](const uint& i) -> double {return M_PI * ((double)i)/((double)valence);};

    // Add outer vertices
    for (uint i = 0u; i < 2u*valence; ++i)
    {
        double alpha = angle(i);
        double r = (i%2==0)? r2 : r1;
        double z = r*std::cos(alpha); if (fabs(z) < 1e-6) z = 0;
        double x = r*std::sin(alpha); if (fabs(x) < 1e-6) x = 0;
        auto vh0 = mesh.add_vertex(Vec3d(x,0,z));
        auto vh1 = mesh.add_vertex(Vec3d(x,h,z));
        mesh.add_edge(vh0, vh1);
    }

    // Add inner vertices
    mesh.add_edge(
        mesh.add_vertex(Vec3d(0,0,0)),
        mesh.add_vertex(Vec3d(0,h,0))
        );

    assert(mesh.n_vertices() == 4*valence+2);

    auto outer_vertex = [&valence](const int& i, const int& o) -> VertexHandle {int m = 4*valence; return VertexHandle((2*i+m+o)%m);};
    auto inner_vertex = [&valence](const int& o) -> VertexHandle {int V=4*valence+2; return VertexHandle(V-2+o);};
    auto rotate = [](const Parameter& p, int r) -> Parameter
    {
        r = r%4;
        switch (r) {
        case 0: return Parameter(+p[0], p[1], +p[2]);
        case 1: return Parameter(+p[2], p[1], -p[0]);
        case 2: return Parameter(-p[0], p[1], -p[2]);
        case 3: return Parameter(-p[2], p[1], +p[0]);
        default: assert(false); break;
        }
        return Parameter();
    };

    // Add Tets and set Parameters

    for (uint i = 0; i < 2*valence; i+=2)
    {
        VertexHandle block_vhs[2][2][2];
        block_vhs[0][0][0] = inner_vertex(0);
        block_vhs[0][1][0] = inner_vertex(1);

        block_vhs[1][0][0] = outer_vertex(i+1, 0);
        block_vhs[1][1][0] = outer_vertex(i+1, 1);

        block_vhs[1][0][1] = outer_vertex(i, 0);
        block_vhs[1][1][1] = outer_vertex(i, 1);

        block_vhs[0][0][1] = outer_vertex(i-1, 0);
        block_vhs[0][1][1] = outer_vertex(i-1, 1);

        for (int j = 0; j < 6; ++j)
        {
            std::vector<VertexHandle> vhs; vhs.reserve(4);
            std::vector<Parameter> params; params.reserve(4);
            for (int k = 0; k < 4; ++k)
            {
                int x = TET_CUBE_6[j][k][0], y = TET_CUBE_6[j][k][1], z = TET_CUBE_6[j][k][2];
                vhs.push_back(block_vhs[1-z][y][x]); // make it tileable for the star mesh
                params.push_back(hexScale * rotate(Parameter(x,y,z), i/2));
            }

            CellHandle ch = mesh.add_cell(vhs);
            for (int k = 0; k < 4; ++k) {
                for (HalfFaceHandle hfh : mesh.cell(ch).halffaces()) {
                    if (mesh.halfface_opposite_vertex(hfh) == vhs[k]) {
                        parameters[hfh] = params[k];
                        break;
                    }
                }
            }
        }
    }
    return result;
}

void randomizeTransitions(TestMesh &tm)
{
#pragma omp parallel
    {
        std::mt19937 rng;
#   pragma omp for
        for (int i = 0; i < static_cast<int>(tm.mesh.n_cells()); ++i)
        {
            CellHandle ch(i);
            auto tau = Transition::random(rng);
            for (HalfFaceHandle hfh : tm.mesh.cell(ch).halffaces()) {
                tm.igm[hfh] = tau.transform_point(tm.igm[hfh]);
            }
        }
    }
}

}
