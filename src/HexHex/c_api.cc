#include <HexHex/c_api.h>
#include <HexHex/HexHex.hh>
#include <HexHex/Utils/Typedefs.hh>

namespace OVM = OpenVolumeMesh;

extern "C"
{

struct HexMesh HEXHEX_EXPORT hexhex_extract_simple(
        size_t n_vertices,
        const double *vpos,           // length: 3 * n_vertices
        size_t n_tets,
        const uint32_t *tet_vertices, // length: 4 * n_cells
        const double *tet_igm         // length: 3 * 4 * n_cells
        )
{
    using namespace HexHex;
    HexHex::TetrahedralMesh mesh;
    mesh.reserve_vertices(n_vertices);
    mesh.reserve_cells(n_tets);
    mesh.reserve_edges(n_vertices);

    for (size_t i = 0; i < n_vertices; ++i) {
        mesh.add_vertex({vpos[3*i], vpos[3*i+1], vpos[3*i+2]});
    }
    HFParam igm = mesh.request_halfface_property<Vec3d>("HexHex::Parametrization");

    using VH = OVM::VH;
    for (size_t tet = 0; tet < n_tets; ++tet) {
        std::array<VH, 4> vhs {
            VH::from_unsigned(tet_vertices[4*tet]),
                VH::from_unsigned(tet_vertices[4*tet+1]),
                VH::from_unsigned(tet_vertices[4*tet+2]),
                VH::from_unsigned(tet_vertices[4*tet+3])};
        auto ch = mesh.add_cell(vhs[0],vhs[1],vhs[2],vhs[3]);
        for (int corner = 0; corner < 4; ++corner) {
            auto idx = (4 * tet + corner)*3;
            auto uv = Vec3d(tet_igm[idx], tet_igm[idx+1], tet_igm[idx+2]);
            igm[mesh.vertex_opposite_halfface(ch, vhs[corner])] = uv;
        }
    }
    auto output = extractHexMesh(mesh, igm);
    auto &ovm_hm = *output.hex_mesh;
    ovm_hm.enable_face_bottom_up_incidences(); // currently required for OVM hex_vertices()

    HexMesh hm;

    hm.n_vertices = ovm_hm.n_vertices();
    hm.vpos = (double*)::malloc(sizeof(double)*3*hm.n_vertices);

    for (size_t v = 0; v < hm.n_vertices; ++v) {
        const auto &p = ovm_hm.vertex(VH(v));
        hm.vpos[3*v+0] = p[0];
        hm.vpos[3*v+1] = p[1];
        hm.vpos[3*v+2] = p[2];
    }

    hm.n_hexes = ovm_hm.n_cells();
    hm.hex_vertices = (uint32_t*)::malloc(sizeof(double)*3*6*hm.n_hexes);

    for (size_t ch = 0; ch < hm.n_hexes; ++ch) {
        int i = 0;
        for (auto vh: ovm_hm.hex_vertices(OVM::CH(ch))) {
            hm.hex_vertices[8*ch+(i++)] = vh.idx();
        }
    }
    return hm;
}

}
