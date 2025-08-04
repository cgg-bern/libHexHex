#include <HexHex/Utils/TetMeshCache.hh>

namespace HexHex
{

TetMeshCache::TetMeshCache(const TetrahedralMesh& tetmesh) :
    tetmesh(tetmesh),

    cellVertices(tetmesh.create_private_cell_property<std::array<VertexHandle,4>>()),
    vertexCells(tetmesh.create_private_vertex_property<std::vector<CellHandle>>()),

    v_owners(tetmesh.create_private_vertex_property<int>("",-1)),
    e_owners(tetmesh.create_private_edge_property<int>("",-1)),
    f_owners(tetmesh.create_private_face_property<int>("",-1))
{}

void TetMeshCache::build()
{
    ScopedStopWatch _{sw::cachePreSanitizeProperties};

    // cache cell vertices and tet owners
#pragma omp parallel for  //default(none) shared(tetmesh, cellVertices, f_owners, e_owners)
    for (int i = 0; i < static_cast<int>(tetmesh.n_cells()); ++i) {
        auto ch = CellHandle(i);

        // Get cell halffaces
        const auto& hfhs = tetmesh.cell(ch).halffaces();

        // Cache the cell vertices fast
        const auto& hfh0_vhs = get_halfface_vertices(hfhs[0]); // vertices of 1st halfface
        const auto& hfh1_vhs = get_halfface_vertices(hfhs[1]); // vertices of 2nd halfface
        cellVertices[ch][0] = hfh0_vhs[0];
        cellVertices[ch][1] = hfh0_vhs[1];
        cellVertices[ch][2] = hfh0_vhs[2];
        for (VertexHandle vh1 : hfh1_vhs) // look for 4th cell vertex
            if (vh1 != hfh0_vhs[0] && vh1 != hfh0_vhs[1] && vh1 != hfh0_vhs[2])
                cellVertices[ch][3] = vh1;

        for (const HalfFaceHandle hfh : hfhs)
        {
            // Set Face Parent
            #pragma omp atomic write
            f_owners[hfh.face_handle()] = i;

            for (const HalfEdgeHandle heh : get_halfface_halfedges(hfh))
            {
                // Set Edge Parent
                #pragma omp atomic write
                e_owners[heh.edge_handle()] = i;
            }
        }
    }

    // Cache Tet cells per vertex and vertex tet parent
#pragma omp parallel for // default(none) shared(tetmesh, vertexCells, v_owners, std::cerr)
    for (int i = 0; i < static_cast<int>(tetmesh.n_vertices()); ++i)
    {
        VertexHandle vh(i);

        std::vector<CellHandle> chs;

        for (auto vc_it = tetmesh.vc_iter(vh); vc_it.is_valid(); ++vc_it)
        {
            chs.push_back(*vc_it);
        }

        if (chs.size()==0) {
            std::cerr << "HexHex: Isolated Vertex" << std::endl;
            continue;
        }

        // Set Vertex Parent
        v_owners[vh] = chs.front().idx();

        // Set Vertex Cells
        vertexCells[vh] = std::move(chs);
    }
}
}

