#pragma once

#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{

/*
 * Input Tet Mesh with some cached info and optimized methods
 */
struct TetMeshCache
{
    TetMeshCache(const TetrahedralMesh& tetmesh);
    void build();
    /*
     * Returns the 1st vertex of the halfedge without copying the edge object
     */
    inline VertexHandle from_vertex_handle(HalfEdgeHandle heh) const
    {
        if (heh.idx() % 2 == 0) return tetmesh.from_vertex_handle(heh);
        return tetmesh.to_vertex_handle(heh.opposite_handle());
    }

    /*
     * Returns the 2nd vertex of the halfedge without copying the edge object
     */
    inline VertexHandle to_vertex_handle(HalfEdgeHandle heh) const
    {
        return from_vertex_handle(heh.opposite_handle());
    }

    /*
     * Returns the 3 halfedges of the halfface without copying the halfface object
     */
    inline std::array<HalfEdgeHandle, 3> get_halfface_halfedges(const HalfFaceHandle hfh) const
    {
        if ((hfh.idx() % 2) == 0)
        {
            const auto& hehs = tetmesh.halfface(hfh).halfedges();
            return {hehs[0],hehs[1],hehs[2]};
        }
        else
        {
            const auto& hehs = tetmesh.halfface(hfh.opposite_handle()).halfedges();
            return {hehs[2].opposite_handle(),hehs[1].opposite_handle(),hehs[0].opposite_handle()};
        }
    }

    /*
     * Returns the 3 vertices of the halfface without copying the halfface object
     */
    inline std::array<VertexHandle, 3> get_halfface_vertices(const HalfFaceHandle hfh) const
    {
        if ((hfh.idx() % 2) == 0)
        {
            const auto& hehs = tetmesh.halfface(hfh).halfedges();
            return {from_vertex_handle(hehs[0]),from_vertex_handle(hehs[1]),from_vertex_handle(hehs[2])};
        }
        else
        {
            const auto& hehs = tetmesh.halfface(hfh.opposite_handle()).halfedges();
            return {to_vertex_handle(hehs[2]),to_vertex_handle(hehs[1]),to_vertex_handle(hehs[0])};
        }
    }

    inline VertexHandle get_nonincident_vertex_in_cell(HalfFaceHandle hfh)
    {
        const auto& cVhs = cellVertices[tetmesh.incident_cell(hfh)];
        const auto& hfVhs = get_halfface_vertices(hfh);
        for (VertexHandle cVh : cVhs) {
            if (std::find(hfVhs.begin(), hfVhs.end(), cVh) == hfVhs.end()) {
                return cVh;
            }
        }
        return VertexHandle(-1);
    }

    /*
     * Returns the 3 halffaces of ch that are incident to the vertex vh
     */
    inline std::vector<HalfFaceHandle> get_incident_halffaces_in_cell(CellHandle ch, VertexHandle vh) const
    {
        const auto& c = tetmesh.cell(ch);
        std::vector<HalfFaceHandle> res;
        res.reserve(3);
        for (const auto& hfh : c.halffaces())
            for (const auto& heh : get_halfface_halfedges(hfh))
                if (from_vertex_handle(heh) == vh) {
                    res.push_back(hfh);
                    if (res.size()==3) return res;
                }
        assert(res.size() == 3);
        return res;
    }

    /*
     * Returns the 2 halffaces of ch that are incident to the edge eh
     */
    inline std::vector<HalfFaceHandle> get_incident_halffaces_in_cell(CellHandle ch, EdgeHandle eh) const
    {
        const auto& c = tetmesh.cell(ch);
        std::vector<HalfFaceHandle> res;
        res.reserve(2);
        for (const auto& hfh : c.halffaces())
            for (const auto& heh : get_halfface_halfedges(hfh))
                if (heh.edge_handle() == eh) {
                    res.push_back(hfh);
                    if (res.size()==2) return res;
                }
        assert(res.size() == 2);
        return res;
    }

    /*
     * Returns the halfface of fh which has ch as the incident cell
     */
    inline HalfFaceHandle get_incident_halfface_in_cell(CellHandle ch, FaceHandle fh) const
    {
        if (tetmesh.incident_cell(fh.halfface_handle(0)) == ch) return fh.halfface_handle(0);
        if (tetmesh.incident_cell(fh.halfface_handle(1)) == ch) return fh.halfface_handle(1);
        return HalfFaceHandle(-1);
    }

    /*
     * Returns the halfface of the same cell as hfh that is incident to the edge but not equal to hfh
     */
    inline HalfFaceHandle get_other_incident_halfface_in_cell(EdgeHandle eh, HalfFaceHandle hfh) const
    {
        auto hfhs = get_incident_halffaces_in_cell(tetmesh.incident_cell(hfh), eh);
        if (hfhs[0] == hfh) return hfhs[1];
        if (hfhs[1] == hfh) return hfhs[0];
        return HalfFaceHandle(-1);
    }

    /*
     * Returns the vertex that is incident to the face but not to the edge
     */
    inline VertexHandle get_nonincident_vertex_in_face(FaceHandle fh, EdgeHandle eh) const
    {
        const auto& e = tetmesh.edge(eh);
        for (const VertexHandle vh : get_halfface_vertices(fh.halfface_handle(0)))
        {
            if (vh != e.from_vertex() && vh != e.to_vertex())
            {
                return vh;
            }
        }
        return VertexHandle(-1);
    }

    /*
     * Returns the 0 (if isolated), 1 (if boundary) or 2 (if inner) incident valid(!) cells to the face
     */
    inline std::vector<CellHandle> get_incident_cells(FaceHandle fh) const
    {
        std::vector<CellHandle> chs;
        CellHandle ch = tetmesh.incident_cell(fh.halfface_handle(0));
        if (ch.is_valid()) chs.push_back(ch);
        ch = tetmesh.incident_cell(fh.halfface_handle(1));
        if (ch.is_valid()) chs.push_back(ch);
        return chs;
    }

    /*
     * Returns whether or not the egde is one of the three edges of the face
     */
    inline bool is_incident(EdgeHandle eh, FaceHandle fh) const
    {
        const auto& e = tetmesh.edge(eh);
        const VertexHandle vh1 = e.from_vertex();
        const VertexHandle vh2 = e.to_vertex();
        const auto& vhs = get_halfface_vertices(fh.halfface_handle(0));
        if (vhs[0] != vh1 && vhs[1] != vh1 && vhs[2] != vh1) return false;
        if (vhs[0] != vh2 && vhs[1] != vh2 && vhs[2] != vh2) return false;
        return true;
    }

    const TetrahedralMesh& tetmesh;

    OVM::CellPropertyT<std::array<VertexHandle, 4>> cellVertices; // tet vertices per tet cell
    OVM::VertexPropertyT<std::vector<CellHandle>> vertexCells; // tet cells per tet vertex

    OVM::VertexPropertyT<int> v_owners; // Parent tet per vertex
    OVM::EdgePropertyT<int> e_owners; // Parent tet per edge
    OVM::FacePropertyT<int> f_owners; // Parent tet per face

};

}

