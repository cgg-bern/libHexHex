#include <HexHex/Commons/MeshElement.hh>

namespace HexHex
{

bool MeshElement::is_incident(const TetMeshCache &cache, const MeshElement& elem) const
{
    assert(this->is_valid() && elem.is_valid());

    if (this->dim() > elem.dim()) return elem.is_incident(cache, *this);
    if (this->dim() == elem.dim()) return false;

    // this.dim < elem.dim
    if (this->is_vertex())
    {
        auto vh = this->vh();
        if (elem._type == ELEMENT_EDGE_TYPE) { // vertex - edge
            auto e = cache.tetmesh.edge(elem.eh());
            return e.from_vertex()==vh || e.to_vertex()==vh;
        } else if (elem._type == ELEMENT_FACE_TYPE) { // vertex - face
            auto vhs = cache.tetmesh.get_cell_vertices(elem.fh().halfface_handle(0));
            for (auto& fvh : vhs) {if (vh==fvh) return true;}
            return false;
        } else { // vertex - cell
            assert(elem._type == ELEMENT_CELL_TYPE);
            auto vhs = cache.cellVertices[elem.ch()];
            for (auto& cvh : vhs) {if (vh==cvh) return true;}
            return false;
        }
    } else if (this->_type == ELEMENT_EDGE_TYPE) {
        auto eh = this->eh();
        if (elem._type == ELEMENT_FACE_TYPE) { // edge - face
            return cache.tetmesh.is_incident(elem.fh(), eh);
        } else { // edge - cell
            assert(elem._type == ELEMENT_CELL_TYPE);
            for (auto ce_iter = cache.tetmesh.ce_iter(elem.ch()); ce_iter.is_valid(); ++ce_iter) {
                if (*ce_iter == eh) {return true;}
            }
            return false;
        }
    } else if (this->_type == ELEMENT_FACE_TYPE) { // face - cell
        assert(elem._type == ELEMENT_CELL_TYPE);
        auto fh = this->fh();
        auto ch = elem.ch();
        return cache.tetmesh.incident_cell(fh.halfface_handle(0))==ch || cache.tetmesh.incident_cell(fh.halfface_handle(1))==ch;
    } else {
        assert(false);
        return false;
    }
}

MeshElement MeshElement::common_smaller_incident_element(const TetMeshCache& cache, const MeshElement& elem) const
{
    const TetrahedralMesh& mesh = cache.tetmesh;

    if (this->dim() != elem.dim() || *this==elem) return MeshElement();
    if (this->is_vertex()) {
        for (auto voh_iter = mesh.voh_iter(this->vh()); voh_iter.is_valid(); ++voh_iter) {
            if (mesh.to_vertex_handle(*voh_iter) == elem.vh()) {
                return MeshElement((*voh_iter).edge_handle());
            }
        }
    } else if (this->is_edge()) {
        auto e1 = mesh.edge(this->eh());
        auto e2 = mesh.edge(elem.eh());
        for (auto vh1 : {e1.from_vertex(), e1.to_vertex()}) {
            for (auto vh2 : {e2.from_vertex(), e2.to_vertex()}) {
                if (vh1 == vh2) {
                    return MeshElement(vh1);
                }
            }
        }
    } else if (this->is_face()) {
        auto f1 = mesh.face(this->fh());
        auto f2 = mesh.face(elem.fh());
        for (auto heh1 : f1.halfedges()) {
            for (auto heh2 : f2.halfedges()) {
                if (heh1.edge_handle() == heh2.edge_handle()) {
                    return MeshElement(heh1.edge_handle());
                }
            }
        }
    } else {
        assert(this->is_cell());
        auto c1 = mesh.cell(this->ch());
        auto c2 = mesh.cell(elem.ch());
        for (auto hfh1 : c1.halffaces()) {
            for (auto hfh2 : c2.halffaces()) {
                if (hfh1.face_handle() == hfh2.face_handle()) {
                    return MeshElement(hfh1.face_handle());
                }
            }
        }
    }
    return MeshElement();
}

std::vector<CellHandle> MeshElement::get_incident_cells(const TetMeshCache &cache) const
{
    assert(is_valid());

    if (is_cell()) return {ch()};

    std::vector<CellHandle> res;
    res.reserve(8);
    if (is_face())
    {
        HalfFaceHandle hf = hfh();
        CellHandle ch = cache.tetmesh.incident_cell(hf);
        if (ch.is_valid()) res.push_back(ch);
        hf = hf.opposite_handle();
        ch = cache.tetmesh.incident_cell(hf);
        if (ch.is_valid()) res.push_back(ch);
        return res;
    }
    else if (is_edge())
    {
        std::vector<CellHandle> res;
        for (auto ec_it = cache.tetmesh.ec_iter(eh()); ec_it.is_valid(); ++ec_it)
        {
            res.push_back(*ec_it);
        }
        return res;
    }

    return cache.vertexCells[vh()];
}

}

HexHex::MeshElement HexHex::getElementBetweenInCell(const TetMeshCache& cache, const CellHandle ch, const MeshElement& m1, const MeshElement& m2)
{
    assert(m1 == ch || m1.is_incident(cache, ch));
    assert(m2 == ch || m2.is_incident(cache, ch));

    // On is a cell? -> In between is the cell
    if (m1.is_cell()) {assert(m1 == ch); return m1;}
    if (m2.is_cell()) {assert(m2 == ch); return m2;}

    // Wlog dim1 <= dim2
    if (m1.dim() > m2.dim()) return getElementBetweenInCell(cache, ch, m2, m1);

    // They are the same element? -> In between is also the same
    if (m1 == m2) {return m1;}

    // Both are vertices? -> Between is edge incident to both
    if (m1.is_vertex() && m2.is_vertex())
    {
        return cache.tetmesh.find_halfedge_in_cell(m1.vh(), m2.vh(), ch).edge_handle();
    }

    // Both are edges? -> Either there is a face incident to both, or otherwise in between is the cell
    if (m1.is_edge() && m2.is_edge())
    {
        for (auto hfh : cache.tetmesh.cell(ch).halffaces())
        {
            bool inc1 = false;
            bool inc2 = false;
            for (auto heh : cache.tetmesh.halfface(hfh).halfedges())
            {
                if (m1.eh() == heh.edge_handle()) {inc1 = true;}
                if (m2.eh() == heh.edge_handle()) {inc2 = true;}
            }
            if (inc1 && inc2)
            {
                return hfh.face_handle();
            }
        }
        return ch;
    }

    // Both are faces? -> In between is cell
    if (m1.is_face() && m2.is_face())
    {
        return ch;
    }

    assert(m1.dim() < m2.dim());

    // If m1 ~ m2 then return m2
    if (m1.is_incident(cache, m2)) {return m2;}

    // Otherwise, return m s.t. m ~ m1, m ~ m2 and dim(m) = dim(m2) + 1
    if (m1.is_vertex())
    {
        if (m2.is_edge())
        {
            // Return face incident to m1, m2
            for (auto hfh : cache.tetmesh.cell(ch).halffaces())
            {
                bool inc1 = false;
                bool inc2 = false;
                for (auto heh : cache.tetmesh.halfface(hfh).halfedges())
                {
                    if (m1.vh() == cache.tetmesh.from_vertex_handle(heh)) {inc1 = true;}
                    if (m2.eh() == heh.edge_handle()) {inc2 = true;}
                }
                if (inc1 && inc2)
                {
                    return hfh.face_handle();
                }
            }
        }
        assert(m2.is_face());
        return ch;
        assert(false);
    }
    assert(m1.is_edge());
    assert(m2.is_face());
    return ch;
}
