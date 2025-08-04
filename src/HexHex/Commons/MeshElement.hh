#pragma once

#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Utils/TetMeshCache.hh>

namespace HexHex
{

enum MeshElementType : char {
    ELEMENT_NULL_TYPE = -1,
    ELEMENT_VERTEX_TYPE = 0,
    ELEMENT_EDGE_TYPE = 1,
    ELEMENT_FACE_TYPE = 2,
    ELEMENT_CELL_TYPE = 3
};

/*
 * Stores either a Vertex, Halfedge, Halfface, or Cell
 */
class MeshElement {
public:
    MeshElement(MeshElementType type, int id) {set(type,id);}
    MeshElement() : MeshElement(ELEMENT_NULL_TYPE, -1) {}
    MeshElement(VertexHandle vh) {set(vh);}
    MeshElement(EdgeHandle eh) {set(eh);}
    MeshElement(FaceHandle fh) {set(fh);}
    MeshElement(CellHandle ch) {set(ch);}
    MeshElement(HalfFaceHandle hfh) {set(hfh);}
    MeshElement(HalfEdgeHandle heh) {set(heh);}

    inline void set(MeshElementType type, int idx) {
        this->_type=type;
        this->_idx=idx;
    }

    inline void set(const VertexHandle& vh) {set(ELEMENT_VERTEX_TYPE,_idx=vh.idx());}
    inline void set(const EdgeHandle& eh) {set(ELEMENT_EDGE_TYPE,eh.halfedge_handle(0).idx());}
    inline void set(const FaceHandle& fh) {set(ELEMENT_FACE_TYPE,fh.halfface_handle(0).idx());}
    inline void set(const CellHandle& ch) {set(ELEMENT_CELL_TYPE,ch.idx());}
    inline void set(const HalfFaceHandle& hfh) {set(ELEMENT_FACE_TYPE,hfh.idx());}
    inline void set(const HalfEdgeHandle& heh) {set(ELEMENT_EDGE_TYPE,heh.idx());}
    inline void set(const MeshElement& elem) {set(elem._type,elem._idx);}

    inline bool is_valid() const {return _type!=ELEMENT_NULL_TYPE&&_idx!=-1;}
    inline bool is_vertex() const {return _type==ELEMENT_VERTEX_TYPE;}
    inline bool is_edge() const {return _type==ELEMENT_EDGE_TYPE;}
    inline bool is_face() const {return _type==ELEMENT_FACE_TYPE;}
    inline bool is_cell() const {return _type==ELEMENT_CELL_TYPE;}

    inline VertexHandle vh() const {return VertexHandle(is_vertex()? _idx : -1);}
    inline HalfEdgeHandle heh() const {return HalfEdgeHandle(is_edge()? _idx : -1);}
    inline HalfFaceHandle hfh() const {return HalfFaceHandle(is_face()? _idx : -1);}
    inline EdgeHandle eh() const {return is_edge()? HalfEdgeHandle(_idx).edge_handle() : EdgeHandle(-1);}
    inline FaceHandle fh() const {return is_face()? HalfFaceHandle(_idx).face_handle() : FaceHandle(-1);}
    inline CellHandle ch() const {return CellHandle(is_cell()? _idx : -1);}

    inline constexpr char dim() const {return _type;}
    inline constexpr int idx() const {return _idx;}
    inline constexpr MeshElementType type() const {return _type;}

    constexpr bool operator<(const MeshElement& elem) const { return (this->dim() < elem.dim()); }
    constexpr bool operator>(const MeshElement& elem) const { return (this->dim() > elem.dim()); }
    constexpr bool operator==(const MeshElement& elem) const { return elem._idx == this->_idx && this->_type == elem._type; }
    constexpr bool operator!=(const MeshElement& elem) const { return elem._idx != this->_idx || this->_type != elem._type;  }

    // Use for Abseil Hash
    template <typename H>
    friend H AbslHashValue(H h, const MeshElement& m)
    {
        return H::combine(std::move(h), (char)(m._type), m._idx);
    }

    // Returns whether this element is part of an elements boundary or vice-versa. Dimensinality must be different.
    bool is_incident(const TetMeshCache& cache, const MeshElement& elem) const;

    MeshElement common_smaller_incident_element(const TetMeshCache& cache, const MeshElement& elem) const;

    // Returns all indicent cells. If the element itself is a cell, returns itself
    std::vector<CellHandle> get_incident_cells(const TetMeshCache& cache) const;

    friend std::ostream& operator<<(std::ostream &os, const MeshElement &elem) {
        switch (elem.type())
        {
        case ELEMENT_CELL_TYPE: os << "Element(CELL"; break;
        case ELEMENT_FACE_TYPE: os << "Element(FACE"; break;
        case ELEMENT_EDGE_TYPE: os << "Element(EDGE"; break;
        case ELEMENT_VERTEX_TYPE: os << "Element(VERTEX"; break;
        default: return os << "Element(NULL, -1)"; break;
        }
        return os << ", " << elem._idx << ")";
    }

private:
    MeshElementType _type = ELEMENT_NULL_TYPE;
    int _idx = -1;

}; // class MeshElement

// Required: m1 ~ ch, m2 ~ ch
MeshElement getElementBetweenInCell(const TetMeshCache& cache, const CellHandle ch, const MeshElement& m1, const MeshElement& m2);

} // namespace HexHex

