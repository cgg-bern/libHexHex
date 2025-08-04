#include <HexHex/Predicates/GeneratorPredicates.hh>
#include <HexHex/Commons/Direction.hh>
#include <HexHex/Commons/MeshElement.hh>
#include <HexHex/Predicates/ExactPredicates.hh>
#include <HexHex/Utils/Utils.hh>

HexHex::MeshElement HexHex::computeVertexGeneratorElement(const TetMeshCache& mesh, const CellHandle ch,
                                          const std::array<VertexHandle,4>& vhs, const std::array<Parameter, 4> &params,
                                                          const Parameter& param)
{
#define VERTEX(i) vhs[i]
#define PARAM(i) params[i]
#define HALFEDGE(i,j) mesh.tetmesh.find_halfedge(VERTEX(i), VERTEX(j))
#define HALFFACE(i,j,k) mesh.tetmesh.find_halfface_in_cell({VERTEX(i), VERTEX(j), VERTEX(k)}, ch)

    assert(HALFEDGE(0,2).is_valid());
    assert(HALFEDGE(0,1).is_valid());
    assert(HALFEDGE(1,2).is_valid());
    assert(HALFEDGE(0,3).is_valid());
    assert(HALFEDGE(2,3).is_valid());
    assert(HALFEDGE(1,3).is_valid());

    assert(ori3d(PARAM(0), PARAM(1), PARAM(2), PARAM(3)) == ORI_ABOVE);

    // Get the 4 vertices of the tet and their parameters in the chart of the tet
    for (char i = 0; i <= 3; ++i) {
        if (PARAM(i) == param) {
            return VERTEX(i);
        }
    }

    auto oris = std::array<ORIENTATION, 4>();
    oris[0] = ori3d(PARAM(0), PARAM(1), PARAM(2), param); if (oris[0] == ORI_BELOW) return MeshElement();
    oris[1] = ori3d(PARAM(0), PARAM(2), PARAM(3), param); if (oris[1] == ORI_BELOW) return MeshElement();
    oris[2] = ori3d(PARAM(0), PARAM(3), PARAM(1), param); if (oris[2] == ORI_BELOW) return MeshElement();
    oris[3] = ori3d(PARAM(1), PARAM(3), PARAM(2), param); if (oris[3] == ORI_BELOW) return MeshElement();

    auto zeros = std::vector<char>();
    zeros.reserve(2);
    for (char i = 0; i <= 3; ++i) {if (oris[i] == ORI_ZERO) {zeros.push_back(i);}}

    if (zeros.size() == 2) {
        if (zeros[0] == 0 && zeros[1] == 1) {return HALFEDGE(0,2);}
        else if (zeros[0] == 0 && zeros[1] == 2) {return HALFEDGE(0,1);}
        else if (zeros[0] == 0 && zeros[1] == 3) {return HALFEDGE(1,2);}
        else if (zeros[0] == 1 && zeros[1] == 2) {return HALFEDGE(0,3);}
        else if (zeros[0] == 1 && zeros[1] == 3) {return HALFEDGE(2,3);}
        else if (zeros[0] == 2 && zeros[1] == 3) {return HALFEDGE(1,3);}
        assert(false);
    }
    if (zeros.size() == 1) {
        if (zeros[0] == 0) {return HALFFACE(0,1,2);}
        else if (zeros[0] == 1) {return HALFFACE(0,2,3);}
        else if (zeros[0] == 2) {return HALFFACE(0,3,1);}
        else if (zeros[0] == 3) {return HALFFACE(1,3,2);}
    }
    if (zeros.size() == 0) {
        return ch;
    }
    return MeshElement();
}

HexHex::MeshElement HexHex::computePropellerHolderElement(
    const TetMeshCache& mesh,
    const CellHandle ch, const std::array<VertexHandle,4>& vhs, const std::array<Parameter,4>& params,

    const MeshElement& generator,
    const Parameter& from, const Direction dir)
{
    assert(computeVertexGeneratorElement(mesh,ch,vhs,params,from).is_valid());

    auto VERTEX_PARAM = [&](VertexHandle vh) -> const Parameter&
    {
        for (int i = 0; i < 4; ++i)
            if (vhs[i] == vh)
                return params[i];
        assert(false);
        throw std::runtime_error("IMPOSSIBLE VERTEX");
    };

    auto HALFFACE_PARAMS = [&](HalfFaceHandle hfh) -> std::array<Parameter,3>
    {
        auto hf_vhs = mesh.get_halfface_vertices(hfh);
        return {VERTEX_PARAM(hf_vhs[0]),VERTEX_PARAM(hf_vhs[1]),VERTEX_PARAM(hf_vhs[2])};
    };

    const Parameter to = from + dir;

    if (generator.is_cell())
    {
        return generator; // from cell into cell
    }
    else if (generator.is_face())
    {
        HalfFaceHandle hfh = mesh.get_incident_halfface_in_cell(ch, generator.fh());
        auto ori = ori3d(HALFFACE_PARAMS(hfh), to);
        if (ori == ORI_ABOVE) {return ch;} // from face into cell
        else if (ori == ORI_ZERO) {return generator;} // from face along face
        else {return MeshElement();}
    }
    else if (generator.is_edge())
    {
        const auto& hfhs = mesh.get_incident_halffaces_in_cell(ch, generator.eh());
        assert(hfhs.size()==2u);
        std::array<ORIENTATION,2> oris;
        oris[0] = (ori3d(HALFFACE_PARAMS(hfhs[0]), to)); if (oris[0] == ORI_BELOW) return MeshElement();
        oris[1] = (ori3d(HALFFACE_PARAMS(hfhs[1]), to)); if (oris[1] == ORI_BELOW) return MeshElement();

        auto zeros = std::vector<char>();
        zeros.reserve(2);
        if (oris[0] == ORI_ZERO) zeros.push_back(0);
        if (oris[1] == ORI_ZERO) zeros.push_back(1);

        if (zeros.size() == 0) {return ch;} // from edge into cell
        if (zeros.size() == 2)
        {
            // from edge along edge (ensure halfedge matches direction)
            auto heh = generator.heh();
            auto hehDir = Direction(VERTEX_PARAM(mesh.tetmesh.to_vertex_handle(heh))
                                    - VERTEX_PARAM(mesh.tetmesh.from_vertex_handle(heh)));
            if (dir != hehDir) {heh = heh.opposite_handle();}
            return heh;
        }
        assert(zeros.size()==1);
        return hfhs[zeros[0]]; // from edge into face
    }
    else
    {
        assert(generator.is_vertex());
        auto hfhs = mesh.get_incident_halffaces_in_cell(ch, generator.vh());
        assert(hfhs.size()==3);
        std::array<ORIENTATION,3> oris;
        oris[0] = (ori3d(HALFFACE_PARAMS(hfhs[0]), to)); if (oris[0] == ORI_BELOW) return MeshElement();
        oris[1] = (ori3d(HALFFACE_PARAMS(hfhs[1]), to)); if (oris[1] == ORI_BELOW) return MeshElement();
        oris[2] = (ori3d(HALFFACE_PARAMS(hfhs[2]), to)); if (oris[2] == ORI_BELOW) return MeshElement();

        auto zeros = std::vector<char>();
        zeros.reserve(2);
        for (char i = 0; i <= 2; ++i) {if (oris[i] == ORI_ZERO) {zeros.push_back(i);}}

        if (zeros.size() == 0) {return ch;} // vertex into cell
        if (zeros.size() == 1) {return hfhs[zeros[0]];} // vertex into face

        // vertex into edge (ensure halfedge matches direction)
        assert(zeros.size()==2);
        MeshElement holder = MeshElement(hfhs[zeros[0]]).common_smaller_incident_element(mesh, hfhs[zeros[1]]);
        if (mesh.tetmesh.from_vertex_handle(holder.heh()) != generator.vh()) holder.set(holder.heh().opposite_handle());
        return holder;
    }
}

HexHex::MeshElement HexHex::computePropellerBladeCasingElement(
    const TetMeshCache& mesh,
    const CellHandle ch, const std::array<VertexHandle,4>& vhs, const std::array<Parameter,4>& params,
    const MeshElement& holder, const Parameter &from,
    const Direction dir1, const Direction dir2)
{

    assert(computePropellerHolderElement(
        mesh, ch, vhs, params,
        computeVertexGeneratorElement(mesh, ch, vhs, params, from),
        from, dir1
    ).is_valid());

    auto VERTEX_PARAM = [&](VertexHandle vh) -> const Parameter&
    {
        for (uint i = 0; i < 4; ++i)
            if (vhs[i] == vh)
                return params[i];
        assert(false);
        throw std::runtime_error("IMPOSSIBLE VERTEX");
    };

    auto HALFFACE_PARAMS = [&](HalfFaceHandle hfh) -> std::array<Parameter,3>
    {
        auto hf_vhs = mesh.get_halfface_vertices(hfh);
        return {VERTEX_PARAM(hf_vhs[0]),VERTEX_PARAM(hf_vhs[1]),VERTEX_PARAM(hf_vhs[2])};
    };

    if (holder.is_cell()) {return holder;}

    const auto& to = from + dir2;

    if (holder.is_face())
    {
        auto ori = ori3d(HALFFACE_PARAMS(mesh.get_incident_halfface_in_cell(ch, holder.fh())), to);
        if (ori == ORI_ABOVE) {return ch;}
        if (ori == ORI_ZERO) {return holder;}
       return MeshElement();
    }
    else
    {
        assert(holder.is_edge());
        const auto& hfhs = mesh.get_incident_halffaces_in_cell(ch, holder.eh());
        assert(hfhs.size()==2);
        auto ori1 = ori3d(HALFFACE_PARAMS(hfhs[0]), to); if (ori1 == ORI_BELOW) return MeshElement();
        auto ori2 = ori3d(HALFFACE_PARAMS(hfhs[1]), to); if (ori2 == ORI_BELOW) return MeshElement();
        assert(ori1 == ORI_ABOVE || ori2 == ORI_ABOVE);
        if (ori1 == ORI_ZERO) {return hfhs[0];}
        if (ori2 == ORI_ZERO) {return hfhs[1];}
        assert(ori1 == ORI_ABOVE && ori2 == ORI_ABOVE);
        return ch;
    }

    assert(false);
    return MeshElement();
}
