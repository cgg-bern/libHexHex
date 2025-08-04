#pragma once

#include <HexHex/Utils/Typedefs.hh>
#include <array>

namespace HexHex
{

/*
 * Stores the 4 Vertices and Parameters of a Tet
 */
class CellIGM
{
public:
    CellIGM()
    {}

    inline Parameter& operator[](VertexHandle vh)
    {
        assert(vh.is_valid());
        int idx = -1;
        for (int i = 0; i < 4; ++i)
            if (vhs[i] == vh)
                return params[i];
            else if (!vhs[i].is_valid())
                idx = i;
        assert(idx >= 0);
        vhs[idx] = vh;
        return params[idx];
    }

    inline const Parameter& operator[](VertexHandle vh) const
    {
        return at(vh);
    }

    inline const Parameter& at(VertexHandle vh) const
    {
        assert(hasValidVertices());
        if (vh == vhs[0]) return params[0];
        if (vh == vhs[1]) return params[1];
        if (vh == vhs[2]) return params[2];
        if (vh == vhs[3]) return params[3];
        throw std::runtime_error("BAD VERTEX: "
            + std::to_string(vh.idx())
            + " for tet with vertices " + std::to_string(vhs[0].idx())
            + ", " + std::to_string(vhs[1].idx())
            + ", " + std::to_string(vhs[2].idx())
            + ", " + std::to_string(vhs[3].idx())
        );
    }

private:
    std::array<VertexHandle, 4> vhs;
    std::array<Parameter, 4> params;

    inline bool hasValidVertices() const
    {
        for (int i = 0; i < 4; ++i) if (!vhs[i].is_valid()) return false;
        return true;
    }
};

using IntegerGridMap = OVM::CellPropertyT<CellIGM>;
using IGM = IntegerGridMap;

}

