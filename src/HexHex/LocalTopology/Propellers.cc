#include <HexHex/LocalTopology/Propellers.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/LocalTopology/HexVertexGenerator.hh>

namespace HexHex
{
    LocalPropeller::LocalPropeller(HexVertexGenerator& generator, MeshElement holder) :
        holder(holder)
    {
    }

    Direction LocalPropeller::getDirection(CellHandle cch)
    {
        assert(images.contains(cch));
        return images.at(cch);
    }

    bool LocalPropeller::hasDirection(CellHandle cch, Direction dir) const
    {
        Direction cDir(0);
        return (images.get(cch, cDir) && cDir==dir);
    }

    void LocalPropeller::initBlade(HexVertexGenerator& generator, MeshElement casing, CellHandle ch, Direction dir)
    {
        blades_casings.push_back(casing);
        FlatMap<CellHandle, Direction> from;
        from.set(ch, dir);

        // For inner face casings, store the blade direction in the 2nd cell chart
        if (casing.is_face() && !generator.hexEx.isBoundary(casing.fh()))
        {
            auto hfh = generator.hexEx.getCache().get_incident_halfface_in_cell(ch, casing.fh());
            auto ch2 = generator.hexEx.getInputMesh().incident_cell(hfh.opposite_handle());
            assert(ch2.is_valid() && ch != ch2);
            from.set(ch2, generator.hexEx.getParametrization().transition(hfh).transform(dir));
        }

        blades_from.push_back(from);

        assert(blades_casings.size() == blades_from.size());
    }

    void LocalPropeller::setBlade(int32 i, LocalPropellerHandle lph, CellHandle to, Transition transition)
    {
        assert(i >= 0);
        assert(lph.is_valid());
        if (i >= (int)blades.size()) {blades.resize(nBlades()); blades_to.resize(nBlades());}
        assert(!bladeIsConnected(i));
        blades[i] = lph;
        blades_to[i] = {to, transition};
    }

    int32 LocalPropeller::findBladeIndex(CellHandle cch, Direction dir) const
    {
        Direction tmp(0);
        for (int j = 0; j < nBlades(); ++j)
        {
            const auto& from = blades_from[j];
            for (size_t l = 0; l < from.size(); ++l)
            {
                if (from.getByIndex(l).first == cch && from.getByIndex(l).second == dir)
                {
                    return j;
                }
            }
        }
        return -1;
    }

    int32 LocalPropeller::findBladeIndex(LocalPropellerHandle lph) const
    {
        for (int j = 0; j < nBlades(); ++j)
            if (blades[j] == lph)
                return j;
        return -1;
    }
}
