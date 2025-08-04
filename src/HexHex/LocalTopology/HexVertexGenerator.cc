#include <HexHex/LocalTopology/HexVertexGenerator.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Preprocessing/Parametrization.hh>
#include <HexHex/LocalTopology/ConstantIndices.hh>
#include <HexHex/Predicates/ExactPredicates.hh>
#include <HexHex/Predicates/TracingPredicates.hh>

namespace HexHex
{

    HexVertexGenerator::HexVertexGenerator(HexExtractor& hexEx, VertexHandle generator, VertexHandle hexVertex) :
        hexEx(hexEx),
        generator(generator),
        vh_range(hexVertex, hexVertex),
        integer_aligned(true),
        global_offsets({0,0,0}),
        type(GENERAL)
    {
    }

    HexVertexGenerator::HexVertexGenerator(HexExtractor& hexEx, EdgeHandle generator, std::pair<VertexHandle, VertexHandle> vh_range) :
        hexEx(hexEx),
        generator(generator),
        vh_range(vh_range),
        integer_aligned(hexEx.edgeIsOnIntegerLattice(generator)),
        global_offsets({0,0,0}),
        type(GENERAL)
    {
    }

    HexVertexGenerator::HexVertexGenerator(HexExtractor& hexEx, HalfFaceHandle generator, std::pair<VertexHandle, VertexHandle> vh_range) :
        hexEx(hexEx),
        generator(generator),
        vh_range(vh_range),
        global_offsets({0,0,0})
    {
        if (hexEx.inputMesh.is_boundary(generator.opposite_handle()))
        {
            assert(hexEx.faceIsOnIntegerLattice(generator.face_handle()));

            type = BOUNDARY_FACE;
            integer_aligned = true;
            const auto& hfh_params = hexEx.parametrization.get(generator);
            for (char i = 0; i < 3; ++i)
            {
                if (hfh_params[0][i] == hfh_params[1][i] && hfh_params[0][i] == hfh_params[2][i])
                {
                    Direction dir(2*i); // positive direction along constant axis
                    if (ori3d(hfh_params, hfh_params[0] + dir) == ORI_BELOW) dir = -dir;
                    boundary_face_direction_into_cell = dir;
#ifndef NDEBUG
                    MeshElement tmp;
                    CellHandle ch = hexEx.inputMesh.incident_cell(generator);
                    assert(hexEx.computePropellerHolderElement(ch, generator, hexEx.getParametrization().getHexVertexParameter(ch,vh_range.first), dir, tmp) && tmp == ch);
#endif
                }
            }
        }
        else
        {
            type = GENERAL;
            integer_aligned = hexEx.faceIsOnIntegerLattice(generator.face_handle());
        }
    }

    HexVertexGenerator::HexVertexGenerator(HexExtractor& hexEx, CellHandle generator, std::pair<VertexHandle, VertexHandle> vh_range) :
        hexEx(hexEx),
        generator(generator),
        vh_range(vh_range),
        integer_aligned(false),
        global_offsets({0,0,0}),
        type(CELL)
    {
    }

    std::array<size_t,3> HexVertexGenerator::extractLocalTopology()
    {
        n_local_blades = n_local_connected_blades = 0;

        // For a cell, we store evrything implicitly
        if (type==CELL)
        {
            localCorners = CELL_GENERATOR_LOCAL_CORNERS;
            //resizeGlobalVectors();
            return {6u * nVertices(), 12u * nVertices(), 8u * nVertices()};
        }

        // For a boundary face (which must be integer aligned) we already know the structure based on the direction into cell
        if (type==BOUNDARY_FACE)
        {
            localCorners = BOUNDARY_FACE_GENERATOR_LOCAL_CORNERS;
            //resizeGlobalVectors();
            return {5u * nVertices(), 8u * nVertices(), 4u * nVertices()};
        }

        // For the other cases (singularities and "rare" intersections), let's do it step by step

        enumerateLocalPropellers();

        enumerateLocalPropellerBlades();

        enumerateLocalCorners();

        if (n_local_blades != n_local_connected_blades)
        {
            hexEx.printErr("HexHex: ", "Generator ", generator, " failed to connect all blades.");
        }
        //resizeGlobalVectors();
        return {(uint)localPropellers.size() * nVertices(), n_local_connected_blades/2 * nVertices(), (uint)localCorners.size() * nVertices()};
    }

    void HexVertexGenerator::enumerateLocalPropellers()
    {
        //ScopedStopWatch _{sw::enumerateLocalPropellers};

        assert(!generator.is_cell());

        auto checkAllDirections = [&](CellHandle ch) -> void
        {
            auto param = hexEx.getParametrization().getHexVertexParameter(ch, vh_range.first);
            for (auto dir : getAll6Directions()) {
                MeshElement holder;
                if (hexEx.computePropellerHolderElement(ch, generator, param, dir, holder) && hexEx.owner(holder)==ch) {
                    addLocalPropeller(holder, ch, dir);
                }
            }
        };

        auto& inputMesh = hexEx.inputMesh;

        if (generator.is_face())
            for (auto ch : hexEx.cache.get_incident_cells(generator.fh()))
                checkAllDirections(ch);
        else if (generator.is_edge())
            for (auto ec_it = inputMesh.ec_iter(generator.eh()); ec_it.is_valid(); ++ec_it)
                checkAllDirections(*ec_it);
        else if (generator.is_vertex())
            for (CellHandle ch : hexEx.cache.vertexCells[generator.vh()])
                checkAllDirections(ch);
        else
            assert(false);
    }

    void HexVertexGenerator::enumerateLocalPropellerBlades()
    {
        //ScopedStopWatch _{sw::enumerateLocalBlades};

        // Directions to blades
        auto checkAllOrthogonalDirections = [&](LocalPropeller& lp, CellHandle ch) -> void
        {
            //auto cch = hexEx.getParametrization().cellClass(ch);
            auto param = hexEx.getParametrization().getHexVertexParameter(ch, vh_range.first);
            auto dir1 = lp.getDirection(ch);
            const auto& dirs2 = getAllOrthogonalDirectionsCCW(dir1);
            std::array<MeshElement, 4> casings;
            bool hasBlade = false;

            // Ensure a ccw ordering of the blades with no gaps
            for (auto j = 0; j < 4; ++j) {
                MeshElement casing;
                if (hexEx.computePropellerBladeCasingElement(ch, lp.holder, param, dir1, dirs2[j], casing) && hexEx.owner(casing)==ch)
                {
                    n_local_blades++;
                    casings[j] = casing;
                    hasBlade = true;
                }
            }
            if (!hasBlade) return;
            char r = 3;
            while (casings[(r+3)%4].is_valid() || !casings[r].is_valid()) r = (r+3)%4;

            for (auto j = r; j < r+4 && casings[j%4].is_valid(); ++j) {
                lp.initBlade(*this, casings[j%4], ch, dirs2[j%4]);
            }
        };

        // Enumerate directions to blades
        for (auto& lp : localPropellers)
        {
            if (lp.holder.is_cell())
            {
                //auto cch = hexEx.getParametrization().cellClass(lp.holder.ch());
                n_local_blades += 4;
                for (const auto& dir2 : getAllOrthogonalDirectionsCCW(lp.getDirection(lp.holder.ch())))
                    lp.initBlade(*this, lp.holder, lp.holder.ch(), dir2);
            }
            else if (lp.holder.is_face())
                for (auto ch : hexEx.cache.get_incident_cells(lp.holder.fh()))
                    checkAllOrthogonalDirections(lp, ch);
            else if (lp.holder.is_edge())
                for (auto hec_it = hexEx.inputMesh.hec_iter(lp.holder.heh()); hec_it.is_valid(); ++hec_it)
                    checkAllOrthogonalDirections(lp, *hec_it);
            else
                assert(false);
        }

        // Rotate to blades
        for (uint i = 0; i < localPropellers.size(); ++i)
        {
            const auto lph1 = LocalPropellerHandle(i);
            auto& lp1 = localPropellers[i];
            for (int j = 0; j < lp1.nBlades(); ++j)
            {
                if (lp1.bladeIsConnected(j)) continue;

                const auto chFrom = lp1.getBladeDirections(j).getByIndex(0).first;
                const auto d2 = lp1.getBladeDirections(j).getByIndex(0).second;
                auto d1 = lp1.getDirection(chFrom);

                const auto& res = findBlade(lp1, chFrom, d1, d2);
                const LocalPropellerHandle lph2 = std::get<0>(res);
                const CellHandle chTo = std::get<1>(res);
                const Transition& transition = std::get<2>(res);
                if (!lph2.is_valid()) continue;
                auto& lp2 = localPropellers[lph2.idx()];

                n_local_connected_blades += 2;

                lp1.setBlade(j, lph2, chTo, transition);
                auto j2 = lp2.findBladeIndex(chTo, transition.transform(d1));
                assert(j2 >= 0);

                lp2.setBlade(j2, lph1, chFrom, transition.inverted());
            }
        }
    }

    std::tuple<LocalPropellerHandle, CellHandle, Transition> HexVertexGenerator::findBlade(const LocalPropeller& lp1, CellHandle ch, Direction d1, Direction d2)
    {
        std::tuple<LocalPropellerHandle, CellHandle, Transition> res;

        //auto& mesh = hexEx.inputMesh;
        auto lph2 = findLocalPropeller(ch, d2);
        auto hfh = HalfFaceHandle();
        auto accumulatedTransition = Transition();

        int i = 0;
        while (!lph2.is_valid()) {
            assert((d1 | d2) == 0);

            // Pick transition halfface
            //auto prevHfh = hfh;
            hfh = findBladeNextHalfFace(lp1, ch, hfh, d1, d2);

            assert(ori3d(hexEx.parametrization.get(hfh), hexEx.getParametrization().getHexVertexParameter(ch, vh_range.first) + d1) == ORI_ABOVE);
            assert(ori3d(hexEx.parametrization.get(hfh), hexEx.getParametrization().getHexVertexParameter(ch, vh_range.first) + d2) == ORI_BELOW);

            // Check bad, bad cases
            if (++i >= 10000) {hexEx.printErr("HexHex: Generator ", generator, " failed to rotate to blade after 10000 iterations."); return res;}
            if (!hfh.is_valid()) {hexEx.printErr("HexHex: Generator", generator, " failed to find next halfface to blade."); return res;}
            if (hexEx.inputMesh.is_boundary(hfh.face_handle())) {hexEx.printErr("HexHex: Generator", generator, " tried to trace through boundary face to blade."); return res;}

            // Update values according to transition into new cell.
            auto& transition = hexEx.parametrization.transition(hfh);
            d1 = transition.transform(d1);
            d2 = transition.transform(d2);
            hfh = hfh.opposite_handle();
            ch = hexEx.inputMesh.incident_cell(hfh);

            lph2 = findLocalPropeller(ch, d2);

            accumulatedTransition = transition * accumulatedTransition;

        }

        std::get<0>(res) = lph2;
        std::get<1>(res) = ch;
        std::get<2>(res) = accumulatedTransition;
        return res;
    }

    HalfFaceHandle HexVertexGenerator::findBladeNextHalfFace(const LocalPropeller& lp, CellHandle ch, HalfFaceHandle hfh, Direction d1, Direction d2)
    {
        const Parameter u = hexEx.getParametrization().getHexVertexParameter(ch, vh_range.first);
        return ::HexHex::findBladeNextHalfFace(
            hexEx.cache, ch,
            hexEx.cache.cellVertices[ch], hexEx.getParametrization().get(ch),
            generator, lp.holder, hfh, u, d1, d2
        );
    }

    void HexVertexGenerator::enumerateLocalCorners()
    {
        //ScopedStopWatch _{sw::enumerateLocalCorners};

        for (uint i1 = 0; i1 < localPropellers.size(); ++i1) {
            const auto& lp1 = localPropellers[i1];
            auto lph1 = LocalPropellerHandle(i1);
            uint m = lp1.nBlades();
            assert(m >= 2);
            bool boundary = !lp1.holder.is_cell() && hexEx.isBoundary(lp1.holder);
            for (uint j = 0; j < m; ++j)
            {
                if (j == m-1 && boundary) continue;
                const auto& lph2 = lp1.getBlade((j+0)%m);
                const auto& lph3 = lp1.getBlade((j+1)%m);
                const auto& lp2 = localPropellers[lph2.idx()];

                if (lp2.getNextBlade(lph3) == lph1)
                {
                    assert(localPropellers[lph3.idx()].getNextBlade(lph1) == lph2);
                    addLocalCorner(lph1, lph2, lph3);
                }
            }
        }
    }

    void HexVertexGenerator::pickTracingStart(const GPH& gph, CellHandle& ch, Direction& d1, Direction& d2, int& bladeIndex1)
    {
        //VertexHandle fromVh = gph.first;
        LPH fromLph = gph.lph;

        bladeIndex1 = 0;

        if (type==CELL)
        {
            ch = generator.ch();
            d1 = Direction(fromLph.idx());
            d2 = Direction::ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW[d1.idx()][bladeIndex1];
        }
        else if (type==BOUNDARY_FACE)
        {
            ch = hexEx.inputMesh.incident_cell(generator.hfh());
            d1 = (fromLph.idx()==0)? boundary_face_direction_into_cell : Direction::ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW[boundary_face_direction_into_cell.idx()][fromLph.idx()-1];
            d2 = (fromLph.idx()==0)? Direction::ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW[d1.idx()][bladeIndex1] : boundary_face_direction_into_cell;
            if (fromLph.idx() != 0) bladeIndex1 = 1;
        }
        else
        {
            const auto& img = localPropellers[fromLph.idx()].getBladeDirections(bladeIndex1).getByIndex(0);
            ch = img.first;
            d1 = localPropellers[fromLph.idx()].getDirection(ch);
            d2 = img.second;
        }
    }

    MeshElement HexVertexGenerator::getHolder(const LPH& lph) const
    {
        if (type == CELL) return generator;
        if (type == BOUNDARY_FACE)
        {
            return lph.idx()==0? hexEx.inputMesh.incident_cell(generator.hfh()) : generator;
        }
        return localPropellers[lph.idx()].holder;
    }

    size_t HexVertexGenerator::getValence(const LPH& lph) const
    {
        if (type == CELL) return 4;
        if (type == BOUNDARY_FACE)
        {
            return lph.idx()==0? 4 : 3;
        }
        return localPropellers[lph.idx()].nBlades();
    }

    EdgeHandle HexVertexGenerator::getHolderEdge(const HalfEdgeHandle& hexHeh) const
    {
        if (!generator.is_vertex() && !generator.is_edge()) return EdgeHandle(-1);
        assert(type == GENERAL);
        assert(hasVertex(hexEx.getOutputMesh().from_vertex_handle(hexHeh)));
        const auto& lph = hexEx.hexHalfEdgeLocalPropellers[hexHeh];
        return localPropellers[lph.idx()].holder.eh();
    }

    FaceHandle HexVertexGenerator::getCasingFace(const HalfEdgeHandle& hexHeh1, const HalfEdgeHandle& hexHeh2) const
    {
        assert(hasVertex(hexEx.getOutputMesh().from_vertex_handle(hexHeh1)));
        if (type == CELL) return FaceHandle(-1);
        const auto& lph1 = hexEx.hexHalfEdgeLocalPropellers[hexHeh1];
        const auto& lph2 = hexEx.hexHalfEdgeLocalPropellers[hexHeh2];
        if (type == BOUNDARY_FACE) return (lph1.idx()!=0&&lph2.idx()!=0)? generator.fh() : FaceHandle(-1);
        assert(type == GENERAL);
        return localPropellers[lph1.idx()].getCasing(lph2).fh();
    }

    LocalPropellerHandle HexVertexGenerator::addLocalPropeller(MeshElement holder, CellHandle ch, Direction dir)
    {
        //if (generator == VertexHandle(43)) std::cout << "NEW Add Propeller " << holder << std::endl;

        auto lp = LocalPropeller(*this, holder);
        lp.setDirection(ch, dir);

        // Add images in other incident cell charts if necessary
        if (holder.is_face())
        {
            // Set 2nd cell chart of face propeller
            auto hfh = holder.hfh();
            auto ch2 = hexEx.inputMesh.incident_cell(hfh.opposite_handle()); assert(ch2 != ch);
            if (ch2.is_valid())
            {
                auto dir2 = hexEx.parametrization.transition(hfh).transform(dir);
                lp.setDirection(ch2, dir2);
            }
        }
        else if (holder.is_edge())
        {
            auto heh = holder.heh();

            // 1st cell chart of halfedge propeller should be 1st incident cell
            auto ch1st = hexEx.inputMesh.incident_cell(*(hexEx.inputMesh.hehf_iter(heh)));
            if (ch != ch1st)
            {
                dir = hexEx.parametrization.getTransition(ch, ch1st, heh.edge_handle()).transform(dir);
                ch = ch1st;
            }

            // Set other cell charts of edge propeller
            auto transition = Transition();
            auto dir2 = dir;
            for (auto hehf_iter = hexEx.inputMesh.hehf_iter(heh); hehf_iter.is_valid(); ++hehf_iter)
            {
                auto hfh2 = *hehf_iter;
                auto ch2 = hexEx.inputMesh.incident_cell(hfh2);
                if (!ch2.is_valid()) break;
                dir2 = transition.transform(dir2);

                {MeshElement tmp; assert(hexEx.computePropellerHolderElement(ch2, generator, hexEx.getParametrization().getHexVertexParameter(ch2, vh_range.first), dir2,  tmp));}
                lp.setDirection(ch2, dir2);

                transition = hexEx.parametrization.transition((*(hehf_iter+1)).opposite_handle());
            }
        }

        localPropellers.push_back(lp);
        return LocalPropellerHandle(localPropellers.size()-1);
    }

    LocalCornerHandle HexVertexGenerator::addLocalCorner(LocalPropellerHandle lph1, LocalPropellerHandle lph2, LocalPropellerHandle lph3)
    {
        // Check if corner already exists
        for (uint i = 0; i < localCorners.size(); ++i)
        {
            const auto& lc = localCorners[i];
                 if (lc[0] == lph1) {if (lc[1] == lph2) {assert(lc[2] == lph3); return i;}}
            else if (lc[1] == lph1) {if (lc[2] == lph2) {assert(lc[0] == lph3); return i;}}
            else if (lc[2] == lph1) {if (lc[0] == lph2) {assert(lc[1] == lph3); return i;}}
        }

        // For consistency: Start with smallest index propeller
        int idx = min3(lph1.idx(), lph2.idx(), lph3.idx());
             if (idx==lph1.idx()) localCorners.push_back({lph1, lph2, lph3});
        else if (idx==lph2.idx()) localCorners.push_back({lph2, lph3, lph1});
        else if (idx==lph3.idx()) localCorners.push_back({lph3, lph1, lph2});
        return LocalCornerHandle(localCorners.size()-1);
    }

    LocalPropellerHandle HexVertexGenerator::findLocalPropeller(const CellHandle &cch, const Direction &dir) const
    {
        if (type==CELL && generator.ch() == cch) return dir.idx();

        if (type==BOUNDARY_FACE && hexEx.inputMesh.incident_cell(generator.hfh()) == cch)
        {
            if (dir == boundary_face_direction_into_cell) return LocalPropellerHandle(0);
            const auto& j = Direction::ORTHOGONAL_INDICES_PER_DIRECTION_CCW[boundary_face_direction_into_cell.idx()][dir.idx()];
            if (j < 0 || j > 3)
            {
                std::cerr << "HexHex: Boundary Face Generator Blade Index Error: " << (int)j << std::endl;
            }
            assert(j >= 0 && j <= 3);
            return LocalPropellerHandle(j + 1);
        }

        for (uint i = 0; i < localPropellers.size(); ++i)
            if (localPropellers[i].hasDirection(cch, dir))
                return LocalPropellerHandle(i);
        return LocalPropellerHandle(-1);
    }

    LocalCornerHandle HexVertexGenerator::findLocalCorner(const LocalPropellerHandle &lph1, const LocalPropellerHandle &lph2) const
    {
        if (type == CELL)
        {
            if (CELL_GENERATORS_LOCAL_CORNERS_FROM_PROPELLER_PAIR.contains({lph1,lph2})) {
                return CELL_GENERATORS_LOCAL_CORNERS_FROM_PROPELLER_PAIR.at({lph1,lph2});
            }
            return LCH(-1);
        }

        for (uint i = 0; i < nLocalCorners(); ++i)
        {
            const auto& lc = localCorners[i];
            for (int c : {0,1,2})
            {
                if (lc[c]==lph1 && lc[(c+1)%3]==lph2)
                {
                    return LCH(i);
                }
            }
        }
        return LCH(-1);
    }

    LocalCornerAxis HexVertexGenerator::findLocalCornerOpposite(const LocalCornerAxis& lca) const
    {
        const LocalCorner& lc = getLocalCorner(lca.lch);
        uint8_t axis = (lc[0]==lca.forward)? 0 : (lc[1]==lca.forward)? 1 : 2;
        assert(lc[axis] == lca.forward);
        const LPH opp_forward = lc[(axis+1)%3];
        const LCH opp_lch = findLocalCorner(opp_forward, lca.forward);
        assert(lca.lch != opp_lch);
        assert(lca.forward != opp_forward);
        return {opp_lch, opp_forward};
    }

    int HexVertexGenerator::nBlades(const LPH& lph) const
    {
        if (type == CELL) return 4;
        else if (type == BOUNDARY_FACE) return (lph.idx()==0)? 4 : 3;
        else
        {
            assert(type==GENERAL);
            return localPropellers[lph.idx()].nBlades();
        }
    }

    void HexVertexGenerator::getDirections(const LPH& lph1, const LPH& lph2, CellHandle& ch, Direction& d1, Direction& d2)
    {
        if (type == CELL)
        {
            ch = generator.ch();
            d1 = lph1.idx();
            d2 = lph2.idx();
        }
        else if (type == BOUNDARY_FACE)
        {
            ch = hexEx.inputMesh.incident_cell(generator.hfh());
            d1 = (lph1.idx()==0)? boundary_face_direction_into_cell :
                     Direction::ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW[boundary_face_direction_into_cell.idx()][findBladeIndex(0, lph1)];
            d2 = (lph2.idx()==0)? boundary_face_direction_into_cell :
                     Direction::ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW[boundary_face_direction_into_cell.idx()][findBladeIndex(0, lph2)];
        }
        else
        {
            const auto& img =localPropellers[lph1.idx()].getBladeDirections(localPropellers[lph1.idx()].findBladeIndex(lph2)).getByIndex(0);
            ch = img.first;
            d1 = localPropellers[lph1.idx()].getDirection(ch);
            d2 = img.second;
        }
    }

    int32 HexVertexGenerator::findBladeIndex(const LPH& lph1, const LPH& lph2) const
    {
        if (type == CELL) return Direction::ORTHOGONAL_INDICES_PER_DIRECTION_CCW[lph1.idx()][lph2.idx()];
        else if (type == BOUNDARY_FACE) {return BOUNDARY_FACE_GENERATOR_BLADE_INDICES_PER_LPH_CCW[lph1.idx()][lph2.idx()];}
        else
        {
            assert(type==GENERAL);
            return localPropellers[lph1.idx()].findBladeIndex(lph2);
        }
    }

    int32 HexVertexGenerator::getBladeIndex(const LPH& lph1, const CellHandle& cch, const Direction& d1, const Direction& d2) const
    {
        assert(findLocalPropeller(cch, d1)==lph1);
        if (type==CELL)
        {
            return Direction::ORTHOGONAL_INDICES_PER_DIRECTION_CCW[lph1.idx()][d2.idx()];
        }
        else if (type==BOUNDARY_FACE)
        {
            if (lph1.idx()==0) return Direction::ORTHOGONAL_INDICES_PER_DIRECTION_CCW[d1.idx()][d2.idx()];
            if (d2 == boundary_face_direction_into_cell) return 1;
            const auto& o2 = Direction::ORTHOGONAL_INDICES_PER_DIRECTION_CCW[d1.idx()][d2.idx()];
            const auto& oR = Direction::ORTHOGONAL_INDICES_PER_DIRECTION_CCW[d1.idx()][boundary_face_direction_into_cell.idx()];
            assert(((4+(o2-oR))%4==1)||((4+o2-oR)%4==-1));
            return 1 + ((4+(o2-oR))%4);
        }
        else
        {
            return localPropellers[lph1.idx()].findBladeIndex(cch, d2);
        }
    }

    LocalPropellerHandle HexVertexGenerator::getBlade(const LPH& lph, const int& j) const
    {
        if (type == CELL) return Direction::ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW[lph.idx()][j];
        else if (type == BOUNDARY_FACE) {return BOUNDARY_FACE_GENERATOR_BLADE_LPHS_PER_LPH_CCW[lph.idx()][j];}
        else
        {
            assert(type==GENERAL);
            return localPropellers[lph.idx()].getBlade(j);
        }
    }

    bool HexVertexGenerator::isBoundary() const {return hexEx.isBoundary(generator);}
}
