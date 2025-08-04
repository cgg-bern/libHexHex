
#include <HexHex/PiecewiseLinear/PLMesh.hh>
#include <HexHex/Predicates/ExactPredicates.hh>
#include <HexHex/Predicates/TracingPredicates.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/Threading.hh>
#include <HexHex/GlobalTopology/HexEdgeExtractor.hh>
#include <HexHex/LocalTopology/HexVertexGenerator.hh>
#include <HexHex/HexExtractor.hh>

namespace HexHex
{

void HexEdgeExtractor::tracePropellers()
{
    ScopedStopWatch _{sw::traceHexEdges};

    global_propellers_visited_in_parallel_edge_extraction = std::vector<std::atomic<bool>>(hexEx.numGlobalPropellers());
    longest_ple_segment = 0;

    // Threaded Edge Tracing
    std::chrono::steady_clock::time_point t1, t2;
    t1 = std::chrono::steady_clock::now();

    // Get Number of Threads
    size_t num_threads = Threading::getNumThreads();
    thread_hex_edges = std::vector<std::vector<HexHalfEdge>>(num_threads);
    size_t num_traced_edges = 0;

    // Each Thread gets its vertices to process
#pragma omp parallel reduction(+:num_traced_edges)
    {
        // Get current thread id
        size_t t_id = Threading::getThreadID();

#pragma omp for nowait
        for (auto i = 0; i < static_cast<int>(hexEx.getOutputMesh().n_vertices()); ++i)
        {
            auto n_local_props = static_cast<int>(hexEx.getHexVertexGenerator(VertexHandle(i)).nLocalPropellers());
            for (int j = 0; j < n_local_props; ++j)
            {
                GlobalPropellerHandle fromGph{VertexHandle(i), LPH(j)};

                if (num_threads >= 2)
                {
                    const int32 from_he_idx = hexEx.toGlobalIndex(fromGph);

                    // Check if propeller has been visited
                    if (global_propellers_visited_in_parallel_edge_extraction[from_he_idx].load(std::memory_order_acquire))
                    {
                        continue;
                    }

                    // Trace Edge
                    bool success = tracePropeller(fromGph, t_id);

                    if (success)
                    {
                        // Get Hex Edge
                        auto& hex_edge = thread_hex_edges[t_id].back();
                        assert(hex_edge.gph1 == fromGph);
                        const int32 to_he_idx = hexEx.toGlobalIndex(hex_edge.gph2);

                        // Mark both propellers as visited
                        bool f = false;
                        if (global_propellers_visited_in_parallel_edge_extraction[from_he_idx]
                                .compare_exchange_strong(f, true, std::memory_order_acq_rel))
                        {
                            global_propellers_visited_in_parallel_edge_extraction[to_he_idx]
                                .store(true, std::memory_order_release);
                        }

                        num_traced_edges++;
                    }
                }
                else
                {
                    // Trace Edge and add edge if using single thread
                    assert(t_id == 0);
                    if (tracePropeller(fromGph, t_id))
                    {
                        num_traced_edges++;
                    }
                }
            }
        }
    }

    t2 = std::chrono::steady_clock::now();
    hexEx.printLog("HexHex: Successfully traced ", num_traced_edges, " Edges in "
                   ,std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count()," ms ("
                   ,(num_traced_edges-hexEx.numGlobalPropellers()/2), " were traced both ways)");


    // If using multiple threads, In a 2nd step, we add the hex edges to the actual mesh and write to the oppsosites vectors
    if (num_threads >= 2)
    {
        t1 = std::chrono::steady_clock::now();
        size_t num_hex_edges_stored = 0;
        size_t num_hex_edges_skipped = 0;

        for (size_t t_id = 0; t_id < num_threads; ++t_id)
        {
            num_hex_edges_stored += thread_hex_edges[t_id].size();

            for (size_t i = 0; i < thread_hex_edges[t_id].size(); ++i)
            {
                if (!hexEx.isGlobalPropellerConnected(thread_hex_edges[t_id][i].gph1))
                {
                    addHexEdge(thread_hex_edges[t_id][i]);
                }
                else
                {
                    num_hex_edges_skipped++;
                }
            }
        }

        t2 = std::chrono::steady_clock::now();
        hexEx.printLog("HexHex: Storing ", num_hex_edges_stored, " Hex Edges, Skipped ", num_hex_edges_skipped);
        hexEx.printLog("HexHex: Combined per Thread Tracing Results in ",
        std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count(), " ms");
    }

    hexEx.printLog("HexHex: Extracted ", hexEx.getOutputMesh().n_edges(), " Edges (",
     100.0*(double)hexEx.getOutputMesh().n_edges()/(hexEx.numGlobalPropellers()/2), " %)"
    );

    if (hexEx.config.extract_piecewise_linear_edges)
    {
        hexEx.printLog("HexHex: Longest Piecewise Linear Edge Segment goes through ", longest_ple_segment, " Tets");
    }
}

bool HexEdgeExtractor::tracePropeller(const GlobalPropellerHandle& fromGph, size_t thread_id)
{
    auto findGlobalPropeller = [&](const CellHandle& ch, const Parameter& u, const Direction& d1) -> GPH
    {
        VertexHandle vh = hexEx.getParametrization().getHexVertexInCellWithParameter(ch, u, hexEx.config);
        if (!vh.is_valid())
        {
            //std::cerr << "Failed to find " << u << " in " << ch << std::endl;
            return {vh, LPH(-1)};
        }
        const auto& lph = hexEx.getHexVertexGenerator(vh).findLocalPropeller(ch, d1);
        assert(lph.is_valid());
        return {vh, lph};
    };

    auto computeTriangleIntersection = [](const Parameter& u, const Direction& d1, const Parameter& param0, const Parameter& param1, const Parameter& param2) -> Parameter
    {
        const auto& v1 = param1 - param0;
        const auto& v2 = param2 - param0;

        const auto& h = d1.vector().cross(v2);
        const auto& a = v1.dot(h);

        const auto& f = 1.0 / a;
        const auto& s = u - param0;

        const auto& q = s.cross(v1);

        const double& t = f * v2.dot(q);
        return u + d1 * t;
    };

    const VertexHandle fromVh = fromGph.vh;
    auto& fromGen = hexEx.getHexVertexGenerator(fromVh);
    //const LPH fromLph = fromGph.lph;

    // If only using a single thread, its safe to check at the start of tracing if the propeller is already connected
    if (Threading::getNumThreads() == 1 && hexEx.isGlobalPropellerConnected(fromGph))
    {
        return false;
    }

    // Start with propeller and a blade in the corresponding cell
    int bladeIndex1 = 0;
    CellHandle ch; Direction d1(0); Direction d2(0);

    fromGen.pickTracingStart(fromGph, ch, d1, d2, bladeIndex1);
    const CellHandle startCh = ch;

    assert(ch.is_valid());
    assert((d1|d2)==0);

    Parameter from = hexEx.getParametrization().getHexVertexParameter(ch, fromVh);
    Transition accumulatedTransition;
    HalfFaceHandle hfh;

    auto toGph = findGlobalPropeller(ch, from + d1, -d1);

    // Piecewise Linear Edge Extraction
    std::vector<PLEPoint> ple_points;

    // Main Tracing Loop - Trace from cell through faces until target is reached
    int i = 0;
    MeshElement previousTraceThrough = fromGen.generator;
    while(!toGph.vh.is_valid())
    {
        if (hexEx.isGlobalPropellerConnected(fromGph))
        {
            return false;
        }
        assert((d1 | d2) == 0);

        // Pick transition halfface
        //auto prevHfh = hfh;
        MeshElement traceThrough; // (<- only matters for piecewise linear extraction)
        hfh = findOppositeNextHalfFace(
            hexEx.cache,
            ch, hexEx.cache.cellVertices[ch], hexEx.parametrization.get(ch),
            hfh, from, d1, d2, traceThrough, hexEx.config.extract_piecewise_linear_edges);

        // Check bad, bad cases
        if (++i >= 10000) {hexEx.printErr("HexHex: Failed to trace to opposite after 10000 iterations from ", fromGen.generator, "."); break;}
        if (!hfh.is_valid()) {hexEx.printErr("HexHex: Failed to find next halfface to opposite from ", fromGen.generator, "."); break;}
        if (hexEx.inputMesh.is_boundary(hfh.face_handle())) {
            {hexEx.printErr("HexHex: Tried to trace through boundary face to opposite from ", fromGen.generator, ".");}
            break;
        }

        assert(ori3d(hexEx.parametrization.get(hfh), from) == ORI_ABOVE);
        assert(ori3d(hexEx.parametrization.get(hfh), from + d1) == ORI_BELOW);

        // Update values according to transition into new cell.
        auto& transition = hexEx.parametrization.transition(hfh);
        d1 = transition.transform(d1);
        d2 = transition.transform(d2);
        hfh = hfh.opposite_handle();
        ch = hexEx.inputMesh.incident_cell(hfh);
        from = transition.transform_point(from);
        toGph = findGlobalPropeller(ch, from + d1, -d1);
        accumulatedTransition = transition * accumulatedTransition;

        // Piecewise Linear Edge Extraction - Add middle Position
        if (hexEx.config.extract_piecewise_linear_edges)
        {
            assert(traceThrough.is_valid());
            if (traceThrough != previousTraceThrough)
            {
                // Trace Through a new element -> Add piecewise linear edge point
    #ifndef NDEBUG
                MeshElement edgeGen;
                if (previousTraceThrough.is_valid())
                {
                    CellHandle previousCh = hexEx.inputMesh.incident_cell(hfh.opposite_handle());
                    edgeGen = getElementBetweenInCell(hexEx.cache, previousCh, previousTraceThrough, traceThrough);
                    assert(edgeGen==previousCh || edgeGen.is_incident(hexEx.cache, previousCh));
                }
    #endif
                const auto& params = hexEx.getParametrization().get(hfh);
                const auto& intersectionParam = computeTriangleIntersection(from, d1, params[0], params[1], params[2]);
                const auto& intersectionPos = hexEx.getParametrization().inverseMapping(hexEx.inputMesh.incident_cell(hfh)).transform_point(intersectionParam);

                // Compute Progress t
                char axis = d1.axis();
                double t = std::clamp(std::abs(intersectionParam[axis] - from[axis]), 0.0, 1.0);

                ple_points.push_back(PLEPoint{.pos = intersectionPos, .t = t, .generator = traceThrough, .ch = ch});
            }
            else
            {
                // Trace through same element as before, update the ch to the latest ch
                ple_points.back().ch = ch;
            }
        }

        // Update Trace Through Element
        previousTraceThrough = traceThrough;
    }

    // Check if we failed to find any opposite
    if (!toGph.vh.is_valid()) {
        {hexEx.addMessage(HHM_FoundNoOppositePropeller());}
        return false;
    }

    // Get opposite
    const VertexHandle toVh = toGph.vh;
    const LPH toLph = toGph.lph;
    const HexVertexGenerator& toGen = hexEx.getHexVertexGenerator(toVh);

    // Compare number of blades
    const size_t fromNumBlades = fromGen.nBlades(fromGph.lph);
    const size_t toNumBlades = toGen.nBlades(toLph);
    if (fromNumBlades != toNumBlades)
    {
        hexEx.addMessage(HHM_PropellerOppositesDifferentBlades(fromNumBlades, toNumBlades));
        return false;
    }

    // Fer Index to compute offset
    int bladeIndex2 = toGen.getBladeIndex(toLph, ch, -d1, d2);

    assert(bladeIndex2 != -1);
    if (bladeIndex2 == -1) {hexEx.addMessage(HHM_FoundNoOppositeBlade()); return false;}

    // Init the Hex Edge
    auto hex_edge = HexHalfEdge{
        .startCh = startCh,
        .gph1 = fromGph,
        .gph2 = toGph,
        .offset = bladeIndex1 + bladeIndex2,
        .ple_points = std::move(ple_points)
    };

    if (Threading::getNumThreads() == 1)
    {
        // If only using a single thread, there is no danger to add the hex edge directly within the loop
        addHexEdge(hex_edge);
        return true;
    }
    else
    {
        // If using multiple threads, we add the information to a per thread vector
        thread_hex_edges[thread_id].emplace_back(hex_edge);
        return true;
    }
}

EdgeHandle HexEdgeExtractor::addHexEdge(HexHalfEdge &hex_edge)
{
    // Add Hex Edge to Mesh, set connectivity opposites

    // Get the two vertices
    const VertexHandle fromVh = hex_edge.gph1.vh;
    const VertexHandle toVh = hex_edge.gph2.vh;

    // Connect the two global propellers
    hexEx.setGlobalPropellerOpposite(hex_edge.gph1, hex_edge.gph2, hex_edge.offset);
    hexEx.setGlobalPropellerOpposite(hex_edge.gph2, hex_edge.gph1, hex_edge.offset);

    // Edge should be the newest edge (did not exist before)
    assert(!hexEx.getOutputMesh().find_halfedge(fromVh, toVh).is_valid());

    // Add the Edge to ovm
    EdgeHandle eh = hexEx.getOutputMesh().add_edge(fromVh, toVh, true);

    // Set propeller - halfedge correspondences
    auto fromHeh = eh.halfedge_handle(0);
    auto toHeh = eh.halfedge_handle(1);

    hexEx.hexHalfEdgeLocalPropellers[fromHeh] = hex_edge.gph1.lph;
    hexEx.hexHalfEdgeLocalPropellers[toHeh] = hex_edge.gph2.lph;

    hexEx.setHexHalfEdge(hex_edge.gph1, fromHeh);
    hexEx.setHexHalfEdge(hex_edge.gph2, toHeh);

    auto &piecewiseLinearMesh = hexEx.getPLMesh();
    // Piecewise Linear Edge Extraction - Add Edge positions
    if (hexEx.config.extract_piecewise_linear_edges || hexEx.config.extract_piecewise_linear_faces)
    {
        // measure longest piecewise segment
        size_t n_points = hex_edge.ple_points.size();
        if (n_points+1 > longest_ple_segment)
        {
            longest_ple_segment = hex_edge.ple_points.size()+1;
        }

        if (n_points==0)
        {
            // No intersections in between - Add straight edge from Start to End
            piecewiseLinearMesh.addPLEdgeSegment(eh, fromVh, toVh, hex_edge.startCh);
        }
        else
        {
            // Add Edge from start to 1st intesection
            VertexHandle ple_vh1 = fromVh;
            VertexHandle ple_vh2 = piecewiseLinearMesh.addPLEVertex(eh, hex_edge.ple_points[0]);
            piecewiseLinearMesh.addPLEdgeSegment(eh, ple_vh1, ple_vh2, hex_edge.startCh);

            // Add Edge between intersections
            ple_vh1 = ple_vh2;
            for (uint32_t i = 1; i < n_points; ++i)
            {
                ple_vh2 = piecewiseLinearMesh.addPLEVertex(eh, hex_edge.ple_points[i]);
                piecewiseLinearMesh.addPLEdgeSegment(eh, ple_vh1, ple_vh2, hex_edge.ple_points[i-1].ch);
                ple_vh1 = ple_vh2;
            }

            // Add Edge from last intersection to end
            ple_vh2 = toVh;

            piecewiseLinearMesh.addPLEdgeSegment(eh, ple_vh1, ple_vh2, hex_edge.ple_points[n_points-1].ch);
        }
    }

    return eh;
}

}
