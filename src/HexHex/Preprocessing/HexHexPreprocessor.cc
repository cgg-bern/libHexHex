#include <HexHex/Preprocessing/HexHexPreprocessor.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Predicates/ExactPredicates.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Preprocessing/Parametrization.hh>
#include <absl/container/flat_hash_set.h>

//#include <HexHex/Utils/TestMeshes.hh>
#if _OPENMP
#  include <omp.h>
#endif

namespace HexHex
{

HexHexPreprocessor::HexHexPreprocessor(HexExtractor& hexEx) :
    hexEx(hexEx),
    vertexHasIncidentNonIdentityTransitions(hexEx.getInputMesh().create_private_vertex_property<bool>()),
    nIncidentSingularEdges(hexEx.getInputMesh().create_private_vertex_property<int>())
{}

void HexHexPreprocessor::preprocess()
{
    ScopedStopWatch _{sw::preprocessSanitization};

    // Note to self: The order of these functions is very important even if it might not seem that way!

    hexEx.cache.build();

    scaleParameters();

    extractTransitionFunctions();

    calculateSingularities();

    if (hexEx.config.sanitize_igm) {
        sanitize();
    }

    if (hexEx.config.assert_igm_validity) {
        checkIGMValidity();
    }
}

void HexHexPreprocessor::scaleParameters()
{
    auto sf = hexEx.config.igm_scaling_factor;
    if (sf > 1) {
        hexEx.printLog("HexHex: Scaling Parametrization by s = ", sf);
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(hexEx.getInputMesh().n_cells()); ++i) {
            CellHandle ch(i);
            for (const auto& vh : hexEx.cache.cellVertices[ch]) {
                hexEx.getParametrization().tetVertexCellParameters[ch][vh] *= sf;

            }
        }
    }
}

void HexHexPreprocessor::extractTransitionFunctions()
{
    ScopedStopWatch _{sw::extractTransitionFunctions};

    // Cache the 24 Rotation transitions
    std::array<Transition, 24> ALL_24_TRANSITIONS;
    #pragma omp parallel for
    for (char i = 0; i < 24; ++i) {ALL_24_TRANSITIONS[i] = Transition(RestrictedRotation(i));}

    auto& inputMesh = hexEx.getInputMesh();
    auto& parametrization = hexEx.getParametrization();
    auto& transitions = parametrization.transitions;

    auto calculateTransitionFunction = [&](const Parameter& u0, const Parameter& v0, const Parameter& w0, const Parameter& u1, const Parameter& v1, const Parameter& w1) -> Transition
    {
        // If each vertex parameter is the same in both charts, the transition is clearly the identity
        if ((u0 == u1) && (v0 == v1) && (w0 == w1))
        {
            return Transition::IDENTITY;
        }

        // Get rotation which fits best (note that most transitions (>95%) are the identity and already detected in the check above)
        auto v0_ = v0 - u0;
        auto w0_ = w0 - u0;

        auto v1_ = v1 - u1;
        auto w1_ = w1 - u1;

        auto min_dist = std::numeric_limits<double>::max();
        auto min_transition = Transition::IDENTITY;

        for (auto i = 0u; i < 24; ++i)
        {
            auto v0_transformed = ALL_24_TRANSITIONS[i].transform_point(v0_);
            auto dist1 = (v1_ - v0_transformed) | (v1_ - v0_transformed);
            auto w0_transformed = ALL_24_TRANSITIONS[i].transform_point(w0_);
            auto dist2 = (w1_ - w0_transformed) | (w1_ - w0_transformed);

            if (dist1 + dist2 < min_dist) {
                min_dist = dist1+dist2;
                min_transition = ALL_24_TRANSITIONS[i];
            }
        }

        // Compute transition
        auto u0_t = min_transition.transform_point(u0);
        auto t = roundVector(u1 - u0_t);
        min_transition.setTranslation(t);

        // If these assertions aren't satisified, our map is nonconforming (IGM1)
        assert((u1 - min_transition.transform_point(u0)).sqrnorm() < hexEx.config.sanitization_epsilon);
        assert((v1 - min_transition.transform_point(v0)).sqrnorm() < hexEx.config.sanitization_epsilon);
        assert((w1 - min_transition.transform_point(w0)).sqrnorm() < hexEx.config.sanitization_epsilon);

        return min_transition;
    };

    auto precomputed_halfface_transitions = inputMesh.get_halfface_property<int>(hexEx.config.hf_transition_prop_name);
    if (precomputed_halfface_transitions) {
        hexEx.printLog("HexHex: Transitions are already computed :)");
    }

    uint nNonIdTransitions = 0;

    auto n_faces = static_cast<int>(hexEx.getInputMesh().n_faces());

    #pragma omp parallel for reduction(+:nNonIdTransitions)
    for (int f_idx = 0; f_idx < n_faces; ++f_idx)
    {
        auto fh = FaceHandle(f_idx);
        auto hfh0 = fh.halfface_handle(0);
        auto hfh1 = fh.halfface_handle(1);

        // For boundary faces set transition to identity
        if (inputMesh.is_boundary(fh)) {
            transitions[hfh0] = Transition::IDENTITY;
            transitions[hfh1] = Transition::IDENTITY;
            continue;
        }

        auto ch0 = inputMesh.incident_cell(hfh0);
        auto ch1 = inputMesh.incident_cell(hfh1);
        auto hehs0 = hexEx.cache.get_halfface_halfedges(hfh0);
        std::array<VertexHandle,3> vhs0 = {
            hexEx.cache.from_vertex_handle(hehs0[0]),
            hexEx.cache.from_vertex_handle(hehs0[1]),
            hexEx.cache.from_vertex_handle(hehs0[2])
        };

        Transition transition0;
        if (precomputed_halfface_transitions) {

            // Take precomputed Rotation
            transition0 = Transition(ALGOHEX_TO_HEXEX_TRANSITION_ROTATION_INDICES[(*precomputed_halfface_transitions)[hfh0]]);

            // Compute Translation
            const auto& u0 = parametrization.get(ch0, vhs0[0]);
            const auto& u1 = parametrization.get(ch1, vhs0[0]);
            const auto& t = roundVector(u1 - transition0.transform_point(u0));
            transition0.setTranslation(t);
        } else {
            const auto& u0 = parametrization.get(ch0, vhs0[0]);
            const auto& v0 = parametrization.get(ch0, vhs0[1]);
            const auto& w0 = parametrization.get(ch0, vhs0[2]);

            const auto& u1 = parametrization.get(ch1, vhs0[0]);
            const auto& v1 = parametrization.get(ch1, vhs0[1]);
            const auto& w1 = parametrization.get(ch1, vhs0[2]);

            // Compute transition from scratch
            transition0 = calculateTransitionFunction(u0,v0,w0,u1,v1,w1);
        }

        if (transition0 == Transition::IDENTITY) {
            transitions[hfh0] = transition0;
            transitions[hfh1] = transition0;
        } else {
            //std::cout << "Found Non Identity Transition" << std::endl;
            transitions[hfh0] = transition0;
            transitions[hfh1] = transition0.inverted();

            nNonIdTransitions += 1;

            // For edges and vertices, remember if entire neighborhood has identity transitions
            #pragma omp critical
            {
                for (HalfEdgeHandle heh : hehs0)
                {
                    parametrization.edgeHasIncidentNonIdentityTransitions[heh.edge_handle()] = true;
                }
                for (VertexHandle vh : vhs0)
                {
                    vertexHasIncidentNonIdentityTransitions[vh] = true;
                }
            }
        }
    }

    hexEx.printLog("HexHex: ", nNonIdTransitions, "/", inputMesh.n_faces(), " Transitions are not the Identity");
}

void HexHexPreprocessor::calculateSingularities()
{
    ScopedStopWatch _{sw::calculateSingularities};

    auto& inputMesh = hexEx.getInputMesh();
    auto& parametrization = hexEx.getParametrization();

    nIncidentSingularEdges.fill(0);

    auto parametrizationAngle = [&](HalfFaceHandle hfh1, HalfFaceHandle hfh2, HalfEdgeHandle heh) -> double
    {
        auto ch = inputMesh.incident_cell(hfh1);

        auto halfedge = inputMesh.halfedge(heh);

        auto nextHe1 = inputMesh.next_halfedge_in_halfface(heh, hfh1);
        auto nextHe2 = inputMesh.next_halfedge_in_halfface(inputMesh.opposite_halfedge_handle(heh), hfh2);

        auto vh0 = halfedge.from_vertex();
        auto vh1 = halfedge.to_vertex();
        auto vh2 = inputMesh.halfedge(nextHe1).to_vertex();
        auto vh3 = inputMesh.halfedge(nextHe2).to_vertex();

        auto u = parametrization.get(ch, vh0);
        auto v = parametrization.get(ch, vh1);
        auto w = parametrization.get(ch, vh2);
        auto t = parametrization.get(ch, vh3);

        auto d1 = v-u;
        auto d2 = w-u;
        auto d3 = d1%d2;
        d2 = d3%d1;

        auto d4 = t-u;
        auto d5 = d1%d4;
        d4 = d5%d1;

        d2.normalize();
        d4.normalize();

        auto alpha = acos(std::max(std::min(d2|d4, 1.0),-1.0));

        return alpha;
    };

    auto calculateAngleAroundEdge = [&](EdgeHandle eh) -> double
    {
        auto heh = inputMesh.halfedge_handle(eh, 0);
        auto hehf_it = inputMesh.hehf_iter(heh,2);

        double alpha = 0.0;
        for (auto i = 0u; i < inputMesh.valence(eh); ++i, ++hehf_it) {
            if (inputMesh.is_boundary(*hehf_it))
                continue;

            auto hf1 = *hehf_it;
            auto hf2 = inputMesh.adjacent_halfface_in_cell(hf1, heh);

            alpha += parametrizationAngle(hf1, hf2, heh);
        }
        return alpha;
    };

    auto precomputed_edge_valence_offsets = inputMesh.get_edge_property<int>(hexEx.config.e_valence_prop_name);
    if (precomputed_edge_valence_offsets) {
        hexEx.printLog("HexHex: Edge Valences are already computed :)");
    }

    uint n_singular_edges = 0;

    #pragma omp parallel for reduction(+:n_singular_edges)
    for (int i = 0; i < static_cast<int>(inputMesh.n_edges()); ++i)
    {
        auto eh = EdgeHandle(i);
        bool isBoundary = inputMesh.is_boundary(eh);

        // Compute valence based on dihedral face angle around edge
        // or just copy it if it is already available
        if (precomputed_edge_valence_offsets)
        {
            // If we have precomputed valences, it is singular if the value (offset from regular) is nonzero
            parametrization.edgeSingularity[eh] = ((*precomputed_edge_valence_offsets)[eh] != 0);
        }
        else
        {
            // Otherwise, compute the angle around the edge
            double alpha = calculateAngleAroundEdge(eh);
            int valence = int(round(alpha / (0.5*M_PI)));
            if (isBoundary)
            {
                parametrization.edgeSingularity[eh] = (valence != 2);
            }
            else
            {
                parametrization.edgeSingularity[eh] = (valence != 4);
            }
        }
        // Edge is singular if angle is not 2pi around inner edge or not pi around boundary edge
        if (parametrization.edgeSingularity[eh])
        {
            n_singular_edges++;
            auto e = inputMesh.edge(eh);
            #pragma omp atomic
            nIncidentSingularEdges[e.from_vertex()] += 1;
            #pragma omp atomic
            nIncidentSingularEdges[e.to_vertex()] += 1;
        }
    }
    hexEx.printLog("HexHex: ", n_singular_edges, "/", inputMesh.n_edges(), " Edges are Singular");
}

void HexHexPreprocessor::sanitize()
{
    hexEx.printLog("HexHex: Sanitize IGM");
    auto& inputMesh = hexEx.getInputMesh();
    auto& parametrization = hexEx.getParametrization();

    auto setTetVertexParameter = [&](const CellHandle& ch, const VertexHandle& vh, const Parameter& param) -> void
    {
        // Check if it is big fix
        auto n = (parametrization.get(ch,vh)-param).sqrnorm();
        if (n > std::pow(hexEx.config.sanitization_epsilon,2)) {
            hexEx.printErr("HexHex: Update igm(", ch, ", ", vh, ") from ", parametrization.get(ch,vh), " to ", param);
            hexEx.addMessage(HHM_LargeSanitizationFix());
        }

        // Update
        for (char i = 0; i < 3; ++i) {
            #pragma omp atomic write
            parametrization.tetVertexCellParameters[ch][vh][i] = param[i];
        }
    };

    auto propagateVertexParameter = [&](const Parameter& param, const VertexHandle& vh, const CellHandle& startCell) -> void
    {
        assert(std::find(hexEx.cache.vertexCells[vh].begin(), hexEx.cache.vertexCells[vh].end(), startCell) != hexEx.cache.vertexCells[vh].end());

        // If we do not have any nonid transitions, just set the parameter for all incident tets
        if (!vertexHasIncidentNonIdentityTransitions[vh]) {
            for (CellHandle ch : hexEx.cache.vertexCells[vh]) {
                setTetVertexParameter(ch, vh, param);
            }
            return;
        }

        // Update the parameter in the start cell
        setTetVertexParameter(startCell, vh, param);

        // Else, move from cell to cell and update the parameter according to the transition

        std::function<void(Parameter,VertexHandle,CellHandle,absl::flat_hash_set<CellHandle>&)> propagateVertexParameterRecursive;
        propagateVertexParameterRecursive = [&](Parameter param, VertexHandle vh, CellHandle ch, absl::flat_hash_set<CellHandle>& toBeProcessed) -> void
        {
            toBeProcessed.erase(ch);

            setTetVertexParameter(ch, vh, param);

            if (toBeProcessed.empty())
                return;

            auto c = inputMesh.cell(ch);
            auto& halffaces = c.halffaces();

            for (const auto& hfh : halffaces)
            {
                auto oppHfh = inputMesh.opposite_halfface_handle(hfh);
                auto oppCh = inputMesh.incident_cell(oppHfh);
                if (!oppCh.is_valid())
                    continue;

                if (toBeProcessed.contains(oppCh)) {

                    auto newParam = parametrization.transition(hfh).transform_point(param);
                    propagateVertexParameterRecursive(newParam, vh, oppCh, toBeProcessed);

                    if (toBeProcessed.empty())
                        return;
                }
            }
        };

        // Vertex has different parameters in different charts
        auto toBeProcessed = absl::flat_hash_set<CellHandle>();
        toBeProcessed.reserve(hexEx.cache.vertexCells[vh].size());
        for (CellHandle ch : hexEx.cache.vertexCells[vh])
        {
            toBeProcessed.insert(ch);
        }
        propagateVertexParameterRecursive(param, vh, startCell, toBeProcessed);
    };

    auto getParameterProjectedToHalfEdge = [&](Parameter param, CellHandle ch, HalfEdgeHandle heh) -> Parameter
    {
        auto he = inputMesh.halfedge(heh);

        const auto& u = parametrization.get(ch, he.from_vertex());
        const auto& v = parametrization.get(ch, he.to_vertex());

        auto maxCoord = argAbsMax(v - u);

        for (auto i = 0u; i < 3u; ++i) {
            if (i != maxCoord) {
                param[i] = round(param[i]);
            }
        }

        return param;
    };

    auto getIncidentSingularEdge = [&](VertexHandle vh) -> HalfEdgeHandle
    {
        for (auto voh_it = inputMesh.voh_iter(vh); voh_it.valid(); ++voh_it) {
            if (parametrization.edgeSingularity[(*voh_it).edge_handle()]) {
                return *voh_it;
            }
        }
        return HalfEdgeHandle();
    };

    //-----------------------------
    // Truncate Precision
    //-----------------------------
    {
    ScopedStopWatch _{sw::sanitizeParametrization};

    uint num_singular_vertices = 0;

#pragma omp parallel for reduction(+:num_singular_vertices)
    for (int vh_idx = 0; vh_idx < static_cast<int>(inputMesh.n_vertices()); ++vh_idx) {

        auto vh = VertexHandle(vh_idx);
        if (hexEx.cache.vertexCells[vh].empty()) {continue;}

        // Compute delta for truncation
        double maxCoord = 0;
        for (CellHandle ch : hexEx.cache.vertexCells[vh])
        {
            const auto& p = parametrization.get(ch, vh);
            for (unsigned char i = 0u; i < 3u; ++i) {
                maxCoord = std::max(std::fabs(p[i]), maxCoord);
            }
        }
        double delta = std::pow(2, std::ceil(std::log2(maxCoord)));

        // Get incident cell (should also be incident to an incident singular edge if one exists)
        auto ch = CellHandle(-1);
        auto hehSingular = HalfEdgeHandle(-1);
        if (nIncidentSingularEdges[vh] > 0) {
            hehSingular = getIncidentSingularEdge(vh);
            ch = *inputMesh.hec_iter(hehSingular);
        } else {
            ch = hexEx.cache.vertexCells[vh][0];
        }
        if (!ch.is_valid()) {hexEx.addMessage(HHM_IsolatedTetVertex()); continue;}
        assert(ch.is_valid());


        // Precision Truncation
        auto param = parametrization.get(ch, vh);
        for (unsigned int i = 0u; i < 3u; ++i)
        {
            int sign = std::signbit(param[i]) ? -1 : 1;
            double tmp = param[i]+sign*delta;
            parametrization.tetVertexCellParameters[ch][vh][i] = tmp - sign*delta;
        }

        // Fix Singularity Point (If we disable boundary alignment, we skip this step if the incident singular edge is boundary)
        param = parametrization.get(ch, vh);

        if (hehSingular.is_valid())
        {
            if (nIncidentSingularEdges[vh] == 2)
            {
                // Vertex is between two singular edges -> Project onto arc
                param = getParameterProjectedToHalfEdge(param, ch, hehSingular);
            } else if (nIncidentSingularEdges[vh] > 0)
            {
                // Vertex is singular -> integer align the little rascal
                param = roundVector(param);
                num_singular_vertices++;
            }
        }

        setTetVertexParameter(ch, vh, param);
        propagateVertexParameter(parametrization.get(ch, vh), vh, ch);
    }

    hexEx.printLog("HexHex: ", num_singular_vertices, "/", inputMesh.n_vertices(), " Vertices are Singular.");

    }

    {
    //------------------------------
    // Project Boundary Faces
    //------------------------------
    ScopedStopWatch _{sw::snapBoundaryFaces};
    auto getParameterNormal = [&](const HalfFaceHandle& hfh) -> Parameter
    {
        auto vhs = inputMesh.get_halfface_vertices(hfh);
        auto ch = inputMesh.incident_cell(hfh);
        if (!ch.is_valid())
            ch = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(hfh));
        const auto& u = parametrization.get(ch, vhs[0]);
        const auto& v = parametrization.get(ch, vhs[1]);
        const auto& w = parametrization.get(ch, vhs[2]);
        return ((v-u)%(w-u)).normalized();
    };

    #pragma omp parallel for
    for (int hfh_idx = 0; hfh_idx < static_cast<int>(inputMesh.n_halffaces()); ++hfh_idx)
    {
        auto hfh = HalfFaceHandle(hfh_idx);
        if (inputMesh.is_boundary(hfh))
        {
            hfh = inputMesh.opposite_halfface_handle(hfh);
            CellHandle ch = inputMesh.incident_cell(hfh);
            if (!ch.is_valid()) {continue;}

            // find projection direction ( = dominant normal direction)
            int projCoord = argAbsMax(getParameterNormal(hfh));

            // find projected parameter (rounded average of all three vertices)
            double projParam = 0.0;
            for (auto hfv_it = inputMesh.hfv_iter(hfh); hfv_it.valid(); ++hfv_it) {
                projParam += parametrization.get(ch, *hfv_it)[projCoord];
            }
            projParam = round(projParam/3.0);

            // set projected parameter for all incident cells
            //bool update = false;
            //auto maxDisplacement = 1e-6;
            for (auto hfv_it = inputMesh.hfv_iter(hfh); hfv_it.valid(); ++hfv_it) {
                auto diff = projParam - parametrization.get(ch, *hfv_it)[projCoord];
                if (diff != 0.) {
                    auto newParam = parametrization.get(ch, *hfv_it);
                    newParam[projCoord] = projParam;
                    propagateVertexParameter(newParam, *hfv_it, ch);
                }
            }
        }
    }
    }
}

void HexHexPreprocessor::checkIGMValidity()
{
    ScopedStopWatch _{sw::checkIGMValidity};

    hexEx.printLog("HexHex: Checking IGM Validity...");

    auto& parametrization = hexEx.getParametrization();

    // Each singular edge must be integer aligned
    uint n_unaligned_s_edges = 0;
#pragma omp parallel for reduction(+:n_unaligned_s_edges)
    for (int i = 0; i < static_cast<int>(hexEx.getInputMesh().n_edges()); ++i)
    {
        if (parametrization.edgeSingularity[EdgeHandle(i)])
        {
            if (!hexEx.edgeIsOnIntegerLattice(EdgeHandle(i)))
            {
                n_unaligned_s_edges++;
            }
        }
    }
    if (n_unaligned_s_edges > 0)
    {
        std::cerr << "HexHex: Found " << n_unaligned_s_edges << " unaligned singular edges" << std::endl;
        assert(false);
    }

    // Each boundary face must be integer aligned if snap boundary is enabled
    uint n_unaligned_b_faces = 0;
    #pragma omp parallel for reduction(+:n_unaligned_b_faces)
    for (int i = 0; i < static_cast<int>(hexEx.getInputMesh().n_faces()); ++i)
    {
        if (hexEx.getInputMesh().is_boundary(FaceHandle(i)))
        {
            if (!hexEx.faceIsOnIntegerLattice(FaceHandle(i)))
            {
                n_unaligned_b_faces++;
            }
        }
    }
    if (n_unaligned_b_faces > 0)
    {
        std::cerr << "HexHex: Found " << n_unaligned_b_faces << " unaligned boundary faces" << std::endl;
        assert(false);
    }


    // For every halfface we must have halfface_transition(halfface_param) = opposite_halfface_param
    uint n_nonmatching_transitions = 0;
    auto& inputMesh = hexEx.getInputMesh();
#pragma omp parallel for reduction(+:n_nonmatching_transitions)
    for (int i = 0; i < static_cast<int>(hexEx.getInputMesh().n_halffaces()); ++i)
    {
        HalfFaceHandle hfh1(i);
        if (hexEx.getInputMesh().is_boundary(hfh1.face_handle())) continue;
        CellHandle ch1 = inputMesh.incident_cell(hfh1);
        HalfFaceHandle hfh2 = hfh1.opposite_handle();
        CellHandle ch2 = inputMesh.incident_cell(hfh2);

        auto vhs = hexEx.getInputMesh().get_halfface_vertices(hfh1);

        std::vector<Parameter> params1 = {parametrization.get(ch1, vhs[0]),parametrization.get(ch1, vhs[1]),parametrization.get(ch1, vhs[2])};
        std::vector<Parameter> params2 = {parametrization.get(ch2, vhs[0]),parametrization.get(ch2, vhs[1]),parametrization.get(ch2, vhs[2])};

        auto trans = parametrization.transitions[hfh1];

        for (int j = 0; j < 3; ++j)
        {
            Parameter u1 = params1[j];
            Parameter u2 = params2[j];
            Parameter uT = trans.transform_point(u1);
            bool match = (uT == u2);
            if (!match)
            {
                //std::cerr << std::hexfloat;
                std::cerr << "Inconsistent Parameters!" << std::endl;
                std::cerr << "vh = " << vhs[j] << ", hf = " << hfh1 << ", ch = " << ch1 << "opp hf = " << hfh2 << ", opp ch = " << ch2 << std::endl;
                std::cerr << "u1 (hf) = " << u1 << std::endl;
                std::cerr << "u2 (opp) = " << u2 << std::endl;
                std::cerr << "trans(u1) = " << uT << std::endl;
                std::cerr << "trans = " << trans.toMatrix() << std::endl;
                std::cerr << "diff = " << (uT - u2).sqrnorm() << std::endl;
                std::cerr << "boundary vertex: " << hexEx.getInputMesh().is_boundary(vhs[j]) << std::endl;
                std::cerr << "n singular edges: " << nIncidentSingularEdges[vhs[j]] << std::endl;
                n_nonmatching_transitions++;
            }
        }
    }
    if (n_nonmatching_transitions > 0)
    {
        std::cerr << "HexHex: Found " << n_nonmatching_transitions << " inconsistent transitions" << std::endl;
        assert(false);
    }

    // Each tet must be ORI_ABOVE oriented
    uint n_flipped = 0;
    uint n_degenerate = 0;
    uint n_proper = 0;
    #pragma omp parallel for reduction(+:n_flipped,n_degenerate,n_proper)
    for (int i = 0; i < static_cast<int>(hexEx.getInputMesh().n_cells()); ++i)
    {
        CellHandle ch = CellHandle(i);

        // Check Tet Volume
        auto ori = ori3d(parametrization.get(ch));
        if (ori == ORI_ABOVE) n_proper++;
        if (ori == ORI_BELOW) n_flipped++;
        if (ori == ORI_ZERO) n_degenerate++;
    }

    if (n_degenerate + n_flipped > 0)
    {
    std::cerr << "HexHex: Detected " << n_proper << " Proper, " << n_flipped << " Flipped and " << n_degenerate << " Degenerate Tets after the Sanitization!" << std::endl;
        hexEx.addMessage(HHM_FlippedOrDegenerateTet());
    }

    if (n_degenerate == 0 && n_flipped == 0 && n_proper == hexEx.getInputMesh().n_cells()
        && n_unaligned_s_edges == 0 && n_unaligned_b_faces == 0 && n_nonmatching_transitions == 0
        )
    {
        hexEx.printLog("HexHex: IGM is valid, Let's GO!");
    }
}

}
