/*
 * Copyright 2023 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */

#include <HexHex/Preprocessing/HexHexPreprocessor.hh>
#include <HexHex/Predicates/ExactPredicates.hh>
#include <HexHex/Utils/Threading.hh>

#include <libTimekeeper/json.hh>
#include <nlohmann/json.hpp>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Config_json.hh>

#include <HexHex/Utils/Utils.hh>
#include <HexHex/Predicates/GeneratorPredicates.hh>
#include <HexHex/Predicates/ExactPredicates.hh>

#include <HexHex/Preprocessing/Parametrization.hh>

#include <HexHex/LocalTopology/HexVertexGenerator.hh>

#include <HexHex/VertexExtraction/HexVertexExtractor.hh>
#include <HexHex/GlobalTopology/HexEdgeExtractor.hh>
#include <HexHex/GlobalTopology/HexCellExtractor.hh>
#include <HexHex/PiecewiseLinear/PLHexFaceExtractor.hh>
#include <HexHex/PiecewiseLinear/PLMesh.hh>

#include <HexHex/Utils/FileAccessor.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/PropertyNames.hh>
#include <HexHex/Commons/Permutation.hh>

#include <HexHex/Utils/MemoryUsage/MemoryUsage.hh>


using namespace OpenVolumeMesh;

namespace HexHex {

HexExtractor::~HexExtractor() = default;

HexExtractor::HexExtractor(const TetrahedralMesh& tetmesh)
    :
    inputMesh(tetmesh),
    cache(inputMesh),
    pwlMesh_{std::make_unique<PLMesh>(cache, getOutputMesh())},
    parametrization(cache, getOutputMesh())
{
    exactinit(); // <- IMPORTANT!!!!!!!!!!!!!!!!!

    config = Config();

    // Output Mesh - only vertex bottom up for faster tracing (keep cell extraction fast without edge bottom up!)
    outputMesh_->enable_bottom_up_incidences(false);
    //outputMesh_->enable_face_bottom_up_incidences(true);
#ifndef NDEBUG
    outputMesh_->enable_vertex_bottom_up_incidences(true);
#endif
}
PolyhedralMesh& HexExtractor::getPiecewiseLinearMesh() {return pwlMesh_->mesh;}

bool HexExtractor::extract(Config config)
{
    MemoryStats::peak_reset();

    this->config = config;
    Threading::setNumThreads(config.num_threads);

#ifndef NDEBUG
    printLog("HexHex: Running in Debug. Assertions are enabled.");
#else
    printLog("HexHex: Running in Release. Assertions are disabled.");
#endif
    hasError = false;
    
    ScopedStopWatch _{sw::extract};
    std::chrono::steady_clock::time_point t1, t2;

    t1 = std::chrono::steady_clock::now();

    printLog("HexHex: Input Tet Mesh has ", inputMesh.n_vertices(), " Vertices, ", inputMesh.n_edges(), " Edges, ", inputMesh.n_faces(), " Faces and ", inputMesh.n_cells(), " Cells");
#if _OPENMP
    printLog("HexHex: Harvest the Power of Parallelization with ", omp_get_max_threads(), "/", omp_get_num_procs(), " available Threads");
#else
    printLog("HexHex: Running single-threaded (built without OpenMP)");
#endif

    // Input Mesh - Ensure bottom up incidences
    if (!inputMesh.has_full_bottom_up_incidences()) {
        std::cerr << "Input Error: The input tet mesh does not have full bottom up incidences!" << std::endl;
        return false;
    }

    preprocessParametrization();
    if (hasError) return false;

    extractHexVertices();
    if (hasError || outputMesh_->n_vertices()==0) return false;

    extractLocalTopology();
    if (hasError) return false;

    extractHexEdges();
    if (hasError) return false;

    extractHexCells();
    if (hasError) return false;

    transferFeatureTags();

    extractPiecewiseLinearHexFaces();

    printLog("HexHex: Output Hex Mesh has ", outputMesh_->n_vertices(), " Vertices, ", outputMesh_->n_edges(), " Edges, ", outputMesh_->n_faces(), " Faces and ", outputMesh_->n_cells(), " Cells");

    t2 = std::chrono::steady_clock::now();
    printLog("HexHex: Finished in ", std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count(), " ms");
    peak_memory_usage_bytes = MemoryStats::peak_get();
    if (peak_memory_usage_bytes) {
        printLog("HexHex: Peak memory was approximately ", peak_memory_usage_bytes, " B (= ", peak_memory_usage_bytes/1e6, " MB)");
    } else {
        printLog("HexHex: Peak memory unknown (compiled without jemalloc)");
    }

    // See if the number of cells matches what we expect from the local topology extraction
    if (inputMesh.n_cells() > 0)
    {
        return numGlobalCorners() == 8*outputMesh_->n_cells();
    }
    return true;
}

void HexExtractor::preprocessParametrization()
{
    printLog("HexHex: Preprocess Parametrization");
    HexHexPreprocessor(*this).preprocess();
}

//=======================
// Predicates
//=======================

bool HexExtractor::computeVertexGeneratorElement(CellHandle ch, const Parameter& param, MeshElement& elem)
{
    // Get the 4 vertices of the tet and their parameters in the chart of the tet
    const auto& vhs = cache.cellVertices[ch];// parametrization.getCellVertices(ch);
    auto params = parametrization.get(ch);

    MeshElement res = ::HexHex::computeVertexGeneratorElement(
        cache,
        ch,
        vhs, params,
        param
    );
    if (res.is_valid()) {elem.set(res); return true;}
    return false;
}

bool HexExtractor::computePropellerHolderElement(const CellHandle& ch, const MeshElement& generator, const Parameter& from, const Direction& dir, MeshElement& holder)
{
    // Get the 4 vertices of the tet and their parameters in the chart of the tet
    const auto& vhs = cache.cellVertices[ch];// parametrization.getCellVertices(ch);
    auto params = parametrization.get(ch);

    MeshElement res = ::HexHex::computePropellerHolderElement(
        cache,
        ch,
        vhs, params, generator,
        from, dir
        );
    if (res.is_valid()) {holder.set(res); return true;}
    return false;
}

bool HexExtractor::computePropellerBladeCasingElement(const CellHandle& ch, const MeshElement& holder, const Parameter &from, const Direction& dir1, const Direction& dir2, MeshElement& casing)
{
    // Get the 4 vertices of the tet and their parameters in the chart of the tet
    const auto& vhs = cache.cellVertices[ch];// parametrization.getCellVertices(ch);
    auto params = parametrization.get(ch);

    MeshElement res = ::HexHex::computePropellerBladeCasingElement(
        cache,
        ch,
        vhs, params, holder,
        from, dir1, dir2
        );
    if (res.is_valid()) {casing.set(res); return true;}
    return false;
}

//========================================
// Geometry Extraction
//========================================

bool HexExtractor::extractHexVertices()
{
    printLog("HexHex: Extract Hex Vertices");

    // Clear the List of Generators
    generators.clear();

    // Clear the Output Mesh completely and request the hex-mesh properties
    outputMesh_->clear();

    hexVertexGenerators.fill({});
    hexVertexTetGeneratorTypes.fill(-1);
    hexVertexTetGeneratorIndices.fill(-1);
    hexHalfEdgeLocalPropellers.fill({});

    return VertexExtractor(*this).extract();
}

void HexExtractor::addHexVertices(const std::vector<Parameter> &params, MeshElement genElem, CellHandle refCh)
{
    auto n = params.size();
    if (n == 0) return;

    auto invMap = parametrization.inverseMapping(refCh);

    // Add Vertices
    VertexHandle vhMin, vhMax;

    vhMin = VertexHandle(outputMesh_->n_vertices());
    auto vh = VertexHandle(vhMin.idx());
    for (auto i = 0u; i < n; ++i) {
        auto param = params[i];
        auto pos = invMap.transform_point(param);
        vh = VertexHandle(vhMin.idx() + i);

        outputMesh_->add_vertex(pos); // add to OVM

        // Add to piecewise linear mesh
        if (config.extract_piecewise_linear_edges || config.extract_piecewise_linear_faces)
        {
            getPLMesh().addVertex(pos, genElem);
        }
    }
    vhMax = vh;

    // Set Hex Vertex Parameters in Reference Cell
    parametrization.setHexVertexParameters(refCh, {vhMin, vhMax}, params);

    // Init Generator
    if (genElem.type() == ELEMENT_VERTEX_TYPE) {
        assert(n == 1);
        // Set Hex Vertex Parameter in other cells around vertex
        for (CellHandle ch : cache.vertexCells[genElem.vh()])
            parametrization.setHexVertexParameters(ch, {vhMin, vhMin}, {parametrization.get(ch, genElem.vh())});

        generators.push_back(HexVertexGenerator(*this, genElem.vh(), vhMin));

    } else if (genElem.type() == ELEMENT_EDGE_TYPE) {
        // Set Hex Vertex Parameters in other cells around edge
        auto ch1st = inputMesh.incident_cell(*inputMesh.hehf_iter(genElem.heh()));
        auto accumulatedTransition = parametrization.getTransition(refCh, ch1st, genElem.eh());

        for (auto hehf_iter = inputMesh.hehf_iter(genElem.heh()); hehf_iter.is_valid(); ++hehf_iter)
        {
            auto ch2 = inputMesh.incident_cell(*hehf_iter);
            if (!ch2.is_valid()) break;

            parametrization.setHexVertexParameters(ch2, {vhMin, vhMax}, accumulatedTransition.transform_points(params));

            auto transition = parametrization.transition((*(hehf_iter+1)).opposite_handle());
            accumulatedTransition = transition * accumulatedTransition;
        }

        // Initialize Edge Generator
        generators.push_back(HexVertexGenerator(*this, genElem.eh(), {vhMin, vhMax}));

    } else if (genElem.type() == ELEMENT_FACE_TYPE) {
        // Set Hex Vertex Parameters in other incident cell to face
        auto hfh = genElem.hfh();
        assert(inputMesh.incident_cell(hfh) == refCh);
        auto ch2 = inputMesh.incident_cell(hfh.opposite_handle());
        if (ch2.is_valid()) parametrization.setHexVertexParameters(ch2, {vhMin, vhMax}, parametrization.transition(hfh).transform_points(params));

        // Initialize Face Generator
        generators.push_back(HexVertexGenerator(*this, hfh, {vhMin, vhMax}));

    } else {
        assert(genElem.type() == ELEMENT_CELL_TYPE);
        generators.push_back(HexVertexGenerator(*this, genElem.ch(), {vhMin, vhMax}));
    }

    // Remember Generator per Vex-Vertex
    for (auto i = vhMin.idx(); i <= vhMax.idx(); ++i)
    {
        VertexHandle vh(i);
        hexVertexGenerators[vh] = generators.size()-1;
        hexVertexTetGeneratorTypes[vh] = genElem.type();

        // Store the index of the whole (not half) entity
        hexVertexTetGeneratorIndices[vh] = genElem.idx();
    }
}

//===================================================
// Topology Extraction
//===================================================

void HexExtractor::extractLocalTopology()
{
    printLog("HexHex: Extract Local Topology");

    ScopedStopWatch _{sw::extractLocalTopology};

    size_t n_edges = 0;
    size_t n_faces = 0;
    size_t n_cells = 0;

    #pragma omp parallel for reduction(+:n_edges,n_faces,n_cells)
    for (int i = 0; i < static_cast<int>(generators.size()); ++i)
    {
        auto n = generators[i].extractLocalTopology();

        n_edges += n[0]; // num propellers (halfedges)
        n_faces += n[1]; // num blades (quarterfaces)
        n_cells += n[2]; // num corners (eights of cells)
    }

    // Set Offsets
    for (int i = 1; i < static_cast<int>(generators.size()); ++i)
    {
        generators[i].global_offsets.hex_vertices_offset =
            generators[i-1].global_offsets.hex_vertices_offset + generators[i-1].nVertices();
        generators[i].global_offsets.global_propellers_offset =
            generators[i-1].global_offsets.global_propellers_offset + generators[i-1].nGlobalPropellers();
        generators[i].global_offsets.global_corners_offset =
            generators[i-1].global_offsets.global_corners_offset + generators[i-1].nGlobalCorners();
    }

    // Resize Global Propeller Vectors
    globalPropellerOpposites.resize(numGlobalPropellers());
    globalPropellerHexHalfEdges.resize(numGlobalPropellers());

    // Get expected number of edges, faces, cells
    n_edges /= 2;
    n_faces /= 4;
    n_cells /= 8;

    outputMesh_->reserve_edges(n_edges);
    outputMesh_->reserve_faces(n_faces);
    outputMesh_->reserve_cells(n_cells);

    printLog("HexHex: Detected ", n_edges, " Edges");
    printLog("HexHex: Detected ", n_faces, " Faces");
    printLog("HexHex: Detected ", n_cells, " Cells");
}

void HexExtractor::extractHexEdges()
{
    printLog("HexHex: Trace ", numGlobalPropellers(), " Propellers");

    // If we extract pl faces we also want pl edges
    if (config.extract_piecewise_linear_faces)
        config.extract_piecewise_linear_edges = true;

    HexEdgeExtractor(*this).tracePropellers();

    if (config.extract_piecewise_linear_edges)
    {
        auto& mesh = getPiecewiseLinearMesh();
        printLog("HexHex: Piecewise Linear Edges Mesh has ", mesh.n_vertices(), " Vertices and ", mesh.n_edges(), " Edges");
    }
}

void HexExtractor::extractHexCells()
{
    printLog("HexHex: Construct Hex Mesh");
    HexCellExtractor(*this).extractHexCells();
}

void HexExtractor::extractPiecewiseLinearHexFaces()
{
    ScopedStopWatch _{sw::extractPiecewiseLinearFaces};
    if (!config.extract_piecewise_linear_faces) return;
    printLog("HexHex: Extract Piecewise Linear Faces");

    std::vector<PLHexFaceExtractor::PolyMesh> face_meshes(outputMesh_->n_faces());

    // Compute the piecewise linear patches
#pragma omp parallel for default(none) schedule(guided) shared(face_meshes)
    for (int i = 0; i < static_cast<int>(outputMesh_->n_faces()); ++i)
    {
        auto pfe = PLHexFaceExtractor(*this);
        pfe.extract(FaceHandle(i));

        face_meshes[i] = pfe.getHexFaceMesh();
    }

    // Append the per Face meshes to the global pl mesh
    for (int f_idx = 0; f_idx < static_cast<int>(face_meshes.size()); ++f_idx)
    {
        const auto& face_mesh = face_meshes[f_idx];

        // Add Vertices
        std::vector<VertexHandle> face_mesh_vhs;
        face_mesh_vhs.reserve(face_mesh.points.size());
        for (uint32 i = 0u; i < face_mesh.points.size(); ++i)
        {
            VertexHandle vh = getPLMesh().addPLFVertex(face_mesh.points[i]);
            face_mesh_vhs.push_back(vh);
        }

        // Add Polygons
        for (uint32 i = 0; i < face_mesh.polygons.size(); ++i)
        {
            std::vector<VertexHandle> poly_vhs;
            poly_vhs.reserve(face_mesh.polygons[i].size());
            for (auto j : face_mesh.polygons[i])
            {
                poly_vhs.push_back(face_mesh_vhs[j]);
            };

            getPLMesh().addPLFacePatch(FaceHandle(f_idx), poly_vhs, face_mesh.tet_per_polygon[i]);
        }
    }

    printLog("HexHex: Piecewise Linear Faces Mesh has ",
             getPLMesh().mesh.n_vertices(), " Vertices, ",
             getPLMesh().mesh.n_edges(), " Edges and ",
             getPLMesh().mesh.n_faces(), " Faces");
}

//=============================
// Saving and Loading
//=============================

void HexExtractor::transferFeatureTags()
{
    printLog("HexHex: Transfer AlgoHex feature tags if they exist");

    ScopedStopWatch _{sw::transferFeatureTags};

    // Vertex Properties
    if (inputMesh.vertex_property_exists<int>(PropertyNames::AlgoHexFeatureVertices)) {
        auto input_vfeature = inputMesh.get_vertex_property<int>(PropertyNames::AlgoHexFeatureVertices).value();
        auto output_vfeature = outputMesh_->request_vertex_property<int>(PropertyNames::AlgoHexFeatureVertices);
        outputMesh_->set_persistent(output_vfeature,true);

        #pragma omp parallel for num_threads(omp_get_max_threads())
        for (int i = 0; i < static_cast<int>(outputMesh_->n_vertices()); ++i) {
            // Generator must be feature vertex
            auto hexVh = VertexHandle(i);
            const auto& tetVh = getHexVertexGenerator(hexVh).generator.vh();
            output_vfeature[hexVh] = (tetVh.is_valid())? input_vfeature[tetVh] : 0;
        }

    }

    // Edge Properties
    if (inputMesh.edge_property_exists<int>(PropertyNames::AlgoHexFeatureEdges)) {
        auto input_efeature = inputMesh.get_edge_property<int>(PropertyNames::AlgoHexFeatureEdges).value();
        auto output_efeature = outputMesh_->request_edge_property<int>(PropertyNames::AlgoHexFeatureEdges);
        outputMesh_->set_persistent(output_efeature,true);

        #pragma omp parallel for num_threads(omp_get_max_threads())
        for (int i = 0; i < static_cast<int>(outputMesh_->n_edges()); ++i) {
            auto hexEh = EdgeHandle(i);
            auto hexHeh = hexEh.halfedge_handle(0);

            // The holder of one of the propellers must be a feature tet edge
            const auto& gen0 = getHexVertexGenerator(outputMesh_->from_vertex_handle(hexHeh));
            auto tetEh = gen0.getHolderEdge(hexHeh);

            output_efeature[hexEh] = (tetEh.is_valid())? input_efeature[tetEh] : 0;
        }
    }

    // Face Properties
    if (inputMesh.face_property_exists<int>(PropertyNames::AlgoHexFeatureFaces)) {
        auto input_ffeature = inputMesh.get_face_property<int>(PropertyNames::AlgoHexFeatureFaces).value();
        auto output_ffeature = outputMesh_->request_face_property<int>(PropertyNames::AlgoHexFeatureFaces);
        outputMesh_->set_persistent(output_ffeature,true);

        #pragma omp parallel for num_threads(omp_get_max_threads())
        for (int i = 0; i < static_cast<int>(outputMesh_->n_faces()); ++i) {
            const auto& hexFh = FaceHandle(i);
            // The casing of one of the propeller face corners must be a feature tet face
            const auto& hexHehs = outputMesh_->face(hexFh).halfedges();
            const auto& heh1 = hexHehs[0];
            const auto& heh2 = hexHehs[3].opposite_handle();
            const auto& hexVh = outputMesh_->from_vertex_handle(heh1);
            assert(hexVh == outputMesh_->from_vertex_handle(heh2));
            const auto& gen = getHexVertexGenerator(hexVh);
            const auto& tetFh = gen.getCasingFace(heh1, heh2);
            output_ffeature[hexFh] = (tetFh.is_valid())? input_ffeature[tetFh] : 0;
        }
    }
}




bool HexExtractor::edgeIsOnIntegerLattice(EdgeHandle eh) const {
    auto ch = *inputMesh.ec_iter(eh);
    auto e = inputMesh.edge(eh);
    const auto& u = parametrization.get(ch, e.from_vertex());
    const auto& v = parametrization.get(ch, e.to_vertex());
    for (auto i = 0u; i < 3; ++i) {
        if (u[i] == v[i] && u[(i+1)%3] == v[(i+1)%3] && u[i] == round(u[i]) && u[(i+1)%3] == round(u[(i+1)%3])) {
            return true;
        }
    }
    return false;
}

bool HexExtractor::faceIsOnIntegerLattice(FaceHandle fh) const {
    auto hfh = fh.halfface_handle(0);
    auto ch = inputMesh.incident_cell(hfh);
    if (!ch.is_valid()) {
        hfh = hfh.opposite_handle();
        ch = inputMesh.incident_cell(hfh);
    }
    auto vhs = inputMesh.get_cell_vertices(hfh);
    const auto& u = parametrization.get(ch, vhs[0]);
    const auto& v = parametrization.get(ch, vhs[1]);
    const auto& w = parametrization.get(ch, vhs[2]);
    for (auto i = 0u; i < 3; ++i) {
        if (u[i] == v[i] && v[i] == w[i] && u[i] == round(u[i])) {
            return true;
        }
    }
    return false;
}

CellHandle HexExtractor::owner(const MeshElement& me) const
{
    if (me.is_cell()) return me.ch();
    if (me.is_face()) return CellHandle(cache.f_owners[me.fh()]);
    if (me.is_edge()) return CellHandle(cache.e_owners[me.eh()]);
    if (me.is_vertex()) return CellHandle(cache.v_owners[me.vh()]);
    return CellHandle(-1);
}
bool HexExtractor::isBoundary(const MeshElement& me) const
{
    if (me.is_cell()) return inputMesh.is_boundary(me.ch());
    if (me.is_face()) return inputMesh.is_boundary(me.fh());
    if (me.is_edge()) return inputMesh.is_boundary(me.eh());
    if (me.is_vertex()) return inputMesh.is_boundary(me.vh());
    return false;
}

GlobalPropellerHandle HexExtractor::getOppositeBlade(const VertexHandle vh, const LPH lph1, const LPH lph2) const
{
    const auto& gen = getHexVertexGenerator(vh);
    const auto& conn = getGlobalPropellerOpposite({.vh = vh, .lph = lph1});
    const int32 i = gen.findBladeIndex(lph1, lph2);
    const int32 o = conn.offset;
    const auto m = gen.nBlades(lph1);

    assert(i >= 0 && i < m);
    assert(isGlobalPropellerConnected({vh, lph1}));

    const LPH ph1Opp = conn.gph.lph;
    const VertexHandle vhOpp = conn.gph.vh;
    const auto& genOpp = getHexVertexGenerator(vhOpp);

    assert(m == genOpp.nBlades(ph1Opp));

    return GPH{.vh = vhOpp, .lph = genOpp.getBlade(ph1Opp, (m + o - i) % m)};
}

uint32 HexExtractor::toGlobalIndex(const GPH& gph) const
{
    const auto& gen = getHexVertexGenerator(gph.vh);
    return gen.global_offsets.global_propellers_offset + gen.toLocalIndex(gph);
}
uint32 HexExtractor::toGlobalIndex(const GCH& gch) const
{
    const auto& gen = getHexVertexGenerator(gch.vh);
    return gen.global_offsets.global_corners_offset + gen.toLocalIndex(gch.vh, gch.lch);
}
size_t HexExtractor::numGlobalPropellers() const
{
    return generators.back().global_offsets.global_propellers_offset + generators.back().nGlobalPropellers();
}

HexVertexGenerator& HexExtractor::getHexVertexGenerator(const VertexHandle& hexVh)
{
    return generators[hexVertexGenerators[hexVh].idx()];
}

const HexVertexGenerator& HexExtractor::getHexVertexGenerator(const VertexHandle& hexVh) const
{
    return generators[hexVertexGenerators[hexVh].idx()];
}

size_t HexExtractor::numGlobalCorners() const
{
    return generators.back().global_offsets.global_corners_offset + generators.back().nGlobalCorners();
}

nlohmann::json HexExtractor::getHexHexReport() const
{
    auto res = Timekeeper::HierarchicalStopWatchResult(HexHex::sw::extract);

    const auto& tetMesh = getInputMesh();
    const auto& hexMesh = getOutputMesh();

    size_t nBoundaryTetVertices = 0;
    size_t nBoundaryTetEdges = 0;
    size_t nBoundaryTetFaces = 0;
    for (auto v_it = tetMesh.vertices_begin(); v_it != tetMesh.vertices_end(); ++v_it) {if (tetMesh.is_boundary(*v_it)) {++nBoundaryTetVertices;}}
    for (auto e_it = tetMesh.edges_begin(); e_it != tetMesh.edges_end(); ++e_it) {if (tetMesh.is_boundary(*e_it)) {++nBoundaryTetEdges;}}
    for (auto f_it = tetMesh.faces_begin(); f_it != tetMesh.faces_end(); ++f_it) {if (tetMesh.is_boundary(*f_it)) {++nBoundaryTetFaces;}}

    nlohmann::json j;
    j["config"] = config;

    j["tet_mesh"] = {
        {"tet_n_vertices", tetMesh.n_vertices()},
        {"tet_n_edges", tetMesh.n_edges()},
        {"tet_n_faces", tetMesh.n_faces()},
        {"tet_n_cells", tetMesh.n_cells()},
        {"n_boundary_vertices", nBoundaryTetVertices},
        {"n_boundary_edges", nBoundaryTetEdges},
        {"n_boundary_faces", nBoundaryTetFaces}
    };
    //j["parametrization"] = { };
    j["stopwatch"] = res;
    j["peak_memory_usage_bytes"] = peak_memory_usage_bytes;

    // Count Propellers and Generators
    size_t nGenerators = generators.size();

    size_t nLocalPropellers = 0;
    size_t nGlobalPropellers = 0;
    size_t nLocalCorners = 0;
    size_t nGlobalCorners = 0;

   for (auto i = 0u; i < generators.size(); ++i)
   {
       const auto& gen = generators[i];
       nLocalPropellers += gen.nLocalPropellers();
       nGlobalPropellers += gen.nLocalPropellers() * gen.nVertices();
       nLocalCorners += gen.nLocalCorners();
       nGlobalCorners += gen.nLocalCorners() * gen.nVertices();
   }

    j["hex_mesh"] = {
        {"hex_n_vertices", hexMesh.n_vertices()},
        {"hex_n_cells", hexMesh.n_cells()},
        {"hex_n_faces", hexMesh.n_faces()},
        {"hex_n_edges", hexMesh.n_edges()}
    };
    j["local_topology"] = {
        {"n_generators", nGenerators},
        {"n_local_propellers", nLocalPropellers},
        {"n_global_propellers", nGlobalPropellers},
        {"n_local_corners", nLocalCorners},
        {"n_global_corners", nGlobalCorners}
    };
    return j;
}

} // namespace HexEx
