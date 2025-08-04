/*
 * Copyright 2025 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */

#pragma once

#include <HexHex/Config/Export.hh>
#include <HexHex/Config.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Utils/HexHexError.hh>
//#include <HexHex/Utils/Utils.hh>
#include <HexHex/LocalTopology/Propellers.hh>
#include <HexHex/Preprocessing/Parametrization.hh>
#include <nlohmann/json_fwd.hpp>

#if _OPENMP
#  include <omp.h>
#endif

namespace HexHex {

class HEXHEX_EXPORT HexExtractor
{
    // Befriend the gang
    //
#if 1
    friend class HexHexPreprocessor;
    friend class Parametrization;
    friend class VertexExtractor;
    friend class HexVertexGenerator;
    friend class HexEdgeExtractor;
    friend class HexCellExtractor;
    friend class PLHexFaceExtractor;
#endif

private:
    HexExtractor(const TetrahedralMesh& tetmesh);


public:
    ~HexExtractor();
    //============================================
    // Constructors
    //============================================

    HexExtractor(const TetrahedralMesh& tetmesh, const OVM::HalfFacePropertyT<Vec3d>& parameters) :
        HexExtractor(tetmesh)
    {
        assert (!tetmesh.needs_garbage_collection());
        parametrization.transferFromHalfFaceProp(tetmesh, parameters);
    }

    template <typename MapT>
    HexExtractor(const TetrahedralMesh& tetmesh, const OVM::CellPropertyT<MapT>& parameters) :
          HexExtractor(tetmesh)
    {
        assert (!tetmesh.needs_garbage_collection());
        parametrization.transferFromCellVertexMapProp(tetmesh, parameters);
    }


    HexExtractor(const HexExtractor& other) = delete;
    HexExtractor(HexExtractor&& other) = delete;
    HexExtractor& operator=(const HexExtractor& other) = delete;
    HexExtractor& operator=(HexExtractor&& other) = delete;

    inline bool hasDetectedError() const {return hasError;}

    //============================================
    // Main functions
    //============================================

    /**
     * The HexHex Pipeline with optional configurations.
     */
    bool extract(Config config = Config());

    void preprocessParametrization();
    bool extractHexVertices();
    void extractLocalTopology();
    void extractHexEdges();
    void extractHexCells();
    void transferFeatureTags();
    void extractPiecewiseLinearHexFaces();

    //============================================
    // mesh getters
    //============================================

    inline const TetrahedralMesh& getInputMesh() const { return inputMesh; }
    inline const HexahedralMesh& getOutputMesh() const {return *outputMesh_;}
    inline std::unique_ptr<HexahedralMesh> takeOutputMesh() {return std::move(outputMesh_);}
    inline HexahedralMesh& getOutputMesh() {return *outputMesh_;}
    PolyhedralMesh& getPiecewiseLinearMesh();
    inline Parametrization& getParametrization() {return parametrization;}
    inline PLMesh & getPLMesh() {return *pwlMesh_;};
    inline PLMesh const & getPLMesh() const {return *pwlMesh_;};
    TetMeshCache const & getCache() const {return cache;}

    template <typename TetMeshT>
    void constructParametrizationMesh(TetMeshT& paramMesh)
    {
        auto tovec = toVec<typename TetMeshT::PointT>;

        paramMesh.clear(false);
        auto trans = paramMesh.template request_face_property<int>("Transitions");
        auto inputVertexPerParamVertex = paramMesh.template request_vertex_property<VertexHandle>("Input Vertices");
        paramMesh.set_persistent(trans);

        for (auto inputCh : inputMesh.cells()) {
            auto paramVhs = std::vector<VertexHandle>();
            for (auto inputVh : cache.cellVertices[inputCh]) {
                auto paramVh = paramMesh.add_vertex(tovec(parametrization.get(inputCh, inputVh)));
                inputVertexPerParamVertex[paramVh] = inputVh;
                paramVhs.push_back(paramVh);
            }

            auto paramCh = paramMesh.add_cell(paramVhs);

            // Store whether transition is identity for visual inspection
            for (auto paramHfh : paramMesh.cell(paramCh).halffaces()) {
                std::vector<VertexHandle> inputVhs; inputVhs.reserve(3);
                for (auto paramVh : paramMesh.get_halfface_vertices(paramHfh)) {
                    inputVhs.push_back(inputVertexPerParamVertex[paramVh]);
                }
                auto inputHfh = inputMesh.find_halfface_in_cell(inputVhs, inputCh);
                if (inputHfh.is_valid())
                    trans[paramHfh.face_handle()] = parametrization.transitions[inputHfh].rotation().index();
            }

        }
    }

    //=================================
    // analysis
    //=================================
    nlohmann::json getHexHexReport() const;

    void printLog() {std::cout << std::endl;}

    template <typename T, typename... Args>
    void printLog(const T& first, const Args&... rest)
    {
        if (!config.verbose) return;
        std::cout << first;
        printLog(rest...);
    }

    void printErr() {std::cerr << std::endl;}

    template <typename T, typename... Args>
    void printErr(const T& first, const Args&... rest)
    {
        std::cerr << first;
        printErr(rest...);
    }

    void addMessage(const HexHexMessage& message)
    {
        if (!message.isOk())
        {
            message.print();
        }
        if (message.isError())
        {
            hasError = true;
        }
    }


    CellHandle owner(const MeshElement& me) const;
    bool isBoundary(const MeshElement& me) const;

private:

    void addHexVertices(const std::vector<Parameter>& params, MeshElement genElem, CellHandle refCh);

    //=======================================
    // predicates
    //=======================================

    // Returns true iff the parameter lies inside the tet or on its boundary. In that case the exact element of intersection (vertex, edge, face or cell) is stored in elem.
    bool computeVertexGeneratorElement(CellHandle ch, const Parameter &param, MeshElement& elem);

    bool computePropellerHolderElement(const CellHandle& ch, const MeshElement& generator, const Parameter& from, const Direction& dir, MeshElement& holder);

    bool computePropellerBladeCasingElement(const CellHandle& ch, const MeshElement& holder, const Parameter& from, const Direction& dir1, const Direction& dir2, MeshElement& casing);

    //=======================================
    // helpers
    //=======================================

    bool edgeIsOnIntegerLattice(EdgeHandle eh) const;
    bool faceIsOnIntegerLattice(FaceHandle fh) const;

    //=============================
    // Simple Getters and Setters
    //=============================
public:

    HexVertexGenerator& getHexVertexGenerator(const VertexHandle& hexVh);

    const HexVertexGenerator& getHexVertexGenerator(const VertexHandle& hexVh) const;

    inline HalfEdgeHandle getHexHalfEdge(const GlobalPropellerHandle& gph) const
    {
        return globalPropellerHexHalfEdges[toGlobalIndex(gph)];
    }

    inline void setHexHalfEdge(const GlobalPropellerHandle& gph, const HalfEdgeHandle heh)
    {
        globalPropellerHexHalfEdges[toGlobalIndex(gph)] = heh;
    }

    uint32 toGlobalIndex(const GPH& gph) const;
    uint32 toGlobalIndex(const GCH& gch) const;

    size_t numGlobalPropellers() const;

    size_t numGlobalCorners() const;

    inline const GlobalPropellerOpposite& getGlobalPropellerOpposite(const GPH& gph) const
    {
        return globalPropellerOpposites[toGlobalIndex(gph)];
    }

    inline bool isGlobalPropellerConnected(const GPH& gph) const
    {
        return getGlobalPropellerOpposite(gph).gph.vh.is_valid();
    }

    inline void setGlobalPropellerOpposite(const GPH& fromGph, const GPH& toGph, const int32 connectionOffset)
    {
        assert(!isGlobalPropellerConnected(fromGph));
        globalPropellerOpposites[toGlobalIndex(fromGph)] = {.gph = toGph, .offset = connectionOffset};
    }

    GlobalPropellerHandle getOppositeBlade(const VertexHandle vh, const LPH lph1, const LPH lph2) const;

public:

    size_t get_peak_memory_usage_bytes() const {return peak_memory_usage_bytes;}

    Config config;

private:

    //=======================================
    // hexhex variables
    //=======================================

    const TetrahedralMesh& inputMesh;
    TetMeshCache cache;

    std::unique_ptr<HexahedralMesh> outputMesh_ = std::make_unique<HexahedralMesh>();
    std::unique_ptr<PLMesh> pwlMesh_;


    bool hasError = false;

    Parametrization parametrization;

    /* Topology */

    std::vector<HexVertexGenerator> generators;
    std::vector<GlobalPropellerOpposite> globalPropellerOpposites;
    std::vector<HalfEdgeHandle> globalPropellerHexHalfEdges;

    /* Output Hex Mesh Properties */
    OVM::VertexPropertyT<GeneratorHandle> hexVertexGenerators = outputMesh_->request_vertex_property<GeneratorHandle>();
    OVM::VertexPropertyT<int> hexVertexTetGeneratorTypes = outputMesh_->create_persistent_vertex_property<int>("v_tet_generator_type", -1).value();
    OVM::VertexPropertyT<int> hexVertexTetGeneratorIndices = outputMesh_->create_persistent_vertex_property<int>("v_tet_generator_index", -1).value();
    OVM::HalfEdgePropertyT<LocalPropellerHandle> hexHalfEdgeLocalPropellers = outputMesh_->request_halfedge_property<LocalPropellerHandle>();

    size_t peak_memory_usage_bytes = 0;
};


} // namespace HexHex
