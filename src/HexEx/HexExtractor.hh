/*
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University
 * Author: Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexEx.
 *
 * HexEx is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * HexEx is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HexEx.  If not, see <http://www.gnu.org/licenses/>.
 */


#pragma once

#include "HPort.hh"
#include "Typedefs.hh"
#include "Dart.hh"
#include "Utils.hh"
#include "MeshConversion.hh"
#include "ExactPredicates.hh"
#include "HexEx/Config/Export.hh"

namespace HexEx {

#ifdef HEXEX_VERBOSE
#define HEXEX_DEBUG_ONLY(x) x
#else
#define HEXEX_DEBUG_ONLY(x)
#endif


class HEXEX_EXPORT HexExtractor
{

    friend Dart;

private:

    enum HVertexType
    {
        VHVertex = 0,
        EHVertex,
        FHVertex,
        CHVertex
    };

    enum CellType
    {
        NotComputed = 0,
        Proper,
        Flipped,
        Degenerate
    };

public:

    HexExtractor();

    HexExtractor(std::string filename);

    template <typename MeshT>
    HexExtractor(const MeshT& tetmesh)
        :
          HexExtractor()
    {
        assert (!tetmesh.needs_garbage_collection());

        convertToHexExTetrahedralMesh(tetmesh, inputMesh);

        // cache cell vertices
        for (auto ch : inputMesh.cells())
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);

    }

    template <typename MeshT, typename ParameterT>
    HexExtractor(const MeshT& tetmesh, PerCellVertexProperty<ParameterT>& parameters)
        :
          HexExtractor(tetmesh)
    {
        // copy parametrization
        for (auto ch: tetmesh.cells())
            for (auto cv_it = tetmesh.cv_iter(ch); cv_it.valid(); ++cv_it)
                vertexParameters[ch][*cv_it] = toVec3d(parameters[ch][*cv_it]);
    }

    ~HexExtractor(){}

    HexExtractor(const HexExtractor& other) = delete;
    HexExtractor(HexExtractor&& other) = delete;
    HexExtractor& operator=(const HexExtractor& other) = delete;
    HexExtractor& operator=(HexExtractor&& other) = delete;

    void extract();

    template <typename HexMeshT>
    void getHexMesh(HexMeshT& hexMesh, bool optimizedMergePosition = true)
    {
        HEXEX_DEBUG_ONLY(std::cout << "converting mesh to hex mesh" << std::endl;)

        hexMesh.clear(false);

        calculateDifferencesInDartTypes();

        mergeEquivalenceClassesOfAllDarts();

        computeEquivalenceClasses();

        addVerticesToHexMesh(hexMesh, optimizedMergePosition);
        addCellsToHexMesh(hexMesh);
    }

    template <typename TetMeshT>
    void getParametrizationMesh(TetMeshT& paramMesh)
    {
        auto tovec = toVec<typename TetMeshT::PointT>;

        paramMesh.clear(false);

        for (auto ch : inputMesh.cells())
        {
            auto vertices = cellVertices[ch];

            auto newVertices = std::vector<VertexHandle>();
            for (auto vh : vertices)
                newVertices.push_back(paramMesh.add_vertex(tovec(parameter(ch, vh))));

            paramMesh.add_cell(newVertices);
        }

    }

    template <typename TetMeshT>
    void getInputMesh(TetMeshT& tetMesh)
    {
        auto tovec = toVec<typename TetMeshT::PointT>;

        tetMesh.clear(false);

        auto parametrization = tetMesh.template request_cell_property<VertexMapProp<Vec3d>>("Parametrization");
        tetMesh.set_persistent(parametrization);

        for (auto v_it = inputMesh.vertices_begin(); v_it != inputMesh.vertices_end(); ++v_it)
            tetMesh.add_vertex(tovec(inputMesh.vertex(*v_it)));

        for (auto c_it = inputMesh.cells_begin(); c_it != inputMesh.cells_end(); ++c_it)
            tetMesh.add_cell(inputMesh.get_cell_vertices(*c_it));;


        for (auto c_it = inputMesh.cells_begin(); c_it != inputMesh.cells_end(); ++c_it)
            for (auto cv_it = inputMesh.cv_iter(*c_it); cv_it.valid(); ++cv_it)
                parametrization[*c_it][*cv_it] = parameter(*c_it, *cv_it);

    }

    void writeToFile(std::string filename);

    template <typename PointT>
    void copySanitizedParametrization(std::map<CellHandle, std::map<VertexHandle, PointT>>& parameters) {
        sanitizeParametrization(false, false);

        // copy parametrization
        for (auto ch: inputMesh.cells())
            for (auto cv_it = inputMesh.cv_iter(ch); cv_it.valid(); ++cv_it) {
                auto pm = vertexParameters[ch][*cv_it];
                parameters[ch][*cv_it] = PointT(pm[0], pm[1], pm[2]);
            }
    }

#ifdef HEXEX_TESTING
public:
#else
private:
#endif
    
    template <typename MeshT>
    void getIrregularEdgeMesh(MeshT& edgeMesh)
    {
        auto tovec = toVec<typename MeshT::PointT>;

        auto edgeValences2 = edgeMesh.template request_edge_property<int>("Edge Valence");
        edgeMesh.set_persistent(edgeValences2);

        calculateValences();

        edgeMesh.clear(false);

        for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
        {
            if (isSingularEdge(*e_it))
            {
                auto ch = *inputMesh.hec_iter(inputMesh.halfedge_handle(*e_it, 0));
                if (!ch.is_valid())
                    continue;


                auto e = inputMesh.edge(*e_it);
                auto vh1 = edgeMesh.add_vertex(tovec(inputMesh.vertex(e.from_vertex())));
                auto vh2 = edgeMesh.add_vertex(tovec(inputMesh.vertex(e.to_vertex())));

                auto eh = edgeMesh.add_edge(vh1, vh2);
                edgeValences2[eh] = edgeValences[*e_it];
            }
        }

    }


    template <typename MeshT>
    void getHPortMesh(MeshT& hportMesh, double scaling = 0.3)
    {
        auto tovec = toVec<typename MeshT::PointT>;

        hportMesh.clear(false);

        for (auto v_it = intermediateHexMesh.vertices_begin(); v_it != intermediateHexMesh.vertices_end(); ++v_it)
        {
            auto& hports = hPortsOnVertex[*v_it];

            if (hports.empty())
                continue;

            auto vertexPos1 = tovec(getPosition(hports.front().parameter(), hports.front().cell()));

            for (auto& p : hports)
            {
                auto vertexParam1 = p.parameter();
                auto vertexParam2 = vertexParam1 + p.dir()*scaling;
                auto vertexPos2 = tovec(getPosition(vertexParam2, p.cell()));
                if ((vertexPos2 - vertexPos1).length() > 10)
                    vertexPos2 = 10.0/(vertexPos2-vertexPos1).length() * (vertexPos2-vertexPos1);

                hportMesh.add_edge(hportMesh.add_vertex(vertexPos1),hportMesh.add_vertex(vertexPos2));
            }
        }
    }

    void extractHVertices();

    void enumerateHPorts();
    void enumerateVertexHPorts(VertexHandle vh);
    void enumerateEdgeHPorts(VertexHandle hexVh);
    void enumerateFaceHPorts(VertexHandle hexVh);
    void enumerateCellHPorts(VertexHandle vh);

    bool isDartInCell(CellHandle ch, HPortHandle port, Direction refDir, Direction normalDir); // assumes port is in cell
    bool isDartInCell(CellHandle ch, VertexHandle hexVh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);

    void enumerateDarts(HPortHandle port, bool secondary = false);
    void enumerateDarts();
    void enumerateSecondaryDarts();

    int getNumDarts(bool exclude_anihillated);

    void traceDarts();
    bool traceDart(Dart& dart);
    bool traceDart2(Dart& dart);

//    bool connectDartToNextDart(Dart& dart);
//    void connectDartsToNextDart();

    bool connectDartToPreviousSecondaryDart(Dart& dart);
    void connectDartsToPreviousSecondaryDart();

//    bool connectDartToNeighborDart(Dart& dart);
//    void connectDartsToNeighborDart();

    bool connectDartToNeighborSecondaryDart2(Dart& dart);
    void connectDartsToNeighborSecondaryDart();

    bool connectDartToOppositeSecondaryDart(Dart& dart);
    void connectDartsToOppositeSecondaryDart();


    bool alpha0FaceTest(HalfFaceHandle hfh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);
    HalfFaceHandle alpha0NextFace(HalfFaceHandle prevFace, CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);

    bool alpha1FaceTest(HalfFaceHandle hfh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);
    HalfFaceHandle alpha1NextFace(HalfFaceHandle prevFace, CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);

    bool alpha2FaceTest(HalfFaceHandle hfh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);
    HalfFaceHandle alpha2NextFace(HalfFaceHandle prevFace, CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);

    bool alpha3FaceTest(HalfFaceHandle hfh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);
    HalfFaceHandle alpha3NextFace(HalfFaceHandle prevFace, CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);


    std::shared_ptr<Dart> getDart(CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);
    std::shared_ptr<Dart> getSecondaryDart(CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir);

    HalfEdgeHandle addEdge(HPortHandle p1, HPortHandle p2);


    void mergeEquivalenceClasses(Dart* dart);
    void mergeEquivalenceClassesOfAllDarts();

    void mergeVerticesOfOppositeFaces(HalfFaceHandle hfh);
    void mergeVerticesOfOppositeFaces();

    // connects two darts if both are not nullptr and both are different from dart and otherDart
    // else they are disconnected
    template <int i>
    void reconnectOrDisconnect(Dart* dart, Dart* otherDart, Dart* d1, Dart* d2)
    {
        if (d2 != dart && d1 != otherDart)
        {
            if (d2 != nullptr && d1 != nullptr)
            {
                d2->connectAlpha<i>(d1);
            }
            else
            {
                if (d1 != nullptr)
                    d1->disconnectAlpha<i>();
                if (d2 != nullptr)
                    d2->disconnectAlpha<i>();
            }
        }
    }

    bool fixProblem(Dart* dart, Dart* otherDart, bool mergeVertices);

    bool fixProblem1(Dart* dart, bool mergeVertices);
    bool fixProblems1(bool mergeVertices);

    bool fixProblem2(Dart* dart, bool mergeVertices);
    bool fixProblems2(bool mergeVertices);

    bool fixProblem3(Dart* dart, bool mergeVertices);
    bool fixProblems3(bool mergeVertices);

    bool fixProblem4(Dart* dart, bool mergeVertices);
    bool fixProblems4(bool mergeVertices);


    void fixProblems(int iterations = 3, bool mergeVertices = true);

//    void updateDarts();

//    void computeLocalUVFromDarts(Dart* dart, CellHandle refCell, Transition transitionToRefCell, std::vector<Dart*>& processedDarts);
//    void computeLocalUVFromDarts(CellHandle ch);
//    void computeLocalUVsFromDarts();

    void computeLocalUVFromSecondaryDarts(Dart* dart, CellHandle refCell, Transition transitionToRefCell, std::vector<Dart*>& processedDarts);
    void computeLocalUVFromSecondaryDarts(CellHandle ch);
    void computeLocalUVsFromSecondaryDarts();

    void computeLocalFaceUVsFromDarts(HalfFaceHandle hexHfh);

    void showLocalUVs(CellHandle ch);

    void computeLocalUVsAlternative(HalfFaceHandle hexHfh);

    HalfFaceHandle extractFaceLazy(HPortHandle& normalPort, HPortHandle& tracePort);

    void deleteHexEdges();
    void deleteHexFaces();
    void deleteHexCells();

//    void extractFaceFromDarts(Dart* dart);
    void extractFaceFromSecondaryDarts(Dart* dart);
    void extractFacesFromDarts();

//    void extractCellFromDarts(HalfFaceHandle hfh);
    void extractCellFromSecondaryDarts(HalfFaceHandle hfh);
    void extractCellsFromDarts();


    std::vector<HalfFaceHandle> getAdjacentHalffaces(HalfFaceHandle hfh);
//    std::vector<HalfFaceHandle> getAdjacentHalffacesFromDarts(HalfFaceHandle hfh);
    std::vector<HalfFaceHandle> getAdjacentHalffacesFromSecondaryDarts(HalfFaceHandle hfh);


    void extractTransitionFunction(FaceHandle fh);
    void extractTransitionFunctions();

    void removeDegeneratedTetrahedra();

    void computeCellTypes();
    CellType computeCellType(CellHandle ch);

    void randomizeParametrization(double offsetSize, double keepBoundary = false);

    void sanitizeParametrization(bool snapBoundary = true, bool extremeTruncation = false);
    void truncatePrecision(bool extremeTruncation = false);

    TetrahedralMesh& getInputMesh() { return inputMesh; }
    PolyhedralMesh& getIntermediateHexMesh()    { return intermediateHexMesh;   }

    void checkThisOneProperty();


    template <typename HexMeshT>
    void addVerticesToHexMesh(HexMeshT& hexMesh, bool optimizedMergePosition)
    {
        auto tovec = toVec<typename HexMeshT::PointT>;

        equivalenceClassVertices.clear();
        equivalenceClassVertices.resize(equivalenceClasses.size());
        for (auto i = 0u; i < equivalenceClasses.size(); ++i)
            if (!equivalenceClasses[i].empty())
            {
                auto pos = tovec(getMergePosition(i, optimizedMergePosition));

                auto vh = hexMesh.add_vertex(pos);
                equivalenceClassVertices[i] = vh;

                if (equivalenceClasses[i].size() > 1)
                {
                    numMerges++;
                    numMergedVertices += equivalenceClasses[i].size();
                }
            }
    }


    template <typename HexMeshT>
    void addCellsToHexMesh(HexMeshT& hexMesh)
    {
        for (auto ch : intermediateHexMesh.cells())
        {
            if (differenceBetweenInvertedAndProperDartsPerCell[ch] == -48)
            {

                auto vertices = std::vector<VertexHandle>();

                auto vertexCube = std::array<std::array<std::array<VertexHandle, 2>,2>,2>();

                for (auto cv_it = intermediateHexMesh.cv_iter(ch); cv_it.valid(); ++cv_it)
                {
                    auto equivalenceClassId = equivalenceClassIds[*cv_it];
                    auto newVh = equivalenceClassVertices[equivalenceClassId];

                    if (localCellUVs[ch].find(*cv_it) == localCellUVs[ch].end())
                        continue;
                    auto localUV = localCellUVs[ch][*cv_it];
                    /*
                    assert((localUV[0] == 0) || (localUV[0] == 1));
                    assert((localUV[1] == 0) || (localUV[1] == 1));
                    assert((localUV[2] == 0) || (localUV[2] == 1));
                    */

                    if (!((localUV[0] == 0) || (localUV[0] == 1)) ||
                            !((localUV[1] == 0) || (localUV[1] == 1)) ||
                            !((localUV[2] == 0) || (localUV[2] == 1)) )
                        HEXEX_DEBUG_ONLY(std::cout << "Error in local UVs" << std::endl);
                    else
                        vertexCube[(int)localUV[0]][(int)localUV[1]][(int)localUV[2]] = newVh;
                }

                vertices.push_back(vertexCube[0][0][1]);
                vertices.push_back(vertexCube[1][0][1]);
                vertices.push_back(vertexCube[1][1][1]);
                vertices.push_back(vertexCube[0][1][1]);
                vertices.push_back(vertexCube[0][0][0]);
                vertices.push_back(vertexCube[0][1][0]);
                vertices.push_back(vertexCube[1][1][0]);
                vertices.push_back(vertexCube[1][0][0]);

                auto vertices2 = vertices;

                std::sort(vertices2.begin(), vertices2.end());
                vertices2.erase(std::unique(vertices2.begin(), vertices2.end()), vertices2.end());

                if (vertices2.size() == 8)
                    if (std::find(vertices.begin(), vertices.end(), VertexHandle()) == vertices.end())
                        hexMesh.add_cell(vertices,false);
            }
        }
    }

    template <typename PolyMeshT>
    void addCellsToHexMeshWithLargeCells(PolyMeshT& polyMesh)
    {
        auto differenceBetweenInvertedAndProperDartsPerHalffacePoly = polyMesh.template request_halfface_property<int>();
        polyMesh.set_persistent(differenceBetweenInvertedAndProperDartsPerHalffacePoly);

        for (auto ch : intermediateHexMesh.cells())
        {
            auto ds = getAllDartsInCell(halffaceDarts[intermediateHexMesh.cell(ch).halffaces().front()].front());

            auto processed_darts = std::set<Dart*>();

            auto halffaces = std::vector<HalfFaceHandle>();

            for (auto d : ds)
                if (d->isPrimary() && std::find(processed_darts.begin(), processed_darts.end(), d) == processed_darts.end())
                {
                    auto halfface_vertices = std::vector<VertexHandle>();

                    auto tmp_d = d;

                    // rotate backwards as far as possible or back to d
                    int steps = 0;
                    do
                    {
                        if (++steps>50)
                            break;

                        if ((tmp_d->getAlpha<1>() == nullptr) || (tmp_d->getAlpha<1>()->getAlpha<0>() == nullptr))
                            break;

                        tmp_d = tmp_d->getAlpha<1>()->getAlpha<0>();
                    }
                    while (tmp_d != d);

                    do
                    {
                        processed_darts.insert(tmp_d);
                        auto equivalenceClassId = equivalenceClassIds[tmp_d->getVertex()];
                        auto newVh = equivalenceClassVertices[equivalenceClassId];
                        halfface_vertices.push_back(newVh);

                        if ((tmp_d->getAlpha<0>() == nullptr) || (tmp_d->getAlpha<0>()->getAlpha<1>() == nullptr))
                            break;

                        tmp_d = tmp_d->getAlpha<0>()->getAlpha<1>();
                    }
                    while (tmp_d != d);

                    if (halfface_vertices.size() > 2)
                    {
                        auto hfh = polyMesh.halfface_extensive(halfface_vertices);
                        if (!hfh.is_valid())
                            hfh = polyMesh.halfface_handle(polyMesh.add_face(halfface_vertices),0);
                        differenceBetweenInvertedAndProperDartsPerHalffacePoly[hfh] = halfface_vertices.size()*2;
                        halffaces.push_back(hfh);
                    }

                }

            if (halffaces.size() > 0)
                polyMesh.add_cell(halffaces, false);
        }
    }


    template <typename PolyMeshT>
    void getHexMeshWithLargeCells(PolyMeshT& polyMesh, bool optimizedMergePosition)
    {
        HEXEX_DEBUG_ONLY(std::cout << "converting mesh to hex mesh" << std::endl;)

        calculateDifferencesInDartTypes();

        mergeEquivalenceClassesOfAllDarts();
        computeEquivalenceClasses();

        addVerticesToHexMesh(polyMesh, optimizedMergePosition);
        addCellsToHexMeshWithLargeCells(polyMesh);
    }

    void computeEquivalenceClasses();
    void degeneracyEquivalenceClassJoin(VertexHandle hexVh);
    void joinEquivalenceClasses(VertexHandle vh1, VertexHandle vh2);

    Position getMergePosition(unsigned int equivalenceClass, bool optimizedMergePosition);


    Position getComplicatedMergePosition(unsigned int equivalenceClass);
    Matrix4x4d getQuadric(Vec3d n, double d);
    std::vector<HalfFaceHandle> getBoundaryHalffacesOfHexVertex(VertexHandle hexVh);
    Vec3d getNormal(HalfFaceHandle hfh);
    double getArea(HalfFaceHandle hfh);
    bool isBoundaryHexVertex(VertexHandle hexVh);

    // predicates

    bool isInCell(CellHandle ch, Parameter param);
    bool isInCellRelaxed(CellHandle ch, Parameter param);
    bool isInFace(HalfFaceHandle hfh, Parameter param);
    bool isInFaceRelaxed(HalfFaceHandle hfh, Parameter param);
    bool isOnEdge(CellHandle ch, EdgeHandle eh, Parameter param);
    bool isOnEdgeRelaxed(CellHandle ch, EdgeHandle eh, Parameter param);

    bool pointsIntoCell(CellHandle ch, VertexHandle vh, Direction dir);
    bool pointsIntoCellRelaxed(CellHandle ch, VertexHandle vh, Direction dir);
    bool pointsIntoCell(CellHandle ch, EdgeHandle eh, Direction dir, Parameter param);
    bool pointsIntoCell(HalfFaceHandle hfh, Direction dir, Parameter param);
    bool pointsIntoCellRelaxed(HalfFaceHandle hfh, Direction dir, Parameter param);
    bool pointsAlongHalfEdge(CellHandle ch, HalfEdgeHandle heh, Direction dir);
    bool pointsAlongHalfEdgeFromVertex(CellHandle ch, HalfEdgeHandle heh, Direction dir);
    bool pointsAlongEdge(CellHandle ch, EdgeHandle eh, Direction dir, Parameter param);
    bool pointsAlongFace(HalfFaceHandle hfh, Direction dir, Parameter param);
    bool pointsIntoFace(HalfFaceHandle hfh, VertexHandle vh, Direction dir);
    bool pointsIntoFace(HalfFaceHandle hfh, EdgeHandle eh, Direction dir, Parameter param);


    bool isCellDegenerate(CellHandle ch);
    bool isCellFlipped(CellHandle ch);
    bool isFaceDegenerate(HalfFaceHandle hfh);
    bool areColinear(CellHandle ch, HalfEdgeHandle heh, Direction dir);

    bool isSingularVertex(VertexHandle vh);
    bool isSingularEdge(EdgeHandle eh);

    bool isFixPointRecursive(Parameter parameter, CellHandle ch, VertexHandle vh, std::set<CellHandle> visited);
    bool isFixPoint(Parameter parameter, CellHandle ch, VertexHandle vh);

    // end predicates

    HalfEdgeHandle getHalfEdgePointingIntoDirection(CellHandle ch, HalfEdgeHandle heh, Direction dir);
    HalfEdgeHandle getHalfEdgePointingIntoDirection(CellHandle ch, EdgeHandle eh, Direction dir);

    bool intersectsFace(HalfFaceHandle hfh, Parameter start, Parameter end);
    bool intersectsFaceRelaxed(HalfFaceHandle hfh, Parameter start, Parameter end);


    HalfFaceHandle rotateAroundHalfedge(CellHandle startCell, HalfEdgeHandle currentEdge, bool ccw = true);

    VertexHandle getVertexWithParam(CellHandle ch, Parameter param);
    EdgeHandle getEdgeWithParam(CellHandle ch, Parameter param);
    HalfFaceHandle getHalffaceWithParam(CellHandle ch, Parameter param);

    int numIncidentSingularEdges(VertexHandle vh);
    HalfEdgeHandle getIncidentSingularEdge(VertexHandle vh);
    HalfEdgeHandle getIncidentSingularEdge(VertexHandle vh, CellHandle ch);
    CellHandle getIncidentCellIncidentToSingularEdge(VertexHandle vh);
    bool isIncident(HalfEdgeHandle heh, CellHandle ch);

    const Position inputPosition(VertexHandle vh) {  return inputMesh.vertex(vh); }

    Parameter& parameter(CellHandle ch, VertexHandle vh) { return vertexParameters[ch][vh]; }
    std::vector<Parameter> getParameters(CellHandle ch, std::vector<VertexHandle> vhs);
    std::vector<Parameter> getParameters(CellHandle ch);
    std::vector<Parameter> getParameters(HalfFaceHandle hfh);
    std::vector<Parameter> getParameters(HalfFaceHandle hfh, HalfEdgeHandle heh);

    Position getPosition(Parameter param, CellHandle ch);
    Parameter getParameter(Position pos, CellHandle ch);

    Parameter getHexVertexParameter(VertexHandle hexVh, CellHandle ch);


    Parameter getParameterNormal(HalfFaceHandle hfh);
    Matrix4x4d getParametrizationMatrix(Position  p, Position  q, Position  r, Position  s,
                                        Parameter u, Parameter v, Parameter w, Parameter t);
    Matrix4x4d getInverseParametrizationMatrix(Position  p, Position  q, Position  r, Position  s,
                                               Parameter u, Parameter v, Parameter w, Parameter t);
    Matrix4x4d getInverseParametrizationMatrix(Position  p, Position  q, Position  r,
                                               Parameter u, Parameter v, Parameter w);

    Matrix4x4d getLocalFrame(Direction dir, Direction refDir, Direction normal, Parameter inputPosition);

    void calculateValences();

    void fixSingularityPoint(VertexHandle vh, CellHandle& ch);
    void projectBoundaryFaces();

    void propagateVertexParameterRecursive(Parameter param, VertexHandle vh, CellHandle ch, std::set<CellHandle>& visited);
    void propagateVertexParameterRecursive2(Parameter param, VertexHandle vh, CellHandle ch, std::set<CellHandle>& toBeProcessed);
    void propagateVertexParameter(Parameter parameter, VertexHandle vh, CellHandle startCell);

    const Transition& getTransitionFunction(HalfFaceHandle hfh);
    void setTransitionFunction(HalfFaceHandle hfh, Transition transitionFunction);

    Transition getTransitionFunctionAroundHalfedge(CellHandle ch, HalfEdgeHandle heh);
    Transition getTransitionFunctionRecursive(CellHandle currentCell, CellHandle toCell, VertexHandle vh, Transition tranFun, std::set<CellHandle>& visited);
    Transition getTransitionFunction(CellHandle fromCell, CellHandle toCell, VertexHandle vh);

    Transition getTransitionFunction(CellHandle fromCell, CellHandle toCell, EdgeHandle eh);

    Matrix4x4d transitionFrame(Parameter u, Parameter v, Parameter w);

    HPortHandle findPort(CellHandle ch, Direction dir, Parameter param);

    Parameter projectedParam(Parameter param, CellHandle ch, HalfEdgeHandle heh);

    double parametrizationAngle(HalfFaceHandle hfh1, HalfFaceHandle hfh2, HalfEdgeHandle heh);

    int edgeValence(EdgeHandle eh);
    void calculateEdgeSingularity(EdgeHandle eh);
    void calculateEdgeSingularities();

    void setTranslation(GridIsomorphism& tranFun, Parameter translation);
    void setTranslation(Matrix4x4d& tranFun, Parameter translation);


    std::vector<Dart*> getDartsBetweenDarts01(Dart* d_s, Dart* d_e);
    std::vector<Dart*> getDartsBetweenDarts12(Dart* d_s, Dart* d_e);
    std::vector<Dart*> getDartsBetweenDarts0121(Dart* d_s, Dart* d_e);

    Transition getTransitionBetweenDarts01(Dart* d_s, Dart* d_e);
    Transition getTransitionBetweenDarts12(Dart* d_s, Dart* d_e);
    Transition getTransitionBetweenDarts0121(Dart* d_s, Dart* d_e);
    Transition getTransitionBetweenDarts01Backward(Dart* d_s, Dart* d_e);
    Transition getTransitionBetweenDarts0121Backward(Dart* d_s, Dart* d_e);

    void reconnectDarts(CellHandle ch, HalfFaceHandle hfh, VertexHandle vh0, VertexHandle vh1);
    void reconnectSecondaryDartsOld(Dart* d1, Dart* d2);

    std::vector<Dart*> getAllDartsInCell(Dart* d);

    void calculateDifferencesInDartTypes(CellHandle ch);
    void calculateDifferencesInDartTypes(HalfFaceHandle hfh);
    void calculateDifferencesInDartTypes();


    void doTransition(HalfFaceHandle hfh, CellHandle& ch);
    void doTransition(HalfFaceHandle hfh, std::vector<Parameter>& params);
    void doTransition(HalfFaceHandle hfh, Parameter& parameter);
    void doTransition(HalfFaceHandle hfh, Direction& dir);
    void doTransition(HalfFaceHandle hfh, Transition& tranFun);

    template <typename T, typename... Rest>
    void doTransition(HalfFaceHandle hfh, T& target, Rest&... rest)
    {
        HEXEX_DEBUG_ONLY(if (isFaceDegenerate(hfh))
            std::cout << "warning: transitioning through a degenerate face. Transition function might be wrong." << std::endl;)
        doTransition(hfh, target);
        doTransition(hfh, rest...);
    }

    TetrahedralMesh inputMesh;
    PolyhedralMesh intermediateHexMesh;

    PerCellVertexProperty<Parameter> vertexParameters;

    CellProperty<std::vector<HPortHandle>> hPortsInCell;

    bool edgeSingularitiesCalculated;
    EdgeProperty<int>  edgeValences;
    EdgeProperty<bool> edgeSingularity;

    HalfFaceProperty<Transition> transitionFunctions;

    VertexProperty<HVertexType> vertexTypes;
    VertexProperty<std::vector<HPortHandle>> hPortsOnVertex;
    VertexProperty<CellHandle> incidentCellInInputMesh;
    VertexProperty<Parameter> hexvertexParameter;
    VertexProperty<VertexHandle> incidentVerticesPerVertex;
    EdgeProperty<std::vector<VertexHandle>> incidentVerticesPerEdge;
    FaceProperty<std::vector<VertexHandle>> incidentVerticesPerFace;
    CellProperty<std::vector<VertexHandle>> incidentVerticesPerCell;
    CellProperty<int> differenceBetweenInvertedAndProperDartsPerCell;
    HalfFaceProperty<int> differenceBetweenInvertedAndProperDartsPerHalfface;
    VertexProperty<int> incidentElementId;
    VertexProperty<std::vector<std::shared_ptr<Dart>>> darts;
    VertexProperty<std::vector<std::shared_ptr<Dart>>> secondaryDarts;
    VertexProperty<Parameter> localUVs;
    PerCellVertexProperty<Parameter> localCellUVs;
    VertexProperty<int> equivalenceClassIds;
    std::vector<std::vector<int>> equivalenceClasses;
    std::vector<VertexHandle> equivalenceClassVertices;

    HalfEdgeProperty<HPortHandle> incidentHPort;

    HalfFaceProperty<std::vector<Dart*>> halffaceDarts;
    HalfFaceProperty<std::vector<Dart*>> halffaceSecondaryDarts;

    static const Transition identity;
    std::vector<Transition> all24Transitions;

    bool transitionFunctionsComputed;
    bool cellTypesComputed;

    int numMerges;
    int numMergedVertices;

    int getNumMerges() { return numMerges; }
    int getNumMergedVertes() { return numMergedVertices; }



    long isCellFlippedCalls;
    long isCellDegenerateCalls;
    long isFaceDegenerateCalls;

    CellProperty<CellType> cellTypes;
    FaceProperty<CellType> faceTypes;
    CellProperty<std::vector<VertexHandle>> cellVertices;

    long numDartsTraced;
    long numDartTraceLoops;

    long callsToFindPort;
    long portsCheckedInFindPort;


};

}

#if (defined(WIN32) || defined(_WIN32)) && (_MSC_VER < 1900)
// specialize the has_input_operator class because the Visual Studio Compiler prior
// to VS 2015 did not support expression SFINAE and thus leads to wrong results.

namespace OpenVolumeMesh
{

#define NO_INPUT_OPERATOR(Type) \
    template <> \
    class has_input_operator<std::istream, Type > \
    { \
    public: \
        enum { bool_value = false }; \
        typedef bool_type<bool_value> type; \
        static type value; \
    };


    NO_INPUT_OPERATOR(HexEx::HexExtractor::HVertexType)
    NO_INPUT_OPERATOR(HexEx::HexExtractor::CellType)
    NO_INPUT_OPERATOR(HexEx::HPortHandle)
    NO_INPUT_OPERATOR(HexEx::GridIsomorphism)
    NO_INPUT_OPERATOR(std::shared_ptr<HexEx::Dart>)
    NO_INPUT_OPERATOR(HexEx::Dart*)

#undef NO_INPUT_OPERATOR

}

#endif
