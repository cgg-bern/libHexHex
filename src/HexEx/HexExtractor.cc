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


#define _USE_MATH_DEFINES
#include <math.h>

#include "HexExtractor.hh"

#include <iomanip>
#include <tuple>
#include <queue>
#include <algorithm>

#include "Utils.hh"
#include "ExactPredicates.hh"

#include "DerivedExactPredicates.hh"
#include "Direction.hh"
#include "FileAccessor.hh"


using namespace OpenVolumeMesh;

namespace HexEx {
    
const Transition HexExtractor::identity = Transition();

HexExtractor::HexExtractor()
    :
      inputMesh(),
      intermediateHexMesh(PolyhedralMesh()),
      vertexParameters(inputMesh.request_cell_property<VertexMapProp<Parameter>>()),
      hPortsInCell(inputMesh.request_cell_property<std::vector<HPortHandle>>()),
      edgeSingularitiesCalculated(false),
      edgeValences(inputMesh.request_edge_property<int>()),
      edgeSingularity(inputMesh.request_edge_property<bool>()),
      transitionFunctions(inputMesh.request_halfface_property<Transition>()),
      vertexTypes(intermediateHexMesh.request_vertex_property<HVertexType>()),
      hPortsOnVertex(intermediateHexMesh.request_vertex_property<std::vector<HPortHandle>>()),
      incidentCellInInputMesh(intermediateHexMesh.request_vertex_property<CellHandle>()),
      hexvertexParameter(intermediateHexMesh.request_vertex_property<Parameter>()),
      incidentVerticesPerVertex(inputMesh.request_vertex_property<VertexHandle>()),
      incidentVerticesPerEdge(inputMesh.request_edge_property<std::vector<VertexHandle>>()),
      incidentVerticesPerFace(inputMesh.request_face_property<std::vector<VertexHandle>>()),
      incidentVerticesPerCell(inputMesh.request_cell_property<std::vector<VertexHandle>>()),
      differenceBetweenInvertedAndProperDartsPerCell(intermediateHexMesh.request_cell_property<int>()),
      differenceBetweenInvertedAndProperDartsPerHalfface(intermediateHexMesh.request_halfface_property<int>()),
      incidentElementId(intermediateHexMesh.request_vertex_property<int>()),
      darts(intermediateHexMesh.request_vertex_property<std::vector<std::shared_ptr<Dart>>>()),
      secondaryDarts(intermediateHexMesh.request_vertex_property<std::vector<std::shared_ptr<Dart>>>()),
      localUVs(intermediateHexMesh.request_vertex_property<Parameter>()),
      localCellUVs(intermediateHexMesh.request_cell_property<VertexMapProp<Parameter>>()),
      equivalenceClassIds(intermediateHexMesh.request_vertex_property<int>()),
      equivalenceClasses(),
      incidentHPort(intermediateHexMesh.request_halfedge_property<HPortHandle>()),
      halffaceDarts(intermediateHexMesh.request_halfface_property<std::vector<Dart*>>()),
      halffaceSecondaryDarts(intermediateHexMesh.request_halfface_property<std::vector<Dart*>>()),
      transitionFunctionsComputed(false),
      cellTypesComputed(false),
      numMerges(0),
      numMergedVertices(0),
      isCellFlippedCalls(0),
      isCellDegenerateCalls(0),
      isFaceDegenerateCalls(0),
      cellTypes(inputMesh.request_cell_property<CellType>("Cell Types", NotComputed)),
      faceTypes(inputMesh.request_face_property<CellType>("Face Types", NotComputed)),
      cellVertices(inputMesh.request_cell_property<std::vector<VertexHandle>>()),
      numDartsTraced(0),
      numDartTraceLoops(0),
      callsToFindPort(0),
      portsCheckedInFindPort(0)
{
    exactinit();

    // initialize the 24 Transitions
    for (auto d1 : getAll6Directions())
        for (auto d2 : getAll6Directions())
            for (auto d3 : getAll6Directions())
                if (d1.vector() % d2.vector() == d3.vector())
                {
                    double entries[] = {d1.vector()[0], d2.vector()[0], d3.vector()[0], 0,
                                        d1.vector()[1], d2.vector()[1], d3.vector()[1], 0,
                                        d1.vector()[2], d2.vector()[2], d3.vector()[2], 0,
                                        0,              0,              0, 1 };
                    all24Transitions.push_back(Transition(Matrix4x4d(entries)));
                }
}

HexExtractor::HexExtractor(std::string filename)
    :
      HexExtractor()
{
    HexEx::readFromFile(filename, inputMesh, vertexParameters);

    // cache cell vertices
    for (auto ch : inputMesh.cells())
        cellVertices[ch] = inputMesh.get_cell_vertices(ch);
}

void HexExtractor::extract()
{
    sanitizeParametrization(true, false);
    extractHVertices();
    enumerateHPorts();
    enumerateDarts();
    traceDarts();
    connectDartsToPreviousSecondaryDart();
    connectDartsToNeighborSecondaryDart();
    connectDartsToOppositeSecondaryDart();
    fixProblems(9999, true);
    extractFacesFromDarts();
    extractCellsFromDarts();
    computeLocalUVsFromSecondaryDarts();
}


void HexExtractor::writeToFile(std::string filename)
{
    HexEx::writeToFile(filename, inputMesh, vertexParameters);
}

void HexExtractor::extractHVertices()
{
    intermediateHexMesh.clear(false);

    for (auto vh : inputMesh.vertices())
        incidentVerticesPerVertex[vh] = VertexHandle();
    for (auto eh : inputMesh.edges())
        incidentVerticesPerEdge[eh].clear();
    for (auto fh : inputMesh.faces())
        incidentVerticesPerFace[fh].clear();
    for (auto ch : inputMesh.cells())
        incidentVerticesPerCell[ch].clear();


    HEXEX_DEBUG_ONLY(std::cout << "Extracting vertices on vertices" << std::endl;)

    for (auto v_it = inputMesh.vertices_begin(); v_it != inputMesh.vertices_end(); ++v_it)
    {
        HEXEX_DEBUG_ONLY(if ((v_it->idx() % 10000) == 0)
          std::cout << "Processing vertex " << v_it->idx() << " of " <<  inputMesh.n_vertices() << std::endl;)

        auto adjacentCell = *inputMesh.vc_iter(*v_it);
        if (!adjacentCell.is_valid())
        {
            HEXEX_DEBUG_ONLY(std::cout << "warning: isolated vertex" << std::endl;)
            continue;
        }
        auto p = inputPosition(*v_it);
//        p = inputMesh.vertex(*v_it);
        auto u = parameter(adjacentCell, *v_it);

        if (u == roundVector(u))
        {
            auto vh = intermediateHexMesh.add_vertex(p);
            vertexTypes[vh] = VHVertex;
            incidentCellInInputMesh[vh] = adjacentCell;
            hexvertexParameter[vh] = u;
            incidentElementId[vh] = v_it->idx();
            incidentVerticesPerVertex[*v_it] = vh;
        }
    }

    HEXEX_DEBUG_ONLY(std::cout << "Extraction vertices on edges" << std::endl;)

    for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
    {

        HEXEX_DEBUG_ONLY(if ((e_it->idx() % 10000) == 0)
          std::cout << "Processing edge " << e_it->idx() << " of " <<  inputMesh.n_edges() << std::endl;)

        auto e = inputMesh.edge(*e_it);
        auto he = inputMesh.halfedge_handle(*e_it,0);
        auto adjacentCell = *inputMesh.hec_iter(he);

        auto p = inputPosition(e.from_vertex());
        auto q = inputPosition(e.to_vertex());

//        p = inputMesh.vertex(e.from_vertex());
//        q = inputMesh.vertex(e.to_vertex());

        auto u = parameter(adjacentCell, e.from_vertex());
        auto v = parameter(adjacentCell, e.to_vertex());
        if (std::isnan(u[0]) ||  std::isnan(u[1]) ||  std::isnan(u[2]))
        {
            std::cerr << "hexex: NaN in parameterisation of vertex "
                << e.from_vertex().idx() << "." << std::endl;
            return; // or completely abort?
        }
        if (std::isnan(v[0]) ||  std::isnan(v[1]) ||  std::isnan(v[2]))
        {
            std::cerr << "hexex: NaN in parameterisation of vertex "
                << e.to_vertex().idx() << "." << std::endl;
            return; // or completely abort?
        }

        double start = -1;
        double end = -1;
        double dist = -1;

        for (int i = 0; i <= 2; ++i)
        if (fabs(u[i] - v[i]) > dist)
        {
            dist = fabs(u[i] - v[i]);
            start = u[i];
            end = v[i];
        }

        if (start > end)
        {
            std::swap(u,v);
            std::swap(p,q);
            std::swap(start,end);
        }

        if (start == end)
        {
            // skip degenerated edge
            continue;
        }

        for (int i = ceil(start); i <= floor(end); ++i)
        {
            double alpha = (i - start) / (end - start);
            if ((0 < alpha) && (alpha < 1))
            {
                auto w = (1-alpha) * u + alpha * v;
                auto roundW = roundVector(w);

//                if (isOnLine(u,v,roundW))
                if (isOnEdge(adjacentCell,*e_it,roundW))
                {
                    auto pos = (1-alpha) * p + alpha * q;
                    auto vh = intermediateHexMesh.add_vertex(pos);
                    vertexTypes[vh] = EHVertex;
                    incidentCellInInputMesh[vh] = adjacentCell;
                    hexvertexParameter[vh] = roundW;
                    incidentElementId[vh] = e_it->idx();
                    incidentVerticesPerEdge[*e_it].push_back(vh);
                }
            }

        }

    }

    HEXEX_DEBUG_ONLY(std::cout << "Extraction vertices on faces" << std::endl;)

    for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
    {
        HEXEX_DEBUG_ONLY(if ((f_it->idx() % 1000) == 0)
          std::cout << "Processing face " << f_it->idx() << " of " <<  inputMesh.n_faces() << std::endl;)


        auto hfh = inputMesh.halfface_handle(*f_it, 0);
        if (inputMesh.is_boundary(hfh))
            hfh = inputMesh.opposite_halfface_handle(hfh);
        auto adjacentCell = inputMesh.incident_cell(hfh);

        auto vertices = inputMesh.get_cell_vertices(hfh);

        auto p = inputPosition(vertices[0]);
        auto q = inputPosition(vertices[1]);
        auto r = inputPosition(vertices[2]);

//        p = inputMesh.vertex(vertices[0]);
//        q = inputMesh.vertex(vertices[1]);
//        r = inputMesh.vertex(vertices[2]);

        auto u = parameter(adjacentCell, vertices[0]);
        auto v = parameter(adjacentCell, vertices[1]);
        auto w = parameter(adjacentCell, vertices[2]);

        double left   = std::min(std::min(u[0], v[0]), w[0]);
        double right  = std::max(std::max(u[0], v[0]), w[0]);
        double bottom = std::min(std::min(u[1], v[1]), w[1]);
        double top    = std::max(std::max(u[1], v[1]), w[1]);
        double front  = std::min(std::min(u[2], v[2]), w[2]);
        double back   = std::max(std::max(u[2], v[2]), w[2]);

        auto invParametrization = getInverseParametrizationMatrix(p,q,r, u,v,w); // todo: only do this if it is needed

        for (int x = ceil(left); x <= floor(right); ++x)
            for (int y = ceil(bottom); y <= floor(top); ++y)
                for (int z = ceil(front); z <= floor(back); ++z) // TODO: this is unnecessary slow
                {
                    auto intPara = Parameter(x,y,z);
//                    if (isInside(u,v,w, intPara))
                    if (isInFace(hfh, intPara))
                    {
                        auto vh = intermediateHexMesh.add_vertex(invParametrization.transform_point(intPara));
                        vertexTypes[vh] = FHVertex;
                        incidentCellInInputMesh[vh] = adjacentCell;
                        hexvertexParameter[vh] = intPara;
                        incidentElementId[vh] = f_it->idx();
                        incidentVerticesPerFace[*f_it].push_back(vh);
                    }
                }
    }

    HEXEX_DEBUG_ONLY(std::cout << "Extraction vertices on cells" << std::endl;)

    for (auto c_it = inputMesh.cells_begin(); c_it != inputMesh.cells_end(); ++c_it)
    {
        HEXEX_DEBUG_ONLY(if ((c_it->idx() % 1000) == 0)
          std::cout << "Processing cell " << c_it->idx() << " of " <<  inputMesh.n_cells() << std::endl;)

        auto vertices = cellVertices[*c_it];

        auto p = inputPosition(vertices[0]);
        auto q = inputPosition(vertices[1]);
        auto r = inputPosition(vertices[2]);
        auto s = inputPosition(vertices[3]);

        auto u = parameter(*c_it, vertices[0]);
        auto v = parameter(*c_it, vertices[1]);
        auto w = parameter(*c_it, vertices[2]);
        auto t = parameter(*c_it, vertices[3]);

        double left   = std::min(std::min(std::min(u[0], v[0]), w[0]), t[0]);
        double right  = std::max(std::max(std::max(u[0], v[0]), w[0]), t[0]);
        double bottom = std::min(std::min(std::min(u[1], v[1]), w[1]), t[1]);
        double top    = std::max(std::max(std::max(u[1], v[1]), w[1]), t[1]);
        double back   = std::min(std::min(std::min(u[2], v[2]), w[2]), t[2]);
        double front  = std::max(std::max(std::max(u[2], v[2]), w[2]), t[2]);

        auto invParametrization = getInverseParametrizationMatrix(p,q,r,s, u,v,w,t);

        for (int x = ceil(left); x <= floor(right); ++x)
            for (int y = ceil(bottom); y <= floor(top); ++y)
                for (int z = ceil(back); z <= floor(front); ++z)
                {
                    auto intPara = Parameter(x,y,z);

//                    if ((!flipped && isInside(u,v,w,t,intPara)) ||
//                        ( flipped && isInside(u,v,t,w,intPara)))
                    if (isInCell(*c_it, intPara))
                    {
                        auto vh = intermediateHexMesh.add_vertex(invParametrization.transform_point(intPara));
                        vertexTypes[vh] = CHVertex;
                        incidentCellInInputMesh[vh] = *c_it;
                        hexvertexParameter[vh] = intPara;
                        incidentElementId[vh] = c_it->idx();
                        incidentVerticesPerCell[*c_it].push_back(vh);
                    }
                }

    }

    HEXEX_DEBUG_ONLY(std::cout << "finished cell extraction..." << std::endl;)
}

void HexExtractor::enumerateHPorts()
{
  for (auto v_it = intermediateHexMesh.vertices_begin(); v_it != intermediateHexMesh.vertices_end(); ++v_it)
    {
        hPortsOnVertex[*v_it].clear();
    }

    for (auto c_it = inputMesh.cells_begin(); c_it != inputMesh.cells_end(); ++c_it)
    {
        hPortsInCell[*c_it].clear();
    }

    HEXEX_DEBUG_ONLY(std::cout << "Enumerating ports" << std::endl;)

    for (auto v_it = intermediateHexMesh.vertices_begin(); v_it != intermediateHexMesh.vertices_end(); ++v_it)
    {
        HEXEX_DEBUG_ONLY(if ((v_it->idx() % 1000) == 0)
          std::cout << "processing vertex " << v_it->idx() << " of " << intermediateHexMesh.n_vertices() << std::endl;)

        if (vertexTypes[*v_it] == VHVertex)
            enumerateVertexHPorts(*v_it);
        else if (vertexTypes[*v_it] == EHVertex)
            enumerateEdgeHPorts(*v_it);
        else if (vertexTypes[*v_it] == FHVertex)
            enumerateFaceHPorts(*v_it);
        else if (vertexTypes[*v_it] == CHVertex)
            enumerateCellHPorts(*v_it);
        else
            assert(false);
    }
}

void HexExtractor::enumerateVertexHPorts(VertexHandle hexVh)
{
    assert(vertexTypes[hexVh] == VHVertex);

    VertexHandle ivh = VertexHandle(incidentElementId[hexVh]);
    assert(ivh.is_valid());

    for (auto vc_it = inputMesh.vc_iter(ivh); vc_it.valid(); ++vc_it)
    {
        auto incidentCell = *vc_it;

        auto param = parameter(incidentCell, ivh);

        // directions along edges
        for (auto heh : getIncidentHalfEdges(inputMesh, incidentCell, ivh))
        {
            for (auto dir : getAll6Directions())
            {
                if (pointsAlongHalfEdgeFromVertex(incidentCell, heh, dir))
                {
                    auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, EdgeHPort, inputMesh.edge_handle(heh).idx()));
                    hPortsOnVertex[hexVh].push_back(port);

                    hPortsInCell[incidentCell].push_back(port);
                }
            }
        }


        // directions into faces
        for (auto hfh : getIncidentHalfFaces(inputMesh, incidentCell, ivh))
        {
            if (isFaceDegenerate(hfh))
                continue;
            for (auto dir : getAll6Directions())
            {
                if (pointsIntoFace(hfh, ivh, dir))
                {
                    auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, FaceHPort, hfh.idx()));
                    hPortsOnVertex[hexVh].push_back(port);
                    hPortsInCell[incidentCell].push_back(port);
                }
            }
        }


        if (!isCellDegenerate(incidentCell))
        {
            // directions into cell
            for (auto dir : getAll6Directions())
                if (pointsIntoCell(incidentCell, ivh, dir))
                {
                    auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, CellHPort, incidentCell.idx()));
                    hPortsOnVertex[hexVh].push_back(port);
                    hPortsInCell[incidentCell].push_back(port);
                }
        }
    }

}

void HexExtractor::enumerateEdgeHPorts(VertexHandle hexVh)
{
    assert(vertexTypes[hexVh] == EHVertex);

    auto eh  = EdgeHandle(incidentElementId[hexVh]);
    auto heh = inputMesh.halfedge_handle(eh,0);

    assert(heh.is_valid());

    // directions pointing into cells
    for (auto hec_it = inputMesh.hec_iter(heh); hec_it.valid(); ++hec_it)
    {
        auto incidentCell = *hec_it;
        if (!incidentCell.is_valid())
            continue;


        auto param = getHexVertexParameter(hexVh, incidentCell);


        // directions pointing along the egde
        for (auto dir : getAll6Directions())
        {
            if (pointsAlongEdge(incidentCell, eh, dir, param))
            {
                auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, EdgeHPort, eh.idx()));
                hPortsOnVertex[hexVh].push_back(port);
                hPortsInCell[incidentCell].push_back(port);
            }
        }


        // directions pointing into faces
        for (auto hfh : getIncidentHalfFaces(inputMesh, incidentCell, eh))
        {
            for (auto dir : getAll6Directions())
            {
                if (pointsIntoFace(hfh, eh, dir, param))
                {
                    auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, FaceHPort, hfh.idx()));
                    hPortsOnVertex[hexVh].push_back(port);
                    hPortsInCell[incidentCell].push_back(port);
                }
            }
        }

        // directions pointing into cell

        if (!isCellDegenerate(incidentCell))
        {
            for (auto dir : getAll6Directions())
            {
                if (pointsIntoCell(incidentCell, eh, dir, param))
                {
                    auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, CellHPort, incidentCell.idx()));
                    hPortsOnVertex[hexVh].push_back(port);
                    hPortsInCell[incidentCell].push_back(port);
                }
            }
        }
    }

}

void HexExtractor::enumerateFaceHPorts(VertexHandle hexVh)
{
    assert(vertexTypes[hexVh] == FHVertex);

    auto fh = FaceHandle(incidentElementId[hexVh]);

    assert(fh.is_valid());

    for (auto hfh : {inputMesh.halfface_handle(fh, 0), inputMesh.halfface_handle(fh, 1)})
    {
        auto incidentCell = inputMesh.incident_cell(hfh);

        if (!incidentCell.is_valid())
            continue;


        auto param = getHexVertexParameter(hexVh, incidentCell);

        // directions along face

        for (auto dir : getAll6Directions())
        {
            if (pointsAlongFace(hfh, dir, param))
            {
                auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, FaceHPort, hfh.idx()));
                hPortsOnVertex[hexVh].push_back(port);
                hPortsInCell[incidentCell].push_back(port);
            }
        }

        // directions into cell

        if (!isCellDegenerate(incidentCell))
        {
            for (auto dir : getAll6Directions())
            {
                if (pointsIntoCell(hfh, dir, param))
                {
                    auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, CellHPort, incidentCell.idx()));
                    hPortsOnVertex[hexVh].push_back(port);
                    hPortsInCell[incidentCell].push_back(port);
                }
            }
        }
    }



}

void HexExtractor::enumerateCellHPorts(VertexHandle hexVh)
{
    assert(vertexTypes[hexVh] == CHVertex);

    auto incidentCell = incidentCellInInputMesh[hexVh];

    if (isCellDegenerate(incidentCell))
        return;

    auto param = getHexVertexParameter(hexVh, incidentCell);
    for (auto dir : getAll6Directions())
    {
        auto port = HPortHandle(std::make_shared<HPort>(dir, param, hexVh, incidentCell, CellHPort, incidentCell.idx()));
        hPortsOnVertex[hexVh].push_back(port);
        hPortsInCell[incidentCell].push_back(port);
    }
}



bool HexExtractor::isDartInCell(CellHandle ch, HPortHandle port, Direction refDir, Direction normalDir)
{

    if (port.type() == CellHPort)
        return true;

    if (port.type() == FaceHPort)
    {
        auto hfh = port.incidentHalfFaceHandle();

        auto vertices = inputMesh.get_cell_vertices(hfh);
        auto params = getParameters(port.cell(), vertices);

        auto isCellRegular = isRegular(params[0], params[1], params[2], params[3]);

        auto cond1 = (isCellRegular && isRegular(params[0], params[1], params[2], port.parameter() + refDir)) ||
                     (!isCellRegular && isFlipped(params[0], params[1], params[2], port.parameter() + refDir));

        if (cond1)
            return true;

        auto cond2 = isDegenerate(params[0], params[1], params[2], port.parameter() + refDir);
        auto cond3 = (isCellRegular && isRegular(params[0], params[1], params[2], port.parameter() + normalDir)) ||
                     (!isCellRegular && isFlipped(params[0], params[1], params[2], port.parameter() + normalDir));

        return cond2 && cond3;

    }

    if (port.type() == EdgeHPort)
    {
        auto eh = port.incidentEdgeHandle();

        auto heh = inputMesh.halfedge_handle(eh, 0);

        auto hfh = getIncidentHalfFaceInCell(inputMesh, ch, heh);

        auto vertices = inputMesh.get_cell_vertices(hfh, heh);
        auto params = getParameters(port.cell(), vertices);


        if ((params[1] | port.dir()) < (params[0] | port.dir()))
            params = getParameters(port.cell(), inputMesh.get_cell_vertices(getIncidentHalfFaceInCell(inputMesh, ch, inputMesh.opposite_halfedge_handle(heh)),inputMesh.opposite_halfedge_handle(heh)));

        if (isCellFlipped(ch))
            std::swap(params[2], params[3]);

        auto cond1 = isRegular(params[0], params[1], params[2], port.parameter() + refDir) && isRegular(params[0], params[3], params[1], port.parameter() + refDir);

        if (cond1) // plane strictly intersects tet
            return true;

        auto cond2a = isRegular(params[0], params[1], params[2], port.parameter() + refDir);// plane intersects on one boundary
        auto cond2b = isDegenerate(params[0], params[3], params[1], port.parameter() + refDir);
        auto cond2c = isRegular(params[0], params[3], params[1], port.parameter() + normalDir);

        auto cond2 = cond2a && cond2b && cond2c;

        if (cond2)
            return true;

        auto cond3 = isDegenerate(params[0], params[1], params[2], port.parameter() + refDir) && // plane intersects on other boundary
                     isRegular(params[0], params[3], params[1], port.parameter() + refDir) &&
                     isRegular(params[0], params[1], params[2], port.parameter() + normalDir);

        return cond3;
    }


    HEXEX_DEBUG_ONLY(std::cout << "Error: reached unreachable code" << std::endl;)
    assert(false);
    return false;

}

bool HexExtractor::isDartInCell(CellHandle ch, VertexHandle hexVh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    if (isCellDegenerate(ch))
        return false;

    auto& ports = hPortsInCell[ch];

    for (auto& p : ports)
        if ((p.vertex() == hexVh) && (p.parameter() == param) && p.dir() == traceDir)
            return isDartInCell(ch, p, refDir, normalDir);

    return false;
}




void HexExtractor::enumerateDarts(HPortHandle port, bool secondary)
{
    auto res = std::vector<Dart>();

    auto incidentCell = port.cell();

    if (isCellDegenerate(incidentCell))
        return;

    auto isFlipped = isCellFlipped(incidentCell);

    for (auto refDir : getAllOrthogonalDirections(port.dir(), false))
    {
        auto normalDir = port.dir() % refDir;
        if (isFlipped)
            normalDir.flip();

        if (secondary)
            normalDir.flip();

        if (isDartInCell(port.cell(), port, refDir, normalDir))
            res.push_back(Dart(*this, port, refDir, normalDir));

    }


    auto ds = &darts[port.vertex()];
    if (secondary)
        ds = &secondaryDarts[port.vertex()];

    for (auto& d : res)
        ds->push_back(std::shared_ptr<Dart>(new Dart(d)));

}


void HexExtractor::enumerateDarts()
{
    auto n = intermediateHexMesh.n_vertices();

    for (auto i = 0u; i < n; ++i)
    {
        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        darts[vh].clear();
    }


    HEXEX_DEBUG_ONLY(std::cout << "Enumerating darts" << std::endl;)

//#pragma omp parallel for
    for (auto i = 0u; i < n; ++i)
    {
        HEXEX_DEBUG_ONLY(if ((i % 1000) == 0)
          std::cout << "Processing vertex " << i << " of " << n << std::endl;)

        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;

        auto& ports = hPortsOnVertex[vh];

        for (auto& port : ports)
        {
            enumerateDarts(port);
        }
    }

    enumerateSecondaryDarts();

}

void HexExtractor::enumerateSecondaryDarts()
{
    auto n = intermediateHexMesh.n_vertices();


    HEXEX_DEBUG_ONLY(std::cout << "Enumerating more darts" << std::endl;)

    for (auto i = 0u; i < n; ++i)
    {

        HEXEX_DEBUG_ONLY(if ((i % 1000) == 0)
          std::cout << "Processing vertex " << i << " of " << n << std::endl;)

        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        secondaryDarts[vh].clear();
    }


//#pragma omp parallel for
    for (auto i = 0u; i < n; ++i)
    {

        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;

        auto& ports = hPortsOnVertex[vh];

        for (auto& port : ports)
        {
            enumerateDarts(port, true);
        }
    }

    HEXEX_DEBUG_ONLY(
      std::cout << "num darts: " << getNumDarts(false) << std::endl;
      std::cout << "size of 1 dart: " << sizeof(Dart) << " byte" << std::endl;
      std::cout << "1 size of HPortHandle: " << sizeof(HPortHandle) << " byte" << std::endl;
      std::cout << "2 size of Direction: " << sizeof(Direction) << " byte" << std::endl;
      std::cout << "4 size of Transition: " << sizeof(Transition) << " byte" << std::endl;
      std::cout << "1 size of Parameter: " << sizeof(Parameter) << " byte" << std::endl;
      std::cout << "1 size of HexExtractor&: " << sizeof(HexExtractor&) << " byte" << std::endl;
      std::cout << "4 size of Dart*: " << sizeof(Dart*) << " byte" << std::endl;
      std::cout << "0 size of bool: " << sizeof(bool) << " byte" << std::endl;
      std::cout << "2 size of HalfEdgeHandle: " << sizeof(HalfEdgeHandle) << " byte" << std::endl;
      std::cout << "size: " << getNumDarts(false) * sizeof(Dart) / 1e9 << " gb" << std::endl;
    )

}

int HexExtractor::getNumDarts(bool exclude_unconnected)
{
    int n_darts = 0;

    auto n = intermediateHexMesh.n_vertices();

    for (auto i = 0u; i < n; ++i)
    {
        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;

        auto& ds = darts[vh];
        auto& ds2 = secondaryDarts[vh];

        if (exclude_unconnected)
        {
            for (auto& d: ds)
                if (d->getAlpha<0>() != nullptr)
                    ++n_darts;
            for (auto& d: ds2)
                if (d->getAlpha<0>() != nullptr)
                    ++n_darts;
        }
        else
        {

            n_darts += ds.size();
            n_darts += ds2.size();
        }
    }

    return n_darts;
}


void HexExtractor::traceDarts()
{
    auto n = (int)intermediateHexMesh.n_vertices();

    HEXEX_DEBUG_ONLY(std::cout << "Tracing darts" << std::endl;)

//#pragma omp parallel for
    for (auto i = 0; i < n; ++i)
    {
        HEXEX_DEBUG_ONLY(if ((i % 1000) == 0)
          std::cout << "Processing vertex " << i << " of " << n << std::endl;)

        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            if (d->getAlpha<0>() == nullptr)
            {
              assert(d->isPrimary());
                traceDart2(*d);
            }
    }

}


bool HexExtractor::traceDart2(Dart& dart)
{
    auto port = dart.getTracePort();
    auto traceDir = port.dir();
    auto refDir = dart.getRefDir();
    auto normalDir = dart.getNormalDir();

    auto nextTraceDir = port.dir().flipped();
    auto nextRefDir = dart.getRefDir();
    auto nextNormalDir = dart.getNormalDir();

    auto currentCell = port.cell();

    if ((traceDir | refDir) != 0)
    {
        HEXEX_DEBUG_ONLY(std::cout << "Kaput" << std::endl;)
        assert(false);
        return false;
    }

    auto currentFlippedness = isCellFlipped(currentCell);

    auto startParam = port.parameter();
    auto endParam = startParam + traceDir;

    auto accumulatedTransitionFunction = identity;

    HEXEX_DEBUG_ONLY(
      std::set<HalfFaceHandle> visited;
    )

    auto i = 0u;

    ++numDartsTraced;

    auto lastFace = HalfFaceHandle();

    // rotate around trace dir from ref dir to normal dir
//    while (!isDartInCell(currentCell, hexVh, currentParameter, nextTraceDir, nextRefDir, nextNormalDir))
    while (getSecondaryDart(currentCell, endParam, nextTraceDir, nextRefDir, nextNormalDir).get() == nullptr)
    {
        HEXEX_DEBUG_ONLY(if (i > 9000)
            std::cout << "ohoh" << std::endl;)



        if (i++ > 10000)
        {
          HEXEX_DEBUG_ONLY(std::cout << "error in next connection" << std::endl;)
          return false;
        }

//        assert(currentPort.type() == EdgeHPort || currentPort.type() == FaceHPort);

        auto transitionFace = alpha0NextFace(lastFace, currentCell, startParam, traceDir, refDir, normalDir);

        if (!transitionFace.is_valid()) {
            HEXEX_DEBUG_ONLY(std::cout << "failed to find alpha0NextFace" << std::endl;)
            //assert(false);
            return false;
        }

        if (isFaceDegenerate(transitionFace))
        {
            HEXEX_DEBUG_ONLY(std::cout << "error in next connection" << std::endl;)
            //assert(false);
            return false;
        }

        if (inputMesh.is_boundary(inputMesh.opposite_halfface_handle(transitionFace)))
        {
            // reached boundary, leave dart without neighbor dart;
            return false;
        }

        doTransition(transitionFace, currentCell, startParam, endParam, traceDir, refDir, normalDir, nextTraceDir, nextRefDir, nextNormalDir, accumulatedTransitionFunction);


        if (currentFlippedness != isCellFlipped(currentCell))
        {
            // reverse search direction
            currentFlippedness = !currentFlippedness;
            std::swap(startParam, endParam);
            std::swap(traceDir, nextTraceDir);
        }

        lastFace = inputMesh.opposite_halfface_handle(transitionFace);
    }

    auto partnerSecondaryDart = getSecondaryDart(currentCell, endParam, nextTraceDir, nextRefDir, nextNormalDir);

    if (partnerSecondaryDart.get() != nullptr)
    {
        dart.connectAlpha<0>(partnerSecondaryDart.get());
//        assert(portFound);
        return true;
    }
    else
    {

        HEXEX_DEBUG_ONLY(std::cout << "could not find partner" << std::endl;)
        return false;

    }


}

bool HexExtractor::connectDartToPreviousSecondaryDart(Dart& dart)
{

    auto traceDir  = dart.getTraceDir();
    auto refDir    = dart.getRefDir();
    auto normalDir = dart.getNormalDir();

    auto nextTraceDir  = refDir;
    auto nextRefDir    = traceDir;
    auto nextNormalDir = normalDir;

    auto currentCell = dart.getCell();
    auto currentParameter = dart.getParameter();

    auto currentFlippedness = isCellFlipped(currentCell);

    auto transitionFunction = identity;

    int i = 0;

    auto lastFace = HalfFaceHandle();

    while (getSecondaryDart(currentCell, currentParameter, nextTraceDir, nextRefDir, nextNormalDir).get() == nullptr)
    {

        if (i++ > 10000)
        {
          HEXEX_DEBUG_ONLY(std::cout << "error in previous connection" << std::endl;)
          return false;
        }

        auto transitionFace = alpha1NextFace(lastFace, currentCell, currentParameter, traceDir, refDir, normalDir);

        if (!transitionFace.is_valid()) {
            connectDartToNeighborSecondaryDart2(dart);
            continue;
        }

        if (isFaceDegenerate(transitionFace))
        {
            assert(false);
            return false;
        }

        if (inputMesh.is_boundary(inputMesh.opposite_halfface_handle(transitionFace)))
        {
            // reached boundary, leave dart without neighbor dart;
            return false;
        }

        doTransition(transitionFace, currentCell, currentParameter, traceDir, refDir, normalDir, nextTraceDir, nextRefDir, nextNormalDir, transitionFunction);


        if (isCellFlipped(currentCell) != currentFlippedness)
        {
            std::swap(nextTraceDir, nextRefDir);
            currentFlippedness = !currentFlippedness;
        }

        lastFace = inputMesh.opposite_halfface_handle(transitionFace);
    }

    auto previousDart = getSecondaryDart(currentCell, currentParameter, nextTraceDir, nextRefDir, nextNormalDir);

    dart.connectAlpha<1>(previousDart.get());
    return true;
}

void HexExtractor::connectDartsToPreviousSecondaryDart()
{
    auto n = (int)intermediateHexMesh.n_vertices();

    HEXEX_DEBUG_ONLY(std::cout << "connecting darts 1 " << std::endl;)

//#pragma omp parallel for
    for (auto i = 0; i < n; ++i)
    {
        HEXEX_DEBUG_ONLY(if ((i % 1000) == 0)
          std::cout << "Processing vertex " << i << " of " << n << std::endl;)

        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            connectDartToPreviousSecondaryDart(*d);
    }

}



bool HexExtractor::connectDartToNeighborSecondaryDart2(Dart& dart)
{

    assert(!isCellDegenerate(dart.getCell()));

    auto traceDir  = dart.getTraceDir();
    auto refDir    = dart.getRefDir();
    auto normalDir = dart.getNormalDir();

    auto nextTraceDir  = traceDir;
    auto nextRefDir    = normalDir;
    auto nextNormalDir = refDir;

    auto currentCell = dart.getCell();
    auto currentParameter = dart.getParameter();

    auto currentFlippedness = isCellFlipped(currentCell);

    auto transitionFunction = identity;

    int i = 0;

    auto lastFace = HalfFaceHandle();

    while (getSecondaryDart(currentCell, currentParameter, nextTraceDir, nextRefDir, nextNormalDir).get() == nullptr)
    {

        if (i++ > 10000)
        {
          HEXEX_DEBUG_ONLY(std::cout << "error in neighbor connection" << std::endl;)
          return false;
        }

        auto transitionFace = alpha2NextFace(lastFace, currentCell, currentParameter, traceDir, refDir, normalDir);

        if (!transitionFace.is_valid())
        {
            HEXEX_DEBUG_ONLY(std::cout << "failed to find alpha2NextFace" << std::endl;)
            //assert(false);
            return false;
        }

        if (isFaceDegenerate(transitionFace))
        {
            //assert(false);
            return false;
        }

        if (inputMesh.is_boundary(inputMesh.opposite_halfface_handle(transitionFace)))
        {
            // reached boundary, leave dart without neighbor dart;
            return false;
        }

        doTransition(transitionFace, currentCell, currentParameter, traceDir, refDir, normalDir, nextTraceDir, nextRefDir, nextNormalDir, transitionFunction);

        if (isCellFlipped(currentCell) != currentFlippedness)
        {
            std::swap(refDir, normalDir);
            std::swap(nextRefDir, nextNormalDir);
            currentFlippedness = !currentFlippedness;
        }

        lastFace = inputMesh.opposite_halfface_handle(transitionFace);
    }



    auto neighborDart = getSecondaryDart(currentCell, currentParameter, nextTraceDir, nextRefDir, nextNormalDir);

    dart.connectAlpha<2>(neighborDart.get());
    return true;

}

void HexExtractor::connectDartsToNeighborSecondaryDart()
{
    auto n = (int)intermediateHexMesh.n_vertices();

    HEXEX_DEBUG_ONLY(std::cout << "connecting darts 3 " << std::endl;)

//#pragma omp parallel for
    for (auto i = 0; i < n; ++i)
    {

        HEXEX_DEBUG_ONLY(if ((i % 1000) == 0)
          std::cout << "Processing vertex " << i << " of " << n << std::endl;)

        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            connectDartToNeighborSecondaryDart2(*d);
    }

}

bool HexExtractor::connectDartToOppositeSecondaryDart(Dart& dart)
{
    auto traceDir  = dart.getTraceDir();
    auto refDir    = dart.getRefDir();
    auto normalDir = dart.getNormalDir();

    auto nextTraceDir  = traceDir;
    auto nextRefDir    = refDir;
    auto nextNormalDir = normalDir.flipped();

    auto currentCell = dart.getCell();
    auto currentParameter = dart.getParameter();

    auto currentFlippedness = isCellFlipped(currentCell);

    auto transitionFunction = identity;


    auto lastFace = HalfFaceHandle();

    int i = 0;
    int j = 0;

    while (getSecondaryDart(currentCell, currentParameter, nextTraceDir, nextRefDir, nextNormalDir).get() == nullptr)
    {

        j++;

        if (!isCellDegenerate(currentCell))
            ++i;

        if (i > 100)
        {
          HEXEX_DEBUG_ONLY(std::cout << "error in opposite dart connection" << std::endl;)
          // leave dart dangling. This should however be fixed, if this still occurs
          return false;
        }

        auto transitionFace = alpha3NextFace(lastFace, currentCell, currentParameter, traceDir, refDir, normalDir);
        if (!transitionFace.is_valid()) {
            std::cerr << "tf invalid" <<lastFace<<" "<<currentCell<<" "<<currentParameter<< std::endl;
            return false;
        }

        if (inputMesh.is_boundary(inputMesh.opposite_halfface_handle(transitionFace)))
        {
            // reached boundary, leave dart without opposite dart;
            return false;
        }

        doTransition(transitionFace, currentCell, currentParameter, traceDir, refDir, normalDir, nextTraceDir, nextRefDir, nextNormalDir, transitionFunction);

        if (isCellFlipped(currentCell) != currentFlippedness)
        {
//            std::swap(nextRefDir, nextNormalDir);
            nextNormalDir.flip();
            currentFlippedness = !currentFlippedness;
        }

        lastFace = inputMesh.opposite_halfface_handle(transitionFace);

    }


    //assert(i < 2);


    auto oppositeDart = getSecondaryDart(currentCell, currentParameter, nextTraceDir, nextRefDir, nextNormalDir);

    if (oppositeDart.get() == nullptr)
        connectDartToOppositeSecondaryDart(dart);

    dart.connectAlpha<3>(oppositeDart.get());
    return true;
}

void HexExtractor::connectDartsToOppositeSecondaryDart()
{
    auto n = (int)intermediateHexMesh.n_vertices();

    HEXEX_DEBUG_ONLY(std::cout << "connecting darts 4 " << std::endl;)

//#pragma omp parallel for
    for (auto i = 0; i < n; ++i)
    {
        HEXEX_DEBUG_ONLY(if ((i % 1000) == 0)
          std::cout << "Processing vertex " << i << " of " << n << std::endl;)

        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            connectDartToOppositeSecondaryDart(*d);
    }

    equivalenceClasses.clear();
    equivalenceClasses.resize(intermediateHexMesh.n_vertices());


    for (auto v_it = intermediateHexMesh.vertices_begin(); v_it != intermediateHexMesh.vertices_end(); ++v_it)
    {
        equivalenceClassIds[*v_it] = v_it->idx();
        equivalenceClasses[v_it->idx()] = {v_it->idx()};
    }

}

bool HexExtractor::alpha0FaceTest(HalfFaceHandle hfh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    auto params = getParameters(hfh);

    // start and end of grid edge must be on opposite sides of the triangles
    auto b1 = (isRegular(params[0], params[1], params[2], param) && isRegular(params[0], params[2], params[1], param + traceDir)) ||
              (isFlipped(params[0], params[1], params[2], param) && isFlipped(params[0], params[2], params[1], param + traceDir));

    if (!b1)
        return false;

    std::array<ORIENTATION, 3> oris;
    oris[0] = sign_orient3d(params[0].data(), params[1].data(), param.data(), (param + traceDir).data());
    oris[1] = sign_orient3d(params[1].data(), params[2].data(), param.data(), (param + traceDir).data());
    oris[2] = sign_orient3d(params[2].data(), params[0].data(), param.data(), (param + traceDir).data());

    if ((oris[0] == oris[1]) && (oris[1] == oris[2]))
    {
        // edge passes through center of triangle
        assert(oris[0] != ORI_ZERO);
        return true;
    }

    bool flipped = std::any_of(oris.begin(), oris.end(), [](ORIENTATION ori){ return ori == ORI_NEGATIVE; });
    bool proper  = std::any_of(oris.begin(), oris.end(), [](ORIENTATION ori){ return ori == ORI_POSITIVE; });

    if (proper && flipped)
        // no intersection
        return false;


    std::vector<int> edgeIntersections;
    for (unsigned int i = 0; i < oris.size(); ++i)
        if (oris[i] == ORI_ZERO)
            edgeIntersections.push_back(i);


    if (edgeIntersections.size() == 1)
    {
        // edge intersection

        // find upper and lower vertex
        auto upperVertex = edgeIntersections[0];
        auto lowerVertex = (edgeIntersections[0]+1)%3;
        auto otherVertex = (edgeIntersections[0]+2)%3;

        if ((params[upperVertex] | normalDir) == (param | normalDir))
        {
            // intersected edge lies in grid face
            assert((params[lowerVertex] | normalDir) == (param | normalDir));

            // check if third point point into normal direction
            return (params[otherVertex] | normalDir) > (param | normalDir);
        }
        else
        {
            if ((params[lowerVertex] | normalDir) > (param | normalDir))
                std::swap(lowerVertex, upperVertex);
            assert((params[lowerVertex] | normalDir) < (param | normalDir));

            // check if triangle is on correct side
            auto p1_2d = Vec2d(params[upperVertex] | refDir, params[upperVertex] | normalDir);
            auto p2_2d = Vec2d(params[lowerVertex] | refDir, params[lowerVertex] | normalDir);
            auto p3_2d = Vec2d(params[otherVertex] | refDir, params[otherVertex] | normalDir);

            auto b3 = isRegular(p1_2d, p2_2d, p3_2d);

            return b3;

        }


    }
    else if (edgeIntersections.size() == 2)
    {
        // vertex intersection
        auto intersectedVertex = 0;
        if (edgeIntersections[0] == 1)
            intersectedVertex = 2;
        else if (edgeIntersections[1] == 1)
            intersectedVertex = 1;

        auto p1 = params[intersectedVertex];
        auto p2 = params[(intersectedVertex+1)%3];
        auto p3 = params[(intersectedVertex+2)%3];

        // check if one point of the triangle is above face and the other below or on the plane of the grid face
        auto b1 = (((p2 | normalDir) > (param | normalDir)) && ((p3 | normalDir) <= (param | normalDir))); // p2 is upper, p3 lower
        auto b2 = (((p3 | normalDir) > (param | normalDir)) && ((p2 | normalDir) <= (param | normalDir))); // p3 is upper, p2 lower


        if (!b1 && !b2)
            return false; // triangle does not intersect cell

        // check if triangle is on correct side
        auto p1_2d = Vec2d(p1 | refDir, p1 | normalDir);
        auto p2_2d = Vec2d(p2 | refDir, p2 | normalDir);
        auto p3_2d = Vec2d(p3 | refDir, p3 | normalDir);

        if (b2)
            std::swap(p2_2d, p3_2d); // now p2 is upper

        auto b3 = isRegular(p1_2d, p3_2d, p2_2d);

        return b3;
    }
    else
    {
        //assert(false);
        return false;
    }




}

HalfFaceHandle HexExtractor::alpha0NextFace(HalfFaceHandle prevFace, CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    for (auto hfh : inputMesh.cell(ch).halffaces())
        if (hfh != prevFace)
            if (!isFaceDegenerate(hfh))
                if (alpha0FaceTest(hfh, param, traceDir, refDir, normalDir))
                    return hfh;

    //assert(false);

    return HalfFaceHandle();
}

bool HexExtractor::alpha1FaceTest(HalfFaceHandle hfh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    auto params = getParameters(hfh);

    if (isInFace(hfh, param))
        return true;

    if (isInFaceRelaxed(hfh, param)) // todo: this is unnecessarily expensive
    {
        // find intersected edge
        int edgeIntersectionId = -1;
        for (unsigned int i = 0; i < 3; ++i)
            if (isOnLineBetween(params[i], params[(i+1)%3], param))
                edgeIntersectionId = i;

        auto tris = std::vector<std::tuple<Parameter, Parameter, Parameter>>(); // stores triangles to be tested, (2 in the case of vertex lying on edge

        if (edgeIntersectionId != -1)
        {

            // check if gridface intersects face

            // pretend we the triangle consists of two triangles, i.e. the original triangle split at the vertex. Perform check for both triangles
            auto tri1 = std::make_tuple(param, params[(edgeIntersectionId+1)%3], params[(edgeIntersectionId+2)%3]);
            auto tri2 = std::make_tuple(param, params[(edgeIntersectionId+2)%3], params[edgeIntersectionId]);

            tris.push_back(tri1);
            tris.push_back(tri2);
        }
        else
        {
            int vertexIntersectionId = -1;
            for (unsigned int i = 0; i < 3; ++i)
                if (params[i] == param)
                    vertexIntersectionId = i;

            if (vertexIntersectionId == -1)
            {
                assert(false);
                return false;
            }

            auto tri1 = std::make_tuple(param, params[(vertexIntersectionId+1)%3], params[(vertexIntersectionId+2)%3]);

            tris.push_back(tri1);

        }


        for (auto tri : tris)
        {
            auto tp1 = std::get<0>(tri);
            auto tp2 = std::get<1>(tri);
            auto tp3 = std::get<2>(tri);

            // check that the two other corners of the grid face are on opposite sides of the triangle
            auto p1 = param + traceDir;
            auto p2 = param + refDir;

            auto b1 = (isRegular(tp1,tp2,tp3,p1) && isRegular(tp1,tp3,tp2,p2));
            auto b2 = (isFlipped(tp1,tp2,tp3,p1) && isFlipped(tp1,tp3,tp2,p2));

            if (!b1 && !b2)
            {
                // face planes do not intersect. No need to test other half triangle
                return false;
            }

            // check that one point is above and the other below or in the plane of the grid face

            auto b3 = (((tp2 | normalDir) > (param | normalDir)) && ((tp3 | normalDir) <= (param | normalDir))); // tp2 is upper, tp3 lower
            auto b4 = (((tp3 | normalDir) > (param | normalDir)) && ((tp2 | normalDir) <= (param | normalDir))); // tp3 is upper, tp2 lower

            // check in 2d that traceDir and refDir point into projected triangle
            auto tp1_2d = Vec2d(tp1 | refDir, tp1 | normalDir);
            auto tp2_2d = Vec2d(tp2 | refDir, tp2 | normalDir);
            auto tp3_2d = Vec2d(tp3 | refDir, tp3 | normalDir);

            auto b5 = (b3 && isRegular(tp1_2d, tp3_2d, tp2_2d)) || (b4 && isRegular(tp1_2d, tp2_2d, tp3_2d));

            if (b5)
                return true;
        }

    }

    return false;

}

HalfFaceHandle HexExtractor::alpha1NextFace(HalfFaceHandle prevFace, CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    for (auto hfh : inputMesh.cell(ch).halffaces())
        if (hfh != prevFace)
            if (!isFaceDegenerate(hfh))
                if (alpha1FaceTest(hfh, param, traceDir, refDir, normalDir))
                    return hfh;

    //assert(false);

    return HalfFaceHandle();
}

bool HexExtractor::alpha2FaceTest(HalfFaceHandle hfh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    auto params = getParameters(hfh);

    if (isInFace(hfh, param))
        return true;

    if (isInFaceRelaxed(hfh, param)) // todo: this is unnecessarily expensive
    {
        // find intersected edge
        int edgeIntersectionId = -1;
        for (unsigned int i = 0; i < 3; ++i)
            if (isOnLineBetween(params[i], params[(i+1)%3], param))
                edgeIntersectionId = i;

        if (edgeIntersectionId != -1)
        {

            // check if trace dir points into face or if it points along boundary or neither
            if (pointsIntoProjectedTriangleFromEdge(params[edgeIntersectionId],params[(edgeIntersectionId+1)%3],params[(edgeIntersectionId+2)%3], param, traceDir))
            {
                return true;
            }
            else if (isOnLine(params[edgeIntersectionId], params[(edgeIntersectionId+1)%3], param + traceDir))
            {
                // trace dir points along edge, test if cell intersects triangle
                auto p = params[(edgeIntersectionId+2)%3];

                return ((p | refDir) > (param | refDir)) && ((p | normalDir) > (param | normalDir));
            }
            else
                return false;

        }
        else
        {
            int vertexIntersectionId = -1;
            for (unsigned int i = 0; i < 3; ++i)
                if (params[i] == param)
                    vertexIntersectionId = i;

            assert(vertexIntersectionId != -1);

            // check if trace dir points into face or if it points along boundary or neither
            if (pointsIntoTriangle(params[vertexIntersectionId],params[(vertexIntersectionId+1)%3],params[(vertexIntersectionId+2)%3], traceDir))
            {
                return true;
            }
            else
            {
                int edgeIntersectionId2 = -1;
                int otherVertex = -1;
                if (isOnLine(params[vertexIntersectionId], params[(vertexIntersectionId+1)%3], param + traceDir) &&
                   ((params[(vertexIntersectionId+1)%3]|traceDir) > (param|traceDir)))
                {
                    edgeIntersectionId2 = (vertexIntersectionId+1)%3;
                    otherVertex = (vertexIntersectionId+2)%3;
                }
                if (isOnLine(params[vertexIntersectionId], params[(vertexIntersectionId+2)%3], param + traceDir) &&
                   ((params[(vertexIntersectionId+2)%3]|traceDir) > (param|traceDir)))
                {
                    edgeIntersectionId2 = (vertexIntersectionId+2)%3;
                    otherVertex = (vertexIntersectionId+1)%3;
                }

                if (edgeIntersectionId2 != -1)
                {

                    // trace dir points along edge, test if cell intersects triangle
                    auto p = params[otherVertex];

                    return ((p | refDir) > (param | refDir)) && ((p | normalDir) > (param | normalDir));

                }
                else
                {
                    // trace dir points neither directly into the face nor along the boundary
                    return false;
                }
            }


        }

    }

    return false;

}

HalfFaceHandle HexExtractor::alpha2NextFace(HalfFaceHandle prevFace, CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    for (auto hfh : inputMesh.cell(ch).halffaces())
        if (hfh != prevFace)
            if (!isFaceDegenerate(hfh))
                if (alpha2FaceTest(hfh, param, traceDir, refDir, normalDir))
                    return hfh;

//    alpha2NextFace(prevFace, ch, param, traceDir, refDir, normalDir);

//    assert(false);
    return HalfFaceHandle();
}

bool HexExtractor::alpha3FaceTest(HalfFaceHandle hfh, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{

    auto params = getParameters(hfh);

    if (isInFace(hfh, param))
        return true;

    if (isInFaceRelaxed(hfh, param)) // todo: this is unnecessarily expensive
    {
        // find intersected edge
        int edgeIntersectionId = -1;
        for (unsigned int i = 0; i < 3; ++i)
            if (isOnLineBetween(params[i], params[(i+1)%3], param))
                edgeIntersectionId = i;

        if (edgeIntersectionId != -1)
        {

            // check if trace dir points into face or if it points along boundary or neither
            if (pointsIntoProjectedTriangleFromEdge(params[edgeIntersectionId],params[(edgeIntersectionId+1)%3],params[(edgeIntersectionId+2)%3], param, traceDir))
            {
                return true;
            }
            else if (isOnLine(params[edgeIntersectionId], params[(edgeIntersectionId+1)%3], param + traceDir))
            {
                // trace dir points along edge, test if face intersects triangle
                auto p = params[(edgeIntersectionId+2)%3];

                // check if triangle and face lie in same plane

                auto faceP = param + refDir;

                if ((faceP | normalDir) == (p | normalDir))
                {
                    // check if triangle lies on correct side
                    if ((p | refDir) > (param | refDir))
                        return true;
                    else
                        return false;
                }
                else
                    return false;
            }
            else
                return false;

        }
        else
        {
            int vertexIntersectionId = -1;
            for (unsigned int i = 0; i < 3; ++i)
                if (params[i] == param)
                    vertexIntersectionId = i;

            assert(vertexIntersectionId != -1);

            // check if trace dir points into face or if it points along boundary or neither
            if (pointsIntoTriangle(params[vertexIntersectionId],params[(vertexIntersectionId+1)%3],params[(vertexIntersectionId+2)%3], traceDir))
            {
                return true;
            }
            else
            {
                int edgeIntersectionId2 = -1;
                int otherVertex = -1;
                if (isOnLine(params[vertexIntersectionId], params[(vertexIntersectionId+1)%3], param + traceDir) &&
                   ((params[(vertexIntersectionId+1)%3]|traceDir) > (param|traceDir)))
                {
                    edgeIntersectionId2 = (vertexIntersectionId+1)%3;
                    otherVertex = (vertexIntersectionId+2)%3;
                }
                if (isOnLine(params[vertexIntersectionId], params[(vertexIntersectionId+2)%3], param + traceDir) &&
                   ((params[(vertexIntersectionId+2)%3]|traceDir) > (param|traceDir)))
                {
                    edgeIntersectionId2 = (vertexIntersectionId+2)%3;
                    otherVertex = (vertexIntersectionId+1)%3;
                }

                if (edgeIntersectionId2 != -1)
                {

                    // trace dir points from vertex along edge, test if face intersects triangle
                    auto p = params[otherVertex];

                    // check if triangle and face lie in same plane

                    auto faceP = param + refDir;

                    if ((faceP | normalDir) == (p | normalDir))
                    {
                        // check if triangle lies on correct side
                        if ((p | refDir) > (param | refDir))
                            return true;
                        else
                            return false;
                    }
                    else
                        return false;

                }
                else
                {
                    // trace dir points neither directly into the face nor along the boundary
                    return false;
                }
            }


        }

    }

    return false;
}

HalfFaceHandle HexExtractor::alpha3NextFace(HalfFaceHandle prevFace, CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{

    for (auto hfh : inputMesh.cell(ch).halffaces())
        if (hfh != prevFace)
            if (!isFaceDegenerate(hfh))
                if (alpha3FaceTest(hfh, param, traceDir, refDir, normalDir))
                    return hfh;

//    alpha2NextFace(prevFace, ch, param, traceDir, refDir, normalDir);

    //assert(false);
    return HalfFaceHandle();
}

std::shared_ptr<Dart> HexExtractor::getDart(CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    int i = 0;

    for (auto& p : hPortsInCell[ch])
        if (p.parameter() == param && p.dir() == traceDir)
            for (auto& d : darts[p.vertex()])
                if ((d->getTraceDir() == traceDir) && (d->getRefDir() == refDir) && (d->getNormalDir() == normalDir) && (d->getTracePort().cell() == ch))
                    ++i;

    HEXEX_DEBUG_ONLY(if (i > 1)
        std::cout << "Bad :(" << std::endl;)



    for (auto& p : hPortsInCell[ch])
        if (p.parameter() == param && p.dir() == traceDir)
            for (auto& d : darts[p.vertex()])
                if ((d->getTraceDir() == traceDir) && (d->getRefDir() == refDir) && (d->getNormalDir() == normalDir) && (d->getTracePort().cell() == ch))
                    return d;



    assert(false);
    return std::shared_ptr<Dart>(nullptr);
}

std::shared_ptr<Dart> HexExtractor::getSecondaryDart(CellHandle ch, Parameter param, Direction traceDir, Direction refDir, Direction normalDir)
{
    int i = 0;

    for (auto& p : hPortsInCell[ch])
        if (p.parameter() == param && p.dir() == traceDir)
            for (auto& d : secondaryDarts[p.vertex()])
                if ((d->getTraceDir() == traceDir) && (d->getRefDir() == refDir) && (d->getNormalDir() == normalDir) && (d->getTracePort().cell() == ch))
                    ++i;

    HEXEX_DEBUG_ONLY(if (i > 1)
        std::cout << "Bad :(" << std::endl;)



    for (auto& p : hPortsInCell[ch])
        if (p.parameter() == param && p.dir() == traceDir)
            for (auto& d : secondaryDarts[p.vertex()])
                if ((d->getTraceDir() == traceDir) && (d->getRefDir() == refDir) && (d->getNormalDir() == normalDir) && (d->getTracePort().cell() == ch))
                    return d;



//    assert(false);
    return std::shared_ptr<Dart>(nullptr);
}


HalfEdgeHandle HexExtractor::addEdge(HPortHandle p1, HPortHandle p2)
{

    if (!p1.isValid() || !p2.isValid())
        return HalfEdgeHandle();

    // check if that edge already exists for the current ports
    auto hehCandidates = getHalfedges(intermediateHexMesh, p1.vertex(), p2.vertex());

    for (auto heh : hehCandidates)
    {
        if (incidentHPort[heh] == p1 && incidentHPort[intermediateHexMesh.opposite_halfedge_handle(heh)] == p2)
        {
            p1.halfedge() = heh;
            p2.halfedge() = intermediateHexMesh.opposite_halfedge_handle(heh);

            return heh;
        }
    }

    // if the edge does not exist, add it

    auto eh = intermediateHexMesh.add_edge(p1.vertex(), p2.vertex(), true);
    p1.halfedge() = intermediateHexMesh.halfedge_handle(eh,0);
    p2.halfedge() = intermediateHexMesh.halfedge_handle(eh, 1);
    incidentHPort[intermediateHexMesh.halfedge_handle(eh,0)] = p1;
    incidentHPort[intermediateHexMesh.halfedge_handle(eh,1)] = p2;

    return intermediateHexMesh.halfedge_handle(eh, 0);

}

void HexExtractor::mergeEquivalenceClasses(Dart* dart)
{

    if (dart->getAlpha<0>() != nullptr)
        if (dart->isAnti() != dart->getAlpha<0>()->isAnti())
        {
            joinEquivalenceClasses(dart->getVertex(), dart->getAlpha<0>()->getVertex());
        }

    if (dart->getAlpha<1>() != nullptr)
        joinEquivalenceClasses(dart->getVertex(), dart->getAlpha<1>()->getVertex());


    if (dart->getAlpha<2>() != nullptr)
        joinEquivalenceClasses(dart->getVertex(), dart->getAlpha<2>()->getVertex());

    if (dart->getAlpha<3>() != nullptr)
        joinEquivalenceClasses(dart->getVertex(), dart->getAlpha<3>()->getVertex());

}

void HexExtractor::mergeEquivalenceClassesOfAllDarts()
{
    auto n = (int)intermediateHexMesh.n_vertices();
    for (auto i = 0; i < n; ++i)
    {
        HEXEX_DEBUG_ONLY(if ((i % 1000) == 0)
          std::cout << "Processing vertex " << i << " of " << n << std::endl;)

        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            mergeEquivalenceClasses(d.get());
        auto& ds2 = secondaryDarts[vh];
        for (auto& d : ds2)
            mergeEquivalenceClasses(d.get());
    }
}

void HexExtractor::mergeVerticesOfOppositeFaces(HalfFaceHandle hfh)
{
    auto& ds = halffaceDarts[hfh];

    if (ds.empty())
        return;

    if (ds.size() != 4)
    {
        HEXEX_DEBUG_ONLY(std::cout << "not merging vertices of opposite faces because number of darts is not 4" << std::endl;)
        return;
    }


    auto getHalffaceDarts = [&](Dart* dart, bool backwards)
    {
        auto res = std::vector<Dart*>();
        auto d = dart;
        do
        {
            res.push_back(d);
            assert(d->getAlpha<0>() != nullptr);
            assert(d->getAlpha<0>()->getAlpha<1>() != nullptr);

            if (backwards)
                d = d->getAlpha<1>()->getAlpha<0>();
            else
                d = d->getAlpha<0>()->getAlpha<1>();
        } while (d != dart);
        return res;
    };

    for (auto d : ds)
    {
        auto oppositeDart = d->getAlpha<3>();
        if (oppositeDart != nullptr)
        {
            if (oppositeDart->getAlpha<1>())
            {
                auto ds1 = getHalffaceDarts(d, false);
                auto ds2 = getHalffaceDarts(oppositeDart->getAlpha<1>(), true);

                if ((ds1.size() != 4) || (ds2.size() != 4))
                {
                    HEXEX_DEBUG_ONLY(std::cout << "not merging vertices of opposite faces because number of darts is not 4" << std::endl;)
                    continue;
                }

                for (auto i = 0u; i < 4; ++i)
                {
                    joinEquivalenceClasses(ds1[i]->getVertex(), ds2[i]->getVertex());
                }
            }

        }
    }
}

void HexExtractor::mergeVerticesOfOppositeFaces()
{

    for (auto hfh : intermediateHexMesh.halffaces())
        mergeVerticesOfOppositeFaces(hfh);

}

bool HexExtractor::fixProblem(Dart* dart, Dart* otherDart, bool mergeVertices)
{
    if (mergeVertices)
    {
        mergeEquivalenceClasses(dart);
        mergeEquivalenceClasses(otherDart);
    }

    assert(dart != nullptr);
    assert(otherDart != nullptr);
    assert(dart->isAnti() != otherDart->isAnti());

    auto partnerDart1 = dart->getAlpha<0>();
    auto previousDart1 = dart->getAlpha<1>();
    auto neighborDart1 = dart->getAlpha<2>();
    auto oppositeDart1 = dart->getAlpha<3>();

    auto partnerDart2 = otherDart->getAlpha<0>();
    auto previousDart2 = otherDart->getAlpha<1>();
    auto neighborDart2 = otherDart->getAlpha<2>();
    auto oppositeDart2 = otherDart->getAlpha<3>();

    reconnectOrDisconnect<0>(dart, otherDart, partnerDart1, partnerDart2);
    reconnectOrDisconnect<1>(dart, otherDart, previousDart1, previousDart2);
    reconnectOrDisconnect<2>(dart, otherDart, neighborDart1, neighborDart2);
    reconnectOrDisconnect<3>(dart, otherDart, oppositeDart1, oppositeDart2);

    dart->disconnectAlpha<0>();
    dart->disconnectAlpha<1>();
    dart->disconnectAlpha<2>();
    dart->disconnectAlpha<3>();

    otherDart->disconnectAlpha<0>();
    otherDart->disconnectAlpha<1>();
    otherDart->disconnectAlpha<2>();
    otherDart->disconnectAlpha<3>();

    return true;
}

bool HexExtractor::fixProblem1(Dart* dart, bool mergeVertices)
{
    assert(dart->getAlpha<0>() != nullptr);
    assert(dart->isAnti() != dart->getAlpha<0>()->isAnti());

    fixProblem(dart, dart->getAlpha<0>(), mergeVertices);

    return true;
}

bool HexExtractor::fixProblems1(bool mergeVertices)
{

    auto n = (int)intermediateHexMesh.n_vertices();

    bool changes = false;

//OMP_PARALLEL_FOR
    for (auto i = 0; i < n; ++i)
    {
        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            if (d->getAlpha<0>() != nullptr)
                if (d->isAnti() != d->getAlpha<0>()->isAnti())
                {
                    fixProblem1(d.get(), mergeVertices);
                    changes = true;
                }
    }

    return changes;
}

bool HexExtractor::fixProblem2(Dart* dart, bool mergeVertices)
{
    assert(dart->getAlpha<1>() != nullptr);
    assert(dart->isAnti() != dart->getAlpha<1>()->isAnti());

    fixProblem(dart, dart->getAlpha<1>(), mergeVertices);

    return true;
}

bool HexExtractor::fixProblems2(bool mergeVertices)
{
    auto n = (int)intermediateHexMesh.n_vertices();

    bool changes = false;

//OMP_PARALLEL_FOR
    for (auto i = 0; i < n; ++i)
    {
        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            if (d->getAlpha<1>() != nullptr)
                if (d->isAnti() != d->getAlpha<1>()->isAnti())
                {
//                    if ((d->getAlpha<0>() != nullptr) && (d->isAnti() != d->getAlpha<0>()->isAnti()))
//                        joinEquivalenceClasses(d->getVertex(), d->getAlpha<0>()->getVertex());
//                    if ((d->getAlpha<1>()->getAlpha<0>() != nullptr) && (d->isAnti() != d->getAlpha<1>()->getAlpha<0>()->isAnti()))
//                        joinEquivalenceClasses(d->getVertex(), d->getAlpha<1>()->getAlpha<0>()->getVertex());

//                    joinEquivalenceClasses(d->getVertex(), d->getAlpha<1>()->getVertex());
                    fixProblem2(d.get(), mergeVertices);
                    changes = true;
                }
    }

    return changes;
}

bool HexExtractor::fixProblem3(Dart* dart, bool mergeVertices)
{
    assert(dart->getAlpha<2>() != nullptr);
    assert(dart->isAnti() != dart->getAlpha<2>()->isAnti());

    fixProblem(dart, dart->getAlpha<2>(), mergeVertices);

    return true;
}

bool HexExtractor::fixProblems3(bool mergeVertices)
{
    auto n = (int)intermediateHexMesh.n_vertices();

    bool changes = false;

//OMP_PARALLEL_FOR
    for (auto i = 0; i < n; ++i)
    {
        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            if (d->getAlpha<2>() != nullptr)
                if (d->isAnti() != d->getAlpha<2>()->isAnti())
                {
//                    joinEquivalenceClasses(d->getVertex(), d->getAlpha<2>()->getVertex());
//                    joinEquivalenceClasses(d->getVertex(), d->getAlpha<2>()->getAlpha<1>()->getVertex());
                    fixProblem3(d.get(), mergeVertices);
                    changes = true;
                }
    }

    return changes;

}

bool HexExtractor::fixProblem4(Dart* dart, bool mergeVertices)
{
    assert(dart->getAlpha<3>() != nullptr);
    assert(dart->isAnti() != dart->getAlpha<3>()->isAnti());

    fixProblem(dart, dart->getAlpha<3>(), mergeVertices);

    return true;

}

bool HexExtractor::fixProblems4(bool mergeVertices)
{
    auto n = (int)intermediateHexMesh.n_vertices();

    bool changes = false;

//OMP_PARALLEL_FOR
    for (auto i = 0; i < n; ++i)
    {
        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            if (d->getAlpha<3>() != nullptr)
                if (d->isAnti() != d->getAlpha<3>()->isAnti())
                {
//                    joinEquivalenceClasses(d->getVertex(), d->getAlpha<3>()->getVertex());
//                    joinEquivalenceClasses(d->getVertex(), d->getAlpha<3>()->getAlpha<1>()->getVertex());
                    fixProblem4(d.get(), mergeVertices);
                    changes = true;
                }
    }

    return changes;

}

void HexExtractor::fixProblems(int iterations, bool mergeVertices)
{

    HEXEX_DEBUG_ONLY(std::cout << "fixing problems" << std::endl;)

    bool changes = true;

    for (int i = 0; i < iterations && changes; ++i)
    {

        HEXEX_DEBUG_ONLY(std::cout << "fixing problems iteration " << i+1 << " of " << iterations << std::endl;)

        changes = true;

        if (fixProblems4(mergeVertices)) continue;

        if (fixProblems3(mergeVertices)) continue;

        if (fixProblems2(mergeVertices)) continue;

        if (fixProblems1(mergeVertices)) continue;

        changes = false;

    }
}


const Transition& HexExtractor::getTransitionFunction(HalfFaceHandle hfh)
{
    if (!transitionFunctionsComputed)
        extractTransitionFunctions();
    return transitionFunctions[hfh];
}

void HexExtractor::setTransitionFunction(HalfFaceHandle hfh, Transition transitionFunction)
{
    transitionFunctions[hfh] = transitionFunction;
}

Transition HexExtractor::getTransitionFunctionAroundHalfedge(CellHandle ch, HalfEdgeHandle heh)
{
    auto tranFun = identity;
    auto valence = inputMesh.valence(inputMesh.edge_handle(heh));
    for (auto i = 0u; i < valence; ++i)
    {
        auto transitionFace =  rotateAroundHalfedge(ch, heh);
        doTransition(transitionFace, ch);
        doTransition(transitionFace, tranFun);
    }

    return tranFun;
}


void HexExtractor::computeLocalUVFromSecondaryDarts(Dart* dart, CellHandle refCell, Transition transitionToRefCell, std::vector<Dart*>& processedDarts)
{
    if (std::find(processedDarts.begin(), processedDarts.end(), dart) != processedDarts.end())
        return;

    processedDarts.push_back(dart);

    auto UV = transitionToRefCell.transform_point(dart->getParameter());

    assert((UV[0] == 0) || (UV[0] == 1));
    assert((UV[1] == 0) || (UV[1] == 1));
    assert((UV[2] == 0) || (UV[2] == 1));

    localCellUVs[refCell][dart->getTracePort().vertex()] = UV;
    dart->getUV() = UV;

    auto neighbor = dart->getAlpha<2>();
    if (neighbor != nullptr)
    {
        auto tranFromNeighbor = dart->getTransitionToAlpha<2>();
        tranFromNeighbor.invert();

        computeLocalUVFromSecondaryDarts(neighbor, refCell, transitionToRefCell*tranFromNeighbor, processedDarts);
    }

    auto previousDart = dart->getAlpha<1>();
    if (previousDart != nullptr)
    {
        auto tranFromPrevious = dart->getTransitionToAlpha<1>();
        tranFromPrevious.invert();

        computeLocalUVFromSecondaryDarts(previousDart, refCell, transitionToRefCell*tranFromPrevious, processedDarts);
    }


    auto nextDart = dart->getAlpha<0>();
    if (nextDart != nullptr)
    {
        auto tranFromPartner = dart->getTransitionToAlpha<0>();
        tranFromPartner.invert();

        computeLocalUVFromSecondaryDarts(nextDart, refCell, transitionToRefCell*tranFromPartner, processedDarts);
    }

}

void HexExtractor::computeLocalUVFromSecondaryDarts(CellHandle ch)
{

    auto cell = intermediateHexMesh.cell(ch);
    auto& halffaces = cell.halffaces();

    auto firstDart = halffaceDarts[halffaces.front()].front();

    auto dir = firstDart->getTraceDir();
    auto refDir = firstDart->getRefDir();
    auto normalDir = firstDart->getNormalDir();
    auto pos = firstDart->getParameter();

    auto transitionToRef = identity;

    transitionToRef = Transition(getLocalFrame(dir, refDir, normalDir, pos));

    if (dir % refDir != normalDir)
        transitionToRef = Transition(getLocalFrame(refDir, dir, normalDir, pos));

    transitionToRef.invert();

    auto processedDarts = std::vector<Dart*>();
    computeLocalUVFromSecondaryDarts(firstDart, ch, transitionToRef, processedDarts);

    auto minUV = Parameter(99999,99999,99999);
    auto maxUV = -1.0 * minUV;

    for (auto vh : intermediateHexMesh.cell_vertices(ch))
    {
        if (localCellUVs[ch].find(vh) == localCellUVs[ch].end())
            continue;
        auto uv = localCellUVs[ch][vh];
        for (auto i = 0u; i < 3u; ++i)
        {
            minUV[i] = std::min(uv[i], minUV[i]);
            maxUV[i] = std::max(uv[i], maxUV[i]);
        }
    }


    if (maxUV - minUV != Parameter(1,1,1))
    {
        for (auto vh : intermediateHexMesh.cell_vertices(ch))
        {
            if (localCellUVs[ch].find(vh) == localCellUVs[ch].end())
                continue;
            auto uv = localCellUVs[ch][vh];
            for (auto i = 0u; i < 3u; ++i)
            {
                minUV[i] = std::min(uv[i], minUV[i]);
                maxUV[i] = std::max(uv[i], maxUV[i]);
            }
        }
    }

    for (auto vh : intermediateHexMesh.cell_vertices(ch))
    {
        localCellUVs[ch][vh] = localCellUVs[ch][vh] - minUV;
    }
}

void HexExtractor::computeLocalUVsFromSecondaryDarts()
{

    for (auto ch : intermediateHexMesh.cells())
        localCellUVs[ch].clear();

    HEXEX_DEBUG_ONLY(std::cout << "computing local uvs" << std::endl;)

    for (auto ch : intermediateHexMesh.cells())
    {
        HEXEX_DEBUG_ONLY(if ((ch.idx() % 10000) == 0)
          std::cout << "processing cell " << ch.idx() << " of " << intermediateHexMesh.n_cells() << std::endl;)

        computeLocalUVFromSecondaryDarts(ch);
    }
}


void HexExtractor::showLocalUVs(CellHandle ch)
{
    // clear all uvs
    for (auto v_it = intermediateHexMesh.vertices_begin(); v_it != intermediateHexMesh.vertices_end(); ++v_it)
        localUVs[*v_it] = Parameter(0.5,0.5,0.5);


    for (auto cv_it = intermediateHexMesh.cv_iter(ch); cv_it.valid(); ++cv_it)
    {
        localCellUVs[ch][*cv_it] = Parameter(0.5,0.5,0.5);
    }

    computeLocalUVFromSecondaryDarts(ch);

    for (auto cv_it = intermediateHexMesh.cv_iter(ch); cv_it.valid(); ++cv_it)
    {
        auto uv = localCellUVs[ch][*cv_it];
        localUVs[*cv_it] = uv;
    }

    auto vs = std::vector<VertexHandle>();
    for (auto cv_it = intermediateHexMesh.cv_iter(ch); cv_it.valid(); ++cv_it)
    {
        vs.push_back(*cv_it);
    }

}

void HexExtractor::deleteHexEdges()
{
    intermediateHexMesh.enable_deferred_deletion(true);

    for (auto vh : intermediateHexMesh.vertices())
    {
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            d->getHalfedge() = HalfEdgeHandle();
        auto& sds = secondaryDarts[vh];
        for (auto& d : sds)
            d->getHalfedge() = HalfEdgeHandle();
    }

    for (auto eh: intermediateHexMesh.edges())
        intermediateHexMesh.delete_edge(eh);

    intermediateHexMesh.collect_garbage();



}

void HexExtractor::deleteHexFaces()
{
    intermediateHexMesh.enable_deferred_deletion(true);

    for (auto vh : intermediateHexMesh.vertices())
    {
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            d->getHalfface() = HalfFaceHandle();
        auto& sds = secondaryDarts[vh];
        for (auto& d : sds)
            d->getHalfface() = HalfFaceHandle();
    }

    for (auto fh: intermediateHexMesh.faces())
        intermediateHexMesh.delete_face(fh);

    deleteHexEdges();

    intermediateHexMesh.collect_garbage();

}

void HexExtractor::deleteHexCells()
{
    intermediateHexMesh.enable_deferred_deletion(true);

    for (auto ch: intermediateHexMesh.cells())
        intermediateHexMesh.delete_cell(ch);

    intermediateHexMesh.collect_garbage();
}

void HexExtractor::extractFaceFromSecondaryDarts(Dart* dart)
{
    auto currentDart = dart;
    auto faceDarts = std::vector<Dart*>();
    auto faceSecondaryDarts = std::vector<Dart*>();
    auto halfedges = std::vector<HalfEdgeHandle>();

    assert(dart->isPrimary());

    auto i = 0;
    do
    {

        if (nullptr == currentDart->getAlpha<0>())
        {
            return;
        }

        assert(currentDart->isPrimary());
        faceDarts.push_back(currentDart);

        if (!currentDart->getHalfedge().is_valid())
        {
            auto heh = addEdge(currentDart->getTracePort(), currentDart->getAlpha<0>()->getTracePort());
            currentDart->getHalfedge() = heh;
        }
        halfedges.push_back(currentDart->getHalfedge());

        currentDart = currentDart->getAlpha<0>();

        assert(currentDart->isSecondary());
        faceSecondaryDarts.push_back(currentDart);

        currentDart = currentDart->getAlpha<1>();

        if (nullptr == currentDart)
            return;

        if (currentDart != dart)
            if (std::find(faceDarts.begin(), faceDarts.end(), currentDart) != faceDarts.end())
            {
                HEXEX_DEBUG_ONLY(std::cout << "abort loop " << "trace port: " << dart->getTracePort() << " ref port: " << dart->getRefDir() << std::endl;
                for (auto d : faceDarts)
                    std::cout << "trace port: " << d->getTracePort() << " ref port: " << d->getRefDir() << "  " << (d == dart ? "==" : "!=") << std::endl;)
                if (i++ > 5)
                    return;
            }

    }while (currentDart != dart);

    auto hfh = getHalfFace(intermediateHexMesh, halfedges);
    if (!hfh.is_valid() && halfedges.size() > 1)
    {
        auto fh = intermediateHexMesh.add_face(halfedges);
        hfh = intermediateHexMesh.halfface_handle(fh, 0);
    }

    if (hfh.is_valid())
    {
        halffaceDarts[hfh] = faceDarts;
        halffaceSecondaryDarts[hfh] = faceSecondaryDarts;

#ifdef DEBUG
        for (auto d : halffaceDarts[hfh])
          assert(d->isPrimary());
        for (auto d : halffaceSecondaryDarts[hfh])
          assert(d->isSecondary());
#endif

        for (auto c : faceDarts)
            c->getHalfface() = hfh;
        for (auto c : faceSecondaryDarts)
            c->getHalfface() = hfh;

    }
}

void HexExtractor::extractFacesFromDarts()
{
    deleteHexFaces();

    HEXEX_DEBUG_ONLY(std::cout << "extracting faces" << std::endl;)

    for (auto v_it = intermediateHexMesh.vertices_begin(); v_it.valid(); ++v_it)
    {
        HEXEX_DEBUG_ONLY(if ((v_it->idx() % 10000) == 0)
          std::cout << "processing vertex " << v_it->idx() << " of " << intermediateHexMesh.n_vertices() << std::endl;)

        auto& ds = darts[*v_it];

        for (auto d: ds)
        {
            if (!d->getHalfface().is_valid())
                extractFaceFromSecondaryDarts(d.get());
        }
    }

    calculateDifferencesInDartTypes();
}

void HexExtractor::extractCellFromSecondaryDarts(HalfFaceHandle hfh)
{

    if (halffaceDarts[hfh].empty())
    {
        // boundary
        return;
    }

    auto halffaceSet = std::set<HalfFaceHandle>();


    auto ds = getAllDartsInCell(halffaceDarts[hfh].front());

    for (auto d: ds)
        if (d->getHalfface().is_valid())
            halffaceSet.insert(d->getHalfface());

    auto halffaceVector = std::vector<HalfFaceHandle>(halffaceSet.begin(), halffaceSet.end());

    for (auto hfh : halffaceVector)
        if (intermediateHexMesh.incident_cell(hfh).is_valid())
        {
            HEXEX_DEBUG_ONLY(std::cout << "Error: Halfface is already incident to a cell" << std::endl);
        }

    intermediateHexMesh.add_cell(halffaceVector, false);
}

void HexExtractor::extractCellsFromDarts()
{
    deleteHexCells();


    HEXEX_DEBUG_ONLY(std::cout << "extracting cells" << std::endl;)

    for (auto hf_it = intermediateHexMesh.halffaces_begin(); hf_it != intermediateHexMesh.halffaces_end(); ++hf_it)
    {
        HEXEX_DEBUG_ONLY(if ((hf_it->idx() % 10000) == 0)
          std::cout << "processing halfface " << hf_it->idx() << " of " << intermediateHexMesh.n_halffaces() << std::endl;)

        auto hfh = *hf_it;
        auto ch = intermediateHexMesh.incident_cell(hfh);
        if (!ch.is_valid()) // ie was not extracted yet
            extractCellFromSecondaryDarts(hfh);
    }

    calculateDifferencesInDartTypes();
}



std::vector<HalfFaceHandle> HexExtractor::getAdjacentHalffacesFromSecondaryDarts(HalfFaceHandle hfh)
{
    auto res = std::vector<HalfFaceHandle>();
    auto& ds = halffaceDarts[hfh];
    for (auto d : ds)
    {
        auto neighbor = d->getAlpha<2>();
        if (neighbor != nullptr)
        {
            auto realNeighbor = neighbor->getAlpha<1>();
            if (realNeighbor != nullptr)
            {
                auto neighborhfh = realNeighbor->getHalfface();
                if (neighborhfh != hfh && neighborhfh.is_valid())
                    res.push_back(neighborhfh);
            }
        }
    }
    std::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());

    return res;
}

void HexExtractor::extractTransitionFunction(FaceHandle fh)
{
    auto hf0 = inputMesh.halfface_handle(fh,0);
    auto hf1 = inputMesh.halfface_handle(fh,1);

    if (inputMesh.is_boundary(fh)) // no transition for boundary faces
    {
        setTransitionFunction(hf0, identity);
        setTransitionFunction(hf1, identity);
    }
    else
    {
        auto c0 = inputMesh.incident_cell(hf0);
        auto c1 = inputMesh.incident_cell(hf1);

        auto vertices0 = inputMesh.get_cell_vertices(hf0);

        auto u0 = parameter(c0, vertices0[0]);
        auto v0 = parameter(c0, vertices0[1]);
        auto w0 = parameter(c0, vertices0[2]);

        auto u1 = parameter(c1, vertices0[0]);
        auto v1 = parameter(c1, vertices0[1]);
        auto w1 = parameter(c1, vertices0[2]);

        // first check if transition function is identity
        if ((u0 == u1) && (v0 == v1) && (w0 == w1))
        {
            setTransitionFunction(hf0, identity);
            setTransitionFunction(hf1, identity);
        }
        else
        {
            auto v0_ = v0 - u0;
            auto w0_ = w0 - u0;

            auto v1_ = v1 - u1;
            auto w1_ = w1 - u1;

            auto min_dist = std::numeric_limits<double>::max();
            auto min_transition = identity;
            for (auto i = 0u; i < all24Transitions.size(); ++i)
            {
                auto v0_transformed = all24Transitions[i].transform_point(v0_);
                auto dist1 = (v1_ - v0_transformed) | (v1_ - v0_transformed);
                auto w0_transformed = all24Transitions[i].transform_point(w0_);
                auto dist2 = (w1_ - w0_transformed) | (w1_ - w0_transformed);

                if (dist1 + dist2 < min_dist)
                {
                    min_dist = dist1+dist2;
                    min_transition = all24Transitions[i];
                }
            }

            auto u0_t = min_transition.transform_point(u0);
            auto t = roundVector(u1 - u0_t);


            setTranslation(min_transition, t);

            setTransitionFunction(hf0, min_transition);
            min_transition.invert();
            setTransitionFunction(hf1, min_transition);


        }
    }
}

void HexExtractor::extractTransitionFunctions()
{

    auto n = inputMesh.n_faces();

//#pragma omp parallel for
    for (auto i = 0u; i < n; ++i)
    {
        auto fh = FaceHandle(i);
        if (inputMesh.is_deleted(fh))
            continue;
        extractTransitionFunction(fh);
    }

transitionFunctionsComputed = true;
}

void HexExtractor::computeCellTypes()
{
    for (auto ch : inputMesh.cells())
        cellTypes[ch] = computeCellType(ch);

    cellTypesComputed = true;
}

HexExtractor::CellType HexExtractor::computeCellType(CellHandle ch)
{
    auto params = getParameters(ch);
    auto sign = sign_orient3d(params[0].data(), params[1].data(), params[2].data(), params[3].data());

    if (sign == ORI_ZERO)
        return Degenerate;
    else if (sign == ORI_BELOW)
        return Flipped;
    else
        return Proper;
}

void HexExtractor::randomizeParametrization(double offsetSize, double keepBoundary)
{
    for (auto vh : inputMesh.vertices())
    {
        if (inputMesh.is_boundary(vh))
        {
            if (keepBoundary)
                continue;

            auto getTwoBoundaryHalffaces = [&](OpenVolumeMesh::EdgeHandle eh)
            {

                auto boundary_faces = std::vector<OpenVolumeMesh::HalfFaceHandle>();
                for (auto hehf_it = inputMesh.hehf_iter(inputMesh.halfedge_handle(eh,0)); hehf_it.valid(); ++hehf_it)
                {
                    if (inputMesh.is_boundary(*hehf_it))
                        boundary_faces.push_back(inputMesh.opposite_halfface_handle(*hehf_it));
                    if (inputMesh.is_boundary(inputMesh.opposite_halfface_handle(*hehf_it)))
                        boundary_faces.push_back(*hehf_it);
                }

                return boundary_faces;
            };


            auto numSingularEdges = numIncidentSingularEdges(vh);

            auto offset = getRandomVector(offsetSize);

            if (numSingularEdges == 0)
            {
                // todo: move in plane
                auto heh = *inputMesh.voh_iter(vh);
                auto iter = inputMesh.voh_iter(vh);
                while (!inputMesh.is_boundary(heh))
                    heh = *++iter;

                auto hfhs = getTwoBoundaryHalffaces(inputMesh.edge_handle(heh));

                auto hfh = HalfFaceHandle();

                if (!isFaceDegenerate(hfhs[0]))
                    hfh = hfhs[0];
                else if (!isFaceDegenerate(hfhs[1]))
                    hfh = hfhs[1];

                if (hfh.is_valid())
                {
                    auto params = getParameters(hfh);

                    auto d1 = (params[1] - params[0]).normalized();
                    auto d2 = (params[2] - params[0]).normalized();
                    auto d3 = d1 % d2;
                    d2 = d1 % d3;

                    offset = offset[0] * d1 + offset[1] * d2;

                    auto incidentCell = inputMesh.incident_cell(hfh);

                    parameter(incidentCell, vh) = parameter(incidentCell, vh) + offset;
                    propagateVertexParameter(parameter(incidentCell, vh), vh, incidentCell);
                }


            }
            else if (numSingularEdges == 1)
            {
                // inner singularity hitting boundary
                // dont do anything
            }
            else if (numSingularEdges == 2)
            {
//                auto heh = getIncidentSingularEdge(vh);

//                auto hfhs = getTwoBoundaryHalffaces(inputMesh.edge_handle(heh));

//                auto hfh = HalfFaceHandle();

//                if (!isFaceDegenerate(hfhs[0]))
//                    hfh = hfhs[0];
//                else if (!isFaceDegenerate(hfhs[1]))
//                    hfh = hfhs[1];

//                if (!inputMesh.next_halfedge_in_halfface(heh, hfh).is_valid())
//                    heh = inputMesh.opposite_halfedge_handle(heh);

//                if (hfh.is_valid())
//                {
//                    auto params = getParameters(hfh, heh);

//                    auto d1 = (params[1] - params[0]).normalized();

//                    offset = offset[0] * d1;

//                    auto incidentCell = inputMesh.incident_cell(hfh);

//                    parameter(incidentCell, vh) = parameter(incidentCell, vh) + offset;
//                    propagateVertexParameter(parameter(incidentCell, vh), vh, incidentCell);
//                }

            }
            else
            {
                // don't offset
            }
        }
        else
        {

            auto numSingularEdges = numIncidentSingularEdges(vh);

            if (numSingularEdges == 0)
            {
                auto offset = getRandomVector(offsetSize);
                auto incidentCell = *inputMesh.vc_iter(vh);
                parameter(incidentCell, vh) = parameter(incidentCell, vh) + offset;
                propagateVertexParameter(parameter(incidentCell, vh), vh, incidentCell);
            }
            else if (numSingularEdges == 1)
            {
                HEXEX_DEBUG_ONLY(std::cout << "error" << std::endl;)
                assert(false);
            }
            else if (numSingularEdges == 2)
            {
                auto incidentHalfedge = getIncidentSingularEdge(vh);
                auto incidentCell = *inputMesh.hec_iter(incidentHalfedge);

                auto tranFun = getTransitionFunctionAroundHalfedge(incidentCell, incidentHalfedge);

                auto n  = 0;
                auto dir = Direction(1,0,0);

                for (auto& dir2 : getAll6Directions())
                {
                    if (tranFun.transform_vector(dir2.vector()) == dir2.vector())
                    {
                        n++;
                        dir = dir2;
                    }
                }

                if (n != 2)
                {
                    // todo: probably error in transition functions
                    // do nothing
                }
                else
                {

                    auto offset = getRandomVector(offsetSize)[0] * dir;

                    parameter(incidentCell, vh) = parameter(incidentCell, vh) + offset;
                    propagateVertexParameter(parameter(incidentCell, vh), vh, incidentCell);
                }

            }
            else
            {
                // don't offset
            }

        }
    }
}


void HexExtractor::sanitizeParametrization(bool snapBoundary, bool extremeTruncation)
{
    computeCellTypes();

    calculateValences();
    calculateEdgeSingularities();

    truncatePrecision(extremeTruncation);

    if (snapBoundary)
        projectBoundaryFaces();

    computeCellTypes();

}

bool HexExtractor::isInCell(CellHandle ch, Parameter param)
{
    auto params = getParameters(ch);

    if (isCellFlipped(ch))
        std::swap(params[2], params[3]);

    return (isInside(params[0], params[1], params[2], params[3], param));
}

bool HexExtractor::isInCellRelaxed(CellHandle ch, Parameter param)
{
    auto params = getParameters(ch);

    for (auto p : params) // TODO: strange bug with fandiskeasy
        if (param == p)
            return true;

    if (isCellFlipped(ch))
        std::swap(params[2], params[3]);

    bool res = isInsideOrOnBoundary(params[0], params[1], params[2], params[3], param);

    return res;
}

bool HexExtractor::isInFace(HalfFaceHandle hfh, Parameter param)
{
    auto params = getParameters(hfh);
    return isInside(params[0], params[1], params[2], param);
}

bool HexExtractor::isInFaceRelaxed(HalfFaceHandle hfh, Parameter param)
{
    Face f = inputMesh.halfface(hfh);
    for (auto heh : f.halfedges())
        if (isOnEdgeRelaxed(inputMesh.incident_cell(hfh), inputMesh.edge_handle(heh), param))
            return true;
    return isInFace(hfh, param);
}

bool HexExtractor::pointsIntoCell(CellHandle ch, VertexHandle vh, Direction dir)
{
    if (isCellDegenerate(ch))
        return false;

    auto vertices = inputMesh.get_cell_vertices(ch, vh);
    auto params = getParameters(ch, vertices);

    if (isCellFlipped(ch))
        std::swap(params[2], params[3]);

    return pointsIntoTet(params[0], params[1], params[2], params[3], dir);
}

bool HexExtractor::pointsIntoCellRelaxed(CellHandle ch, VertexHandle vh, Direction dir)
{
    if (isCellDegenerate(ch))
        return false;

    auto vertices = inputMesh.get_cell_vertices(ch, vh);
    auto params = getParameters(ch, vertices);

    if (isCellFlipped(ch))
        std::swap(params[2], params[3]);

    return pointsIntoTetRelaxed(params[0], params[1], params[2], params[3], dir);
}

bool HexExtractor::pointsIntoCell(CellHandle ch, EdgeHandle eh, Direction dir, Parameter param)
{
    if (isCellDegenerate(ch))
        return false;

    HalfEdgeHandle heh = inputMesh.halfedge_handle(eh, 0);
    HalfFaceHandle hfh = getIncidentHalfFaceInCell(inputMesh, ch, heh);
    auto vertices = inputMesh.get_cell_vertices(hfh, heh);
    auto params = getParameters(ch, vertices);


    if (isCellFlipped(ch))
        std::swap(params[2], params[3]);

    return isRegular(params[0], params[1], params[2], param+dir) &&
            isRegular(params[0], params[3], params[1], param+dir);
}

bool HexExtractor::pointsIntoCell(HalfFaceHandle hfh, Direction dir, Parameter param)
{
    if (isCellDegenerate(inputMesh.incident_cell(hfh)))
        return false;
    auto vertices = inputMesh.get_cell_vertices(hfh);
    auto params = getParameters(inputMesh.incident_cell(hfh), vertices);

    if (isCellFlipped(inputMesh.incident_cell(hfh)))
        std::swap(params[1], params[2]);

    return isRegular(params[0], params[1], params[2], param+dir) &&
           isInside(params[0], params[1], params[2], param);

}

bool HexExtractor::pointsIntoCellRelaxed(HalfFaceHandle hfh, Direction dir, Parameter param)
{
    if (isCellDegenerate(inputMesh.incident_cell(hfh)))
        return false;
    auto vertices = inputMesh.get_cell_vertices(hfh);
    auto params = getParameters(inputMesh.incident_cell(hfh), vertices);

    if (isCellFlipped(inputMesh.incident_cell(hfh)))
        std::swap(params[1], params[2]);

    return !isFlipped(params[0], params[1], params[2], param+dir) &&
           isInsideOrOnBoundary(params[0], params[1], params[2], param);
}

bool HexExtractor::pointsAlongHalfEdge(CellHandle ch, HalfEdgeHandle heh, Direction dir)
{
    auto e = inputMesh.halfedge(heh);
    auto param1 = parameter(ch, e.from_vertex());
    auto param2 = parameter(ch, e.to_vertex());
    if (param1 == param2)
        return false;
    return areCoLinear(param1, param2, dir);
}

bool HexExtractor::pointsAlongHalfEdgeFromVertex(CellHandle ch, HalfEdgeHandle heh, Direction dir)
{
    auto e = inputMesh.halfedge(heh);
    auto u = parameter(ch, e.from_vertex());
    auto v = parameter(ch, e.to_vertex());
    if (u == v)
      return false;
    return (isOnLine(u, v, u+dir) && ((v|dir) > (u|dir)));
}

bool HexExtractor::pointsAlongEdge(CellHandle ch, EdgeHandle eh, Direction dir, Parameter param)
{
    auto e = inputMesh.edge(eh);
    auto param1 = parameter(ch, e.from_vertex());
    auto param2 = parameter(ch, e.to_vertex());
    if (param1 == param2)
        return false;
    return isOnLineBetween(param1, param2, param) &&
            (areCoLinear(param1, param2, dir) || areCoLinear(param2, param1, dir));
}

bool HexExtractor::pointsAlongFace(HalfFaceHandle hfh, Direction dir, Parameter param)
{
    if (isFaceDegenerate(hfh))
        return false;
    auto params = getParameters(hfh);
    return isInPlane(params[0], params[1], params[2], param) &&
           isInPlane(params[0], params[1], params[2], param+dir);
}

bool HexExtractor::pointsIntoFace(HalfFaceHandle hfh, VertexHandle vh, Direction dir)
{
    if (isFaceDegenerate(hfh))
        return false;
    auto vertices = inputMesh.get_halfface_vertices(hfh, vh);
    auto params   = getParameters(inputMesh.incident_cell(hfh), vertices);
    return pointsIntoTriangle(params[0], params[1], params[2], dir);
}

bool HexExtractor::pointsIntoFace(HalfFaceHandle hfh, EdgeHandle eh, Direction dir, Parameter param)
{
    if (isFaceDegenerate(hfh))
        return false;
    auto f = inputMesh.halfface(hfh);
    auto heh = inputMesh.halfedge_handle(eh, 0);
    if (std::find(f.halfedges().begin(), f.halfedges().end(), heh) == f.halfedges().end())
        heh = inputMesh.opposite_halfedge_handle(heh);
    assert(std::find(f.halfedges().begin(), f.halfedges().end(), heh) != f.halfedges().end());

    auto vertices = inputMesh.get_halfface_vertices(hfh, heh);
    auto params = getParameters(inputMesh.incident_cell(hfh), vertices);
    return isInPlane(params[0], params[1], params[2], param) &&
           isInPlane(params[0], params[1], params[2], param + dir) &&
           pointsIntoProjectedTriangleFromEdge(params[0], params[1], params[2], param, dir);
}


bool HexExtractor::isOnEdge(CellHandle ch, EdgeHandle eh, Parameter param)
{
    auto e = inputMesh.edge(eh);
    auto param1 = parameter(ch, e.from_vertex());
    auto param2 = parameter(ch, e.to_vertex());
    return isOnLineBetween(param1, param2, param);
}

bool HexExtractor::isOnEdgeRelaxed(CellHandle ch, EdgeHandle eh, Parameter param)
{
    auto e = inputMesh.edge(eh);
    auto param1 = parameter(ch, e.from_vertex());
    auto param2 = parameter(ch, e.to_vertex());
    return (param1 == param) || (param2 == param) || isOnLineBetween(param1, param2, param);
}

VertexHandle HexExtractor::getVertexWithParam(CellHandle ch, Parameter param)
{
    for (auto vh : cellVertices[ch])
        if (parameter(ch, vh) == param)
            return vh;
    return VertexHandle();
}

EdgeHandle HexExtractor::getEdgeWithParam(CellHandle ch, Parameter param)
{
    auto c = inputMesh.cell(ch);
    for (auto hfh : c.halffaces())
    {
        auto f = inputMesh.halfface(hfh);
        for (auto heh : f.halfedges())
        {
            auto e = inputMesh.halfedge(heh);
            auto param1 = parameter(ch, e.from_vertex());
            auto param2 = parameter(ch, e.to_vertex());
            if (isOnLineBetween(param1, param2, param))
                return inputMesh.edge_handle(heh);
        }
    }
    return EdgeHandle();
}

HalfFaceHandle HexExtractor::getHalffaceWithParam(CellHandle ch, Parameter param)
{
    auto c = inputMesh.cell(ch);
    for (auto hfh : c.halffaces())
    {
        auto params = getParameters(hfh);
        if (isInside(params[0], params[1], params[2], param))
            return hfh;
    }
    return HalfFaceHandle();
}

int HexExtractor::numIncidentSingularEdges(VertexHandle vh)
{
    int res = 0;
    for (auto voh_it = inputMesh.voh_iter(vh); voh_it.valid(); ++voh_it)
        if (isSingularEdge(inputMesh.edge_handle(*voh_it)))
            ++res;
    return res;
}

HalfEdgeHandle HexExtractor::getIncidentSingularEdge(VertexHandle vh)
{
    for (auto voh_it = inputMesh.voh_iter(vh); voh_it.valid(); ++voh_it)
        if (isSingularEdge(inputMesh.edge_handle(*voh_it)))
            return *voh_it;

    return HalfEdgeHandle();
}

HalfEdgeHandle HexExtractor::getIncidentSingularEdge(VertexHandle vh, CellHandle ch)
{
    for (auto voh_it = inputMesh.voh_iter(vh); voh_it.valid(); ++voh_it)
        if (isSingularEdge(inputMesh.edge_handle(*voh_it)))
            if (isIncident(*voh_it, ch))
                return *voh_it;

    return HalfEdgeHandle();
}

CellHandle HexExtractor::getIncidentCellIncidentToSingularEdge(VertexHandle vh)
{
    auto heh = getIncidentSingularEdge(vh);
    return *inputMesh.hec_iter(heh);
}

bool HexExtractor::isIncident(HalfEdgeHandle heh, CellHandle ch)
{
    for (auto hec_it = inputMesh.hec_iter(heh); hec_it.valid(); ++hec_it)
        if (*hec_it == ch)
            return true;
    return false;
}

std::vector<Parameter> HexExtractor::getParameters(CellHandle ch, std::vector<VertexHandle> vhs)
{
    auto res = std::vector<Parameter>();
    for (auto vh : vhs)
        res.push_back(vertexParameters[ch][vh]);
    return res;
}

std::vector<Parameter> HexExtractor::getParameters(CellHandle ch)
{
    auto vertices = cellVertices[ch];
    return getParameters(ch, vertices);
}

std::vector<Parameter> HexExtractor::getParameters(HalfFaceHandle hfh)
{
    auto vertices = inputMesh.get_halfface_vertices(hfh);
    return getParameters(inputMesh.incident_cell(hfh), vertices);
}

std::vector<Parameter> HexExtractor::getParameters(HalfFaceHandle hfh, HalfEdgeHandle heh)
{
    auto vertices = inputMesh.get_halfface_vertices(hfh, heh);
    return getParameters(inputMesh.incident_cell(hfh), vertices);
}

Position HexExtractor::getPosition(Parameter param, CellHandle ch)
{
    auto vertices = cellVertices[ch];
    auto positions = std::vector<Position>();
    auto parameters = std::vector<Parameter>();
    for (auto vh : vertices)
    {
        positions.push_back(inputPosition(vh));
        parameters.push_back(parameter(ch, vh));
    }
    auto matrix = getInverseParametrizationMatrix(positions[0],positions[1],positions[2],positions[3], parameters[0],parameters[1],parameters[2],parameters[3]);
    return matrix.transform_point(param);
}

Parameter HexExtractor::getParameter(Position pos, CellHandle ch)
{
    auto vertices = cellVertices[ch];
    auto positions = std::vector<Position>();
    auto parameters = std::vector<Parameter>();
    for (auto vh : vertices)
    {
        positions.push_back(inputPosition(vh));
        parameters.push_back(parameter(ch, vh));
    }
    auto matrix = getParametrizationMatrix(positions[0],positions[1],positions[2],positions[3], parameters[0],parameters[1],parameters[2],parameters[3]);
    return matrix.transform_point(pos);
}

Parameter HexExtractor::getHexVertexParameter(VertexHandle hexVh, CellHandle ch)
{
    auto param = hexvertexParameter[hexVh];
    auto incidentCell = incidentCellInInputMesh[hexVh];
    auto type = vertexTypes[hexVh];
    if (type == VHVertex)
    {
        auto incidentVertex = VertexHandle(incidentElementId[hexVh]);
        /*auto tranFun = getTransitionFunction(incidentCell, ch, incidentVertex);
        param = tranFun.transform_point(param);*/
        param = parameter(ch, incidentVertex);
    }
    else if (type == EHVertex)
    {
        auto incidentEdge = EdgeHandle(incidentElementId[hexVh]);
        auto tranFun = getTransitionFunction(incidentCell, ch, incidentEdge);
        param = tranFun.transform_point(param);
    }
    else if (type == FHVertex)
    {
        auto tranFun = identity;
        if (incidentCell != ch)
        {
            auto incidentHalfFace = getIncidentHalfFaceInCell(inputMesh, incidentCell, FaceHandle(incidentElementId[hexVh]));
            tranFun = getTransitionFunction(incidentHalfFace);
            assert(inputMesh.incident_cell(inputMesh.opposite_halfface_handle(incidentHalfFace)) == ch);
        }
        param = tranFun.transform_point(param);
    }
    else if (type == CHVertex)
    {
        // nothing todo
        assert(incidentCell == ch);
    }
    else
    {
        // there should be no other types
        assert(false);
    }

    return param;

}

Matrix4x4d HexExtractor::getParametrizationMatrix(Position p, Position q, Position r, Position s, Parameter u, Parameter v, Parameter w, Parameter t)
{
    Matrix4x4d mat1;
    mat1(0,0) = q[0]-p[0]; mat1(0,1) = r[0]-p[0]; mat1(0,2) = s[0]-p[0]; mat1(0,3) = p[0];
    mat1(1,0) = q[1]-p[1]; mat1(1,1) = r[1]-p[1]; mat1(1,2) = s[1]-p[1]; mat1(1,3) = p[1];
    mat1(2,0) = q[2]-p[2]; mat1(2,1) = r[2]-p[2]; mat1(2,2) = s[2]-p[2]; mat1(2,3) = p[2];
    mat1(3,0) = 0;         mat1(3,1) = 0;         mat1(3,2) = 0;         mat1(3,3) = 1;

    Matrix4x4d mat2;
    mat2(0,0) = v[0]-u[0]; mat2(0,1) = w[0]-u[0]; mat2(0,2) = t[0]-u[0]; mat2(0,3) = u[0];
    mat2(1,0) = v[1]-u[1]; mat2(1,1) = w[1]-u[1]; mat2(1,2) = t[1]-u[1]; mat2(1,3) = u[1];
    mat2(2,0) = v[2]-u[2]; mat2(2,1) = w[2]-u[2]; mat2(2,2) = t[2]-u[2]; mat2(2,3) = u[2];
    mat2(3,0) = 0;         mat2(3,1) = 0;         mat2(3,2) = 0;         mat2(3,3) = 1;

    mat1.invert();

    return mat2*mat1;
}

Matrix4x4dd HexExtractor::getInverseParametrizationMatrix(Position p, Position q, Position r, Parameter u, Parameter v, Parameter w)
{
    auto s = p + ((q-p)%(r-p)).normalized();
    auto t = u + ((v-u)%(w-u)).normalized();
    return getInverseParametrizationMatrix(p,q,r,s, u,v,w,t);
}

Matrix4x4dd HexExtractor::getLocalFrame(Direction dir, Direction refDir, Direction normal, Parameter position)
{
    auto res = Matrix4x4d();
    res(0,0) = dir.vector()[0]; res(0,1) = refDir.vector()[0]; res(0,2) = normal.vector()[0]; res(0,3) = position[0];
    res(1,0) = dir.vector()[1]; res(1,1) = refDir.vector()[1]; res(1,2) = normal.vector()[1]; res(1,3) = position[1];
    res(2,0) = dir.vector()[2]; res(2,1) = refDir.vector()[2]; res(2,2) = normal.vector()[2]; res(2,3) = position[2];
    res(3,0) =               0; res(3,1) =                  0; res(3,2) =                  0; res(3,3) =           1;

    return res;
}

Matrix4x4dd HexExtractor::getInverseParametrizationMatrix(Position p, Position q, Position r, Position s, Parameter u, Parameter v, Parameter w, Parameter t)
{
    Matrix4x4d mat1;
    mat1(0,0) = p[0]; mat1(0,1) = q[0]; mat1(0,2) = r[0]; mat1(0,3) = s[0];
    mat1(1,0) = p[1]; mat1(1,1) = q[1]; mat1(1,2) = r[1]; mat1(1,3) = s[1];
    mat1(2,0) = p[2]; mat1(2,1) = q[2]; mat1(2,2) = r[2]; mat1(2,3) = s[2];
    mat1(3,0) = 1;    mat1(3,1) = 1;    mat1(3,2) = 1;    mat1(3,3) = 1;

    Matrix4x4d mat2;
    mat2(0,0) = u[0]; mat2(0,1) = v[0]; mat2(0,2) = w[0]; mat2(0,3) = t[0];
    mat2(1,0) = u[1]; mat2(1,1) = v[1]; mat2(1,2) = w[1]; mat2(1,3) = t[1];
    mat2(2,0) = u[2]; mat2(2,1) = v[2]; mat2(2,2) = w[2]; mat2(2,3) = t[2];
    mat2(3,0) = 1;    mat2(3,1) = 1;    mat2(3,2) = 1;    mat2(3,3) = 1;

    mat2.invert();

    return mat1*mat2;
}

Parameter HexExtractor::getParameterNormal(HalfFaceHandle hfh)
{
    auto vertices = inputMesh.get_halfface_vertices(hfh);
    auto ch = inputMesh.incident_cell(hfh);
    if (!ch.is_valid())
        ch = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(hfh));
    auto u = parameter(ch, vertices[0]);
    auto v = parameter(ch, vertices[1]);
    auto w = parameter(ch, vertices[2]);

    auto n = ((v-u)%(w-u));

    if (n.length() < 1e-6)
    {
        // dont trust n
        if (isDegenerate(u,v,w,u+n))
            return -1.0*n.normalized();
        else
            return n.normalized();
    }
    else
        return n.normalized();

}

void HexExtractor::calculateValences()
{
    extractTransitionFunctions();
    for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
    {
        edgeValences[*e_it] = edgeValence(*e_it);
    }
}

void HexExtractor::truncatePrecision(bool extremeTruncation)
{
    for (auto v_it = inputMesh.vertices_begin(); v_it != inputMesh.vertices_end(); ++v_it)
    {
        double maxCoord = 0;
        for (auto vc_it = inputMesh.vc_iter(*v_it); vc_it.valid(); ++vc_it)
        {
            auto p = parameter(*vc_it, *v_it);
            for (unsigned int i = 0; i < 3; ++i)
                maxCoord = std::max(std::fabs(p[i]), maxCoord);
        }

        double delta = std::pow(2, std::ceil(std::log(maxCoord)/std::log(2)));

        auto ch = CellHandle();
        if (isSingularVertex(*v_it))
            ch = getIncidentCellIncidentToSingularEdge(*v_it);
        else
            ch = *inputMesh.vc_iter(*v_it);

        {
            if (!ch.is_valid())
            {
                auto ch = CellHandle();
                if (isSingularVertex(*v_it))
                    ch = getIncidentCellIncidentToSingularEdge(*v_it);
                else
                    ch = *inputMesh.vc_iter(*v_it);
            }
        }

        if(!ch.is_valid()) // tolerate isolated vertices but output warning in debug mode
        {
            HEXEX_DEBUG_ONLY( std::cerr << "Warning: vertex skipped during sanitization since invalid cell handle found for vertex with valence " << inputMesh.valence(*v_it)  << std::endl;)
            continue;
        }
        // assert(ch.is_valid());


        // precision truncation
        auto p = parameter(ch, *v_it);
        for (unsigned int i = 0; i < 3; ++i)
        {
            int sign = std::signbit(p[i]) ? -1 : 1;
            volatile double tmp = p[i]+sign*delta;
            parameter(ch, *v_it)[i] = tmp - sign*delta;
        }


        // TODO: remove this hotfix
        for (unsigned int i = 0; i < 3; ++i)
        {
            volatile double tmp = p[i]+pow(2,10);
            parameter(ch, *v_it)[i] = tmp - pow(2,10);
        }

// extreme truncation
        if (extremeTruncation)
        {
            for (unsigned int i = 0; i < 3; ++i)
            {
    //          if (std::abs(parameter(ch, *v_it)[i] - (int)round(parameter(ch, *v_it)[i])) < 0.0001)
                parameter(ch, *v_it)[i] = (int)round(parameter(ch, *v_it)[i]);
            }
        }

        if (isSingularVertex(*v_it))
            fixSingularityPoint(*v_it, ch);

        propagateVertexParameter(parameter(ch, *v_it), *v_it, ch);
    }
}

void HexExtractor::checkThisOneProperty()
{
    auto n = (int)intermediateHexMesh.n_vertices();

//#pragma omp parallel for
    for (auto i = 0; i < n; ++i)
    {
        auto vh = VertexHandle(i);
        if (intermediateHexMesh.is_deleted(vh))
            continue;
        auto& ds = darts[vh];
        for (auto& d : ds)
            if (d->isPrimary())
            {
                auto ds2 = getDartsBetweenDarts12(d.get(), d.get());

                auto posDarts = 0;
                auto negDarts = 0;

                for (auto& d2: ds2)
                {
                    if (d2->isAnti())
                        negDarts++;
                    else
                        posDarts++;
                }

                auto diff = posDarts - negDarts;

                if ((diff != 6) && (diff != 0) && (diff != -6))
                {
                    HEXEX_DEBUG_ONLY(std::cout << "alarmalarmalarm1 " << ds2.size() << "diff: " << diff << " vertex: " << d->getVertex();)
                    HEXEX_DEBUG_ONLY(std::cout << " " <<  d->getAlpha<0>()->getVertex()  << std::endl;)

                }
            }
    }
}


void HexExtractor::computeEquivalenceClasses()
{
    for (auto hexVh : intermediateHexMesh.vertices())
        degeneracyEquivalenceClassJoin(hexVh);
}

void HexExtractor::degeneracyEquivalenceClassJoin(VertexHandle hexVh)
{
    if (vertexTypes[hexVh] == FHVertex)
    {
        auto fh = FaceHandle(incidentElementId[hexVh]);

        auto hfhs = std::vector<HalfFaceHandle>();

        auto tmphfh = inputMesh.halfface_handle(fh,0);
        if (!inputMesh.is_boundary(tmphfh) && isCellDegenerate(inputMesh.incident_cell(tmphfh)))
            hfhs.push_back(tmphfh);

        tmphfh = inputMesh.halfface_handle(fh,1);
        if (!inputMesh.is_boundary(tmphfh) && isCellDegenerate(inputMesh.incident_cell(tmphfh)))
            hfhs.push_back(tmphfh);


        for (auto hfh : hfhs)
        {
            auto ch = inputMesh.incident_cell(hfh);
            auto param = getHexVertexParameter(hexVh, ch);

            auto oppVh = inputMesh.get_cell_vertices(hfh)[3];
            if (parameter(ch, oppVh) == param)
            {
                auto incidentHexh = incidentVerticesPerVertex[oppVh];
                if (incidentHexh.is_valid())
                    joinEquivalenceClasses(incidentHexh, hexVh);
            }

            for (auto heh : getIncidentHalfEdges(inputMesh, ch, oppVh))
            {
                auto eh = inputMesh.edge_handle(heh);
                if (isOnEdge(ch, eh, param))
                {
                    for (auto incidentHexh : incidentVerticesPerEdge[eh])
                        if (param == getHexVertexParameter(incidentHexh, ch))
                            joinEquivalenceClasses(incidentHexh, hexVh);
                }
            }


            for (auto hfh : getIncidentHalfFaces(inputMesh, ch, oppVh))
            {
                if (isInFace(hfh, param))
                {
                    auto fh = inputMesh.face_handle(hfh);
                    for (auto incidentHexVh : incidentVerticesPerFace[fh])
                        if (param == getHexVertexParameter(incidentHexVh, ch))
                            joinEquivalenceClasses(incidentHexVh, hexVh);
                }
            }
        }
    }
    else if (vertexTypes[hexVh] == EHVertex)
    {
        auto eh = EdgeHandle(incidentElementId[hexVh]);


        auto heh = inputMesh.halfedge_handle(eh, 0);

        for (auto hec_it = inputMesh.hec_iter(heh); hec_it.valid(); ++hec_it)
        {
            auto ch = *hec_it;
            auto param = getHexVertexParameter(hexVh, ch);

            if (isCellDegenerate(ch))
            {
                for (auto eh2 : getCellEdges(inputMesh, ch))
                {
                    if (eh2 != eh)
                    {
                        for (auto incidentHexVh : incidentVerticesPerEdge[eh2])
                            if (param == getHexVertexParameter(incidentHexVh, ch))
                                joinEquivalenceClasses(incidentHexVh, hexVh);
                    }
                }


                for (auto vh : cellVertices[ch])
                {
                    auto incidentHexVh = incidentVerticesPerVertex[vh];
                    if (incidentHexVh.is_valid())
                        if (param == getHexVertexParameter(incidentHexVh, ch))
                            joinEquivalenceClasses(incidentHexVh, hexVh);
                }

            }
        }
    }
    else if (vertexTypes[hexVh] == VHVertex)
    {
        auto vh = VertexHandle(incidentElementId[hexVh]);

        for (auto vc_it = inputMesh.vc_iter(vh); vc_it.valid(); ++vc_it)
        {
            auto ch = *vc_it;
            if (isCellDegenerate(ch))
            {
                auto param = getHexVertexParameter(hexVh, ch);

                for (auto vh2 : cellVertices[ch])
                {
                    if (vh2 != vh)
                    {
                        auto incidentHexVh = incidentVerticesPerVertex[vh2];
                        if (incidentHexVh.is_valid())
                            if (param == getHexVertexParameter(incidentHexVh, ch))
                                joinEquivalenceClasses(incidentHexVh, hexVh);
                    }
                }
            }
        }
    }

}

void HexExtractor::joinEquivalenceClasses(VertexHandle vh1, VertexHandle vh2)
{

    auto id1 = equivalenceClassIds[vh1];
    auto id2 = equivalenceClassIds[vh2];

    if (id1 == id2)
        return;

    if (equivalenceClasses[id1].size() < equivalenceClasses[id2].size())
        std::swap(id1,id2);

    for (auto i : equivalenceClasses[id2])
    {
        equivalenceClassIds[VertexHandle(i)] = id1;
        equivalenceClasses[id1].push_back(i);
    }
    equivalenceClasses[id2].clear();

}


void HexExtractor::fixSingularityPoint(VertexHandle vh, CellHandle& ch)
{
    auto param = parameter(ch, vh);

    if (2 == numIncidentSingularEdges(vh))
    {
        auto heh = getIncidentSingularEdge(vh, ch);
        param = projectedParam(param, ch, heh);
    }
    else
    {
        param = roundVector(param);
    }

    auto maxDisplacement = 1e-4;

    if ((parameter(ch, vh) - param).length() < maxDisplacement)
        parameter(ch,vh) = param;

}

void HexExtractor::projectBoundaryFaces()
{

    for (HalfFaceIter hf_it = inputMesh.halffaces_begin(); hf_it != inputMesh.halffaces_end(); ++hf_it)
    {
        if (inputMesh.is_boundary(*hf_it))
        {
            HalfFaceHandle hfh = inputMesh.opposite_halfface_handle(*hf_it);
            CellHandle ch = inputMesh.incident_cell(hfh);
            if (!ch.is_valid())
            {
                continue;
            }

            // find projection direction ( = dominant normal direction)
            //TODO: this is a bad idea for degenerate boundary triangles
//            auto n = getParameterNormal(hfh);
//            std::cout << "normal: " << n << std::endl;
//            int projectionCoord = getMaxCoord(n);

            // for now try to snap in all direction
            for (unsigned int i = 0; i < 3; ++i)
            {

                int projectionCoord = i;

                // find projected parameter (rounded average of all three vertices)
                double tmpProjectedParameter = 0.0;
                for (HalfFaceVertexIter hfv_it = inputMesh.hfv_iter(hfh); hfv_it.valid(); ++hfv_it)
                    tmpProjectedParameter += parameter(ch, *hfv_it)[projectionCoord];
                int projectedParameter = round(tmpProjectedParameter/3.0);

                bool update = false;
                for (HalfFaceVertexIter hfv_it = inputMesh.hfv_iter(hfh); hfv_it.valid(); ++hfv_it)
                    if (projectedParameter != parameter(ch, *hfv_it)[projectionCoord])
                        update = true;

                auto maxDisplacement = 1e-6;
                if (update)
                    for (HalfFaceVertexIter hfv_it = inputMesh.hfv_iter(hfh); hfv_it.valid(); ++hfv_it)
                    {
                        if (std::abs(projectedParameter - parameter(ch, *hfv_it)[projectionCoord]) > maxDisplacement)
                        {
                            update = false;
                        }
                    }



                // set projected parameter for all incident cells
                if (update)
                    for (HalfFaceVertexIter hfv_it = inputMesh.hfv_iter(hfh); hfv_it.valid(); ++hfv_it)
                    {
                        auto newParameter = parameter(ch,*hfv_it);
                        newParameter[projectionCoord] = projectedParameter;
                        propagateVertexParameter(newParameter, *hfv_it, ch);
                    }

            }

        }
    }
}


/// finds a port that could be stored in ch but might be stored in a neighboring cell
HPortHandle HexExtractor::findPort(CellHandle ch, Direction dir, Parameter param)
{
    assert(isInCellRelaxed(ch, param));

    ++callsToFindPort;

    // port must be stored in this cell
    auto& ports = hPortsInCell[ch];
    for (auto p : ports)
    {
        ++portsCheckedInFindPort;
        if ((p.dir() == dir) && (p.parameter() == param))
            return p;
    }

    // not found
    return HPortHandle();
}


Parameter HexExtractor::projectedParam(Parameter param, CellHandle ch, HalfEdgeHandle heh)
{
    auto he = inputMesh.halfedge(heh);
    auto incidentCell = *inputMesh.hec_iter(heh);

    auto param1 = parameter(incidentCell, he.from_vertex());
    auto param2 = parameter(incidentCell, he.to_vertex());
    auto dir = param2-param1;

    auto tranFun = getTransitionFunction(incidentCell, ch, he.from_vertex());

    dir = tranFun.transform_vector(dir);

    auto maxCoord = getMaxCoord(dir); // the maximal coordinate is the one that is not projected

    for (auto i = 0u; i < 3u; ++i)
    {
        if (i != maxCoord)
            param[i] = round(param[i]);
    }

    return param;

}

void HexExtractor::propagateVertexParameterRecursive(Parameter param, VertexHandle vh, CellHandle ch, std::set<CellHandle>& visited)
{
    visited.insert(ch);

    parameter(ch, vh) = param;

    auto halffaces = getIncidentHalfFaces(inputMesh, ch, vh);
    for (auto hfh : halffaces)
    {
        auto oppHfh = inputMesh.opposite_halfface_handle(hfh);
        auto oppCh = inputMesh.incident_cell(oppHfh);
        if (!oppCh.is_valid())
            continue;
        if (visited.count(oppCh) < 1)
        {
            auto newParameter = getTransitionFunction(hfh).transform_point(param);
            propagateVertexParameterRecursive(newParameter, vh, oppCh, visited);
        }
    }
}

void HexExtractor::propagateVertexParameterRecursive2(Parameter param, VertexHandle vh, CellHandle ch, std::set<CellHandle>& toBeProcessed)
{
    toBeProcessed.erase(ch);

    parameter(ch, vh) = param;

    if (toBeProcessed.empty())
        return;

    auto c = inputMesh.cell(ch);
    auto& halffaces = c.halffaces();

    for (auto hfh : halffaces)
    {
        auto oppHfh = inputMesh.opposite_halfface_handle(hfh);
        auto oppCh = inputMesh.incident_cell(oppHfh);
        if (!oppCh.is_valid())
            continue;
        if (toBeProcessed.count(oppCh) == 1)
        {
            auto newParameter = getTransitionFunction(hfh).transform_point(param);
            propagateVertexParameterRecursive2(newParameter, vh, oppCh, toBeProcessed);

            if (toBeProcessed.empty())
                return;
        }
    }
}

void HexExtractor::propagateVertexParameter(Parameter parameter, VertexHandle vh, CellHandle startCell)
{
    auto toBeProcessed = std::set<CellHandle>();
    for (auto vc_it = inputMesh.vc_iter(vh); vc_it.valid(); ++vc_it)
        toBeProcessed.insert(*vc_it);

    propagateVertexParameterRecursive2(parameter, vh, startCell, toBeProcessed);
}

Transition HexExtractor::getTransitionFunctionRecursive(CellHandle currentCell, CellHandle toCell, VertexHandle vh, Transition tranFun, std::set<CellHandle>& visited)
{
    visited.insert(currentCell);

    if (currentCell == toCell)
        return tranFun;

    auto halffaces = getIncidentHalfFaces(inputMesh, currentCell, vh);
    for (auto hfh : halffaces)
    {
        auto oppHfh = inputMesh.opposite_halfface_handle(hfh);
        auto oppCh = inputMesh.incident_cell(oppHfh);
        if (!oppCh.is_valid())
            continue;
        if (visited.count(oppCh) < 1)
        {
            auto newTranFun = tranFun;
            doTransition(hfh, newTranFun);

            auto finalTranFun = getTransitionFunctionRecursive(oppCh, toCell, vh, newTranFun, visited);
            if (visited.count(toCell) >= 1)
                return finalTranFun;
        }
    }
    return identity;
}

Transition HexExtractor::getTransitionFunction(CellHandle fromCell, CellHandle toCell, VertexHandle vh)
{
    auto visited = std::set<CellHandle>();
    return getTransitionFunctionRecursive(fromCell, toCell, vh, identity, visited);
}

Transition HexExtractor::getTransitionFunction(CellHandle fromCell, CellHandle toCell, EdgeHandle eh)
{
    auto heh = inputMesh.halfedge_handle(eh, 0);
    if (inputMesh.is_boundary(eh))
    {
        for (auto ccw : {true, false})
        {
            auto tranFun = identity;
            auto currentCell = fromCell;
            while ((currentCell != toCell) && currentCell.is_valid())
            {
                auto transitionFace =  rotateAroundHalfedge(currentCell, heh, ccw);
                doTransition(transitionFace, currentCell, tranFun);
            }
            if (currentCell.is_valid())
                return tranFun;
        }

        assert(false);
        return identity;
    }
    else
    {
        auto tranFun = identity;
        auto currentCell = fromCell;
        while ((currentCell != toCell))
        {
            auto transitionFace =  rotateAroundHalfedge(currentCell, heh, true);
            doTransition(transitionFace, currentCell, tranFun);
        }
        return tranFun;
    }
}

double HexExtractor::parametrizationAngle(HalfFaceHandle hfh1, HalfFaceHandle hfh2, HalfEdgeHandle heh)
{
    auto ch = inputMesh.incident_cell(hfh1);

    auto halfedge = inputMesh.halfedge(heh);

    auto nextHe1 = inputMesh.next_halfedge_in_halfface(heh, hfh1);
    auto nextHe2 = inputMesh.next_halfedge_in_halfface(inputMesh.opposite_halfedge_handle(heh), hfh2);

    auto vh0 = halfedge.from_vertex();
    auto vh1 = halfedge.to_vertex();
    auto vh2 = inputMesh.halfedge(nextHe1).to_vertex();
    auto vh3 = inputMesh.halfedge(nextHe2).to_vertex();

    auto u = parameter(ch, vh0);
    auto v = parameter(ch, vh1);
    auto w = parameter(ch, vh2);
    auto t = parameter(ch, vh3);

    auto d1 = v-u;
    auto d2 = w-u;
    auto d3 = d1%d2;
    d2 = d3%d1;

    auto d4 = t-u;
    auto d5 = d1%d4;
    d4 = d5%d1;

    if ((d2.length() == 0) || (d4.length() == 0))
    {
      HEXEX_DEBUG_ONLY(std::cerr << "cannot compute dihedral angle for degenerate triangle" << std::endl);
      return 0;
    }

    d2.normalize();
    d4.normalize();

    auto alpha = acos(std::max(std::min(d2|d4, 1.0),-1.0));

    return alpha;

}


int HexExtractor::edgeValence(EdgeHandle eh)
{
    auto angleSum = 0.0;

    auto heh = inputMesh.halfedge_handle(eh, 0);

    auto hehf_it = inputMesh.hehf_iter(heh,2);

    for (auto i = 0u; i < inputMesh.valence(eh); ++i, ++hehf_it)
    {
        if (inputMesh.is_boundary(*hehf_it))
            continue;

        auto hf1 = *hehf_it;
        auto hf2 = inputMesh.adjacent_halfface_in_cell(hf1, heh);

        auto ch = inputMesh.incident_cell(hf1);

        auto alpha = parametrizationAngle(hf1, hf2, heh);

        if (isCellFlipped(ch))
            angleSum -= alpha;
        else
            angleSum += alpha;

    }

    return round(angleSum / (M_PI/2.0));
}

void HexExtractor::calculateEdgeSingularity(EdgeHandle eh)
{
    if (inputMesh.is_boundary(eh))
    {

        // for safety, all boundary edges incident to a degenerate cell are considered singular
        for (auto hec_it = inputMesh.hec_iter(inputMesh.halfedge_handle(eh, 0)); hec_it.valid(); ++hec_it)
            if (isCellDegenerate(*hec_it))
            {
                edgeSingularity[eh] = true;
                return;
            }


        // get first boundary halfface
        auto heh = inputMesh.halfedge_handle(eh, 0);
        auto boundaryHalfFace1 = *inputMesh.hehf_iter(heh);
        while (!inputMesh.is_boundary(inputMesh.opposite_halfface_handle(boundaryHalfFace1)))
        {
            auto currentCell = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(boundaryHalfFace1));
            auto transitionFace =  rotateAroundHalfedge(currentCell, heh, false);
            boundaryHalfFace1 = transitionFace;
        }

        // find other boundary face and transition
        auto tranFun = identity;
        auto boundaryHalfFace2 = inputMesh.adjacent_halfface_in_cell(boundaryHalfFace1,heh);

        for (auto i = 0; i < (int)inputMesh.valence(eh)-2; ++i)
        {
            auto transitionFace = boundaryHalfFace2;
            if (isFaceDegenerate(transitionFace)) // todo: get rotational part of transition as input
            {
                edgeSingularity[eh] = true;
                return;
            }
            doTransition(transitionFace, tranFun);
            boundaryHalfFace2 = inputMesh.adjacent_halfface_in_cell(inputMesh.opposite_halfface_handle(boundaryHalfFace2), heh);
        }

        auto n0 = getParameterNormal(boundaryHalfFace1);
        auto n1 = getParameterNormal(boundaryHalfFace2);
        n0 = tranFun.transform_vector(n0);

        auto res = (n0|n1) < 0.5;

        if ((res == true) && (edgeValences[eh] == 2))
        {
            edgeSingularity[eh] = false;
            return;
        }

        edgeSingularity[eh] = res;
        return;

    }
    else
    {
        auto tranFun = identity;
        auto heh = inputMesh.halfedge_handle(eh, 0);
        auto currentCell = *inputMesh.hec_iter(heh);
        for (auto i = 0u; i < inputMesh.valence(eh); ++i)
        {
            auto transitionFace =  rotateAroundHalfedge(currentCell, heh);
            if (isFaceDegenerate(transitionFace)) // todo: get rotational part of transition as input
            {
                edgeSingularity[eh] = true;
                return;
            }
            doTransition(transitionFace, currentCell, tranFun);
        }

        if (tranFun == identity)
        {
            edgeSingularity[eh] = edgeValences[eh] > 6;
            if (edgeValences[eh] > 6)
            {
                HEXEX_DEBUG_ONLY(std::cout << "wow" << std::endl;)
                // i don't believe in valence 8 singularities
                edgeSingularity[eh] = false;
            }
            return;
        }
        else
        {
            if (edgeValences[eh] == 4)
            {
                HEXEX_DEBUG_ONLY(
                std::cout << "warning: problem, edge valence 4 but not identity transition function" << std::endl;
                std::cout << tranFun << std::endl;
                auto heh = inputMesh.halfedge_handle(eh, 0);
                auto currentCell = *inputMesh.hec_iter(heh);

                auto accTranFun = identity;


                for (auto i = 0u; i < inputMesh.valence(eh); ++i)
                {
                    auto transitionFace =  rotateAroundHalfedge(currentCell, heh);
                    doTransition(transitionFace, accTranFun, currentCell);

                    auto vertices = inputMesh.get_halfface_vertices(transitionFace);
                    std::cout << "tran fun " << std::endl << getTransitionFunction(transitionFace) << std::endl;
                    std::cout << "accumulated tran fun " << std::endl << accTranFun << std::endl;
                    for (auto vh : vertices)
                    {
                        std::cout << "face 1: " << parameter(inputMesh.incident_cell(transitionFace), vh)  << (isCellFlipped(inputMesh.incident_cell(transitionFace)) ? " inverted" : " not inverted") << std::endl;
                        std::cout << "face 2: " << parameter(currentCell, vh)  << (isCellFlipped(currentCell) ? " inverted" : " not inverted") << std::endl;
                    }


                }
                )
            }

            edgeSingularity[eh] = true;
            return;
        }
    }


}

void HexExtractor::calculateEdgeSingularities()
{
    for (auto eh : inputMesh.edges())
    {
        calculateEdgeSingularity(eh);
    }
    edgeSingularitiesCalculated = true;
}

void HexExtractor::setTranslation(GridIsomorphism& tranFun, Parameter translation)
{
    tranFun.setTranslation(translation);
}

void HexExtractor::setTranslation(Matrix4x4d& tranFun, Parameter translation)
{
    tranFun(0,3) = translation[0];
    tranFun(1,3) = translation[1];
    tranFun(2,3) = translation[2];
}

std::vector<Dart*> HexExtractor::getDartsBetweenDarts01(Dart* d_s, Dart* d_e)
{
    auto res = std::vector<Dart*>();

    auto d_c = d_s;

    do
    {
        res.push_back(d_c);
        d_c = d_c->getAlpha<0>();
        if (d_c == d_e || d_c == nullptr)
            break;
        res.push_back(d_c);
        d_c = d_c->getAlpha<1>();
    }
    while (d_c != d_e || d_c == nullptr);

    if (d_c == nullptr)
      res.clear();

    return res;
}

std::vector<Dart*> HexExtractor::getDartsBetweenDarts12(Dart* d_s, Dart* d_e)
{
  auto res = std::vector<Dart*>();

  auto d_c = d_s;

  do
  {
      res.push_back(d_c);
      d_c = d_c->getAlpha<1>();
      if (d_c == d_e || d_c == nullptr)
          break;
      res.push_back(d_c);
      d_c = d_c->getAlpha<2>();
  }
  while (d_c != d_e && d_c != nullptr);

  if (d_c == nullptr)
    res.clear();

  return res;
}

std::vector<Dart*> HexExtractor::getDartsBetweenDarts0121(Dart* d_s, Dart* d_e)
{
    auto res = std::vector<Dart*>();

    auto d_c = d_s;

    do
    {
        res.push_back(d_c);
        d_c = d_c->getAlpha<0>();
        if (d_c == d_e || d_c == nullptr)
            break;
        res.push_back(d_c);
        d_c = d_c->getAlpha<1>();
        if (d_c == d_e || d_c == nullptr)
            break;
        res.push_back(d_c);
        d_c = d_c->getAlpha<2>();
        if (d_c == d_e || d_c == nullptr)
            break;
        res.push_back(d_c);
        d_c = d_c->getAlpha<1>();
    }
    while (d_c != d_e && d_c != d_s && d_c != nullptr);

    if (d_c == d_s && d_s != d_e)
    {
        HEXEX_DEBUG_ONLY(std::cout << "error: not connected" << std::endl);
    }


    if (d_c == nullptr)
      res.clear();

    return res;
}

Transition HexExtractor::getTransitionBetweenDarts01(Dart* d_s, Dart* d_e)
{
    auto res = identity;

    auto d_c = d_s;

    do
    {
        res = d_c->getTransitionToAlpha<0>() * res;
        d_c = d_c->getAlpha<0>();
        if (d_c == d_e)
            break;
        res = d_c->getTransitionToAlpha<1>() * res;
        d_c = d_c->getAlpha<1>();
    }
    while (d_c != d_e);

    return res;
}

Transition HexExtractor::getTransitionBetweenDarts12(Dart* d_s, Dart* d_e)
{
  auto res = identity;

  auto d_c = d_s;

  do
  {
      res = d_c->getTransitionToAlpha<1>() * res;
      d_c = d_c->getAlpha<0>();
      if (d_c == d_e)
          break;
      res = d_c->getTransitionToAlpha<2>() * res;
      d_c = d_c->getAlpha<1>();
  }
  while (d_c != d_e);

  return res;

}

Transition HexExtractor::getTransitionBetweenDarts0121(Dart* d_s, Dart* d_e)
{
    auto res = identity;

    auto d_c = d_s;

    do
    {
        res = d_c->getTransitionToAlpha<0>() * res;
        d_c = d_c->getAlpha<0>();
        if (d_c == d_e)
            break;
        res = d_c->getTransitionToAlpha<1>() * res;
        d_c = d_c->getAlpha<1>();
        if (d_c == d_e)
            break;
        res = d_c->getTransitionToAlpha<2>() * res;
        d_c = d_c->getAlpha<2>();
        if (d_c == d_e)
            break;
        res = d_c->getTransitionToAlpha<1>() * res;
        d_c = d_c->getAlpha<1>();
    }
    while (d_c != d_e);

    return res;
}

Transition HexExtractor::getTransitionBetweenDarts01Backward(Dart* d_s, Dart* d_e)
{
    auto res = identity;

    auto d_c = d_s;

    do
    {
        res = d_c->getTransitionToAlpha<1>() * res;
        d_c = d_c->getAlpha<1>();
        if (d_c == d_e)
            break;
        res = d_c->getTransitionToAlpha<0>() * res;
        d_c = d_c->getAlpha<0>();
    }
    while (d_c != d_e);

    return res;
}

Transition HexExtractor::getTransitionBetweenDarts0121Backward(Dart* d_s, Dart* d_e)
{
    auto res = identity;

    auto d_c = d_s;

    do
    {
        res = d_c->getTransitionToAlpha<2>() * res;
        d_c = d_c->getAlpha<2>();
        if (d_c == d_e)
            break;
        res = d_c->getTransitionToAlpha<0>() * res;
        d_c = d_c->getAlpha<0>();
        if (d_c == d_e)
            break;
        res = d_c->getTransitionToAlpha<0>() * res;
        d_c = d_c->getAlpha<0>();
        if (d_c == d_e)
            break;
        res = d_c->getTransitionToAlpha<0>() * res;
        d_c = d_c->getAlpha<0>();
    }
    while (d_c != d_e);

    return res;
}

void HexExtractor::reconnectDarts(CellHandle ch, HalfFaceHandle hfh, VertexHandle vh0, VertexHandle vh1)
{
    if (vh0 == vh1)
        return;

    auto& ds = halffaceDarts[hfh];

    bool done = false;
    for (auto d0 : ds)
    {
        if (d0->getVertex() == vh0)
        {
            for (auto d3 : ds)
            {
                if (d3->getVertex() == vh1)
                {
                    reconnectSecondaryDartsOld(d0, d3);
                    done = true;
                }
                if (done)
                    break;
            }
        }
        if (done)
            break;
    }


    for (auto d : ds)
        d->getHalfface() = HalfFaceHandle();

    for (auto d : ds)
        if (d->getHalfface() == HalfFaceHandle())
            extractFaceFromSecondaryDarts(d);


    auto darts = std::vector<Dart*>();

    auto c = intermediateHexMesh.cell(ch);

    for (auto hfh : c.halffaces())
        darts.insert(darts.end(), halffaceDarts[hfh].begin(), halffaceDarts[hfh].end());

    intermediateHexMesh.delete_cell(ch);

    if (!intermediateHexMesh.incident_cell(intermediateHexMesh.opposite_halfface_handle(hfh)).is_valid())
        intermediateHexMesh.delete_face(intermediateHexMesh.face_handle(hfh));

    for (auto d : darts)
        if (d->getHalfface().is_valid())
            if (!intermediateHexMesh.incident_cell(d->getHalfface()).is_valid())
                extractCellFromSecondaryDarts(d->getHalfface());

}


void HexExtractor::reconnectSecondaryDartsOld(Dart* d1, Dart* d2)
{
    auto d1prev = d1->getAlpha<1>();
    auto d2prev = d2->getAlpha<1>();

    auto d1neighbor = d1->getAlpha<2>()->getAlpha<1>();
    auto d2neighbor = d2->getAlpha<2>()->getAlpha<1>();

    auto d1prevNeighbor = d1prev->getAlpha<2>();
    auto d2prevNeighbor = d2prev->getAlpha<2>();

    auto d1prevNeighborPrev = d1prevNeighbor->getAlpha<1>();
    auto d2prevNeighborPrev = d2prevNeighbor->getAlpha<1>();

    d1->connectAlpha<1>(d2prev);
    d2->connectAlpha<1>(d1prev);

    d1neighbor->connectAlpha<2>(d2prevNeighborPrev);
    d2neighbor->connectAlpha<2>(d1prevNeighborPrev);

}


std::vector<Dart*> HexExtractor::getAllDartsInCell(Dart* d)
{
    auto dartsToBeProcessed = std::vector<Dart*>();
    auto darts = std::set<Dart*>();
    dartsToBeProcessed.push_back(d);
    while (!dartsToBeProcessed.empty())
    {
        Dart* d = dartsToBeProcessed.back();
        dartsToBeProcessed.resize(dartsToBeProcessed.size()-1);

        darts.insert(d);

        Dart* d2 = d->getAlpha<0>();
        if (d2 != nullptr)
            if (darts.find(d2) == darts.end())
                dartsToBeProcessed.push_back(d2);

        d2 = d->getAlpha<1>();
        if (d2 != nullptr)
            if (darts.find(d2) == darts.end())
                dartsToBeProcessed.push_back(d2);

        d2 = d->getAlpha<2>();
        if (d2 != nullptr)
            if (darts.find(d2) == darts.end())
                dartsToBeProcessed.push_back(d2);

    }

    auto res = std::vector<Dart*>();

    res.insert(res.end(), darts.begin(), darts.end());

    return res;
}

void HexExtractor::calculateDifferencesInDartTypes(CellHandle ch)
{

    auto cell = intermediateHexMesh.cell(ch);
    auto& halffaces = cell.halffaces();
    if (halffaces.empty())
    {
        HEXEX_DEBUG_ONLY(std::cout << "error: face empty" << std::endl;)
        return;
    }

    auto& facedarts = halffaceDarts[halffaces.front()];
    if (facedarts.empty())
    {
        HEXEX_DEBUG_ONLY(std::cout << "error: no darts in face" << std::endl;)
        return;
    }

    auto darts = getAllDartsInCell(facedarts.front());

    auto num_inverted = 0;
    auto num_proper   = 0;

    for (auto d : darts)
    {
        if (d->getHalfface().is_valid())
            if (intermediateHexMesh.incident_cell(d->getHalfface()) != ch)
            {
                HEXEX_DEBUG_ONLY(std::cout << "ERROR :( cell id: " << ch << " incident cell: " << intermediateHexMesh.incident_cell(d->getHalfface()).idx() << std::endl);
            }

        if (d->isAnti())
          ++num_inverted;
        else
          ++num_proper;
    }


  differenceBetweenInvertedAndProperDartsPerCell[ch] = num_inverted - num_proper;

}

void HexExtractor::calculateDifferencesInDartTypes(HalfFaceHandle hfh)
{

    auto dartsToBeProcessed = std::vector<Dart*>();
    auto darts = std::set<Dart*>();

    if(!halffaceDarts[hfh].empty())
        dartsToBeProcessed.push_back(halffaceDarts[hfh].front());

    while (!dartsToBeProcessed.empty())
    {
        Dart* d = dartsToBeProcessed.back();
        dartsToBeProcessed.resize(dartsToBeProcessed.size()-1);

        darts.insert(d);

        Dart* d2 = d->getAlpha<0>();
        if (d2 != nullptr)
            if (darts.find(d2) == darts.end())
                dartsToBeProcessed.push_back(d2);

        d2 = d->getAlpha<1>();
        if (d2 != nullptr)
            if (darts.find(d2) == darts.end())
                dartsToBeProcessed.push_back(d2);

    }

    auto num_inverted = 0;
    auto num_proper   = 0;

    for (auto d : darts)
        if (d->isAnti())
          ++num_inverted;
        else
          ++num_proper;


    differenceBetweenInvertedAndProperDartsPerHalfface[hfh] = num_inverted - num_proper;

}

void HexExtractor::calculateDifferencesInDartTypes()
{
  for (auto ch : intermediateHexMesh.cells())
  {
    calculateDifferencesInDartTypes(ch);
  }

  for (auto hfh : intermediateHexMesh.halffaces())
  {
    calculateDifferencesInDartTypes(hfh);
  }

}


bool HexExtractor::isCellDegenerate(CellHandle ch)
{
    if (!cellTypesComputed)
        computeCellTypes();

    ++isCellDegenerateCalls;

    return cellTypes[ch] == Degenerate;
}

bool HexExtractor::isCellFlipped(CellHandle ch)
{

    if (!cellTypesComputed)
        computeCellTypes();

    ++isCellFlippedCalls;

    return cellTypes[ch] == Flipped;
}

bool HexExtractor::isFaceDegenerate(HalfFaceHandle hfh)
{
    if (faceTypes[inputMesh.face_handle(hfh)] == NotComputed)
    {
        ++isFaceDegenerateCalls;

        if (!inputMesh.incident_cell(hfh).is_valid())
            return false;
        auto params = getParameters(hfh);
        if (isOnLine(params[0], params[1], params[2]))
            faceTypes[inputMesh.face_handle(hfh)] = Degenerate;
        else
            faceTypes[inputMesh.face_handle(hfh)] = Proper;


    }

    return faceTypes[inputMesh.face_handle(hfh)] == Degenerate;
}

bool HexExtractor::areColinear(CellHandle ch, HalfEdgeHandle heh, Direction dir)
{
    auto e = inputMesh.halfedge(heh);
    auto param1 = parameter(ch, e.from_vertex());
    auto param2 = parameter(ch, e.to_vertex());
    return areCoLinear(param1, param2, dir);
}

bool HexExtractor::isSingularVertex(VertexHandle vh)
{
    for (auto voh_it = inputMesh.voh_iter(vh); voh_it.valid(); ++voh_it)
    {
        auto eh = inputMesh.edge_handle(*voh_it);
        if (isSingularEdge(eh))
            return true;
    }
    return false;
}

bool HexExtractor::isSingularEdge(EdgeHandle eh)
{
    if (!edgeSingularitiesCalculated)
        calculateEdgeSingularities();
    return edgeSingularity[eh];
}

bool HexExtractor::isFixPointRecursive(Parameter param, CellHandle ch, VertexHandle vh, std::set<CellHandle> visited)
{
    if ((visited.count(ch) > 0) || visited.size() > 15)
        return param == parameter(ch, vh);

    visited.insert(ch);

    auto halffaces = getIncidentHalfFaces(inputMesh, ch, vh);
    for (auto hfh : halffaces)
    {
        auto oppHfh = inputMesh.opposite_halfface_handle(hfh);
        auto oppCh = inputMesh.incident_cell(oppHfh);
        if (!oppCh.is_valid())
            continue;

        auto newParameter = getTransitionFunction(hfh).transform_point(param);
        if (!isFixPointRecursive(newParameter, oppCh, vh, visited))
            return false;
    }

    return true;
}

bool HexExtractor::isFixPoint(Parameter parameter, CellHandle ch, VertexHandle vh)
{
    return isFixPointRecursive(parameter, ch, vh, std::set<CellHandle>());
}

HalfFaceHandle HexExtractor::rotateAroundHalfedge(CellHandle startCell, HalfEdgeHandle currentEdge, bool ccw)
{
    if (ccw)
        currentEdge = inputMesh.opposite_halfedge_handle(currentEdge);

    for (auto hehf_it = inputMesh.hehf_iter(currentEdge); hehf_it.valid(); ++hehf_it)
        if (startCell == inputMesh.incident_cell(*hehf_it))
            return *hehf_it;

    assert(false);
    return HalfFaceHandle();
}

void HexExtractor::doTransition(HalfFaceHandle hfh, CellHandle& ch)
{
    ch = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(hfh));
}

void HexExtractor::doTransition(HalfFaceHandle hfh, std::vector<Parameter>& params)
{
    auto ch = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(hfh));
    if (ch.is_valid())
        params = getParameters(ch);
}

void HexExtractor::doTransition(HalfFaceHandle hfh, Parameter& parameter)
{
    auto tranFun = getTransitionFunction(hfh);
    parameter = tranFun.transform_point(parameter);
}

void HexExtractor::doTransition(HalfFaceHandle hfh, Direction& dir)
{
    auto tranFun = getTransitionFunction(hfh);
    dir.transform(tranFun);
}

void HexExtractor::doTransition(HalfFaceHandle hfh, Transition& tranFun)
{
    auto additionalTranFun = getTransitionFunction(hfh);
    tranFun = additionalTranFun * tranFun;
}

Matrix4x4dd HexExtractor::transitionFrame(Parameter u, Parameter v, Parameter w)
{
    auto d0 = (v-u).normalized();
    auto d1 = (w-u).normalized();
    auto d2 = (d0%d1).normalized();
    d1 = (d2%d0).normalized();

    Matrix4x4dd frame;
    frame(0,0) = d0[0]; frame(0,1) = d1[0]; frame(0,2) = d2[0]; frame(0,3) = u[0];
    frame(1,0) = d0[1]; frame(1,1) = d1[1]; frame(1,2) = d2[1]; frame(1,3) = u[1];
    frame(2,0) = d0[2]; frame(2,1) = d1[2]; frame(2,2) = d2[2]; frame(2,3) = u[2];
    frame(3,0) =     0; frame(3,1) =     0; frame(3,2) =     0; frame(3,3) =    1;

    return frame;
}

} // namespace HexEx
