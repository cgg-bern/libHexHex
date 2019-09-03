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


#include "Dart.hh"
#include "HexExtractor.hh"


namespace HexEx
{


Dart::Dart(HexExtractor& hexExtractor, HPortHandle tracePort, Direction refDir, Direction normalDir)
    :
      he(hexExtractor),
      tracePort(tracePort),
      refDir(refDir),
      normalDir(normalDir),
      alphas({nullptr, nullptr, nullptr, nullptr}),
      halfedge(HalfEdgeHandle()),
      halfface(HalfFaceHandle()),
      uv(Parameter(0,0,0))
{

  static int global_id = 0;
  id = global_id++;
    assert(hexExtractor.isDartInCell(tracePort.cell(), tracePort, refDir,normalDir));
}

CellHandle Dart::getHexCell()
{
    if (!halfface.is_valid())
        return CellHandle();

    return he.intermediateHexMesh.incident_cell(halfface);
}

bool Dart::operator==(const Dart& other)
{
    return tracePort == other.tracePort &&
           refDir == other.refDir &&
           normalDir == other.normalDir;
}


template <>
Transition Dart::getTransitionToAlpha<0>()
{
    auto d1 = getTraceDir();
    auto d2 = getRefDir();
    auto d3 = getNormalDir();
    auto p1 = getParameter();

    auto other = getAlpha<0>();

    auto d4 = other->getTraceDir().flipped();
    auto d5 = other->getRefDir();
    auto d6 = other->getNormalDir();
    auto p2 = other->getParameter() - d4;

    if (isAnti() != other->isAnti())
    {
        d4 = other->getTraceDir();
        p2 = other->getParameter();
    }

    return getTransitionBetweenFrames(d1,d2,d3,p1,d4,d5,d6,p2);
}

template <>
Transition Dart::getTransitionToAlpha<1>()
{
    auto d1 = getTraceDir();
    auto d2 = getRefDir();
    auto d3 = getNormalDir();
    auto p1 = getParameter();

    auto other = getAlpha<1>();

    auto d4 = other->getRefDir();
    auto d5 = other->getTraceDir();
    auto d6 = other->getNormalDir();
    auto p2 = other->getParameter();

    if (isAnti() != other->isAnti())
        std::swap(d4,d5);

    return getTransitionBetweenFrames(d1,d2,d3,p1,d4,d5,d6,p2);
}

template <>
Transition Dart::getTransitionToAlpha<2>()
{
    auto d1 = getTraceDir();
    auto d2 = getRefDir();
    auto d3 = getNormalDir();
    auto p1 = getParameter();

    auto other = getAlpha<2>();

    auto d4 = other->getTraceDir();
    auto d5 = other->getNormalDir();
    auto d6 = other->getRefDir();
    auto p2 = other->getParameter();

    if (isAnti() != other->isAnti())
        std::swap(d5,d6);

    return getTransitionBetweenFrames(d1,d2,d3,p1,d4,d5,d6,p2);
}

template <>
Transition Dart::getTransitionToAlpha<3>()
{
    auto d1 = getTraceDir();
    auto d2 = getRefDir();
    auto d3 = getNormalDir();
    auto p1 = getParameter();

    auto other = getAlpha<3>();

    auto d4 = other->getTraceDir();
    auto d5 = other->getRefDir();
    auto d6 = other->getNormalDir().flipped();
    auto p2 = other->getParameter();

    if (isAnti() != other->isAnti())
        d6.flip();

    return getTransitionBetweenFrames(d1,d2,d3,p1,d4,d5,d6,p2);
}

bool Dart::isAnti()
{
    return he.isCellFlipped(tracePort.cell());
}

Transition Dart::getTransitionBetweenFrames(Direction d1, Direction d2, Direction d3, Parameter p1,
                                            Direction d4, Direction d5, Direction d6, Parameter p2)
{
    double entries[] = {d1.vector()[0], d1.vector()[1], d1.vector()[2], 0,
                        d2.vector()[0], d2.vector()[1], d2.vector()[2], 0,
                        d3.vector()[0], d3.vector()[1], d3.vector()[2], 0,
                                 p1[0],          p1[1],          p1[2], 1 };
    auto t1 = Transition(Matrix4x4d(entries));

    double entries2[] = {d4.vector()[0], d4.vector()[1], d4.vector()[2], 0,
                         d5.vector()[0], d5.vector()[1], d5.vector()[2], 0,
                         d6.vector()[0], d6.vector()[1], d6.vector()[2], 0,
                                  p2[0],          p2[1],          p2[2], 1 };
    auto t2 = Transition(Matrix4x4d(entries2));


    t1.invert();

    return t2*t1;


}

std::istream& deserialize(std::istream& _istr, Dart&)
{
  return _istr;
}


}
