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

namespace HexEx
{

class HexExtractor;

class Dart
{
public:
    Dart(HexExtractor& hexExtractor, HPortHandle tracePort, Direction refDir, Direction normalDir);

    HPortHandle getTracePort()  { return tracePort;  }

    Direction getTraceDir()  const { return tracePort.dir(); }
    Direction getRefDir()    const  { return refDir; }
    Direction getNormalDir() const  { return normalDir; }

    CellHandle getCell() const { return tracePort.cell(); }
    CellHandle getHexCell();

    Parameter getParameter() const { return tracePort.parameter(); }

    bool operator==(const Dart& other);

    template <int i>
    void connectAlpha(Dart* dart)
    {
        static_assert(i >= 0 && i <= 3, "wrong i");

        assert(this->isPrimary());
        assert(dart->isSecondary());
        assert(dart != nullptr);

        alphas[i] = dart;
        dart->alphas[i] = this;
    }

    template <int i>
    void disconnectAlpha()
    {
        static_assert(i >= 0 && i <= 3, "wrong i");

        if (alphas[i] != nullptr)
            if (alphas[i]->alphas[i] == this)
                alphas[i]->alphas[i] = nullptr;

        alphas[i] = nullptr;
    }

    template <int i>
    Transition getTransitionToAlpha();


    template <int i>
    Dart* getAlpha()
    {
        static_assert(i >= 0 && i <= 3, "wrong i");
        return alphas[i];
    }


    VertexHandle& getVertex() { return tracePort.vertex(); }
    HalfEdgeHandle& getHalfedge()  { return halfedge; }
    HalfFaceHandle& getHalfface() { return halfface; }

    bool isRighthanded() { return getTraceDir() % getRefDir() == getNormalDir(); }
    bool isPrimary() { return !isAnti() ? isRighthanded() : !isRighthanded(); }
    bool isSecondary() { return !isPrimary(); }
    bool isAnti();

    void setNeighbor(Dart* neighborDart, Transition transitionToNeighbor);

    Parameter& getUV() { return uv; }

private:

    Transition getTransitionBetweenFrames(Direction d1, Direction d2, Direction d3, Parameter p1,
                                          Direction d4, Direction d5, Direction d6, Parameter p2);

    HexExtractor& he;

    HPortHandle tracePort;
    Direction refDir;
    Direction normalDir;

    std::array<Dart*, 4> alphas;

    HalfEdgeHandle halfedge;
    HalfFaceHandle halfface;

    Parameter uv;
public:
    int id;

};

}
