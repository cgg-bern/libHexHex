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

#include <vector>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

#include "Direction.hh"
#include "Typedefs.hh"

namespace HexEx
{

enum HPortType
{
    CellHPort,
    FaceHPort,
    EdgeHPort
};

std::ostream& operator<<(std::ostream& os, HPortType t);

class HPort;

class HPortHandle
{
public:
    explicit HPortHandle(std::shared_ptr<HPort> hport)
        :
          hport(hport)
    {}

    HPortHandle()
        :
          hport(nullptr)
    {}

    Direction& dir();
    const Direction& dir() const;
    Parameter& parameter();
    const Parameter& parameter() const;
    VertexHandle& vertex();
    const VertexHandle& vertex() const;
    CellHandle& cell();
    const CellHandle& cell() const;
    HPortType& type();

    EdgeHandle incidentEdgeHandle();
    HalfFaceHandle incidentHalfFaceHandle();
    CellHandle incidentCellHandle();

    HalfEdgeHandle& halfedge();

    bool isValid() const { return hport.get() != nullptr; }

    bool operator==(const HPortHandle& other) const;
    bool operator!=(const HPortHandle& other) const;

#ifdef DEBUG
public:
#else
private:
#endif
    std::shared_ptr<HPort> hport;
};

std::ostream& operator<<(std::ostream& os, const HPortHandle& hph);

class HPort
{
    friend class HPortHandle;
public:
    HPort(Direction dir, Parameter parameter, VertexHandle vertex, CellHandle cell, HPortType type, int incidentElementID);

private:

    Direction dir;
    Parameter parameter;
    VertexHandle vertex;
    CellHandle cell;
    HPortType type;
    int incidentElementID;

    HalfEdgeHandle halfedge;

    bool operator==(const HPort& other) const;

};


}
