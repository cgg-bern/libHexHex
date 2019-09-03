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


#include "HPort.hh"
#include "Utils.hh"

using namespace HexEx;
using namespace OpenVolumeMesh;

namespace HexEx
{

    std::ostream& operator<<(std::ostream& os, HPortType t)
    {
        switch (t)
        {
        case(HPortType::CellHPort):
            os << "CellHPort";
            break;
        case(HPortType::FaceHPort):
            os << "FaceHPort";
            break;
        case(HPortType::EdgeHPort):
            os << "EdgeHPort";
            break;
        }
        return os;
    }


HPort::HPort(Direction dir, Parameter parameter, VertexHandle vertex, CellHandle cell, HPortType type, int incidentElementID)
    :
      dir(dir),
      parameter(parameter),
      vertex(vertex),
      cell(cell),
      type(type),
      incidentElementID(incidentElementID)
{
    assert(parameter == roundVector(parameter));
#ifdef DEBUG
    if (dir.vector().length() < 0.5)
        std::cout << "warning: bad direction in " << __FUNCTION__ << std::endl;
#endif

}

bool HPort::operator==(const HPort& other) const
{
    return dir == other.dir &&
           parameter == other.parameter &&
           vertex == other.vertex &&
           cell == other.cell &&
           incidentElementID == other.incidentElementID;
}

Direction& HPortHandle::dir()
{
    return hport->dir;
}

const Direction& HPortHandle::dir() const
{
    return hport->dir;
}

Parameter& HPortHandle::parameter()
{
    return hport->parameter;
}

const Parameter&HPortHandle::parameter() const
{
    return hport->parameter;
}

VertexHandle& HPortHandle::vertex()
{
    return hport->vertex;
}

const VertexHandle& HPortHandle::vertex() const
{
    return hport->vertex;
}

CellHandle& HPortHandle::cell()
{
    return hport->cell;
}

const CellHandle& HPortHandle::cell() const
{
    return hport->cell;
}

HPortType& HPortHandle::type()
{
    return hport->type;
}

EdgeHandle HPortHandle::incidentEdgeHandle()
{
    assert(type() == EdgeHPort);
    return EdgeHandle(hport->incidentElementID);
}

HalfFaceHandle HPortHandle::incidentHalfFaceHandle()
{
    assert(type() == FaceHPort);
    return HalfFaceHandle(hport->incidentElementID);
}

CellHandle HPortHandle::incidentCellHandle()
{
    assert(type() == CellHPort);
    return CellHandle(hport->incidentElementID);
}

HalfEdgeHandle& HPortHandle::halfedge()
{
    return hport->halfedge;
}

bool HPortHandle::operator ==(const HPortHandle& other) const
{
    if (!isValid() || !other.isValid())
        return isValid() == other.isValid();
    return *hport == *other.hport;
}

bool HPortHandle::operator!=(const HPortHandle& other) const
{
    return !operator==(other);
}

std::ostream& operator<<(std::ostream& os, const HPortHandle& hph)
{
    if (!hph.isValid())
        return os << "HPort(invalid)";
    else
        return os << "HPort(dir: " << hph.dir() << ", param: " << hph.parameter() << ", " << hph.cell() << ", " << hph.vertex() << ")";
}

}

