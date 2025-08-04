#pragma once

#include <HexHex/Commons/CellIGM.hh>
#include <HexHex/Utils/Utils.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Commons/Matrix4x4T.hh>

bool check_that_all_halffaces_referenced_by_cells_exists(HexHex::TetrahedralMesh& mesh);

HexHex::Vec3d getRandomVector(int absmax = 10);

HexHex::Vec3d getRandomVectorDouble();

HexHex::Matrix4x4d getRandomMatrix();

inline void createCube(HexHex::TetrahedralMesh& mesh, double parametrizationSideLength)
{
    using VertexHandle = OpenVolumeMesh::VertexHandle;
    using CellHandle   = OpenVolumeMesh::CellHandle;

    std::vector<VertexHandle> vhs;
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0, 1.0)));

    std::vector<CellHandle> chs;
    std::vector<VertexHandle> c;

    c.clear();
    c.push_back(vhs[0]);
    c.push_back(vhs[3]);
    c.push_back(vhs[1]);
    c.push_back(vhs[4]);
    chs.push_back(mesh.add_cell(c));

    c.clear();
    c.push_back(vhs[4]);
    c.push_back(vhs[6]);
    c.push_back(vhs[7]);
    c.push_back(vhs[3]);
    chs.push_back(mesh.add_cell(c));

    c.clear();
    c.push_back(vhs[4]);
    c.push_back(vhs[5]);
    c.push_back(vhs[6]);
    c.push_back(vhs[1]);
    chs.push_back(mesh.add_cell(c));

    c.clear();
    c.push_back(vhs[1]);
    c.push_back(vhs[3]);
    c.push_back(vhs[2]);
    c.push_back(vhs[6]);
    chs.push_back(mesh.add_cell(c));

    c.clear();
    c.push_back(vhs[3]);
    c.push_back(vhs[4]);
    c.push_back(vhs[6]);
    c.push_back(vhs[1]);
    chs.push_back(mesh.add_cell(c));

    HexHex::IGM parametrization = mesh.request_cell_property<HexHex::CellIGM>("Parametrization");
    mesh.set_persistent(parametrization);

    HexHex::Vec3d offset = HexHex::Vec3d(0.0,0.0,0.0);
    double scale = parametrizationSideLength;

    for (auto ch : chs)
    {
        auto vertices = mesh.get_cell_vertices(ch);
        for (auto vh : vertices)
        {
            auto pos = mesh.vertex(vh);
            HexHex::Vec3d pos2 = HexHex::Vec3d(pos[0], pos[1], pos[2]);
            parametrization[ch][vh] = pos2*scale+offset;
        }
    }
}


inline void createCube(HexHex::TetrahedralMesh& mesh)
{
    using VertexHandle = OpenVolumeMesh::VertexHandle;
    using HalfFaceHandle = OpenVolumeMesh::HalfFaceHandle;

    std::vector<VertexHandle> vhs;
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0, 1.0)));

    std::vector<HalfFaceHandle> hfhs;

    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[0],vhs[1],vhs[2],vhs[3]}), 0));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[7],vhs[6],vhs[5],vhs[4]}), 0));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[1],vhs[0],vhs[4],vhs[5]}), 0));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[2],vhs[1],vhs[5],vhs[6]}), 0));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[3],vhs[2],vhs[6],vhs[7]}), 0));
    hfhs.push_back(mesh.halfface_handle(mesh.add_face({vhs[0],vhs[3],vhs[7],vhs[4]}), 0));

    mesh.add_cell(hfhs);

}


inline void createCubeWithTransition(HexHex::TetrahedralMesh& mesh, double parametrizationSideLength)
{
    using VertexHandle = OpenVolumeMesh::VertexHandle;
    using CellHandle   = OpenVolumeMesh::CellHandle;

    std::vector<VertexHandle> vhs;
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0, 0.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0, 1.0)));
    vhs.push_back(mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0, 1.0)));

    std::vector<CellHandle> chs;
    std::vector<VertexHandle> c;

    c.clear();
    c.push_back(vhs[0]);
    c.push_back(vhs[3]);
    c.push_back(vhs[1]);
    c.push_back(vhs[4]);
    chs.push_back(mesh.add_cell(c));

    c.clear();
    c.push_back(vhs[4]);
    c.push_back(vhs[6]);
    c.push_back(vhs[7]);
    c.push_back(vhs[3]);
    chs.push_back(mesh.add_cell(c));

    c.clear();
    c.push_back(vhs[4]);
    c.push_back(vhs[5]);
    c.push_back(vhs[6]);
    c.push_back(vhs[1]);
    chs.push_back(mesh.add_cell(c));

    c.clear();
    c.push_back(vhs[1]);
    c.push_back(vhs[3]);
    c.push_back(vhs[2]);
    c.push_back(vhs[6]);
    chs.push_back(mesh.add_cell(c));

    c.clear();
    c.push_back(vhs[3]);
    c.push_back(vhs[4]);
    c.push_back(vhs[6]);
    c.push_back(vhs[1]);
    chs.push_back(mesh.add_cell(c));

    HexHex::IGM parametrization = mesh.request_cell_property<HexHex::CellIGM>("Parametrization");
    mesh.set_persistent(parametrization);

    HexHex::Vec3d offset = HexHex::Vec3d(2.0,0.0,0.0);
    double scale = parametrizationSideLength;

    double heightThreshold = 0.5;
    HexHex::Matrix4x4d transform = HexHex::Matrix4x4d();
    transform(0,0) = 0; transform(0,1) = 1; transform(0,2) = 0; transform(0,3) = 3;
    transform(1,0) = 0; transform(1,1) = 0; transform(1,2) = 1; transform(1,3) = 3;
    transform(2,0) = 1; transform(2,1) = 0; transform(2,2) = 0; transform(2,3) = 3;
    transform(3,0) = 0; transform(3,1) = 0; transform(3,2) = 0; transform(3,3) = 1;

    for (auto ch : chs)
    {
        auto pos = mesh.barycenter(ch);
        if (pos[1] > heightThreshold)
        {
            auto vertices = mesh.get_cell_vertices(ch);
            for (auto vh : vertices)
            {
                auto pos = mesh.vertex(vh);
                HexHex::Vec3d pos2 = HexHex::Vec3d(pos[0], pos[1], pos[2]);
                parametrization[ch][vh] = transform.transform_point(pos2*scale+offset);
            }
        }
        else
        {
            auto vertices = mesh.get_cell_vertices(ch);
            for (auto vh : vertices)
            {
                auto pos = mesh.vertex(vh);
                HexHex::Vec3d pos2 = HexHex::Vec3d(pos[0], pos[1], pos[2]);
                parametrization[ch][vh] = pos2*scale+offset;
            }
        }
    }

}

inline OpenVolumeMesh::CellHandle addCube(OpenVolumeMesh::VertexHandle vh0, OpenVolumeMesh::VertexHandle vh1,
                                   OpenVolumeMesh::VertexHandle vh2, OpenVolumeMesh::VertexHandle vh3,
                                   OpenVolumeMesh::VertexHandle vh4, OpenVolumeMesh::VertexHandle vh5,
                                   OpenVolumeMesh::VertexHandle vh6, OpenVolumeMesh::VertexHandle vh7,
                                   HexHex::TetrahedralMesh& mesh, bool other)
{
    using CellHandle   = OpenVolumeMesh::CellHandle;
    if (!other)
    {
        CellHandle ret = mesh.add_cell(vh0, vh1, vh3, vh4);
        mesh.add_cell(vh1, vh2, vh3, vh6);
        mesh.add_cell(vh1, vh6, vh3, vh4);
        mesh.add_cell(vh1, vh6, vh4, vh5);
        mesh.add_cell(vh4, vh7, vh6, vh3);
        return ret;
    }
    else
    {
        CellHandle ret = mesh.add_cell(vh0, vh1, vh2, vh5);
        mesh.add_cell(vh0, vh7, vh4, vh5);
        mesh.add_cell(vh0, vh5, vh2, vh7);
        mesh.add_cell(vh0, vh3, vh7, vh2);
        mesh.add_cell(vh5, vh7, vh6, vh2);
        return ret;
    }
}


inline void createMasterVertexMesh(HexHex::TetrahedralMesh& mesh, bool centerOffset = false)
{
    using VertexHandle = OpenVolumeMesh::VertexHandle;

    VertexHandle luh = mesh.add_vertex(HexHex::Vec3d(-1.0,-1.5,-1.0));
    VertexHandle ruh = mesh.add_vertex(HexHex::Vec3d( 1.0,-1.0,-1.5));
    VertexHandle loh = mesh.add_vertex(HexHex::Vec3d(-1.0, 1.5,-1.0));
    VertexHandle roh = mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0,-1.5));
    VertexHandle luv = mesh.add_vertex(HexHex::Vec3d(-1.0,-1.5, 1.0));
    VertexHandle ruv = mesh.add_vertex(HexHex::Vec3d( 1.0,-1.0, 1.5));
    VertexHandle lov = mesh.add_vertex(HexHex::Vec3d(-1.0, 1.5, 1.0));
    VertexHandle rov = mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0, 1.5));

    VertexHandle lmhu = mesh.add_vertex(HexHex::Vec3d(-1.0,-0.5,-1.0));
    VertexHandle rmh  = mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0,-1.5));
    VertexHandle lmvu = mesh.add_vertex(HexHex::Vec3d(-1.0,-0.5, 1.0));
    VertexHandle rmv  = mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0, 1.5));

    VertexHandle lmho = mesh.add_vertex(HexHex::Vec3d(-1.0, 0.5,-1.0));
    VertexHandle lmvo = mesh.add_vertex(HexHex::Vec3d(-1.0, 0.5, 1.0));

    VertexHandle muh = mesh.add_vertex(HexHex::Vec3d( 0.0,-1.0,-1.0));
    VertexHandle moh = mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0,-1.0));
    VertexHandle muv = mesh.add_vertex(HexHex::Vec3d( 0.0,-1.0, 1.0));
    VertexHandle mov = mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0, 1.0));

    VertexHandle lum  = mesh.add_vertex(HexHex::Vec3d(-1.0,-1.5, 0.0));
    VertexHandle rumh = mesh.add_vertex(HexHex::Vec3d( 1.0,-1.0,-0.5));
    VertexHandle rumv = mesh.add_vertex(HexHex::Vec3d( 1.0,-1.0, 0.5));
    VertexHandle lom  = mesh.add_vertex(HexHex::Vec3d(-1.0, 1.5, 0.0));
    VertexHandle romh = mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0,-0.5));
    VertexHandle romv = mesh.add_vertex(HexHex::Vec3d( 1.0, 1.0, 0.5));

    VertexHandle mmh = mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0,-1.0));
    VertexHandle mmv = mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0, 1.0));

    VertexHandle lmmu = mesh.add_vertex(HexHex::Vec3d(-1.0,-0.5, 0.0));
    VertexHandle lmmo = mesh.add_vertex(HexHex::Vec3d(-1.0, 0.5, 0.0));
    VertexHandle rmmh = mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0,-0.5));
    VertexHandle rmmv = mesh.add_vertex(HexHex::Vec3d( 1.0, 0.0, 0.5));

    VertexHandle mum = mesh.add_vertex(HexHex::Vec3d( 0.0,-1.0, 0.0));
    VertexHandle mom = mesh.add_vertex(HexHex::Vec3d( 0.0, 1.0, 0.0));

    VertexHandle mmm = mesh.add_vertex(HexHex::Vec3d( 0.0, 0.0, 0.0));

    addCube(luv, muv, mum, lum, lmvu, mmv, mmm, lmmu, mesh, false);
    addCube(muv, ruv, rumv, mum, mmv, rmv, rmmv, mmm, mesh, true);
    addCube(lum, mum, muh, luh, lmmu, mmm, mmh, lmhu, mesh, true);
    addCube(mum, rumh, ruh, muh, mmm, rmmh, rmh, mmh, mesh, false);
    addCube(lmvo, mmv, mmm, lmmo, lov, mov, mom, lom, mesh, true);
    addCube(mmv, rmv, rmmv, mmm, mov, rov, romv, mom, mesh, false);
    addCube(lmmo, mmm, mmh, lmho, lom, mom, moh, loh, mesh, false);
    addCube(mmm, rmmh, rmh, mmh, mom, romh, roh, moh, mesh, true);

    HexHex::IGM parametrization = mesh.request_cell_property<HexHex::CellIGM>("Parametrization");
    mesh.set_persistent(parametrization);

    HexHex::Vec3d offset = HexHex::Vec3d(0.0,0.0,0.0);

    HexHex::Matrix4x4d transformleftshiftup = HexHex::Matrix4x4d();
    transformleftshiftup(0,0) = 1; transformleftshiftup(0,1) = 0; transformleftshiftup(0,2) = 0; transformleftshiftup(0,3) = 0;
    transformleftshiftup(1,0) =-0.5; transformleftshiftup(1,1) = 1; transformleftshiftup(1,2) = 0; transformleftshiftup(1,3) = 0;
    transformleftshiftup(2,0) = 0; transformleftshiftup(2,1) = 0; transformleftshiftup(2,2) = 1; transformleftshiftup(2,3) = 0;
    transformleftshiftup(3,0) = 0; transformleftshiftup(3,1) = 0; transformleftshiftup(3,2) = 0; transformleftshiftup(3,3) = 1;

    HexHex::Matrix4x4d transformleftshiftdown = HexHex::Matrix4x4d();
    transformleftshiftdown(0,0) = 1; transformleftshiftdown(0,1) = 0; transformleftshiftdown(0,2) = 0; transformleftshiftdown(0,3) = 0;
    transformleftshiftdown(1,0) =0.5; transformleftshiftdown(1,1) = 1; transformleftshiftdown(1,2) = 0; transformleftshiftdown(1,3) = 0;
    transformleftshiftdown(2,0) = 0; transformleftshiftdown(2,1) = 0; transformleftshiftdown(2,2) = 1; transformleftshiftdown(2,3) = 0;
    transformleftshiftdown(3,0) = 0; transformleftshiftdown(3,1) = 0; transformleftshiftdown(3,2) = 0; transformleftshiftdown(3,3) = 1;

    HexHex::Matrix4x4d transformrightshiftback = HexHex::Matrix4x4d();
    transformrightshiftback(0,0) = 1; transformrightshiftback(0,1) = 0; transformrightshiftback(0,2) = 0; transformrightshiftback(0,3) = 0;
    transformrightshiftback(1,0) = 0; transformrightshiftback(1,1) = 1; transformrightshiftback(1,2) = 0; transformrightshiftback(1,3) = 0;
    transformrightshiftback(2,0) =-0.5; transformrightshiftback(2,1) = 0; transformrightshiftback(2,2) = 1; transformrightshiftback(2,3) = 0;
    transformrightshiftback(3,0) = 0; transformrightshiftback(3,1) = 0; transformrightshiftback(3,2) = 0; transformrightshiftback(3,3) = 1;

    HexHex::Matrix4x4d transformrightshiftfront = HexHex::Matrix4x4d();
    transformrightshiftfront(0,0) = 1; transformrightshiftfront(0,1) = 0; transformrightshiftfront(0,2) = 0; transformrightshiftfront(0,3) = 0;
    transformrightshiftfront(1,0) = 0; transformrightshiftfront(1,1) = 1; transformrightshiftfront(1,2) = 0; transformrightshiftfront(1,3) = 0;
    transformrightshiftfront(2,0) = 0.5; transformrightshiftfront(2,1) = 0; transformrightshiftfront(2,2) = 1; transformrightshiftfront(2,3) = 0;
    transformrightshiftfront(3,0) = 0; transformrightshiftfront(3,1) = 0; transformrightshiftfront(3,2) = 0; transformrightshiftfront(3,3) = 1;

    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
    {
        auto ch = *c_it;
        HexHex::Matrix4x4d transform;
        if (ch.idx() < 5)
            transform = transformleftshiftup;
        if ((ch.idx() > 9) && (ch.idx() < 15))
            transform = transformleftshiftup;
        if ((ch.idx() > 19) && (ch.idx() < 25))
            transform = transformleftshiftdown;
        if ((ch.idx() > 29) && (ch.idx() < 35))
            transform = transformleftshiftdown;
        if ((ch.idx() > 4) && (ch.idx() < 10))
            transform = transformrightshiftback;
        if ((ch.idx() > 14) && (ch.idx() < 20))
            transform = transformrightshiftfront;
        if ((ch.idx() > 24) && (ch.idx() < 30))
            transform = transformrightshiftback;
        if ((ch.idx() > 34) && (ch.idx() < 40))
            transform = transformrightshiftfront;

        auto vertices = mesh.get_cell_vertices(ch);
        for (auto vh : vertices)
        {
            auto pos = mesh.vertex(vh);
            HexHex::Vec3d pos2 = HexHex::Vec3d(pos[0], pos[1], pos[2]);
            parametrization[ch][vh] = transform.transform_point(pos2)+offset;
            if (centerOffset && (vh == mmm))
                parametrization[ch][vh] += 0.01*getRandomVector();
        }

    }


}

inline void createPrism(HexHex::TetrahedralMesh& mesh, bool centeroffset = false, bool smaller = false)
{
    using VertexHandle = OpenVolumeMesh::VertexHandle;

    auto vertices = std::vector<VertexHandle>();
    vertices.resize(38);

    for (auto i = 0u; i < 38; ++i)
        vertices[i] = mesh.add_vertex(HexHex::Vec3d(0,0,0));

    mesh.set_vertex(VertexHandle(0), HexHex::Vec3d(-1,-1,1));
    mesh.set_vertex(VertexHandle(18), HexHex::Vec3d(1,-1,1));
    mesh.set_vertex(VertexHandle(12), HexHex::Vec3d(0,-1,-1));
    mesh.set_vertex(VertexHandle(2), 0.5*(mesh.vertex(vertices[0]) + mesh.vertex(vertices[18])));
    mesh.set_vertex(VertexHandle(6), 0.5*(mesh.vertex(vertices[0]) + mesh.vertex(vertices[12])));
    mesh.set_vertex(VertexHandle(14), 0.5*(mesh.vertex(vertices[12]) + mesh.vertex(vertices[18])));
    mesh.set_vertex(VertexHandle(8), (1.0/3.0)*(mesh.vertex(vertices[6]) + mesh.vertex(vertices[14]) + mesh.vertex(vertices[2])));

    mesh.set_vertex(VertexHandle(3), 0.5*(mesh.vertex(vertices[0]) + mesh.vertex(vertices[6])));
    mesh.set_vertex(VertexHandle(9), 0.5*(mesh.vertex(vertices[6]) + mesh.vertex(vertices[12])));
    mesh.set_vertex(VertexHandle(13), 0.5*(mesh.vertex(vertices[12]) + mesh.vertex(vertices[14])));
    mesh.set_vertex(VertexHandle(15), 0.5*(mesh.vertex(vertices[14]) + mesh.vertex(vertices[18])));
    mesh.set_vertex(VertexHandle(17), 0.5*(mesh.vertex(vertices[18]) + mesh.vertex(vertices[2])));
    mesh.set_vertex(VertexHandle(1), 0.5*(mesh.vertex(vertices[2]) + mesh.vertex(vertices[0])));

    mesh.set_vertex(VertexHandle(5), 0.5*(mesh.vertex(vertices[2]) + mesh.vertex(vertices[8])));
    mesh.set_vertex(VertexHandle(7), 0.5*(mesh.vertex(vertices[6]) + mesh.vertex(vertices[8])));
    mesh.set_vertex(VertexHandle(11), 0.5*(mesh.vertex(vertices[14]) + mesh.vertex(vertices[8])));

    mesh.set_vertex(VertexHandle(4), 0.5*(mesh.vertex(vertices[3]) + mesh.vertex(vertices[5])));
    mesh.set_vertex(VertexHandle(10), 0.5*(mesh.vertex(vertices[9]) + mesh.vertex(vertices[11])));
    mesh.set_vertex(VertexHandle(16), 0.5*(mesh.vertex(vertices[15]) + mesh.vertex(vertices[5])));

    for (auto i = 0u; i < 19u; ++i)
        mesh.set_vertex(VertexHandle(i+19), mesh.vertex(vertices[i])+HexHex::Vec3d(0,2,0));

    int a,b,c,d;
    for (auto i : {0u,6u})
    {
        a = 0+i; b = 1+i; c = 4+i; d = 3+i;
        addCube(vertices[a], vertices[b], vertices[c], vertices[d], vertices[a+19], vertices[b+19], vertices[c+19], vertices[d+19], mesh, false);
        a = 1+i; b = 2+i; c = 5+i; d = 4+i;
        addCube(vertices[a], vertices[b], vertices[c], vertices[d], vertices[a+19], vertices[b+19], vertices[c+19], vertices[d+19], mesh, true);
        a = 3+i; b = 4+i; c = 7+i; d = 6+i;
        addCube(vertices[a], vertices[b], vertices[c], vertices[d], vertices[a+19], vertices[b+19], vertices[c+19], vertices[d+19], mesh, true);
        a = 4+i; b = 5+i; c = 8+i; d = 7+i;
        addCube(vertices[a], vertices[b], vertices[c], vertices[d], vertices[a+19], vertices[b+19], vertices[c+19], vertices[d+19], mesh, false);
    }

    a = 2; b = 17; c = 16; d = 5;
    addCube(vertices[a], vertices[b], vertices[c], vertices[d], vertices[a+19], vertices[b+19], vertices[c+19], vertices[d+19], mesh, false);
    a = 17; b = 18; c = 15; d = 16;
    addCube(vertices[a], vertices[b], vertices[c], vertices[d], vertices[a+19], vertices[b+19], vertices[c+19], vertices[d+19], mesh, true);
    a = 5; b = 16; c = 11; d = 8;
    addCube(vertices[a], vertices[b], vertices[c], vertices[d], vertices[a+19], vertices[b+19], vertices[c+19], vertices[d+19], mesh, true);
    a = 16; b = 15; c = 14; d = 11;
    addCube(vertices[a], vertices[b], vertices[c], vertices[d], vertices[a+19], vertices[b+19], vertices[c+19], vertices[d+19], mesh, false);


    auto params = std::vector<HexHex::Vec3d>(38,HexHex::Vec3d(-10,0,0));
//    params.resize(38);

    for (auto i = 0u; i < 15; ++i)
    {
        params[i] =    HexHex::Vec3d(-1+0.5*(i%3),-1,1-0.5*(i/3));
        params[i+19] = HexHex::Vec3d(-1+0.5*(i%3),1,1-0.5*(i/3));
    }

    auto params2 = std::vector<HexHex::Vec3d>(38,HexHex::Vec3d(-10,0,0));
//    params2.resize(38);
    params2[8] = HexHex::Vec3d(0,-1,0);
    params2[5] = HexHex::Vec3d(0.5,-1,0);
    params2[2] = HexHex::Vec3d(1,-1,0);
    params2[11] = HexHex::Vec3d(0,-1,-0.5);
    params2[16] = HexHex::Vec3d(0.5,-1,-0.5);
    params2[17] = HexHex::Vec3d(1,-1,-0.5);
    params2[14] = HexHex::Vec3d(0,-1,-1);
    params2[15] = HexHex::Vec3d(0.5,-1,-1);
    params2[18] = HexHex::Vec3d(1,-1,-1);

    for (auto i : {2,8,5,11,16,17,14,15,18})
        params2[i+19] = params2[i]+HexHex::Vec3d(0,2,0);


    if (smaller)
    {
        for (auto& p : params)
            if (p[1] > 0)
                p += HexHex::Vec3d(0,-1,0);
        for (auto& p : params2)
            if (p[1] > 0)
                p += HexHex::Vec3d(0,-1,0);

        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
          if (mesh.vertex(*v_it)[1] > 0)
            mesh.set_vertex(*v_it, mesh.vertex(*v_it) + HexHex::Vec3d(0,-1,0));
    }


    HexHex::IGM parameters = mesh.request_cell_property<HexHex::CellIGM>("Parametrization");
    mesh.set_persistent(parameters);

    HexHex::Vec3d offset = HexHex::Vec3d(0.0,0.0,0.0);

    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
    {
        auto ch = *c_it;
        auto vertices = mesh.get_cell_vertices(ch);
        for (auto vh : vertices)
        {
            if (ch.idx() < 8*5)
                parameters[ch][vh] = params[vh.idx()];
            else
                parameters[ch][vh] = params2[vh.idx()];
            if (centeroffset && ((vh.idx() == 8u) || (vh.idx() == 8u+19u)))
                parameters[ch][vh] += 0.01*getRandomVector();
        }

    }


}

inline void createCrazyCorner(HexHex::TetrahedralMesh& mesh, bool offsets)
{
    using VertexHandle = OpenVolumeMesh::VertexHandle;
    using CellHandle   = OpenVolumeMesh::CellHandle;

    auto vertices = std::vector<VertexHandle>();

    auto scale = 0.5;
    vertices.push_back(mesh.add_vertex(scale * HexHex::Vec3d(   0,   0,   0)));
    vertices.push_back(mesh.add_vertex(scale * HexHex::Vec3d(   0,   1,   1)));
    vertices.push_back(mesh.add_vertex(scale * HexHex::Vec3d(  -1,   0,   0)));
    vertices.push_back(mesh.add_vertex(scale * HexHex::Vec3d(-0.5,   0,   1)));
    vertices.push_back(mesh.add_vertex(scale * HexHex::Vec3d( 0.5,   0,   1)));
    vertices.push_back(mesh.add_vertex(scale * HexHex::Vec3d(   1,   0,   0)));
    vertices.push_back(mesh.add_vertex(scale * HexHex::Vec3d( 0.5,   1,   0)));
    vertices.push_back(mesh.add_vertex(scale * HexHex::Vec3d(-0.5,   1,   0)));

    auto cells = std::vector<CellHandle>();

    cells.push_back(mesh.add_cell(vertices[0], vertices[2], vertices[3], vertices[1]));
    cells.push_back(mesh.add_cell(vertices[0], vertices[3], vertices[4], vertices[1]));
    cells.push_back(mesh.add_cell(vertices[0], vertices[4], vertices[5], vertices[1]));
    cells.push_back(mesh.add_cell(vertices[0], vertices[5], vertices[6], vertices[1]));
    cells.push_back(mesh.add_cell(vertices[0], vertices[6], vertices[7], vertices[1]));
    cells.push_back(mesh.add_cell(vertices[0], vertices[7], vertices[2], vertices[1]));

    HexHex::IGM parametrization = mesh.request_cell_property<HexHex::CellIGM>("Parametrization");
    mesh.set_persistent(parametrization);

    for (auto ch : cells)
        for (auto vh : mesh.get_cell_vertices(ch))
            if (offsets)
            {
                if (vh == VertexHandle(3) || vh == VertexHandle(4))
                {
                    auto pos = mesh.vertex(VertexHandle(3+(1-(vh.idx()-3))));
                    HexHex::Vec3d pos2 = HexHex::Vec3d(pos[0], pos[1], pos[2]);
                    parametrization[ch][vh] = pos2;
                }
                else
                {
                    auto pos = mesh.vertex(vh);
                    HexHex::Vec3d pos2 = HexHex::Vec3d(pos[0], pos[1], pos[2]);
                    parametrization[ch][vh] = pos2;
                }
            }
            else
            {
                auto pos = mesh.vertex(vh);
                HexHex::Vec3d pos2 = HexHex::Vec3d(pos[0], pos[1], pos[2]);
                parametrization[ch][vh] = pos2;
            }


}


