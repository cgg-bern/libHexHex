
#include <HexHex/GlobalTopology/cell_connectivity.hh>
#include <HexHex/Predicates/ExactPredicates.hh>
#include <gtest/gtest.h>
#include <HexHex/HexExtractor.hh>
#include <random>

namespace HexHex
{

TEST(MeshElementTest, tetIncidencyTest)
{
    TetrahedralMesh mesh;
    VertexHandle vh0 = mesh.add_vertex(Vec3d(0,0,0));
    VertexHandle vh1 = mesh.add_vertex(Vec3d(0,0,1));
    VertexHandle vh2 = mesh.add_vertex(Vec3d(1,0,0));
    VertexHandle vh3 = mesh.add_vertex(Vec3d(0,1,0));
    CellHandle ch = mesh.add_cell({vh0,vh1,vh2,vh3});

    TetMeshCache cache(mesh);
    cache.build();

    ASSERT_EQ(MeshElement(ch), ch);

    for (auto cv_it = mesh.cv_iter(ch); cv_it.is_valid(); ++cv_it)
    {
        ASSERT_TRUE(MeshElement(*cv_it).is_incident(cache, ch));
    }

    for (auto ce_it = mesh.ce_iter(ch); ce_it.is_valid(); ++ce_it)
    {
        ASSERT_TRUE(MeshElement(*ce_it).is_incident(cache, ch));
    }

    for (auto cf_it = mesh.cf_iter(ch); cf_it.is_valid(); ++cf_it)
    {
        ASSERT_TRUE(MeshElement(*cf_it).is_incident(cache, ch));
    }
}

TEST(HashTest, HexHalfFaceEqualTest)
{
    HexHalfFace hf1 = create_canonical_halfface(
        VertexHandle(0), LCH(0), LPH(0),
        VertexHandle(1), LCH(0), LPH(0),
        VertexHandle(2), LCH(0), LPH(0),
        VertexHandle(3), LCH(0), LPH(0)
        );
    HexHalfFace hf2 = create_canonical_halfface(
        VertexHandle(2), LCH(0), LPH(0),
        VertexHandle(3), LCH(0), LPH(0),
        VertexHandle(0), LCH(0), LPH(0),
        VertexHandle(1), LCH(0), LPH(0)
        );
    assert(hf1 == hf2);
}

}

int main(int _argc, char** _argv) {

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();
}
 
