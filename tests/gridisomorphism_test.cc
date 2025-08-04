/*
 * Copyright 2025 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */

#include <chrono>
#include <gtest/gtest.h>
#include "common.hh"
#include <HexHex/Utils/Utils.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Commons/Transition.hh>

using namespace HexHex;

TEST(GridIsomorphismTest, matrixTest) {

    double entries[] = {1, 0, 0, 0,
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        3, 2, 0, 1 };
    auto t1 = Matrix4x4d(entries);

    Vec3d test = Vec3d(1,2,3);

    auto test2 = test;

    test2 = t1.transform_point(test2);

    EXPECT_EQ(test + Vec3d(3,2,0), test2);
}


TEST(GridIsomorphismTest, coversionTest) {

    for (auto i = 0; i < 24; i++)
    {
        auto rr = HexHex::RestrictedRotation(i);
        auto m = rr.toMatrix();

        auto rr2 = HexHex::RestrictedRotation(m);

        EXPECT_EQ(rr, rr2);
    }

}

TEST(GridIsomorphismTest, transformationTest) {

    for (auto i = 0; i < 24; i++)
    {
        auto rr = HexHex::RestrictedRotation(i);
        auto m = rr.toMatrix();

        for (auto j = 0; j < 100; ++j)
        {
            auto randomVec = getRandomVector(10);
            EXPECT_EQ(m.transform_point(randomVec), rr.transform(randomVec)) << i;
        }
    }

}

TEST(GridIsomorphismTest, inversionTest) {

    for (char i = 0; i < 24; i++)
    {
        auto rr = HexHex::RestrictedRotation(i);
        auto m = rr.toMatrix();
        m.invert();

        auto rr2 = rr;
        rr2.invert();
        auto m2 = rr2.toMatrix();

        EXPECT_EQ(m2, m) << i;

    }

}

TEST(GridIsomorphismTest, multiplicationTest) {

    for (auto i = 0u; i < 24; i++)
        for (auto j = 0u; j < 24; j++)
        {
            auto rr1 = HexHex::RestrictedRotation(i);
            auto rr2 = HexHex::RestrictedRotation(j);

            auto res = rr1 * rr2;


            EXPECT_EQ(res.toMatrix(), rr1.toMatrix()*rr2.toMatrix()) << i << " " << j;

        }

}

TEST(GridIsomorphismTest, isomorphismInversionTest) {

    for (auto j = 0; j < 3; ++j)
        for (auto k = 0; k < 3; ++k)
            for (auto l = 0; l < 3; ++l)
            {
                auto t = HexHex::Transition::Translation(j,k,l);
                for (auto i = 0; i < 24; i++)
                {
                    auto gi = HexHex::Transition(i, t);
                    auto m = gi.toMatrix();
                    m.invert();

                    auto gi2 = gi.inverted();
                    auto m2 = gi2.toMatrix();

                    EXPECT_EQ(m2, m) << i;

                }
            }

}

//TEST(GridIsomorphismTest, performanceMulitplicitionTest) {

//    using namespace std::chrono;

//    auto N = 1000000;
//#ifdef DEBUG
//    N /= 100;
//#endif
//    auto gis = std::vector<HexHex::GridIsomorphism>();
//    for (auto i = 0; i < 24; ++i)
//        gis.push_back(HexHex::GridIsomorphism(i));

//    auto matrices = std::vector<Matrix4x4dd>();
//    for (auto i = 0; i < 24; ++i)
//        matrices.push_back(gis[i].toMatrix());

//    auto startGI = steady_clock::now();

//    auto res = HexHex::GridIsomorphism(0);
//    for (auto n = 0; n < N; ++n)
//    {
//        for (auto i = 0u; i < 24; i++)
//            for (auto j = 0u; j < 24; j++)
//            {
//                res = gis[i] * gis[j];
//            }
//    }

//    auto stopGI = steady_clock::now();

//    auto durationGI = duration_cast<duration<double>>(stopGI - startGI);

//    std::cout << res.transform_point(Vec3d(0,0,0)) << std::endl;

//    auto startMatrices = steady_clock::now();

//    auto res2 = Matrix4x4dd();
//    for (auto n = 0; n < N; ++n)
//    {
//        for (auto i = 0u; i < 24; i++)
//            for (auto j = 0u; j < 24; j++)
//            {
//                res2 = matrices[i] * matrices[j];
//            }
//    }

//    auto stopMatrices = steady_clock::now();

//    auto durationMatrices = duration_cast<duration<double>>(stopMatrices - startMatrices);

//    std::cout << res2.transform_point(Vec3d(0,0,0)) << std::endl;

//    std::cout << "Duration GridIsomorphism: " << durationGI.count() << std::endl;
//    std::cout << "Duration Matrices       : " << durationMatrices.count() << std::endl;
//    std::cout << "Ratio                   : " << static_cast<double>(durationGI.count())/durationMatrices.count() << std::endl;

//}




//TEST(GridIsomorphismTest, performanceInversionTest) {

//    using namespace std::chrono;

//    auto N = 10000000;
//#ifdef DEBUG
//    N /= 100;
//#endif
//    auto gis = std::vector<HexHex::GridIsomorphism>();
//    for (auto i = 0; i < 24; ++i)
//        gis.push_back(HexHex::GridIsomorphism(i));

//    auto matrices = std::vector<Matrix4x4dd>();
//    for (size_t i = 0; i < 24; ++i)
//        matrices.push_back(gis[i].toMatrix());

//    auto startGI = steady_clock::now();

//    for (auto n = 0; n < N; ++n)
//    {
//        for (auto i = 0u; i < 24; i++)
//            gis[i].invert();
//    }

//    auto stopGI = steady_clock::now();

//    auto durationGI = duration_cast<duration<double>>(stopGI - startGI);

//    auto startMatrices = steady_clock::now();

//    for (auto n = 0; n < N; ++n)
//    {
//        for (auto i = 0u; i < 24; i++)
//            matrices[i].invert();
//    }

//    auto stopMatrices = steady_clock::now();

//    auto durationMatrices = duration_cast<duration<double>>(stopMatrices - startMatrices);

//    std::cout << "Duration GridIsomorphism: " << durationGI.count() << std::endl;
//    std::cout << "Duration Matrices       : " << durationMatrices.count() << std::endl;
//    std::cout << "Ratio                   : " << (double)durationGI.count()/durationMatrices.count() << std::endl;

//}


//TEST(GridIsomorphismTest, performanceTransformationTest) {

//    using namespace std::chrono;

//    auto N = 10000000;
//#ifdef DEBUG
//    N /= 100;
//#endif
//    auto gis = std::vector<HexHex::GridIsomorphism>();
//    for (auto i = 0; i < 24; ++i)
//        gis.push_back(HexHex::GridIsomorphism(i));

//    auto matrices = std::vector<Matrix4x4dd>();
//    for (auto i = 0; i < 24; ++i)
//        matrices.push_back(gis[i].toMatrix());

//    auto startGI = steady_clock::now();

//    auto res = Vec3d(3,3,3);
//    for (auto n = 0; n < N; ++n)
//    {
//        for (auto i = 0u; i < 24; i++)
//            res += gis[i].transform_point(res);
//    }

//    auto stopGI = steady_clock::now();

//    auto durationGI = duration_cast<duration<double>>(stopGI - startGI);
//    std::cout << res << std::endl;

//    auto startMatrices = steady_clock::now();

//    auto res2 = Vec3d(3,3,3);
//    for (auto n = 0; n < N; ++n)
//    {
//        for (auto i = 0u; i < 24; i++)
//            res2 += matrices[i].transform_point(res2);
//    }

//    auto stopMatrices = steady_clock::now();

//    auto durationMatrices = duration_cast<duration<double>>(stopMatrices - startMatrices);

//    std::cout << res2 << std::endl;


//    std::cout << "Duration GridIsomorphism: " << durationGI.count() << std::endl;
//    std::cout << "Duration Matrices       : " << durationMatrices.count() << std::endl;
//    std::cout << "Ratio                   : " << (double)durationGI.count()/durationMatrices.count() << std::endl;

//}
