#!/bin/bash

COMPILER=$1
LANGUAGE=$2

# Exit script on any error
set -e

OPTIONS=""
MAKE_OPTIONS=""
BUILDPATH=""

if [ "$COMPILER" == "gcc" ]; then
  echo "Building with GCC";
  BUILDPATH="gcc"

  # without icecc: no options required
  OPTIONS="-DCMAKE_CXX_COMPILER=/usr/lib/icecc/bin/g++ -DCMAKE_C_COMPILER=/usr/lib/icecc/bin/gcc"
  MAKE_OPTIONS="-j16"
  export ICECC_CXX=/usr/bin/g++ ; export ICECC_CC=/usr/bin/gcc

elif [ "$COMPILER" == "clang" ]; then

  OPTIONS="$OPTIONS -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
  echo "Building with CLANG";
  BUILDPATH="clang"
fi

if [ "$LANGUAGE" == "C++98" ]; then
  echo "Building with C++98";
  BUILDPATH="$BUILDPATH-cpp98"
elif [ "$LANGUAGE" == "C++11" ]; then
  echo "Building with C++11";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++11' "
  BUILDPATH="$BUILDPATH-cpp11"
elif [ "$LANGUAGE" == "C++14" ]; then
  echo "Building with C++14";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++14' "
  BUILDPATH="$BUILDPATH-cpp14"
fi

#=====================================
# Color Settings:
#=====================================
NC='\033[0m'
OUTPUT='\033[0;32m'
WARNING='\033[0;93m'

#clone OpenVolumeMesh
if [ ! -d OpenVolueMesh ]; then
  git clone https://graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh --branch cmake-overhaul-3
else
  cd OpenVolumeMesh
  git fetch
  git reset --hard origin/cmake-overhaul-3
  cd ..
fi


echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Basic configuration details:"
echo "======================================================================"
echo -e "${NC}"

echo "Compiler:     $COMPILER"
echo "Options:      $OPTIONS"
echo "Language:     $LANGUAGE"
echo "Make Options: $OPTIONS"
echo "BuildPath:    $BUILDPATH"
echo "Path:         $PATH"
echo "Language:     $LANGUAGE"

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Building Release version"
echo "======================================================================"
echo -e "${NC}"


if [ ! -d build-release-$BUILDPATH ]; then
  mkdir build-release-$BUILDPATH
fi

cd build-release-$BUILDPATH

cmake -DCMAKE_BUILD_TYPE=Release -DHEXEX_BUILD_UNIT_TESTS=TRUE $OPTIONS ../

#build it
make $MAKE_OPTIONS

#build the unit tests
make  $MAKE_OPTIONS unittests

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Running unittests Release version"
echo "======================================================================"
echo -e "${NC}"

cd Build

#execute tests
./unittests --gtest_color=yes --gtest_output=xml

cd ..
cd ..

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Building Debug version"
echo "======================================================================"
echo -e "${NC}"


if [ ! -d build-debug-$BUILDPATH ]; then
  mkdir build-debug-$BUILDPATH
fi

cd build-debug-$BUILDPATH

cmake -DCMAKE_BUILD_TYPE=Debug -DHEXEX_BUILD_UNIT_TESTS=TRUE $OPTIONS ../

#build it
make $MAKE_OPTIONS

#build the unit tests
make  $MAKE_OPTIONS unittests

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Running unittests Debug version"
echo "======================================================================"
echo -e "${NC}"


cd Build

#execute tests
./unittests --gtest_color=yes --gtest_output=xml

cd ..
cd ..
