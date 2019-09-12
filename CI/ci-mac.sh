#!/bin/bash

#Exit on any error
set -e

LANGUAGE=$1


PATH=$PATH:/opt/local/bin
export PATH

export GTEST_ROOT="${HOME}/sw/gtest-1.7.0"

OPTIONS=""

if [ "$LANGUAGE" == "C++98" ]; then
  echo "Building with C++98";
  BUILDPATH="cpp98"
elif [ "$LANGUAGE" == "C++11" ]; then
  echo "Building with C++11";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++11' "
  BUILDPATH="cpp11"
elif [ "$LANGUAGE" == "C++14" ]; then
  echo "Building with C++14";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++14' "
  BUILDPATH="cpp14"
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

echo "Options:    $OPTIONS"
echo "BuildPath:  $BUILDPATH"
echo "Path:       $PATH"
echo "Language:   $LANGUAGE"

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Building Release version with vectorchecks enabled"
echo "======================================================================"
echo -e "${NC}"


if [ ! -d build-release-$BUILDPATH ]; then
  mkdir build-release-$BUILDPATH
fi

cd build-release-$BUILDPATH

cmake -DCMAKE_BUILD_TYPE=Release -DHEXEX_BUILD_UNIT_TESTS=TRUE $OPTIONS ../

#build it
make

#build the unit tests
make unittests

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


if [ ! -d build-debug-$BUILDPATH-Vector-Checks ]; then
  mkdir build-debug-$BUILDPATH-Vector-Checks
fi

cd build-debug-$BUILDPATH-Vector-Checks

cmake -DCMAKE_BUILD_TYPE=Debug -DHEXEX_BUILD_UNIT_TESTS=TRUE $OPTIONS ../

#build it
make

#build the unit tests
make unittests

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
