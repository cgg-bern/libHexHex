#include <iostream>
#include <fstream>

#include <HexEx.hh>
#include <HexEx/Stopwatches.hh>
#include <libTimekeeper/StopWatchPrinting.hh>

bool isInFileGood(const std::string& filename)
{
  std::ifstream is(filename.c_str());
  return is.good();
}

bool isOutFileGood(const std::string& filename)
{
  std::ofstream os(filename.c_str());
  return os.good();
}

void printUsage(const std::string& progname)
{
  std::cout << "usage: " << progname << " <inFile> <outFile>" << std::endl
            << std::endl
            << "Reads input tet mesh with parametrization from <inFile> " << std::endl
            << "and writes the resulting hex mesh to <outFile> in ovm format." << std::endl;
}

int main(int argc, const char* argv[])
{

  if (argc != 3)
  {
    printUsage(argv[0]);
    return 0;
  }

  auto inFilename  = std::string(argv[1]);
  auto outFilename = std::string(argv[2]);

  if (!isInFileGood(inFilename))
  {
    std::cout << "Could not open input File " << inFilename << std::endl;
    return 1;
  }

  if (!isOutFileGood(outFilename))
  {
    std::cout << "Could not open output File " << outFilename << std::endl;
    return 1;
  }

  if (isInFileGood(inFilename) && isOutFileGood(outFilename))
    HexEx::extractHexMesh(inFilename, outFilename);

  std::cout << HexEx::sw::root << std::endl;
  return 0;
}
