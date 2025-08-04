#include <HexHex/Config.hh>
#include <HexHex/Config_json.hh>

#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>

namespace HexHex {

bool loadConfig(const std::filesystem::path& filename, Config& config)
{
    std::ifstream file(filename);
    if (!file) return false;
    try
    {
        nlohmann::json j;
        file >> j;
        config = j.get<Config>();
        return true;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading config: " << e.what() << std::endl;
        return false;
    }
}

} // namespace HexHex
