#include <HexHex/Utils/FileAccessor.hh>
#include <HexHex/Utils/Stopwatches.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/PropertyCodecsT_impl.hh>
#include <OpenVolumeMesh/IO/ovmb_read.hh>
#include <HexHex/Utils/Utils.hh>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>

namespace HexHex
{

namespace FileExtensions
{
enum Type : unsigned char {INVALID, HEXEX, OVM, OVMB, MESH};
inline Type get(const std::filesystem::path& filename)
{
    std::string ext = filename.extension().string();
    if (ext == ".hexex") return HEXEX;
    if (ext == ".ovm") return OVM;
    if (ext == ".ovmb") return OVMB;
    if (ext == ".mesh") return MESH;
    return INVALID;
}
}

static bool loadInputFromOVM(const std::filesystem::path& filename, TetrahedralMesh& tetmesh, OVM::HalfFacePropertyT<Vec3d>& igm)
{
    OpenVolumeMesh::IO::FileManager fm;
    if (!fm.readFile(filename.string(), tetmesh)) {return false;}
    auto prop = tetmesh.get_halfface_property<Vec3d>("HexHex::Parametrization");
    if (!prop.has_value()) {
        std::cerr << "Input Error: Tet mesh must have the property OVM::HalfFacePropertyT<Vec3d>(\"HexHex::Parametrization\")" << std::endl;
        return false;
    }
    igm = prop.value();
    return true;
}

static bool loadInputFromOVMB(const std::filesystem::path& filename, TetrahedralMesh& tetmesh, OVM::HalfFacePropertyT<Vec3d>& igm)
{
    auto codecs = OVM::IO::g_default_property_codecs;
    auto res = OVM::IO::ovmb_read(filename, tetmesh, OVM::IO::ReadOptions(), codecs);
    if (res != OVM::IO::ReadResult::Ok) {
        std::cerr << OVM::IO::to_string(res) << std::endl;
        return false;
    }
    auto prop = tetmesh.get_halfface_property<Vec3d>("HexHex::Parametrization");
    if (!prop.has_value()) {
        std::cerr << "Input Error: Tet mesh must have the property OVM::HalfFacePropertyT<Vec3d>(\"HexHex::Parametrization\")" << std::endl;
        return false;
    }
    igm = prop.value();
    return true;
}

static bool loadInputFromHEXEX(const std::filesystem::path& filename, TetrahedralMesh& tetmesh, OVM::HalfFacePropertyT<Vec3d>& igm)
{
    std::ifstream is(filename, std::ifstream::in);
    if (!is.is_open()) {
        std::cout << "could not open " << filename << std::endl;
        return false;
    }

    /*
     * .hexex File Format:
     *
     * n_vertices
     * x1 y1 z1
     * ...
     * xnv ynv znv
     * n_cells
     * v11 v12 v13 v14 px11 py11 pz11 px12 py12 pz12 px13 py13 pz13 px14 py14 pz14
     * ...
     * vnc1 vnc2 vnc3 vnc4 ...
     */

    auto n_vertices = 0u;
    is >> n_vertices;

    auto read_dbl = [](std::istream &is) -> double {
        // "is >> dbl" does not work with subnormal numbers (at least on MacOS libc++),
        // use this as workaround.
        std::string s;
        is >> s;
        return std::atof(s.c_str());
    };
    auto read_vec3d = [&read_dbl](std::istream &is) -> Vec3d {
        auto x = read_dbl(is);
        auto y = read_dbl(is);
        auto z = read_dbl(is);
        return Vec3d(x, y, z);
    };
    for (unsigned int i = 0; i < n_vertices; ++i) {
        Position pos;
        tetmesh.add_vertex(read_vec3d(is));
        if (is.fail()) {
            std::cout << "reading .hexex failed after reading vertex " << i << std::endl;
            return false;
        }
    }

    auto n_cells = 0u;
    is >> n_cells;

    for (auto i = 0u; i < n_cells; ++i)
    {
        // Read the 4 vertices
        VertexHandle vh0, vh1, vh2, vh3;
        is >> vh0;
        is >> vh1;
        is >> vh2;
        is >> vh3;

        // Create the 4 halffaces
        std::vector<HalfFaceHandle> hfhs;
        hfhs.push_back(tetmesh.add_halfface(vh0, vh1, vh2)); // hf0 - v3
        hfhs.push_back(tetmesh.add_halfface(vh0, vh2, vh3)); // hf1 - v1
        hfhs.push_back(tetmesh.add_halfface(vh0, vh3, vh1)); // hf2 - v2
        hfhs.push_back(tetmesh.add_halfface(vh1, vh3, vh2)); // hf3 - v0

        // Add the tet cell
        tetmesh.add_cell(hfhs);

        // Read the parameters, match to opposite halfface
        auto param = Parameter();
        is >> param;
        igm[hfhs[3]] = param;
        is >> param;
        igm[hfhs[1]] = param;
        is >> param;
        igm[hfhs[2]] = param;
        is >> param;
        igm[hfhs[0]] = param;
    }
    if (is.fail()) {
        std::cout << ".hexex loading: iostream failbit is set!" << std::endl;
        return false;
    }
    return true;
}

bool loadInputFromFile(const std::filesystem::path& filename, TetrahedralMesh& tetmesh, OVM::HalfFacePropertyT<Vec3d>& igm)
{
    auto ext = FileExtensions::get(filename);

    switch (ext) {
    case FileExtensions::OVM:  return loadInputFromOVM(filename, tetmesh, igm);
    case FileExtensions::OVMB:  return loadInputFromOVMB(filename, tetmesh, igm);
    case FileExtensions::HEXEX: return loadInputFromHEXEX(filename, tetmesh, igm);
    default:
        std::cerr << "Unknown Tet Mesh File Format!" << std::endl;
        return false;
    }
}

bool saveInputToHEXEX(const std::filesystem::path& filename, const TetrahedralMesh& tetmesh, const OVM::HalfFacePropertyT<Vec3d>& igm)
{
    std::ofstream os(filename, std::ofstream::out);
    os << std::setprecision(100);

    os << tetmesh.n_vertices() << std::endl;

    for (auto v_it = tetmesh.v_iter(); v_it.is_valid(); ++v_it) {
        os << tetmesh.vertex(*v_it) << std::endl;
    }

    os << tetmesh.n_cells() << std::endl;

    std::vector<Parameter> params;
    for (auto c_it = tetmesh.c_iter(); c_it.is_valid(); ++c_it)
    {
        const auto& vhs = tetmesh.get_cell_vertices(*c_it);
        params.clear();

        // TODO: Use vertex_opposite_halfface of new OVM version
        for (VertexHandle vh : vhs) {
            for (HalfFaceHandle hfh : tetmesh.cell(*c_it).halffaces()) {
                if (tetmesh.halfface_opposite_vertex(hfh) == vh) {
                    params.push_back(igm[hfh]);
                    break;
                }
            }
        }
        if (params.size() != 4) {return false;}

        for (const auto v : vhs) {os << v << " ";}
        for (const auto& param : params) {os << param << " ";}
        os << std::endl;
    }

    return true;
}

static bool saveOutputToOVMB(const std::filesystem::path& filename, const HexahedralMesh& mesh)
{
    auto res = OpenVolumeMesh::IO::ovmb_write(filename, mesh);
    if (res != OVM::IO::WriteResult::Ok) {
        std::cerr << OVM::IO::to_string(res) << std::endl;
        return false;
    }
    return true;
}

static bool saveOutputToOVM(const std::filesystem::path& filename, const HexahedralMesh& mesh)
{
    OpenVolumeMesh::IO::FileManager fm;
    if (!fm.writeFile(filename.string(), mesh)) {return false;}
    return true;
}

static bool saveOutputToMESH(const std::filesystem::path& filename, const HexahedralMesh& mesh)
{
    std::ofstream os(filename, std::ofstream::out);
    //os << std::setprecision(100);

    os << "MeshVersionFormatted 1" << std::endl << "Dimension 3" << std::endl;

    // Vertices
    os << "Vertices " << mesh.n_vertices() << std::endl;
    for (auto v_it = mesh.v_iter(); v_it.is_valid(); ++v_it) {
        os << mesh.vertex(*v_it) << " 0" << std::endl;
    }

    // Hexahedra
    os << "Hexahedra " << mesh.n_cells() << std::endl;
    for (auto c_it = mesh.c_iter(); c_it.is_valid(); ++c_it)
    {
        // Get the four vertices of a halfface
        HalfFaceHandle hfh1 = mesh.cell(*c_it).halffaces()[0];
        const auto& vhs1 = mesh.get_halfface_vertices(hfh1);

        // For each vertex, find its opposite in opposing halfface
        HalfFaceHandle hfh2 = mesh.opposite_halfface_handle_in_cell(hfh1, *c_it);
        std::vector<VertexHandle> vhs2;
        for (VertexHandle vh1 : vhs1) {
            for (VertexHandle vh2 : mesh.get_halfface_vertices(hfh2)) {
                for (auto che_it = mesh.che_iter(*c_it); che_it.is_valid(); ++che_it) {
                    if ((mesh.from_vertex_handle(*che_it) == vh1 && mesh.to_vertex_handle(*che_it) == vh2)) {
                        vhs2.push_back(vh2);
                        break;
                    }
                }
            }
        }
        if (vhs2.size() != 4) {return false;}

        // The order is slighlty different for OVM and MESH, also .MESH is 1-based!
        os << (vhs1[0].idx()+1) << " ";
        os << (vhs1[1].idx()+1) << " ";
        os << (vhs1[2].idx()+1) << " ";
        os << (vhs1[3].idx()+1) << " ";

        os << (vhs2[0].idx()+1) << " ";
        os << (vhs2[1].idx()+1) << " ";
        os << (vhs2[2].idx()+1) << " ";
        os << (vhs2[3].idx()+1) << " ";

        // End of Hex
        os << "0" << std::endl;
    }
    os << "End" << std::endl;

    return true;
}

bool saveOutputToFile(const std::filesystem::path& filename, const HexahedralMesh& mesh)
{
    auto ext = FileExtensions::get(filename);

    switch (ext) {
    case FileExtensions::OVM:  return saveOutputToOVM(filename, mesh);
    case FileExtensions::OVMB:  return saveOutputToOVMB(filename, mesh);
    case FileExtensions::MESH: return saveOutputToMESH(filename, mesh);
    default:
        std::cerr << "Unknown Hex Mesh File Format!" << std::endl;
        return false;
    }
}


} // namespace HexHex
