#include <iostream>

#include <HexEx.hh>
#include <HexEx/Stopwatches.hh>
#include <libTimekeeper/StopWatchPrinting.hh>

int main()
{
  using namespace HexEx;
  using namespace OpenVolumeMesh;
  using Vec3d = HexEx::Vec3d;

  // create a tet mesh representing a cube
  TetrahedralMesh tetMesh;

  // add the eight corners
  auto corners = std::vector<OpenVolumeMesh::VertexHandle>();
  corners.push_back(tetMesh.add_vertex(Vec3d(0,0,0)));
  corners.push_back(tetMesh.add_vertex(Vec3d(0,0,1)));
  corners.push_back(tetMesh.add_vertex(Vec3d(1,0,1)));
  corners.push_back(tetMesh.add_vertex(Vec3d(1,0,0)));
  corners.push_back(tetMesh.add_vertex(Vec3d(0,1,0)));
  corners.push_back(tetMesh.add_vertex(Vec3d(0,1,1)));
  corners.push_back(tetMesh.add_vertex(Vec3d(1,1,1)));
  corners.push_back(tetMesh.add_vertex(Vec3d(1,1,0)));

  // add the five tets
  tetMesh.add_cell(corners[0], corners[1], corners[2], corners[5]);
  tetMesh.add_cell(corners[2], corners[3], corners[0], corners[7]);
  tetMesh.add_cell(corners[5], corners[4], corners[7], corners[0]);
  tetMesh.add_cell(corners[7], corners[6], corners[5], corners[2]);
  tetMesh.add_cell(corners[0], corners[5], corners[2], corners[7]);

  // create a property that can store the parameter location for each of the four vertices of each tet
  auto parametrization = tetMesh.request_cell_property<std::map<OpenVolumeMesh::VertexHandle, Vec3d> >("Parametrization");

  // optionally set parametrization persistent so the property
  // stays on tetMesh even after parametrization goes out of scope
  tetMesh.set_persistent(parametrization);

  // fill parametrization property
  // iterate over all cells
  for (auto c_it = tetMesh.cells_begin(); c_it  != tetMesh.cells_end(); ++c_it)
    // iterate over all vertices of that cell
    for (auto cv_it = tetMesh.cv_iter(*c_it); cv_it.valid(); ++cv_it)
      // set parameter to 2 times its position
      parametrization[*c_it][*cv_it] = 2 * tetMesh.vertex(*cv_it);

  // create a hex mesh
  HexahedralMesh hexMesh;

  // extract the hex mesh
  extractHexMesh(tetMesh, parametrization, hexMesh);

  // print some interesting information
  std::cout << "The extracted hex mesh has " << hexMesh.n_cells() << " cells, "
            << hexMesh.n_faces() << " faces, " << hexMesh.n_edges() << " edges, and "
            << hexMesh.n_vertices() << " vertices." << std::endl;

  std::cout << HexEx::sw::root << std::endl;
}
