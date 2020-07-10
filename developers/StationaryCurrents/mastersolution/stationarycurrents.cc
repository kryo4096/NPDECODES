/** @file
 * @brief Solving a second-order elliptic BVP with special mixed boundary
 * conditions
 * @author Ralf Hiptmair
 * @date July 2020
 * @copyright MIT License
 */

#include "stationarycurrents.h"

namespace dmxbc {

std::pair<std::shared_ptr<const lf::mesh::Mesh>,
          lf::mesh::utils::CodimMeshDataSet<int>>
readMeshWithTags(std::string filename) {
  std::cout << "Reading mesh from " << filename << std::endl;
  // Total number of contacts
  const int NPhysGrp = 2;
  // physical names for contacts
  const std::array<std::string, NPhysGrp> contactnames{"Contact0", "Contact1"};
  // load the mesh from a .msh file produced by Gmsh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), filename);

  // Obtain list of physical groups for edges (codim = 1 entities)
  std::vector<std::pair<lf::base::size_type, std::string>> phys_ent_list{
      reader.PhysicalEntities(1)};
  // Search for names and find associated ids
  std::array<lf::base::size_type, NPhysGrp> ids{0};
  for (int j = 0; j < NPhysGrp; j++) {
    auto found = std::find_if(phys_ent_list.begin(), phys_ent_list.end(),
                              [j, &contactnames](auto phys) -> bool {
                                return (phys.second == contactnames[j]);
                              });
    if (found == phys_ent_list.end()) {
      // A physical name was not found in the edge groups
      std::cout << contactnames[j] << " not found" << std::endl;
    } else {
      // Physical name found
      ids[j] = found->first;
      std::cout << contactnames[j] << " <-> id = " << ids[j] << std::endl;
    }
  }
  // Obtain pointer to mesh object
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{reader.mesh()};
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Output information on the mesh
  std::cout << mesh << std::endl;
  // A set of integers associated with edges of the mesh (codim = 1 entities)
  lf::mesh::utils::CodimMeshDataSet<int> edgeids{mesh_p, 1, -1};
  // Counter for nodes on a particular part of the boundary
  std::array<int, NPhysGrp> edcnt{0};
  // Loop over edges, check their physical groups, and mark their endpoints
  for (const lf::mesh::Entity* edge : mesh.Entities(1)) {
    LF_ASSERT_MSG(edge->RefEl() == lf::base::RefEl::kSegment(),
                  " edge must be a SEGMENT!");
    for (int j = 0; j < NPhysGrp; ++j) {
      if (reader.IsPhysicalEntity(*edge, ids[j])) {
        // The edge belongs to the physical group with j-th name
        edcnt[j]++;
        edgeids(*edge) = j;
      }
    }
  }
  // Print number of edges  on contacts
  for (int j = 0; j < NPhysGrp; ++j) {
    std::cout << edcnt[j] << " edges in group " << contactnames[j] << std::endl;
  }
  return {mesh_p, edgeids};
}  // function readMeshWithTags()

// Spreading tags to nodes
lf::mesh::utils::CodimMeshDataSet<int> tagNodes(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<int> edgeids) {
  // Current mesh object
  const lf::mesh::Mesh& mesh{*mesh_p};
  // A set of integer ids associated with nodes of the mesh (codim = 2 entities)
  lf::mesh::utils::CodimMeshDataSet<int> nodeids{mesh_p, 2, -1};

  // Loop over edges and spread their ids, if non-negative
  for (const lf::mesh::Entity* edge : mesh.Entities(1)) {
    LF_ASSERT_MSG(edgeids.DefinedOn(*edge),
                  "No flag available for edge " << *edge);
    const int id = edgeids(*edge);
    if (id >= 0) {
      // Obtain iterator over set of endpoints
      nonstd::span<const lf::mesh::Entity* const> sub_ent_range{
          edge->SubEntities(1)};
      LF_ASSERT_MSG(sub_ent_range.size() == 2, " Edge with #endpoints != 2!");
      // Access endpoints
      const lf::mesh::Entity& ep0{*sub_ent_range[0]};
      const lf::mesh::Entity& ep1{*sub_ent_range[1]};
      // Set physical group ids for the endpoints
      nodeids(ep0) = id;
      nodeids(ep1) = id;
    }
  }
  return nodeids;
}  // end tagNodes

// Debugging function
void printNodeTags(const lf::mesh::Mesh& mesh,
                   lf::mesh::utils::CodimMeshDataSet<int>& nodeids) {
  // Container for counters
  std::map<int, unsigned int> counters;
  // Loop over nodes of the mesh and counter occurrence of ids
  for (const lf::mesh::Entity* node : mesh.Entities(2)) {
    if (nodeids.DefinedOn(*node)) {
      counters[nodeids(*node)]++;
    } else {
      std::cout << " Node " << *node << " has no id!" << std::endl;
    }
  }
  for (auto& cnt : counters) {
    std::cout << "id = " << cnt.first << ": " << cnt.second << " nodes"
              << std::endl;
  }
}  // end printNodeTags

std::tuple<Eigen::Matrix<double, 2, 3>, Eigen::Matrix<double, 2, 3>, double>
getTriangleGradLambdaNormals(Eigen::Matrix<double, 2, 3> vertices) {
  // Compute gradients of barycentric coordinate functions for a flat triangle,
  // whose vertex coordinates are passed in the columns of the argument matrix
  // The algorithm is explained in Remark 2.4.5.9
  Eigen::Matrix<double, 3, 3> X;
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  // Compute the gradients of the barycentric coordinate functions
  // and store them in the columns of a 2x3 matrix grad\_bary\_coords
  const auto grad_bary_coords{X.inverse().block<2, 3>(1, 0)};
  // The determinant of the auxiliary matrix also supplies the area of the
  // triangle
  const double twicearea = std::abs(X.determinant());
  return {grad_bary_coords,
          -twicearea * grad_bary_coords *
              Eigen::PermutationMatrix<3>(Eigen::Vector3i(2, 0, 1)),
          twicearea / 2};
}  // end getTriangleGradLambdaNormals

Eigen::Matrix<double, 2, 3> exteriorTriangleNormals(
    Eigen::Matrix<double, 2, 3> vertices) {
  // Compute gradients of barycentric coordinate functions for a flat triangle,
  // whose vertex coordinates are passed in the columns of the argument matrix
  // The algorithm is explained in Remark 2.4.5.9
  Eigen::Matrix<double, 3, 3> X;
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  // Compute the gradients of the barycentric coordinate functions
  // and store them in the columns of a 2x3 matrix grad\_bary\_coords
  const auto grad_bary_coords{X.inverse().block<2, 3>(1, 0)};
  // The determinant of the auxiliary matrix also supplies the area of the
  // triangle
  const double twicearea = std::abs(X.determinant());
  return -twicearea * grad_bary_coords *
         Eigen::PermutationMatrix<3>(Eigen::Vector3i(2, 0, 1));
}  // end exteriorTriangleNormals

Eigen::MatrixXd exteriorCellNormals(const Eigen::MatrixXd& corners) {
  LF_ASSERT_MSG(corners.rows() == 2,
                "Columns must contains coordinates of 2D points");
  // Number of vertices
  unsigned int m = corners.cols();
  // Center of gravity
  const Eigen::Vector2d c{corners.rowwise().sum() / m};
  // Matrix in whose columns we return the normals
  Eigen::MatrixXd normals(corners.rows(), corners.cols());
  for (int j = 0; j < m; ++j) {
    // Edge vector
    const Eigen::Vector2d edge_vec{corners.col((j + 1) % m) - corners.col(j)};
    // Vector from barycenter to endpoint of the edge
    const Eigen::Vector2d a{corners.col(j) - c};
    // Normal vector
    normals.col(j) = Eigen::Vector2d(-edge_vec[1], edge_vec[0]);
    // Adjust direction
    if (a.dot(normals.col(j)) < 0) {
      normals.col(j) *= -1;
    }
  }  // End: loop over straight edges
  return normals;
}  // end: exteriorCellNormals

lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> exteriorEdgeWeightedNormals(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Array indexed by edges for returning exterior normals
  lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> extnormals(
      mesh_p, 1, Eigen::Vector2d(0.0, 0.0));
  // Find edges on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Loop over all cells
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    // Make sure the cell is of triangular shape
    const lf::base::RefEl ref_el_type{cell->RefEl()};
    LF_ASSERT_MSG(ref_el_type == lf::base::RefEl::kTria(),
                  "Only implemented for triangular cells");
    const auto vertices{lf::geometry::Corners(*(cell->Geometry()))};
    const Eigen::Matrix<double, 2, 3> normals{
        exteriorTriangleNormals(vertices)};
    // Obtain array of edge pointers (relative co-dimension = 1)
    nonstd::span<const lf::mesh::Entity* const> sub_ent_range{
        cell->SubEntities(1)};
    // loop over the edges and check whether they belong to the boundary
    for (lf::base::sub_idx_t j = 0; j < ref_el_type.NumSubEntities(1); ++j) {
      const lf::mesh::Entity& edge{*sub_ent_range[j]};
      if (bd_flags(edge)) {
        // Found edge on the boundary. Set normal vector
        extnormals(edge) = normals.col(j);
      }
    }
  }
  return extnormals;
}  // end exteriorEdgeWeightedNormals

// A debugging function
bool validateNormals(const lf::mesh::Mesh& mesh) {
  // Run through cells and compute the normals
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    const lf::base::RefEl ref_el_type{cell->RefEl()};
    LF_ASSERT_MSG(ref_el_type == lf::base::RefEl::kTria(),
                  "implemented for triangles only");
    const lf::geometry::Geometry& geo{*(cell->Geometry())};
    const auto vertices{lf::geometry::Corners(geo)};

    auto etn{exteriorTriangleNormals(vertices)};
    auto ecn{exteriorCellNormals(vertices)};
    if ((etn - ecn).norm() > 1E-8) {
      std::cout << "Mismatch of exterior normals for " << std::endl
                << vertices << std::endl;
      std::cout << "etn = \n" << etn << std::endl;
      std::cout << "ecn = \n" << ecn << std::endl;
      return false;
    }
  }  // end loop over edges
  std::cout << "Normals ok" << std::endl;
  return true;
}  // end validate normals

double computeMeshwidth(const lf::mesh::Mesh& mesh) {
  double h = 0.0;
  for (const lf::mesh::Entity* edge : mesh.Entities(1)) {
    h = std::max(h, lf::geometry::Volume(*(edge->Geometry())));
  }
  return h;
}

// The weight function for the volume-based flux formula
static auto psi = [](Eigen::Vector2d x) -> double {
  if (x[0] > 3)
    return 1.0;
  else if (x[0] < 2)
    return 0.0;
  const double p = std::cos((3.0 - x[0]) * M_PI / 2);
  return (p * p);
};

// Gradient of the weight function
static auto psi_grad = [](Eigen::Vector2d x) -> Eigen::Vector2d {
  double gx = 0.0;
  if ((x[0] > 2.0) && (x[0] < 3.0)) {
    const double a = (3.0 - x[0]) * M_PI / 2;
    gx = M_PI * std::cos(a) * std::sin(a);
  }
  return Eigen::Vector2d(gx, 0.0);
};

std::tuple<double, double, double> computePotential(std::string basename) {
  // Path to mesh .msh file
  std::string mesh_path = CURRENT_SOURCE_DIR "/../meshes/" + basename + ".msh";
  // Read mesh and label nodes
  auto [mesh_p, edgeids] = readMeshWithTags(mesh_path);
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Distribute tags to nodes
  auto nodeids{tagNodes(mesh_p, edgeids)};
  // Print diagnostic information
  printNodeTags(mesh, nodeids);
  // Set up global FE space; lowest order Lagrangian finite elements
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Coefficient function for conductivity
  auto sigma = [](Eigen::Vector2d x) -> double {
    return std::max(1.0, x[0] - 1.0);
  };
  // Solve discrete boundary value problem. Potential values at contacts are
  // passed as a list of values corresponding to "nodal ids" 0,1,... stored in
  // the nodeids array.
  Eigen::VectorXd sol_vec{solveMixedBVP(fe_space, nodeids, {0.0, 1.0}, sigma)};

  // Output solution into a VTK file
  std::cout << "VTK output for " << basename << std::endl;
  auto mds = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2, 0.0);
  for (const lf::mesh::Entity* node : mesh.Entities(2)) {
    auto gdof_idx{(fe_space->LocGlobMap()).GlobalDofIndices(*node)};
    LF_ASSERT_MSG(gdof_idx.size() == 1, " A node can hold only one dof!");
    const lf::assemble::gdof_idx_t nd_idx = gdof_idx[0];
    (*mds)(*node) = sol_vec[nd_idx];
  }
  // Write VTK output
  {
    lf::io::VtkWriter vtk_writer(mesh_p, (basename + ".vtk").c_str());
    vtk_writer.WritePointData("potential", *mds);
  }
  // Compute the contact flux
  std::cout << basename << ", compute contact flux: boundary formula"
            << std::endl;
  const double contact_flux = contactFlux(fe_space, sol_vec, sigma, edgeids, 1);
  const double contact_flux_mf =
      contactFluxMF(fe_space, sol_vec, sigma, edgeids, 1);
  std::cout << "Contact flux 1 = " << contact_flux << " = " << contact_flux_mf
            << std::endl;
  std::cout << basename << ", compute contact flux: volume formula"
            << std::endl;
  const double stab_flux = stabFlux(fe_space, sol_vec, sigma, psi_grad);
  const double stab_flux_mf = stabFluxMF(fe_space, sol_vec, sigma, psi_grad);
  const double stab_flux_mpr = stabFluxMPR(fe_space, sol_vec, sigma, psi_grad);
  std::cout << "Stabilized flux = " << std::setprecision(16) << stab_flux
            << " = " << stab_flux_mf << ", by MPR = " << stab_flux_mpr
            << std::endl;
  return {computeMeshwidth(mesh),contact_flux, stab_flux_mpr};
}  // end computePotential

}  // namespace dmxbc
