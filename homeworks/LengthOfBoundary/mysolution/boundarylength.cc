/**
 * @ file boundarylength.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "boundarylength.h"

#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

namespace LengthOfBoundary {

/* SAM_LISTING_BEGIN_1 */
double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double volume = 0.0;

  for(const auto& triangle : mesh_p->Entities(0)) {
    volume += lf::geometry::Volume(*triangle->Geometry());
  }

  return volume;
}
/* SAM_LISTING_END_1 */

using lf::mesh::utils::CodimMeshDataSet;

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double length = 0.0;
  //====================
  // Your code goes here
  //====================

  CodimMeshDataSet<lf::base::size_type> adj_cell_count = lf::mesh::utils::CountNumSuperEntities(mesh_p, 1, 1);

  for(auto& edge : mesh_p->Entities(1)) {
    if(adj_cell_count(*edge) == 1) {
      auto vs = lf::geometry::Corners(*edge->Geometry());
      length += (vs.col(0) - vs.col(1)).norm();
    }
  }

  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(const std::string& filename) {

  auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(factory), filename);

  std::shared_ptr<lf::mesh::Mesh> mesh = reader.mesh();


  return {volumeOfDomain(mesh), lengthOfBoundary(mesh)};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
