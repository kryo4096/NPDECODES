/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearfeelementmatrix.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 4, 4> MyLinearFEElementMatrix::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());
  // Matrix for returning element matrix

  auto lapl_mat = lf::uscalfe::LinearFELaplaceElementMatrix().Eval(cell);

  if (ref_el == lf::base::RefElType::kQuad) {
    Eigen::Matrix<double, 4, 4> elem_mat;
    elem_mat << 4,2,1,2,
                2,4,2,1,
                1,2,4,2,
                2,1,2,4;

    double x0 = vertices.row(0).minCoeff();
    double x1 = vertices.row(0).maxCoeff();
    double y0 = vertices.row(1).minCoeff();
    double y1 = vertices.row(1).maxCoeff();

    double area = (x1 - x0) * (y1 - y0);

    return area * elem_mat / 36. + lapl_mat;
  }

  else if(ref_el == lf::base::RefElType::kTria) {
      Eigen::Matrix<double, 4, 4> elem_mat;
      elem_mat << 2,1,1,0,
                  1,2,1,0,
                  1,1,2,0,
                  0,0,0,0;

      Eigen::Matrix2d side_lengths;

      side_lengths.col(0) = vertices.col(1) - vertices.col(0);
      side_lengths.col(1) = vertices.col(2) - vertices.col(0);

      double area = 0.5 * side_lengths.determinant();

      return area * elem_mat / 12. + lapl_mat;

  }
  else {
    assert(0 && "unsupported cell type");
  }

}
/* SAM_LISTING_END_1 */
}  // namespace ElementMatrixComputation
