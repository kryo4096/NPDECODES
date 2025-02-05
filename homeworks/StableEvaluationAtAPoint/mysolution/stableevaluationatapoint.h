/**
 * @file stableevaluationatapoint.h
 * @brief NPDE homework StableEvaluationAtAPoint
 * @author Amélie Loher
 * @date 22/04/2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <complex>

namespace StableEvaluationAtAPoint {

/* Returns the mesh size for the given mesh. */
double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh> &mesh_p) {
  double mesh_size = 0.0;
  // Find maximal edge length
  double edge_length;
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    // Compute the length of the edge
    auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
    edge_length = (endpoints.col(0) - endpoints.col(1)).norm();
    if (mesh_size < edge_length) {
      mesh_size = edge_length;
    }
  }
  return mesh_size;
}

/* Returns G(x,y). */
double G(Eigen::Vector2d x, Eigen::Vector2d y) {
  double res;
  LF_ASSERT_MSG(x != y, "G not defined for these coordinates!");
  // Straightforward implementation
  res = (-1.0 / (2.0 * M_PI)) * std::log((x - y).norm());
  return res;
}

/* Returns the gradient of G(x,y). */
Eigen::Vector2d gradG(Eigen::Vector2d x, Eigen::Vector2d y) {
  Eigen::Vector2d res;
  LF_ASSERT_MSG(x != y, "G not defined for these coordinates!");
  // Straightforward implementation
  res = (x - y) / (2.0 * M_PI * (x - y).squaredNorm());
  return res;
}

/* Evaluates the Integral P_SL using the local midpoint rule
 * on the partitioning of the boundary of Omega induced by the mesh.
 * The supplied meshes are unitary squares.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
double PSL(std::shared_ptr<const lf::mesh::Mesh> mesh, FUNCTOR &&v,
           const Eigen::Vector2d x) {
  double PSLval = 0.0;
  //====================
  // Your code goes here
  //====================
  return PSLval;
}

/* SAM_LISTING_END_1 */

/* Evaluates the Integral P_DL using the local midpoint rule
 * on the partitioning of the boundary of Omega induced by the mesh.
 * The supplied meshes are unitary squares.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
double PDL(std::shared_ptr<const lf::mesh::Mesh> mesh, FUNCTOR &&v,
           const Eigen::Vector2d x) {
  double PDLval = 0.0;
  //====================
  // Your code goes here
  //====================
  return PDLval;
}

/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
/* This function computes u(x) = P_SL(grad u * n) - P_DL(u).
 * For u(x) = log( (x + (1, 0)^T).norm() ) and x = (0.3, 0.4)^T,
 * it computes the difference between the analytical and numerical
 * evaluation of u. The mesh is supposed to cover the unit square.
 */
double pointEval(std::shared_ptr<const lf::mesh::Mesh> mesh) {
  double error = 0.0;
  //====================
  // Your code goes here
  //====================
  return error;
}

/* SAM_LISTING_END_3 */

/* Computes Psi_x(y). */
double Psi(const Eigen::Vector2d y) {
  double Psi_xy;
  const Eigen::Vector2d half(0.5, 0.5);
  const double constant = M_PI / (0.5 * std::sqrt(2) - 1.0);
  const double dist = (y - half).norm();

  if (dist <= 0.25 * std::sqrt(2)) {
    Psi_xy = 0.0;
  } else if (dist >= 0.5) {
    Psi_xy = 1.0;
  } else {
    Psi_xy = std::pow(std::cos(constant * (dist - 0.5)), 2);
  }

  return Psi_xy;
}

/* Computes grad(Psi_x(y)). */
Eigen::Vector2d gradPsi(const Eigen::Vector2d y) {
  Eigen::Vector2d gradPsi_xy;

  Eigen::Vector2d half(0.5, 0.5);
  double constant = M_PI / (0.5 * std::sqrt(2) - 1.0);
  double dist = (y - half).norm();

  if (dist <= 0.25 * std::sqrt(2)) {
    gradPsi_xy(0) = 0.0;
    gradPsi_xy(1) = 0.0;

  } else if (dist >= 0.5) {
    gradPsi_xy(0) = 0.0;
    gradPsi_xy(1) = 0.0;

  } else {
    gradPsi_xy = -2.0 * std::cos(constant * (dist - 0.5)) *
                 std::sin(constant * (dist - 0.5)) * (constant / dist) *
                 (y - half);
  }

  return gradPsi_xy;
}

/* Computes Laplacian of Psi_x(y). */
double laplPsi(const Eigen::Vector2d y) {
  double laplPsi_xy;
  Eigen::Vector2d half(0.5, 0.5);
  double constant = M_PI / (0.5 * std::sqrt(2) - 1.0);
  double dist = (y - half).norm();

  if (dist <= 0.25 * std::sqrt(2)) {
    laplPsi_xy = 0.0;
  } else if (dist >= 0.5) {
    laplPsi_xy = 0.0;
  } else {
    laplPsi_xy =
        (2 * std::pow(constant, 2) / (y - half).squaredNorm()) *
            (y - half).dot(y - half) *
            (std::pow(std::sin(constant * (dist - 0.5)), 2) -
             std::pow(std::cos(constant * (dist - 0.5)), 2)) -
        (2 * constant / dist) * std::cos(constant * (dist - 0.5)) *
            std::sin(constant * (dist - 0.5)) *
            (1.0 - ((y - half).dot(y - half) / (y - half).squaredNorm()));
  }
  return laplPsi_xy;
}

/* Computes Jstar
 * fe_space: finite element space defined on a triangular mesh of the square
 * domain u: Function handle for u x: Coordinate vector for x
 */
/* SAM_LISTING_BEGIN_4 */
template <typename FUNCTOR>
double Jstar(std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
             FUNCTOR &&u, const Eigen::Vector2d x) {
  double val = 0.0;
  //====================
  // Your code goes here
  //====================
  return val;
}

/* SAM_LISTING_END_4 */

/* Evaluates u(x) according to (3.11.14).
 * u: Function Handle for u
 * x: Coordinate vector for x
 */
template <typename FUNCTOR>
double stab_pointEval(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    FUNCTOR &&u, const Eigen::Vector2d x) {
  double res = 0.0;

//====================
// Your code goes here
//====================
  return res;
}

} /* namespace StableEvaluationAtAPoint */
