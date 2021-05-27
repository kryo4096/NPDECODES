/**
 * @file matode_main.cc
 * @brief NPDE homework MatODE code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "matode.h"

int main(int /*argc*/, char** /*argv*/) {
  /* SAM_LISTING_BEGIN_6 */
  double h = 0.01;  // stepsize
  Eigen::Vector3d norms;
  // Build M
  Eigen::Matrix3d M;
  M << 8, 1, 6, 3, 5, 7, 9, 9, 2;
  // Build A
  Eigen::Matrix3d A;
  A << 0, 1, 1, -1, 0, 1, -1, -1, 0;
  Eigen::MatrixXd I = Eigen::Matrix3d::Identity();

  auto QR = M.householderQr();

  Eigen::Matrix3d Q = QR.householderQ() * Eigen::Matrix3d::Identity();

  Eigen::Matrix3d explicit_euler, implicit_euler, implicit_midpoint;

  explicit_euler = Q;
  implicit_euler = Q;
  implicit_midpoint = Q;

	std::cout << (explicit_euler * explicit_euler.transpose() - I).norm() << "\t"
	          <<  (implicit_euler * implicit_euler.transpose() - I).norm() << "\t"
	          <<  (implicit_midpoint * implicit_midpoint.transpose() - I).norm() << std::endl;

	std::cout << 0 << "\t" <<(explicit_euler * explicit_euler.transpose() - I).norm() << "\t"
	          <<  (implicit_euler * implicit_euler.transpose() - I).norm() << "\t"
	          <<  (implicit_midpoint * implicit_midpoint.transpose() - I).norm() << std::endl;

  for(int i = 0; i < 20; i++) {
	explicit_euler = MatODE::eeulstep(A, explicit_euler, h);
	implicit_euler = MatODE::ieulstep(A, implicit_euler, h);
	implicit_midpoint = MatODE::impstep(A, implicit_midpoint, h);

	std::cout << i + 1 << "\t" <<(explicit_euler * explicit_euler.transpose() - I).norm() << "\t"
				<<  (implicit_euler * implicit_euler.transpose() - I).norm() << "\t"
				<<  (implicit_midpoint * implicit_midpoint.transpose() - I).norm() << std::endl;


  }





  /* SAM_LISTING_END_6 */
  return 0;
}
