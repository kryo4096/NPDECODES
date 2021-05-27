#ifndef InitCondLV_CC_
#define InitCondLV_CC_
/**
 * @file initcondlv.cc
 * @brief NPDE homework InitCondLV code
 * @author lfilippo, tille, jgacon, dcasati
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>
#include <utility>

#include "../../../lecturecodes/Ode45/ode45.h"

namespace InitCondLV {

/* Compute the maps Phi(t,y0) and W(t,y0) at final time T.
 * Use initial data given by u0 and v0. */
/* SAM_LISTING_BEGIN_1 */
std::pair<Eigen::Vector2d, Eigen::Matrix2d> PhiAndW(double u0, double v0,
                                                    double T) {
  // Save the values of Phi and W at time T in PaW.first and PaW.second resp.
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW;

  auto y_dot = [&](Eigen::Vector2d y) {
  	return Eigen::Vector2d((2 - y[1]) * y[0], (y[0] - 1) * y[1]);
  };

  auto Df = [&](Eigen::Vector2d y, Eigen::Matrix2d W) {

  	Eigen::Matrix2d Df;

  	Df <<   2 - y[1], -y[0],
  	        y[1],  y[0] - 1;

  	return Df;
  };

  auto yW_dot = [&](Eigen::Matrix<double, 2, 3> state) {

  	auto y = state.col(0);
  	auto W = state.rightCols<2>();

  	return (Eigen::Matrix<double, 2, 3>() << y_dot(y), Df(y, W) * W).finished();
  };

  auto integrator = Ode45<Eigen::Matrix<double, 2, 3>>(yW_dot);

  integrator.options.atol = 1e-12;
  integrator.options.rtol = 1e-14;

  Eigen::Matrix<double, 2, 3> yW_0;
  yW_0 << Eigen::Vector2d(u0,v0), Eigen::Matrix2d::Identity();

  auto [y, _t] = integrator.solve(yW_0, T).back();

  return std::pair(y.col(0), y.rightCols<2>());
}
/* SAM_LISTING_END_1 */

}  // namespace InitCondLV

#endif  // #define InitCondLV_CC_
