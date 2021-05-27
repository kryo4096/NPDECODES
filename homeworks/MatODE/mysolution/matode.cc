/**
 * @file matode.cc
 * @brief NPDE homework MatODE code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>

namespace MatODE {

/* SAM_LISTING_BEGIN_3 */
	Eigen::MatrixXd eeulstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
	                         double h) {
		return Y0 + A * Y0 * h;
	}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
	Eigen::MatrixXd ieulstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
	                         double h) {

		Eigen::MatrixXd Y1 = (Eigen::MatrixXd::Identity(A.rows(), A.cols()) - h * A).lu().solve(Y0);

		return Y1;
	}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
	Eigen::MatrixXd impstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
	                        double h) {

		Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(), A.cols());

		return 	(I - 0.5 * h * A).lu().solve((I + 0.5 * h * A) * Y0);
	}
/* SAM_LISTING_END_5 */

}  // namespace MatODE
