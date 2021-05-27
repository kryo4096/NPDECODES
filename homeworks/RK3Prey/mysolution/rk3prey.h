/**
 * @file rk3prey_main.cc
 * @brief NPDE homework RK3Prey code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <utility>
#include <vector>

namespace RK3Prey {

// Butcher Tableau based Runge-Kutta explicit solver for autonomous ODEs
/* SAM_LISTING_BEGIN_0 */
class RKIntegrator {
 public:
  RKIntegrator(Eigen::MatrixXd A, Eigen::VectorXd b) : A(std::move(A)), b(std::move(b)) {

  }

  // Explicit Runge-Kutta numerical integrator
  template <class Function>
  std::vector<Eigen::VectorXd> solve(Function &&f, double T,
                                     const Eigen::VectorXd &y0, int M) const;

 private:
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
};
/* SAM_LISTING_END_0 */

/* Solves an autonomous ODE y' = f(y), y(0) = y0, using a
 * RK scheme from the Butcher tableau provided by the
 * constructor. Performs N equidistant steps up to time T */
/* SAM_LISTING_BEGIN_1 */
template <typename Function>
std::vector<Eigen::VectorXd> RKIntegrator::solve(Function &&f, double T,
                                                 const Eigen::VectorXd &y0,
                                                 int M) const {
  long dim = y0.size();  // dimension
  double h = T / M;     // step size
  std::vector<Eigen::VectorXd> sol(M+1);

  sol[0] = y0;

  for (long n = 1; n <= M; n++) {

	  long s = A.rows();

	  Eigen::MatrixXd K(dim, s);

	  for (long i = 0; i < s; i++) {
		  K.col(i) = f(sol[n-1] + h * K.leftCols(std::fmax(1, i)) * A.row(i).transpose().head(std::fmax(1, i)));
	  }

	  sol[n] = sol[n-1] + h * K * b;
  }

  return sol;
}
/* SAM_LISTING_END_1 */

}  // namespace RK3Prey
