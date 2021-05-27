/**
 * @file rk3prey_main.cc
 * @brief NPDE homework RK3Prey code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <limits>

#include "rk3prey.h"

typedef std::numeric_limits<double> dbl;

int main() {
	/* SAM_LISTING_BEGIN_0 */
	// Data for the Butcher scheme describing the explicit 3-stage Runge-Kutta
	// method
	Eigen::Matrix3d A;
	Eigen::Vector3d b;
	A << 0, 0, 0, 1.0 / 3.0, 0, 0, 0, 2.0 / 3.0, 0;
	b << 0.25, 0, 0.75;

	auto lv = [](Eigen::Vector2d y) {
		return Eigen::Vector2d((3 - 0.1 * y[1]) * y[0], (0.1 * y[0] - 2) * y[1]);
	};

	// Final time for model
	double T = 10.;

	// Initial value for model
	Eigen::Vector2d y0;
	y0 << 100, 5;

	// Array of number of steps (for convergence study)
	int M[11] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384, 2 * 16384, 4 * 16384, 8 * 16384};

	// Reference "exact" value y(10) at final time T = 10 (approximated)
	Eigen::Vector2d y_ref;
	y_ref << 0.319465882659820, 9.730809352326228;

	// Initialize RK with Butcher scheme
	RK3Prey::RKIntegrator RK(A, b);

	// Start convergence study
	std::cout << std::setw(15) << "M" << std::setw(15) << "error" << std::setw(15)
	          << "rate" << std::endl;

	std::vector<double> errors(std::size(M) + 1);


	for (int i = 0; i < std::size(M); i++) {
		auto sol = RK.solve(lv, T, y0, M[i]);

		errors[i] = (sol[M[i]] - y_ref).norm();

		std::cout << std::setw(15) << M[i] << std::setw(15) << errors[i] << std::setw(15);

		if (i > 0) {
			double rate = (std::log(errors[i - 1]) - std::log(errors[i])) / (std::log(M[i - 1]) - std::log(M[i]));
			std::cout << rate << std::endl;
		} else {
			std::cout << std::endl;
		}
	}


	/* SAM_LISTING_END_0 */
}
