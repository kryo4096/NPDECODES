/**
 * @file conslawwithsource.cc
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 20.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "conslawwithsource.h"

#include <cmath>

namespace ConsLawWithSource {

/* SAM_LISTING_BEGIN_1 */
double godnfn(double v, double w) {
  auto f = [](double u) { return std::exp(u) - u; };

  //====================
  // Your code goes here
  // Replace the dummy return value below

  return 0.0;
  //====================
}
/* SAM_LISTING_END_1 */

}  // namespace ConsLawWithSource
