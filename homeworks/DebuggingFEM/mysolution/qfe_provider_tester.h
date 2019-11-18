/**
 * @file qfe_provider_tester.h
 * @brief NPDE homework DebuggingFEM code
 * @author Oliver Rietmann
 * @date 03/04/2019
 * @copyright Developed at ETH Zurich
 */

#ifndef QFE_PROV
#define QFE_PROV

#include "qfe_interpolator.h"

#include <lf/assemble/assemble.h>

namespace DebuggingFEM {

/* SAM_LISTING_BEGIN_1 */
template <typename ENTITY_MATRIX_PROVIDER>
class QFEProviderTester {
 public:
  /**
   * @brief Sets up and stores the Galerkin matrix
   * @param dofh dof handler
   * @param element_matrix_provider an instance of LocalLaplaceQFEX, X=1,2,3
   */
  QFEProviderTester(lf::assemble::DofHandler &dofh,
                    ENTITY_MATRIX_PROVIDER &element_matrix_provider);

  /**
   * @brief Computes the energy (H_1-seminorm) of u
   * @param u function of type double(Eigen::Vector2d)
   */
  template <typename FUNCTOR>
  double energyOfInterpolant(FUNCTOR &&u) const;

 private:
  lf::assemble::DofHandler &dofh_;
  ENTITY_MATRIX_PROVIDER &element_matrix_provider_;
  Eigen::SparseMatrix<double> A_;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename ENTITY_MATRIX_PROVIDER>
QFEProviderTester<ENTITY_MATRIX_PROVIDER>::QFEProviderTester(
    lf::assemble::DofHandler &dofh,
    ENTITY_MATRIX_PROVIDER &element_matrix_provider)
    : dofh_(dofh), element_matrix_provider_(element_matrix_provider) {
  //====================
  // Your code goes here
  //====================
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <typename ENTITY_MATRIX_PROVIDER>
template <typename FUNCTOR>
double QFEProviderTester<ENTITY_MATRIX_PROVIDER>::energyOfInterpolant(
    FUNCTOR &&u) const {
  double energy;
  //====================
  // Your code goes here
  //====================
  return energy;
}

/* SAM_LISTING_END_3 */

}  // namespace DebuggingFEM

#endif
