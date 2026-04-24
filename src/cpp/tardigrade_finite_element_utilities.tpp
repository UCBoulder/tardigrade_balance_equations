/**
 ******************************************************************************
 * \file tardigrade_finite_element_utilities.tpp
 ******************************************************************************
 * The source file for utilities which can assist using the balance equations
 * in finite element codes. We here assume that unknown quantities are
 * computed at the evaluation point using interpolation functions (here called
 * interp) and are projected to the nodes using test functions (here called
 * test).
 ******************************************************************************
 */

#include <Eigen/Dense>

namespace tardigradeBalanceEquations {

    namespace finiteElement {

        /*!
         * Compute the derivative of the spatial gradient of a quantity a
         * in the current configuration w.r.t. the spatial degrees of freedom
         *
         * \f$ \frac{D}{Du_a} \left( a_{i,j} \right) \f$
         *
         * \param &grad_a_start: An iterator representing the start of the gradient of the quantity
         * \param &grad_a_size: The size of the gradient of a
         * \param &grad_interp: The gradient of the interpolation function
         * \param &index: The index of the spatial degree of freedom to compute the derivative for (0, 1 or 2)
         * \param &dgrad_adui_start: An iterator representing the start of the output
         */
        template <typename grad_iterator, typename output_iterator>
        void computeGradientSpatialJacobian(const grad_iterator &grad_a_start, const unsigned int grad_a_size,
                                            floatVector grad_interp, const unsigned int index,
                                            output_iterator dgrad_adui_start) {
            TARDIGRADE_ERROR_TOOLS_CHECK((grad_a_size % dim) == 0, "The incoming spatial gradient has a dimension of " +
                                                                       std::to_string(grad_a_size) +
                                                                       " which is not a multiple of " +
                                                                       std::to_string(dim));

            const unsigned int a_size = grad_a_size / dim;

            std::fill(dgrad_adui_start, dgrad_adui_start + a_size, 0);

            for (unsigned int i = 0; i < a_size; ++i) {
                for (unsigned int j = 0; j < dim; ++j) {
                    *(dgrad_adui_start + dim * i + j) -= *(grad_a_start + dim * i + index) * grad_interp[j];
                }
            }
        }

        /*!
         * Compute the current rate of change of a degree of freedom vector
         * by assuming a generalized trapezoidal rule integration scheme
         *
         * \param &dt: The change in time
         * \param &v_t_begin: The starting iterator of the previous dof vector
         * \param &v_t_end: The stopping iterator of the previous dof vector
         * \param &v_tp1_begin: The starting iterator of the current dof vector
         * \param &v_tp1_end: The stopping iterator of the current dof vector
         * \param &vDot_t_begin: The starting iterator of the previous dof rate of change vector
         * \param &vDot_t_end: The stopping iterator of the previous dof rate of change vector
         * \param &alpha: The integration parameter (0 is explicit, 1 is implicit)
         * \param &vDot_tp1_begin: The starting iterator of the current dof rate of change vector
         * \param &vDot_tp1_end: The stopping iterator of the current dof rate of change vector
         * \param &dVDotdV: The derivative of the current dof rate of change vector with respect
         *     to the dof vector
         */
        template <typename dt_type, class v_t_in, class v_tp1_in, class vDot_t_in, typename alpha_type,
                  class vDot_tp1_out, typename dVDotdV_type>
        void compute_current_rate_of_change(const dt_type &dt, const v_t_in &v_t_begin, const v_t_in &v_t_end,
                                            const v_tp1_in &v_tp1_begin, const v_tp1_in &v_tp1_end,
                                            const vDot_t_in &vDot_t_begin, const vDot_t_in &vDot_t_end,
                                            const alpha_type &alpha, vDot_tp1_out vDot_tp1_begin,
                                            vDot_tp1_out vDot_tp1_end, dVDotdV_type &dVDotdV) {
            dVDotdV = 1. / (alpha * dt);

            for (unsigned int i = 0; i < (unsigned int)(v_t_end - v_t_begin); ++i) {
                *(vDot_tp1_begin + i) = ((*(v_tp1_begin + i)) - (*(v_t_begin + i))) / (alpha * dt) -
                                        ((1 - alpha) / alpha) * (*(vDot_t_begin + i));
            }
        }

        /*!
         * Compute the current acceleration (second time derivative) of a degree of freedom vector
         * by assuming a generalized trapezoidal rule integration scheme
         *
         * \param &dt: The change in time
         * \param &v_t_begin: The starting iterator of the previous dof vector
         * \param &v_t_end: The stopping iterator of the previous dof vector
         * \param &v_tp1_begin: The starting iterator of the current dof vector
         * \param &v_tp1_end: The stopping iterator of the current dof vector
         * \param &vDot_t_begin: The starting iterator of the previous dof rate of change vector
         * \param &vDot_t_end: The stopping iterator of the previous dof rate of change vector
         * \param &vDDot_t_begin: The starting iterator of the previous dof acceleration (second time derivative) vector
         * \param &vDDot_t_end: The stopping iterator of the previous dof acceleration (second time derivative) vector
         * \param &alpha: The integration parameter for the velocity (0 is explicit, 1 is implicit)
         * \param &beta: The integration parameter for the acceleration (second time derivative) (0 is explicit, 1 is
         * implicit)
         * \param &vDDot_tp1_begin: The starting iterator of the current dof acceleration (second time derivative)
         * vector
         * \param &vDDot_tp1_end: The stopping iterator of the current dof acceleration (second time derivative) vector
         * \param &dVDDotdV: The derivative of the acceleration (second time derivative) dof vector with respect
         *     to the dof vector
         */
        template <typename dt_type, class v_t_in, class v_tp1_in, class vDot_t_in, class vDDot_t_in,
                  typename alpha_type, typename beta_type, class vDDot_tp1_out, typename dVDDotdV_type>
        void compute_current_acceleration(const dt_type &dt, const v_t_in &v_t_begin, const v_t_in &v_t_end,
                                          const v_tp1_in &v_tp1_begin, const v_tp1_in &v_tp1_end,
                                          const vDot_t_in &vDot_t_begin, const vDot_t_in &vDot_t_end,
                                          const vDDot_t_in &vDDot_t_begin, const vDDot_t_in &vDDot_t_end,
                                          const alpha_type &alpha, const beta_type &beta, vDDot_tp1_out vDDot_tp1_begin,
                                          vDDot_tp1_out vDDot_tp1_end, dVDDotdV_type &dVDDotdV) {
            dVDDotdV = 1. / (dt * dt * alpha * beta);

            for (unsigned int i = 0; i < (unsigned int)(v_t_end - v_t_begin); ++i) {
                *(vDDot_tp1_begin + i) = ((*(v_tp1_begin + i)) - (*(v_t_begin + i)) - dt * (*(vDot_t_begin + i)) -
                                          (dt * dt) * alpha * (1 - beta) * (*(vDDot_t_begin + i))) *
                                         dVDDotdV;
            }
        }

    }  // namespace finiteElement

}  // namespace tardigradeBalanceEquations
