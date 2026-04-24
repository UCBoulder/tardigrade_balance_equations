/**
 ******************************************************************************
 * \file tardigrade_finite_element_utilities.h
 ******************************************************************************
 * The header file for utilities which can assist using the balance equations
 * in finite element codes. We here assume that unknown quantities are
 * computed at the evaluation point using interpolation functions (here called
 * interp) and are projected to the nodes using test functions (here called
 * test).
 ******************************************************************************
 */

#ifndef TARDIGRADE_FINITE_ELEMENT_UTILITIES_H
#define TARDIGRADE_FINITE_ELEMENT_UTILITIES_H

#include <array>

#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations {

    namespace finiteElement {

        typedef unsigned int size_type;  //!< Define unsigned int as the default size type

        constexpr unsigned int dim = 3;  //!< Set the dimension as 3D by default

        constexpr unsigned int sot_dim = dim * dim;  //!< Set the dimensions of a standard second order tensor

        typedef double floatType;  //!< Define the float type as a double

        typedef std::array<floatType, dim> floatVector;  //!< Define a standard vector

        typedef std::array<floatType, sot_dim> secondOrderTensor;  //!< Define a standard second-order tensor

        template <typename grad_iterator, typename output_iterator>
        void computeGradientSpatialJacobian(const grad_iterator &grad_a_start, const unsigned int grad_a_size,
                                            floatVector grad_interp, const unsigned int index,
                                            output_iterator dgrad_adui_start);

        template <typename dt_type, class v_t_in, class v_tp1_in, class vDot_t_in, typename alpha_type, class vDot_tp1_out,
                  typename dVDotdV_type>
        void compute_current_rate_of_change(const dt_type &dt, const v_t_in &v_t_begin, const v_t_in &v_t_end,
                                            const v_tp1_in &v_tp1_begin, const v_tp1_in &v_tp1_end,
                                            const vDot_t_in &vDot_t_begin, const vDot_t_in &vDot_t_end, const alpha_type &alpha,
                                            vDot_tp1_out vDot_tp1_begin, vDot_tp1_out vDot_tp1_end, dVDotdV_type &dVDotdV);
        
        template <typename dt_type, class v_t_in, class v_tp1_in, class vDot_t_in, class vDDot_t_in, typename alpha_type,
                  typename beta_type, class vDDot_tp1_out, typename dVDDotdV_type>
        void compute_current_acceleration(const dt_type &dt, const v_t_in &v_t_begin, const v_t_in &v_t_end,
                                          const v_tp1_in &v_tp1_begin, const v_tp1_in &v_tp1_end, const vDot_t_in &vDot_t_begin,
                                          const vDot_t_in &vDot_t_end, const vDDot_t_in &vDDot_t_begin,
                                          const vDDot_t_in &vDDot_t_end, const alpha_type &alpha, const beta_type &beta,
                                          vDDot_tp1_out vDDot_tp1_begin, vDDot_tp1_out vDDot_tp1_end, dVDDotdV_type &dVDDotdV);

    }  // namespace finiteElement

}  // namespace tardigradeBalanceEquations

#include "tardigrade_finite_element_utilities.tpp"

#endif
