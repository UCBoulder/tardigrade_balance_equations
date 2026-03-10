/**
 ******************************************************************************
 * \file tardigrade_QuadraticHex.h
 ******************************************************************************
 * The header file for the Quadratic Hexahedral element
 ******************************************************************************
 */

#ifndef TARDIGRADE_QUADRATICHEX_H
#define TARDIGRADE_QUADRATICHEX_H

#include "tardigrade_FiniteElementBase.h"

namespace tardigradeBalanceEquations {

    namespace finiteElement {

        //! An implementation of a quadratic hexahedral element
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out, class local_point_out,
                  typename weight_type>
        class QuadraticHex : public FiniteElementBase<3, 3, 20, node_in, typename std::array<T, 3 * 8>::const_iterator,
                                                      local_point_in, shape_functions_out, grad_shape_functions_out,
                                                      local_point_out, weight_type> {
           public:
            //! The local nodes for an isoparametric quadratic hex element
            constexpr static std::array<T, 3 * 20> local_nodes = {
                -1, -1, -1, 1, -1, -1, 1,  1, -1, -1, 1, -1, -1, -1, 1, 1,  -1, 1,  1, 1,
                1,  -1, 1,  1, 0,  -1, -1, 1, 0,  -1, 0, 1,  -1, -1, 0, -1, 0,  -1, 1, 1,
                0,  1,  0,  1, 1,  -1, 0,  1, -1, -1, 0, 1,  -1, 0,  1, 1,  0,  -1, 1, 0};

            //! The integration points for a fully integrated isoparametric linear hex element
            constexpr static std::array<T, 3 * 8> integration_points = {
                -0.57735027, -0.57735027, -0.57735027, 0.57735027,  -0.57735027, -0.57735027, 0.57735027, 0.57735027,
                -0.57735027, -0.57735027, 0.57735027,  -0.57735027, -0.57735027, -0.57735027, 0.57735027, 0.57735027,
                -0.57735027, 0.57735027,  0.57735027,  0.57735027,  0.57735027,  -0.57735027, 0.57735027, 0.57735027};

            //! The integration weights for a fully integrated isoparametric linear hex element
            constexpr static std::array<T, 8> integration_weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

            QuadraticHex(const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin,
                         const node_in &_X_end);

            using FiniteElementBase<3, 3, 20, node_in, typename std::array<T, 3 * 8>::const_iterator, local_point_in,
                                    shape_functions_out, grad_shape_functions_out, local_point_out,
                                    weight_type>::FiniteElementBase;

            virtual void GetShapeFunctions(const local_point_in &xi_begin, const local_point_in &xi_end,
                                           shape_functions_out N_begin, shape_functions_out N_end) override;

            virtual void GetLocalShapeFunctionGradients(const local_point_in &xi_begin, const local_point_in &xi_end,
                                                        grad_shape_functions_out dNdxi_begin,
                                                        grad_shape_functions_out dNdxi_end) override;

            virtual void GetGlobalShapeFunctionGradients(const local_point_in &xi_begin, const local_point_in &xi_end,
                                                         const node_in           &node_positions_begin,
                                                         const node_in           &node_positions_end,
                                                         grad_shape_functions_out value_begin,
                                                         grad_shape_functions_out value_end) override;

            virtual void GetVolumeIntegralJacobianOfTransformation(
                const local_point_in &xi_begin, const local_point_in &xi_end,
                typename std::iterator_traits<node_in>::value_type &value, const bool configuration = 1) override;

            virtual void GetVolumeIntegrationPointData(const unsigned int i, local_point_out xi_begin,
                                                       local_point_out xi_end, weight_type &weight) override;
        };

    }  // namespace finiteElement

}  // namespace tardigradeBalanceEquations

#include "tardigrade_QuadraticHex.tpp"

#endif
