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

        /*!
         * The configuration of the linear hex element
         */
        class QuadraticHexConfiguration : public FiniteElementConfigurationBase<3, 3, 20, double>{

            public:

                //! The number of integration points
                constexpr static unsigned int num_volume_integration_points = 8;

                //! The type for the element local node coordinates
                using node_value_type = double;

                //! The type for the volume integration point weights
                using volume_integration_point_weight_value_type = double;
        };

        //! An implementation of a quadratic hexahedral element
        template <class element_configuration, typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out, class local_point_out, typename weight_type>
        class QuadraticHex
            : public FiniteElementBase<element_configuration, node_in, local_point_in,
                                       shape_functions_out, grad_shape_functions_out, local_point_out, weight_type> {
           public:
            //! The local nodes for an isoparametric quadratic hex element
            constexpr static std::array<typename element_configuration::local_node_value_type, element_configuration::local_dim * element_configuration::node_count> local_nodes = {
                -1, -1, -1, 1, -1, -1, 1,  1, -1, -1, 1, -1, -1, -1, 1, 1,  -1, 1,  1, 1,
                1,  -1, 1,  1, 0,  -1, -1, 1, 0,  -1, 0, 1,  -1, -1, 0, -1, 0,  -1, 1, 1,
                0,  1,  0,  1, 1,  -1, 0,  1, -1, -1, 0, 1,  -1, 0,  1, 1,  0,  -1, 1, 0};

            //! The integration points for a fully integrated isoparametric linear hex element
            constexpr static std::array<typename element_configuration::local_node_value_type, element_configuration::local_dim * element_configuration::num_volume_integration_points> volume_integration_points = {
                -0.57735027, -0.57735027, -0.57735027, 0.57735027,  -0.57735027, -0.57735027, 0.57735027, 0.57735027,
                -0.57735027, -0.57735027, 0.57735027,  -0.57735027, -0.57735027, -0.57735027, 0.57735027, 0.57735027,
                -0.57735027, 0.57735027,  0.57735027,  0.57735027,  0.57735027,  -0.57735027, 0.57735027, 0.57735027};

            //! The integration weights for a fully integrated isoparametric linear hex element
            constexpr static std::array<typename element_configuration::volume_integration_point_weight_value_type, element_configuration::num_volume_integration_points> volume_integration_weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

            //! The surface integration points for a fully integrated isoparametric linear hex element
            constexpr static std::array<T, 6 * 4 * 3> surface_integration_points = {
                -1,          -0.57735027, -0.57735027, -1,         0.57735027,  -0.57735027, -1,          0.57735027,
                0.57735027,  -1,          -0.57735027, 0.57735027, 1,           -0.57735027, -0.57735027, 1,
                0.57735027,  -0.57735027, 1,           0.57735027, 0.57735027,  1,           -0.57735027, 0.57735027,
                -0.57735027, -1,          -0.57735027, 0.57735027, -1,          -0.57735027, 0.57735027,  -1,
                0.57735027,  -0.57735027, -1,          0.57735027, -0.57735027, 1,           -0.57735027, 0.57735027,
                1,           -0.57735027, 0.57735027,  1,          0.57735027,  -0.57735027, 1,           0.57735027,
                -0.57735027, -0.57735027, -1,          0.57735027, -0.57735027, -1,          0.57735027,  0.57735027,
                -1,          -0.57735027, 0.57735027,  -1,         -0.57735027, -0.57735027, 1,           0.57735027,
                -0.57735027, 1,           0.57735027,  0.57735027, 1,           -0.57735027, 0.57735027,  1};

            //! The surface integration weights for a fully integrated isoparametric linear hex element
            constexpr static std::array<T, 4 * 6> surface_integration_weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

            //! The surface normals
            constexpr static std::array<T, 3 * 6> surface_normals = {-1, 0, 0, 1, 0, 0,  0, -1, 0,
                                                                     0,  1, 0, 0, 0, -1, 0, 0,  1};

            QuadraticHex(const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin,
                         const node_in &_X_end);

            using FiniteElementBase<element_configuration, node_in, local_point_in,
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

            virtual void GetSurfaceIntegralJacobianOfTransformation(
                const unsigned int s, const local_point_in &xi_begin, const local_point_in &xi_end,
                typename std::iterator_traits<node_in>::value_type &value, const bool configuration = 1) override;

            virtual void GetVolumeIntegrationPointData(const unsigned int i, local_point_out xi_begin,
                                                       local_point_out xi_end, weight_type &weight) override;

            virtual void GetSurfaceIntegrationPointData(const unsigned int s, const unsigned int i,
                                                        local_point_out xi_begin, local_point_out xi_end,
                                                        weight_type &weight) override;
        };

    }  // namespace finiteElement

}  // namespace tardigradeBalanceEquations

#include "tardigrade_QuadraticHex.tpp"

#endif
