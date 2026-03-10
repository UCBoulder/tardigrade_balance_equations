/**
 ******************************************************************************
 * \file tardigrade_FiniteElementBase.h
 ******************************************************************************
 * The header file for the FiniteElementBase class
 ******************************************************************************
 */

#ifndef TARDIGRADE_FINITEELEMENTBASE_H
#define TARDIGRADE_FINITEELEMENTBASE_H

#include "tardigrade_finite_element_utilities.h"

namespace tardigradeBalanceEquations {

    namespace finiteElement {

        //! A base class for a simple finite element formulation useful for testing
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out, class local_point_out,
                  typename weight_type>
        class FiniteElementBase {
           public:
            FiniteElementBase(const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin,
                              const node_in &_X_end, const local_node_in &_local_node_xi_begin,
                              const local_node_in &_local_node_xi_end);

            virtual void GetShapeFunctions(const local_point_in &xi_begin, const local_point_in &xi_end,
                                           shape_functions_out N_begin, shape_functions_out N_end);

            virtual void GetLocalShapeFunctionGradients(const local_point_in &xi_begin, const local_point_in &xi_end,
                                                        grad_shape_functions_out dNdxi_begin,
                                                        grad_shape_functions_out dNdxi_end);

            virtual void GetGlobalShapeFunctionGradients(const local_point_in &xi_begin, const local_point_in &xi_end,
                                                         const node_in           &node_positions_begin,
                                                         const node_in           &node_positions_end,
                                                         grad_shape_functions_out value_begin,
                                                         grad_shape_functions_out value_end);

            virtual void GetVolumeIntegralJacobianOfTransformation(
                const local_point_in &xi_begin, const local_point_in &xi_end,
                typename std::iterator_traits<node_in>::value_type &value, const bool configuration = 1);

            virtual void GetSurfaceIntegralJacobianOfTransformation(
                const unsigned int s,
                const local_point_in &xi_begin, const local_point_in &xi_end,
                typename std::iterator_traits<node_in>::value_type &value, const bool configuration = 1);

            virtual void GetVolumeIntegrationPointData(const unsigned int i, local_point_out xi_begin,
                                                       local_point_out xi_end, weight_type &weight);

            virtual void GetSurfaceIntegrationPointData(const unsigned int s, const unsigned int i,
                                                        local_point_out xi_begin, local_point_out xi_end,
                                                        weight_type &weight);

            template <class quantity_in, class quantity_out>
            void InterpolateQuantity(const local_point_in &xi_begin, const local_point_in &xi_end,
                                     const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                     quantity_out value_begin, quantity_out value_end);

            template <class quantity_in, class quantity_gradient_out>
            void GetLocalQuantityGradient(const local_point_in &xi_begin, const local_point_in &xi_end,
                                          const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                          quantity_gradient_out value_begin, quantity_gradient_out value_end);

            template <class quantity_in, class quantity_gradient_out>
            void GetGlobalQuantityGradient(const local_point_in &xi_begin, const local_point_in &xi_end,
                                           const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                           quantity_gradient_out value_begin, quantity_gradient_out value_end,
                                           const bool configuration = true);

           protected:
            const node_in x_begin;  //!< Starting iterator for the current position of the nodal coordinates
            const node_in x_end;    //!< Stopping iterator for the current position of the nodal coordinates

            const node_in X_begin;  //!< Starting iterator for the reference position of the nodal coordinates
            const node_in X_end;    //!< Stopping iterator for the reference position of the nodal coordinates

            const local_node_in local_node_xi_begin;  //!< Starting iterator for the local nodal coordinates
            const local_node_in local_node_xi_end;    //!< Stopping iterator for the local nodal coordinates

            std::array<typename std::iterator_traits<shape_functions_out>::value_type, node_count>
                _shapefunctions;  //!< A temporary storage container for shapefunction values

            std::array<typename std::iterator_traits<shape_functions_out>::value_type, local_dim * node_count>
                _local_gradshapefunctions;  //!< A temporary storage container for local grad shapefunction values

            std::array<typename std::iterator_traits<shape_functions_out>::value_type, dim * node_count>
                _global_gradshapefunctions;  //!< A temporary storage container for global grad shapefunction values
        };

    }  // namespace finiteElement

}  // namespace tardigradeBalanceEquations

#include "tardigrade_FiniteElementBase.tpp"

#endif
