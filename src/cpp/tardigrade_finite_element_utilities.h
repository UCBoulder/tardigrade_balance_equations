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

    namespace finiteElementUtilities {

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

        //! A base class for a simple finite element formulation useful for testing
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
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

    }  // namespace finiteElementUtilities

}  // namespace tardigradeBalanceEquations

#include "tardigrade_finite_element_utilities.tpp"

#endif
