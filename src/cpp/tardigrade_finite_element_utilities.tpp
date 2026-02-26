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

    namespace finiteElementUtilities {

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
         * Constructor for the base finite element class
         *
         * \param &_x_begin: The starting iterator for the current node positions
         * \param &_x_end: The stopping iterator for the current node positions
         * \param &_X_begin: The starting iterator for the reference node positions
         * \param &_X_end: The stopping iterator for the reference node positions
         * \param &_local_node_xi_begin: The starting iterator for the local node positions
         * \param &_local_node_xi_end: The stopping iterator for the local node positions
         */
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
        FiniteElementBase<dim, local_dim, node_count, node_in, local_node_in, local_point_in, shape_functions_out,
                          grad_shape_functions_out>::FiniteElementBase(const node_in &_x_begin, const node_in &_x_end,
                                                                       const node_in &_X_begin, const node_in &_X_end,
                                                                       const local_node_in &_local_node_xi_begin,
                                                                       const local_node_in &_local_node_xi_end)
            : x_begin(_x_begin),
              x_end(_x_end),
              X_begin(_X_begin),
              X_end(_X_end),
              local_node_xi_begin(_local_node_xi_begin),
              local_node_xi_end(_local_node_xi_end) {
        }

        /*!
         * Get the values of the shape functions
         *
         * \param &xi_begin: The starting iterator of the local point
         * \param &xi_end: The stopping iterator of the local point
         * \param N_begin: The starting iterator of the shape functions
         * \param N_end: The stopping iterator of the shape functions
         */
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
        void FiniteElementBase<dim, local_dim, node_count, node_in, local_node_in, local_point_in,shape_functions_out,grad_shape_functions_out>::GetShapeFunctions(const local_point_in &xi_begin, const local_point_in &xi_end,
                                       shape_functions_out N_begin, shape_functions_out N_end) {

            throw std::logic_error("Function not implemented");
        }

        /*!
         * Get the values of the gradients of the shape functions w.r.t. the local coordinates
         *
         * \param &xi_begin: The starting iterator of the local point
         * \param &xi_end: The stopping iterator of the local point
         * \param dNdxi_begin: The starting iterator of the shape function gradients
         * \param dNdxi_end: The stopping iterator of the shape function gradients
         */
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
        void FiniteElementBase<dim, local_dim, node_count, node_in, local_node_in, local_point_in,shape_functions_out,grad_shape_functions_out>::GetLocalShapeFunctionGradients(const local_point_in &xi_begin, const local_point_in &xi_end,
                                                    grad_shape_functions_out dNdxi_begin,
                                                    grad_shape_functions_out dNdxi_end) {

            throw std::logic_error("Function not implemented");
        }

        /*!
         * Compute the global gradient of the shape functions
         *
         * \param &xi_begin: The starting iterator of the local point
         * \param &xi_end: The stopping iterator of the local point
         * \param &node_positions_begin: The starting iterator of the nodal positions (row major)
         * \param &node_positions_end: The stopping iterator of the nodal positions (row major)
         * \param &value_begin: The starting iterator of the shape function global gradient (row major)
         * \param &value_end: The stopping iterator of the shape function global gradient (row major)
         */
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
        void FiniteElementBase<dim, local_dim, node_count, node_in, local_node_in, local_point_in,shape_functions_out,grad_shape_functions_out>::GetGlobalShapeFunctionGradients(const local_point_in &xi_begin, const local_point_in &xi_end,
                                                     const node_in           &node_positions_begin,
                                                     const node_in           &node_positions_end,
                                                     grad_shape_functions_out value_begin,
                                                     grad_shape_functions_out value_end) {

            throw std::logic_error("Function not implemented");
        }

        /*!
         * Compute the value of the Jacobian of transformation from the local coordinates to the configuration
         * for volume integrals
         *
         * \param &xi_begin: The starting iterator of the local coordinates
         * \param &xi_end: The stopping iterator of the local coordinates
         * \param value: The Jacobian of transformation going from the local coordinates to the indicated
         * configuration
         * \param configuration: Compute the gradient w.r.t. the current configuration ( true ) or reference
         * configuration ( false )
         */
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
        void FiniteElementBase<dim, local_dim, node_count, node_in, local_node_in, local_point_in,shape_functions_out,grad_shape_functions_out>::GetVolumeIntegralJacobianOfTransformation(
            const local_point_in &xi_begin, const local_point_in &xi_end,
            typename std::iterator_traits<node_in>::value_type &value, const bool configuration) {

            throw std::logic_error("Function not implemented");
        }

        /*!
         * Interpolate the provided quantity to the local point
         *
         * \param &xi_begin: The starting iterator of the local coordinates
         * \param &xi_end: The stopping iterator of the local coordinates
         * \param &quantity_begin: The starting iterator of the quantity at the nodes (row-major)
         * \param &quantity_end: The stopping iterator of the quantity at the nodes (row-major)
         * \param value_begin: The starting iterator for the interpolated value
         * \param value_end: The stopping iterator for the interpolated value
         */
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
        template <class quantity_in, class quantity_out>
        void FiniteElementBase<dim, local_dim, node_count, node_in, local_node_in, local_point_in,shape_functions_out,grad_shape_functions_out>::InterpolateQuantity(const local_point_in &xi_begin, const local_point_in &xi_end,
                                 const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                 quantity_out value_begin, quantity_out value_end) {

            const size_type quantity_dim = (size_type)(value_end - value_begin);

            TARDIGRADE_ERROR_TOOLS_CHECK(quantity_dim * node_count == (size_type)(quantity_end - quantity_begin),
                                         "The returned value size (" + std::to_string(quantity_dim) +
                                             ") and the quantity dimension (" +
                                             std::to_string((size_type)(quantity_end - quantity_begin)) +
                                             ") are inconsistent with the node count (" +
                                             std::to_string(node_count) + ")");

            GetShapeFunctions(xi_begin, xi_end, std::begin(_shapefunctions), std::end(_shapefunctions));

            std::fill(value_begin, value_end, 0);

            for (auto N = std::begin(_shapefunctions); N != std::end(_shapefunctions); ++N) {
                unsigned int node = (size_type)(N - std::begin(_shapefunctions));

                for (auto v = value_begin; v != value_end; ++v) {
                    *v += (*N) * (*(quantity_begin + quantity_dim * node + (size_type)(v - value_begin)));
                }
            }
        }

        /*!
         * Compute the gradient of the quantity at a local point
         *
         * \param &xi_begin: The starting iterator of the local coordinates
         * \param &xi_end: The stopping iterator of the local coordinates
         * \param &quantity_begin: The starting iterator of the quantity at the nodes (row-major)
         * \param &quantity_end: The stopping iterator of the quantity at the nodes (row-major)
         * \param value_begin: The starting iterator for the computed local gradient of the quantity
         * \param value_end: The stopping iterator for the computed local gradient of the quantity in row-major
         * form
         */
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
        template <class quantity_in, class quantity_gradient_out>
        void FiniteElementBase<dim, local_dim, node_count, node_in, local_node_in, local_point_in,shape_functions_out,grad_shape_functions_out>::GetLocalQuantityGradient(const local_point_in &xi_begin, const local_point_in &xi_end,
                                      const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                      quantity_gradient_out value_begin, quantity_gradient_out value_end) {

            const size_type quantity_dim = (size_type)(value_end - value_begin) / local_dim;

            TARDIGRADE_ERROR_TOOLS_CHECK(quantity_dim * node_count == (size_type)(quantity_end - quantity_begin),
                                         "The returned value size (" + std::to_string(quantity_dim) +
                                             ") and the quantity dimension (" +
                                             std::to_string((size_type)(quantity_end - quantity_begin)) +
                                             ") are inconsistent with the node count (" +
                                             std::to_string(node_count) + ")");

            TARDIGRADE_ERROR_TOOLS_CATCH(GetLocalShapeFunctionGradients(xi_begin, xi_end,
                                                                        std::begin(_local_gradshapefunctions),
                                                                        std::end(_local_gradshapefunctions)));

            std::fill(value_begin, value_end, 0);

            for (unsigned int node = 0; node < node_count; ++node) {
                for (unsigned int row = 0; row < quantity_dim; ++row) {
                    for (unsigned int col = 0; col < local_dim; ++col) {
                        *(value_begin + local_dim * row + col) +=
                            _local_gradshapefunctions[local_dim * node + col] *
                            (*(quantity_begin + quantity_dim * node + row));
                    }
                }
            }
        }

        /*!
         * Compute the global gradient of the quantity at a local point
         *
         * \param &xi_begin: The starting iterator of the local coordinates
         * \param &xi_end: The stopping iterator of the local coordinates
         * \param &quantity_begin: The starting iterator of the quantity at the nodes (row-major)
         * \param &quantity_end: The stopping iterator of the quantity at the nodes (row-major)
         * \param value_begin: The starting iterator for the computed global gradient of the quantity
         * \param value_end: The stopping iterator for the computed global gradient of the quantity in row-major
         * form
         * \param configuration: Compute the gradient w.r.t. the current configuration ( true ) or reference
         * configuration ( false )
         */
        template <int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in,
                  class shape_functions_out, class grad_shape_functions_out>
        template <class quantity_in, class quantity_gradient_out>
        void FiniteElementBase<dim, local_dim, node_count, node_in, local_node_in, local_point_in,shape_functions_out,grad_shape_functions_out>::GetGlobalQuantityGradient(const local_point_in &xi_begin, const local_point_in &xi_end,
                                       const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                       quantity_gradient_out value_begin, quantity_gradient_out value_end,
                                       const bool configuration) {

            const size_type quantity_dim = (size_type)(value_end - value_begin) / dim;

            TARDIGRADE_ERROR_TOOLS_CHECK(quantity_dim * node_count == (size_type)(quantity_end - quantity_begin),
                                         "The returned value size (" + std::to_string(quantity_dim) +
                                             ") and the quantity dimension (" +
                                             std::to_string((size_type)(quantity_end - quantity_begin)) +
                                             ") are inconsistent with the node count (" +
                                             std::to_string(node_count) + ")");

            if (configuration) {
                TARDIGRADE_ERROR_TOOLS_CATCH(GetGlobalShapeFunctionGradients(xi_begin, xi_end, x_begin, x_end,
                                                                             std::begin(_global_gradshapefunctions),
                                                                             std::end(_global_gradshapefunctions)));

            } else {
                TARDIGRADE_ERROR_TOOLS_CATCH(GetGlobalShapeFunctionGradients(xi_begin, xi_end, X_begin, X_end,
                                                                             std::begin(_global_gradshapefunctions),
                                                                             std::end(_global_gradshapefunctions)));
            }

            std::fill(value_begin, value_end, 0);

            for (unsigned int node = 0; node < node_count; ++node) {
                for (unsigned int row = 0; row < quantity_dim; ++row) {
                    for (unsigned int col = 0; col < local_dim; ++col) {
                        *(value_begin + local_dim * row + col) +=
                            _global_gradshapefunctions[local_dim * node + col] *
                            (*(quantity_begin + quantity_dim * node + row));
                    }
                }
            }
        }

    }  // namespace finiteElementUtilities

}  // namespace tardigradeBalanceEquations
