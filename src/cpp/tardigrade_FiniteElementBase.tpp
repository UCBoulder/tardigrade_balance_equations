/**
 ******************************************************************************
 * \file tardigrade_LinearHex.tpp
 ******************************************************************************
 * The template file for the FiniteElementBase class
 ******************************************************************************
 */

namespace tardigradeBalanceEquations {

    namespace finiteElement {

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
        template <class element_configuration>
        FiniteElementBase<element_configuration
                          >::FiniteElementBase(const typename element_configuration::node_in &_x_begin, const typename element_configuration::node_in &_x_end,
                                                                       const typename element_configuration::node_in &_X_begin, const typename element_configuration::node_in &_X_end,
                                                                       const typename element_configuration::local_node_in &_local_node_xi_begin,
                                                                       const typename element_configuration::local_node_in &_local_node_xi_end)
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
        template <class element_configuration>
        void FiniteElementBase<element_configuration>::GetShapeFunctions(const typename element_configuration::local_point_in &xi_begin, const typename element_configuration::local_point_in &xi_end,
                                       typename element_configuration::shape_functions_out N_begin, typename element_configuration::shape_functions_out N_end) {

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
        template <class element_configuration>
        void FiniteElementBase<element_configuration>::GetLocalShapeFunctionGradients(const typename element_configuration::local_point_in &xi_begin, const typename element_configuration::local_point_in &xi_end,
                                                    typename element_configuration::grad_shape_functions_out dNdxi_begin,
                                                    typename element_configuration::grad_shape_functions_out dNdxi_end) {

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
        template <class element_configuration>
        void FiniteElementBase<element_configuration>::GetGlobalShapeFunctionGradients(const typename element_configuration::local_point_in &xi_begin, const typename element_configuration::local_point_in &xi_end,
                                                     const typename element_configuration::node_in           &node_positions_begin,
                                                     const typename element_configuration::node_in           &node_positions_end,
                                                     typename element_configuration::grad_shape_functions_out value_begin,
                                                     typename element_configuration::grad_shape_functions_out value_end) {

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
        template <class element_configuration>
        void FiniteElementBase<element_configuration>::GetVolumeIntegralJacobianOfTransformation(
            const typename element_configuration::local_point_in &xi_begin, const typename element_configuration::local_point_in &xi_end,
            typename std::iterator_traits<typename element_configuration::node_in>::value_type &value, const bool configuration) {

            throw std::logic_error("Function not implemented");
        }

        /*!
         * Compute the value of the Jacobian of transformation from the local coordinates to the configuration
         * for surface integrals
         *
         * \param s: The surface to compute the Jacobian of transformation for
         * \param &xi_begin: The starting iterator of the local coordinates
         * \param &xi_end: The stopping iterator of the local coordinates
         * \param value: The Jacobian of transformation going from the local coordinates to the indicated
         * configuration
         * \param configuration: Compute the gradient w.r.t. the current configuration ( true ) or reference
         * configuration ( false )
         */
        template <class element_configuration>
        void FiniteElementBase<element_configuration>::GetSurfaceIntegralJacobianOfTransformation(
                const unsigned int s,
                const typename element_configuration::local_point_in &xi_begin, const typename element_configuration::local_point_in &xi_end,
                typename std::iterator_traits<typename element_configuration::node_in>::value_type &value, const bool configuration) {

            throw std::logic_error("Function not implemented");
        }

        /*!
         * Get the integration point information for a volumetric integral
         *
         * \param i: The index of the integration point
         * \param xi_begin: The starting iterator of the local coordinates of the integration point
         * \param xi_end: The stopping iterator of the local coordinates of the integration point
         * \param &weight: The weight to be applied to the integration point
         */
        template <class element_configuration>
        void FiniteElementBase<element_configuration>::GetVolumeIntegrationPointData(const unsigned int i, typename element_configuration::local_point_out xi_begin, typename element_configuration::local_point_out xi_end, typename element_configuration::volume_integration_point_weight_value_type &weight){

            throw std::logic_error("Function not implemented");
        }

        /*!
         * Get the integration point information for a surface
         *
         * \param s: The index of the surface
         * \param i: The index of the integration point
         * \param xi_begin: The starting iterator of the local coordinates of the integration point
         * \param xi_end: The stopping iterator of the local coordinates of the integration point
         * \param &weight: The weight to be applied to the integration point
         */
        template <class element_configuration>
        void FiniteElementBase<element_configuration>::GetSurfaceIntegrationPointData(const unsigned int s, const unsigned int i, typename element_configuration::local_point_out xi_begin,
                                                    typename element_configuration::local_point_out xi_end, typename element_configuration::surface_integration_point_weight_value_type &weight){

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
        template <class element_configuration>
        template <class quantity_in, class quantity_out>
        void FiniteElementBase<element_configuration>::InterpolateQuantity(const typename element_configuration::local_point_in &xi_begin, const typename element_configuration::local_point_in &xi_end,
                                 const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                 quantity_out value_begin, quantity_out value_end) {

            const size_type quantity_dim = (size_type)(value_end - value_begin);

            TARDIGRADE_ERROR_TOOLS_CHECK(quantity_dim * element_configuration::node_count == (size_type)(quantity_end - quantity_begin),
                                         "The returned value size (" + std::to_string(quantity_dim) +
                                             ") and the quantity dimension (" +
                                             std::to_string((size_type)(quantity_end - quantity_begin)) +
                                             ") are inconsistent with the node count (" +
                                             std::to_string(element_configuration::node_count) + ")");

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
        template <class element_configuration>
        template <class quantity_in, class quantity_gradient_out>
        void FiniteElementBase<element_configuration>::GetLocalQuantityGradient(const typename element_configuration::local_point_in &xi_begin, const typename element_configuration::local_point_in &xi_end,
                                      const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                      quantity_gradient_out value_begin, quantity_gradient_out value_end) {

            const size_type quantity_dim = (size_type)(value_end - value_begin) / element_configuration::local_dim;

            TARDIGRADE_ERROR_TOOLS_CHECK(quantity_dim * element_configuration::node_count == (size_type)(quantity_end - quantity_begin),
                                         "The returned value size (" + std::to_string(quantity_dim) +
                                             ") and the quantity dimension (" +
                                             std::to_string((size_type)(quantity_end - quantity_begin)) +
                                             ") are inconsistent with the node count (" +
                                             std::to_string(element_configuration::node_count) + ")");

            TARDIGRADE_ERROR_TOOLS_CATCH(GetLocalShapeFunctionGradients(xi_begin, xi_end,
                                                                        std::begin(_local_gradshapefunctions),
                                                                        std::end(_local_gradshapefunctions)));

            std::fill(value_begin, value_end, 0);

            for (unsigned int node = 0; node < element_configuration::node_count; ++node) {
                for (unsigned int row = 0; row < quantity_dim; ++row) {
                    for (unsigned int col = 0; col < element_configuration::local_dim; ++col) {
                        *(value_begin + element_configuration::local_dim * row + col) +=
                            _local_gradshapefunctions[element_configuration::local_dim * node + col] *
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
        template <class element_configuration>
        template <class quantity_in, class quantity_gradient_out>
        void FiniteElementBase<element_configuration>::GetGlobalQuantityGradient(const typename element_configuration::local_point_in &xi_begin, const typename element_configuration::local_point_in &xi_end,
                                       const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                       quantity_gradient_out value_begin, quantity_gradient_out value_end,
                                       const bool configuration) {

            const size_type quantity_dim = (size_type)(value_end - value_begin) / element_configuration::dim;

            TARDIGRADE_ERROR_TOOLS_CHECK(quantity_dim * element_configuration::node_count == (size_type)(quantity_end - quantity_begin),
                                         "The returned value size (" + std::to_string(quantity_dim) +
                                             ") and the quantity dimension (" +
                                             std::to_string((size_type)(quantity_end - quantity_begin)) +
                                             ") are inconsistent with the node count (" +
                                             std::to_string(element_configuration::node_count) + ")");

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

            for (unsigned int node = 0; node < element_configuration::node_count; ++node) {
                for (unsigned int row = 0; row < quantity_dim; ++row) {
                    for (unsigned int col = 0; col < element_configuration::local_dim; ++col) {
                        *(value_begin + element_configuration::local_dim * row + col) +=
                            _global_gradshapefunctions[element_configuration::local_dim * node + col] *
                            (*(quantity_begin + quantity_dim * node + row));
                    }
                }
            }
        }

    }

}
