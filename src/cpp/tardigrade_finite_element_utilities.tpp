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

        /*!
         * Constructor for the linear 8 node hex finite element class
         *
         * \param &_x_begin: The starting iterator for the current node positions
         * \param &_x_end: The stopping iterator for the current node positions
         * \param &_X_begin: The starting iterator for the reference node positions
         * \param &_X_end: The stopping iterator for the reference node positions
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        LinearHex<T, node_in, local_point_in, shape_functions_out, grad_shape_functions_out>::LinearHex(
            const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin, const node_in &_X_end)
            : FiniteElementBase<3, 3, 8, node_in, typename std::array<T, 3 * 8>::const_iterator, local_point_in,
                                shape_functions_out, grad_shape_functions_out>(_x_begin, _x_end, _X_begin, _X_end,
                                                                               std::cbegin(local_nodes),
                                                                               std::cend(local_nodes)) {
        }

        /*!
         * Get the shape functions of the linear hexahedral element
         *
         * \param &xi_begin: The starting iterator of the local coordinates (must have dimension 3)
         * \param &xi_end: The stopping iterator of the local coordinates (must have dimension 3)
         * \param N_begin: The starting iterator of the shape functions (must have dimension 8)
         * \param N_end: The stopping iterator of the shape functions (must have dimension 8)
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        void LinearHex<T, node_in, local_point_in, shape_functions_out, grad_shape_functions_out>::GetShapeFunctions(
            const local_point_in &xi_begin, const local_point_in &xi_end, shape_functions_out N_begin,
            shape_functions_out N_end) {

            TARDIGRADE_ERROR_TOOLS_CHECK((size_type)(N_end - N_begin) == 8,
                                         "The dimension of the shape-function iterator is " +
                                             std::to_string((size_type)(N_end - N_begin)));

            for (std::pair<size_type, shape_functions_out> i(0, N_begin); i.second != N_end; ++i.first, ++i.second) {
                *i.second = (1 + (*(this->local_node_xi_begin + 3 * i.first + 0)) * (*(xi_begin + 0))) *
                            (1 + (*(this->local_node_xi_begin + 3 * i.first + 1)) * (*(xi_begin + 1))) *
                            (1 + (*(this->local_node_xi_begin + 3 * i.first + 2)) * (*(xi_begin + 2))) / 8;
            }
        }

        /*!
         * Get the shape functions of the linear hexahedral element
         *
         * \param &xi_begin: The starting iterator of the local coordinates (must have dimension 3)
         * \param &xi_end: The stopping iterator of the local coordinates (must have dimension 3)
         * \param dNdxi_begin: The starting iterator of the local gradient of the shape functions (must have
         * dimension 24)
         * \param dNdxi_end: The stopping iterator of the local gradient of the shape functions (must have dimension
         * 24)
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        void LinearHex<T, node_in, local_point_in, shape_functions_out,
                       grad_shape_functions_out>::GetLocalShapeFunctionGradients(const local_point_in    &xi_begin,
                                                                                 const local_point_in    &xi_end,
                                                                                 grad_shape_functions_out dNdxi_begin,
                                                                                 grad_shape_functions_out dNdxi_end) {

            TARDIGRADE_ERROR_TOOLS_CHECK((size_type)(dNdxi_end - dNdxi_begin) == 24,
                                         "The dimension of the shape-function iterator is " +
                                             std::to_string((size_type)(dNdxi_end - dNdxi_begin)));

            for (size_type i = 0; i < 8; ++i) {
                *(dNdxi_begin + 3 * i + 0) = (*(this->local_node_xi_begin + 3 * i + 0)) *
                                             (1 + (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1))) *
                                             (1 + (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) / 8;
                *(dNdxi_begin + 3 * i + 1) = (1 + (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0))) *
                                             (*(this->local_node_xi_begin + 3 * i + 1)) *
                                             (1 + (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) / 8;
                *(dNdxi_begin + 3 * i + 2) = (1 + (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0))) *
                                             (1 + (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1))) *
                                             (*(this->local_node_xi_begin + 3 * i + 2)) / 8;
            }
        }

        /*!
         * Compute the global gradient of the quantity in the global coordinates
         *
         * \param &xi_begin: The starting iterator of the local point
         * \param &xi_end: The stopping iterator of the local point
         * \param &node_positions_begin: The starting iterator of the nodal positions (row major)
         * \param &node_positions_end: The stopping iterator of the nodal positions (row major)
         * \param &value_begin: The starting iterator of the shape function global gradient (row major)
         * \param &value_end: The stopping iterator of the shape function global gradient (row major)
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        void LinearHex<T, node_in, local_point_in, shape_functions_out,
                       grad_shape_functions_out>::GetGlobalShapeFunctionGradients(const local_point_in &xi_begin,
                                                                                  const local_point_in &xi_end,
                                                                                  const node_in &node_positions_begin,
                                                                                  const node_in &node_positions_end,
                                                                                  grad_shape_functions_out value_begin,
                                                                                  grad_shape_functions_out value_end) {

            TARDIGRADE_ERROR_TOOLS_CHECK((size_type)(value_end - value_begin) == 24,
                                         "The shape function global gradient must have a size of 24");

            std::array<typename std::iterator_traits<node_in>::value_type, 9> dxdxi;
            std::array<typename std::iterator_traits<node_in>::value_type, 9> dxidx;

            TARDIGRADE_ERROR_TOOLS_CATCH(this->GetLocalQuantityGradient(xi_begin, xi_end, node_positions_begin,
                                                                        node_positions_end, std::begin(dxdxi),
                                                                        std::end(dxdxi)));

            Eigen::Map<const Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            Eigen::Map<Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _dxidx(dxidx.data());

            _dxidx = (_dxdxi.inverse()).eval();

            std::fill(value_begin, value_end, 0);

            for (unsigned int node = 0; node < 8; ++node) {
                for (unsigned int inner = 0; inner < 3; ++inner) {
                    for (unsigned int outer = 0; outer < 3; ++outer) {
                        *(value_begin + 3 * node + outer) +=
                            this->_local_gradshapefunctions[3 * node + inner] * dxidx[3 * inner + outer];
                    }
                }
            }
        }

        /*!
         * Compute the value of the Jacobian of transformation from the local coordinates to the configuration for
         * volume integrals
         *
         * \param &xi_begin: The starting iterator of the local coordinates
         * \param &xi_end: The stopping iterator of the local coordinates
         * \param value: The Jacobian of transformation going from the local coordinates to the indicated
         * configuration
         * \param configuration: Compute the gradient w.r.t. the current configuration ( true ) or reference
         * configuration ( false )
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        void LinearHex<T, node_in, local_point_in, shape_functions_out, grad_shape_functions_out>::
            GetVolumeIntegralJacobianOfTransformation(const local_point_in &xi_begin, const local_point_in &xi_end,
                                                      typename std::iterator_traits<node_in>::value_type &value,
                                                      const bool configuration) {

            std::array<typename std::iterator_traits<node_in>::value_type, 9> dxdxi;

            if (configuration) {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->x_begin, this->x_end, std::begin(dxdxi),
                                               std::end(dxdxi));

            } else {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->X_begin, this->X_end, std::begin(dxdxi),
                                               std::end(dxdxi));
            }

            Eigen::Map<const Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            value = _dxdxi.determinant();
        }

        /*!
         * Constructor for the quadratic 20 node hex element class
         *
         * \param &_x_begin: The starting iterator for the current node positions
         * \param &_x_end: The stopping iterator for the current node positions
         * \param &_X_begin: The starting iterator for the reference node positions
         * \param &_X_end: The stopping iterator for the reference node positions
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        QuadraticHex<T, node_in, local_point_in, shape_functions_out, grad_shape_functions_out>::QuadraticHex(
            const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin, const node_in &_X_end)
            : FiniteElementBase<3, 3, 20, node_in, typename std::array<T, 3 * 20>::const_iterator, local_point_in,
                                shape_functions_out, grad_shape_functions_out>(_x_begin, _x_end, _X_begin, _X_end,
                                                                               std::cbegin(local_nodes),
                                                                               std::cend(local_nodes)) {
        }

        /*!
         * Get the shape functions of the quadratic 20 node hexahedral element
         *
         * \param &xi_begin: The starting iterator of the local coordinates (must have dimension 3)
         * \param &xi_end: The stopping iterator of the local coordinates (must have dimension 3)
         * \param N_begin: The starting iterator of the shape functions (must have dimension 8)
         * \param N_end: The stopping iterator of the shape functions (must have dimension 8)
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        void QuadraticHex<T, node_in, local_point_in, shape_functions_out, grad_shape_functions_out>::GetShapeFunctions(
            const local_point_in &xi_begin, const local_point_in &xi_end, shape_functions_out N_begin,
            shape_functions_out N_end) {

            TARDIGRADE_ERROR_TOOLS_CHECK((size_type)(N_end - N_begin) == 20,
                                         "The dimension of the shape-function iterator is " +
                                             std::to_string((size_type)(N_end - N_begin)) + " but should be 20");

            // Set the first 8 nodes
            for (unsigned int i = 0; i < 8; ++i) {
                *(N_begin + i) = -0.125 * (1 + (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0))) *
                                 (1 + (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1))) *
                                 (1 + (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) *
                                 (2 - (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0)) -
                                  (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1)) -
                                  (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2)));
            }

            // Set the next 12 nodes
            for (unsigned int i = 0; i < 4; ++i) {
                *(N_begin + 8 + 2 * i + 0) =
                    0.250 * ((1 - (*(xi_begin + 0))) * (1 + (*(xi_begin + 0)))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 0) + 1)) * (*(xi_begin + 1))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 0) + 2)) * (*(xi_begin + 2)));

                *(N_begin + 8 + 2 * i + 1) =
                    0.250 * ((1 - (*(xi_begin + 1))) * (1 + (*(xi_begin + 1)))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 1) + 0)) * (*(xi_begin + 0))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 1) + 2)) * (*(xi_begin + 2)));

                *(N_begin + 16 + 1 * i + 0) =
                    0.250 * ((1 - (*(xi_begin + 2))) * (1 + (*(xi_begin + 2)))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (16 + 1 * i + 0) + 0)) * (*(xi_begin + 0))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (16 + 1 * i + 0) + 1)) * (*(xi_begin + 1)));
            }
        }

        /*!
         * Get the shape functions of the quadratic hexahedral element
         *
         * \param &xi_begin: The starting iterator of the local coordinates (must have dimension 3)
         * \param &xi_end: The stopping iterator of the local coordinates (must have dimension 3)
         * \param dNdxi_begin: The starting iterator of the local gradient of the shape functions (must have
         * dimension 24)
         * \param dNdxi_end: The stopping iterator of the local gradient of the shape functions (must have dimension
         * 24)
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        void QuadraticHex<T, node_in, local_point_in, shape_functions_out, grad_shape_functions_out>::
            GetLocalShapeFunctionGradients(const local_point_in &xi_begin, const local_point_in &xi_end,
                                           grad_shape_functions_out dNdxi_begin, grad_shape_functions_out dNdxi_end) {

            TARDIGRADE_ERROR_TOOLS_CHECK((size_type)(dNdxi_end - dNdxi_begin) == 60,
                                         "The dimension of the shape-function iterator is " +
                                             std::to_string((size_type)(dNdxi_end - dNdxi_begin)) +
                                             " but should be 60");

            // Set the first 8 nodes
            for (unsigned int i = 0; i < 8; ++i) {
                *(dNdxi_begin + 3 * i + 0) =
                    -0.125 * ((*(this->local_node_xi_begin + 3 * i + 0)) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1))) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) *
                                  (2 - (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0)) -
                                   (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1)) -
                                   (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) -
                              (1 + (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0))) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1))) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) *
                                  (*(this->local_node_xi_begin + 3 * i + 0)));

                *(dNdxi_begin + 3 * i + 1) =
                    -0.125 * ((1 + (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0))) *
                                  (*(this->local_node_xi_begin + 3 * i + 1)) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) *
                                  (2 - (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0)) -
                                   (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1)) -
                                   (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) -
                              (1 + (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0))) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1))) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) *
                                  (*(this->local_node_xi_begin + 3 * i + 1)));

                *(dNdxi_begin + 3 * i + 2) =
                    -0.125 * ((1 + (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0))) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1))) *
                                  (*(this->local_node_xi_begin + 3 * i + 2)) *
                                  (2 - (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0)) -
                                   (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1)) -
                                   (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) -
                              (1 + (*(this->local_node_xi_begin + 3 * i + 0)) * (*(xi_begin + 0))) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 1)) * (*(xi_begin + 1))) *
                                  (1 + (*(this->local_node_xi_begin + 3 * i + 2)) * (*(xi_begin + 2))) *
                                  (*(this->local_node_xi_begin + 3 * i + 2)));
            }

            // Set the next 12 nodes
            for (unsigned int i = 0; i < 4; ++i) {
                *(dNdxi_begin + 3 * (8 + 2 * i + 0) + 0) =
                    0.500 * (-(*(xi_begin + 0))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 0) + 1)) * (*(xi_begin + 1))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 0) + 2)) * (*(xi_begin + 2)));

                *(dNdxi_begin + 3 * (8 + 2 * i + 0) + 1) =
                    0.250 * ((1 - (*(xi_begin + 0))) * (1 + (*(xi_begin + 0)))) *
                    (*(this->local_node_xi_begin + 3 * (2 * i + 0) + 1)) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 0) + 2)) * (*(xi_begin + 2)));

                *(dNdxi_begin + 3 * (8 + 2 * i + 0) + 2) =
                    0.250 * ((1 - (*(xi_begin + 0))) * (1 + (*(xi_begin + 0)))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 0) + 1)) * (*(xi_begin + 1))) *
                    (*(this->local_node_xi_begin + 3 * (2 * i + 0) + 2));

                *(dNdxi_begin + 3 * (8 + 2 * i + 1) + 0) =
                    0.250 * ((1 - (*(xi_begin + 1))) * (1 + (*(xi_begin + 1)))) *
                    (*(this->local_node_xi_begin + 3 * (2 * i + 1) + 0)) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 1) + 2)) * (*(xi_begin + 2)));

                *(dNdxi_begin + 3 * (8 + 2 * i + 1) + 1) =
                    0.500 * (-(*(xi_begin + 1))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 1) + 0)) * (*(xi_begin + 0))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 1) + 2)) * (*(xi_begin + 2)));

                *(dNdxi_begin + 3 * (8 + 2 * i + 1) + 2) =
                    0.250 * ((1 - (*(xi_begin + 1))) * (1 + (*(xi_begin + 1)))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (2 * i + 1) + 0)) * (*(xi_begin + 0))) *
                    (*(this->local_node_xi_begin + 3 * (2 * i + 1) + 2));

                *(dNdxi_begin + 3 * (16 + 1 * i + 0) + 0) =
                    0.250 * ((1 - (*(xi_begin + 2))) * (1 + (*(xi_begin + 2)))) *
                    (*(this->local_node_xi_begin + 3 * (16 + 1 * i + 0) + 0)) *
                    (1 + (*(this->local_node_xi_begin + 3 * (16 + 1 * i + 0) + 1)) * (*(xi_begin + 1)));

                *(dNdxi_begin + 3 * (16 + 1 * i + 0) + 1) =
                    0.250 * ((1 - (*(xi_begin + 2))) * (1 + (*(xi_begin + 2)))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (16 + 1 * i + 0) + 0)) * (*(xi_begin + 0))) *
                    (*(this->local_node_xi_begin + 3 * (16 + 1 * i + 0) + 1));

                *(dNdxi_begin + 3 * (16 + 1 * i + 0) + 2) =
                    0.500 * (-(*(xi_begin + 2))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (16 + 1 * i + 0) + 0)) * (*(xi_begin + 0))) *
                    (1 + (*(this->local_node_xi_begin + 3 * (16 + 1 * i + 0) + 1)) * (*(xi_begin + 1)));
            }
        }

        /*!
         * Compute the global gradient of the quantity in the global coordinates
         *
         * \param &xi_begin: The starting iterator of the local point
         * \param &xi_end: The stopping iterator of the local point
         * \param &node_positions_begin: The starting iterator of the nodal positions (row major)
         * \param &node_positions_end: The stopping iterator of the nodal positions (row major)
         * \param &value_begin: The starting iterator of the shape function global gradient (row major)
         * \param &value_end: The stopping iterator of the shape function global gradient (row major)
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        void QuadraticHex<T, node_in, local_point_in, shape_functions_out, grad_shape_functions_out>::
            GetGlobalShapeFunctionGradients(const local_point_in &xi_begin, const local_point_in &xi_end,
                                            const node_in &node_positions_begin, const node_in &node_positions_end,
                                            grad_shape_functions_out value_begin, grad_shape_functions_out value_end) {

            TARDIGRADE_ERROR_TOOLS_CHECK((size_type)(value_end - value_begin) == 60,
                                         "The shape function global gradient has a size of " +
                                             std::to_string((unsigned int)(value_end - value_begin)) +
                                             " but must have a size of 60");

            std::array<typename std::iterator_traits<node_in>::value_type, 9> dxdxi;
            std::array<typename std::iterator_traits<node_in>::value_type, 9> dxidx;

            TARDIGRADE_ERROR_TOOLS_CATCH(this->GetLocalQuantityGradient(xi_begin, xi_end, node_positions_begin,
                                                                        node_positions_end, std::begin(dxdxi),
                                                                        std::end(dxdxi)));

            Eigen::Map<const Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            Eigen::Map<Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _dxidx(dxidx.data());

            _dxidx = (_dxdxi.inverse()).eval();

            std::fill(value_begin, value_end, 0);

            for (unsigned int node = 0; node < 20; ++node) {
                for (unsigned int inner = 0; inner < 3; ++inner) {
                    for (unsigned int outer = 0; outer < 3; ++outer) {
                        *(value_begin + 3 * node + outer) +=
                            this->_local_gradshapefunctions[3 * node + inner] * dxidx[3 * inner + outer];
                    }
                }
            }
        }

        /*!
         * Compute the value of the Jacobian of transformation from the local coordinates to the configuration for
         * volume integrals
         *
         * \param &xi_begin: The starting iterator of the local coordinates
         * \param &xi_end: The stopping iterator of the local coordinates
         * \param value: The Jacobian of transformation going from the local coordinates to the indicated
         * configuration
         * \param configuration: Compute the gradient w.r.t. the current configuration ( true ) or reference
         * configuration ( false )
         */
        template <typename T, class node_in, class local_point_in, class shape_functions_out,
                  class grad_shape_functions_out>
        void QuadraticHex<T, node_in, local_point_in, shape_functions_out, grad_shape_functions_out>::
            GetVolumeIntegralJacobianOfTransformation(const local_point_in &xi_begin, const local_point_in &xi_end,
                                                      typename std::iterator_traits<node_in>::value_type &value,
                                                      const bool configuration) {

            std::array<typename std::iterator_traits<node_in>::value_type, 9> dxdxi;

            if (configuration) {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->x_begin, this->x_end, std::begin(dxdxi),
                                               std::end(dxdxi));

            } else {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->X_begin, this->X_end, std::begin(dxdxi),
                                               std::end(dxdxi));
            }

            Eigen::Map<const Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            value = _dxdxi.determinant();
        }

    }  // namespace finiteElementUtilities

}  // namespace tardigradeBalanceEquations
