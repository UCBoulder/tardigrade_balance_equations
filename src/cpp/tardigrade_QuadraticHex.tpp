/**
 ******************************************************************************
 * \file tardigrade_QuadraticHex.tpp
 ******************************************************************************
 * The template file for the Quadratic Hexahedral element
 ******************************************************************************
 */

namespace tardigradeBalanceEquations {

    namespace finiteElement {

        /*!
         * Constructor for the quadratic 20 node hex element class
         *
         * \param &_x_begin: The starting iterator for the current node positions
         * \param &_x_end: The stopping iterator for the current node positions
         * \param &_X_begin: The starting iterator for the reference node positions
         * \param &_X_end: The stopping iterator for the reference node positions
         */
        template <class element_configuration>
        QuadraticHex<element_configuration>::QuadraticHex(const typename element_configuration::node_in &_x_begin,
                                                          const typename element_configuration::node_in &_x_end,
                                                          const typename element_configuration::node_in &_X_begin,
                                                          const typename element_configuration::node_in &_X_end)
            : FiniteElementBase<QuadraticHexConfiguration>(_x_begin, _x_end, _X_begin, _X_end, std::cbegin(local_nodes),
                                                           std::cend(local_nodes)) {}

        /*!
         * Get the shape functions of the quadratic 20 node hexahedral element
         *
         * \param &xi_begin: The starting iterator of the local coordinates (must have dimension 3)
         * \param &xi_end: The stopping iterator of the local coordinates (must have dimension 3)
         * \param N_begin: The starting iterator of the shape functions (must have dimension 8)
         * \param N_end: The stopping iterator of the shape functions (must have dimension 8)
         */
        template <class element_configuration>
        void QuadraticHex<element_configuration>::GetShapeFunctions(
            const typename element_configuration::local_point_in &xi_begin,
            const typename element_configuration::local_point_in &xi_end,
            typename element_configuration::shape_functions_out   N_begin,
            typename element_configuration::shape_functions_out   N_end) {
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
        template <class element_configuration>
        void QuadraticHex<element_configuration>::GetLocalShapeFunctionGradients(
            const typename element_configuration::local_point_in    &xi_begin,
            const typename element_configuration::local_point_in    &xi_end,
            typename element_configuration::grad_shape_functions_out dNdxi_begin,
            typename element_configuration::grad_shape_functions_out dNdxi_end) {
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
        template <class element_configuration>
        void QuadraticHex<element_configuration>::GetGlobalShapeFunctionGradients(
            const typename element_configuration::local_point_in    &xi_begin,
            const typename element_configuration::local_point_in    &xi_end,
            const typename element_configuration::node_in           &node_positions_begin,
            const typename element_configuration::node_in           &node_positions_end,
            typename element_configuration::grad_shape_functions_out value_begin,
            typename element_configuration::grad_shape_functions_out value_end) {
            TARDIGRADE_ERROR_TOOLS_CHECK((size_type)(value_end - value_begin) == 60,
                                         "The shape function global gradient has a size of " +
                                             std::to_string((unsigned int)(value_end - value_begin)) +
                                             " but must have a size of 60");

            std::array<typename std::iterator_traits<typename element_configuration::node_in>::value_type, 9> dxdxi;
            std::array<typename std::iterator_traits<typename element_configuration::node_in>::value_type, 9> dxidx;

            TARDIGRADE_ERROR_TOOLS_CATCH(this->GetLocalQuantityGradient(xi_begin, xi_end, node_positions_begin,
                                                                        node_positions_end, std::begin(dxdxi),
                                                                        std::end(dxdxi)));

            Eigen::Map<
                const Eigen::Matrix<typename std::iterator_traits<typename element_configuration::node_in>::value_type,
                                    3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            Eigen::Map<Eigen::Matrix<typename std::iterator_traits<typename element_configuration::node_in>::value_type,
                                     3, 3, Eigen::RowMajor> >
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
        template <class element_configuration>
        void QuadraticHex<element_configuration>::GetVolumeIntegralJacobianOfTransformation(
            const typename element_configuration::local_point_in                               &xi_begin,
            const typename element_configuration::local_point_in                               &xi_end,
            typename std::iterator_traits<typename element_configuration::node_in>::value_type &value,
            const bool                                                                          configuration) {
            std::array<typename std::iterator_traits<typename element_configuration::node_in>::value_type, 9> dxdxi;

            if (configuration) {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->x_begin, this->x_end, std::begin(dxdxi),
                                               std::end(dxdxi));

            } else {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->X_begin, this->X_end, std::begin(dxdxi),
                                               std::end(dxdxi));
            }

            Eigen::Map<
                const Eigen::Matrix<typename std::iterator_traits<typename element_configuration::node_in>::value_type,
                                    3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            value = _dxdxi.determinant();
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
        void QuadraticHex<element_configuration>::GetSurfaceIntegralJacobianOfTransformation(
            const unsigned int s, const typename element_configuration::local_point_in &xi_begin,
            const typename element_configuration::local_point_in                               &xi_end,
            typename std::iterator_traits<typename element_configuration::node_in>::value_type &value,
            const bool                                                                          configuration) {
            using dxdxi_type = typename std::iterator_traits<typename element_configuration::node_in>::value_type;

            std::array<dxdxi_type, 9> dxdxi;
            std::array<dxdxi_type, 9> A, Ainv;

            if (configuration) {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->x_begin, this->x_end, std::begin(dxdxi),
                                               std::end(dxdxi));

            } else {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->X_begin, this->X_end, std::begin(dxdxi),
                                               std::end(dxdxi));
            }

            Eigen::Map<
                const Eigen::Matrix<typename std::iterator_traits<typename element_configuration::node_in>::value_type,
                                    3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            auto Jvol = _dxdxi.determinant();

            std::fill(std::begin(A), std::end(A), dxdxi_type());
            for (unsigned int i = 0; i < 3; ++i) {
                for (unsigned int I = 0; I < 3; ++I) {
                    for (unsigned int J = 0; J < 3; ++J) {
                        A[3 * I + J] += dxdxi[dim * i + I] * dxdxi[dim * i + J];
                    }
                }
            }

            Eigen::Map<
                const Eigen::Matrix<typename std::iterator_traits<typename element_configuration::node_in>::value_type,
                                    3, 3, Eigen::RowMajor> >
                _A(A.data());
            Eigen::Map<Eigen::Matrix<typename std::iterator_traits<typename element_configuration::node_in>::value_type,
                                     3, 3, Eigen::RowMajor> >
                _Ainv(Ainv.data());

            _Ainv = _A.inverse().eval();

            value = dxdxi_type();

            for (unsigned int I = 0; I < 3; ++I) {
                for (unsigned int J = 0; J < 3; ++J) {
                    value += surface_normals[3 * s + I] * Ainv[3 * I + J] * surface_normals[3 * s + J];
                }
            }

            value = std::sqrt(value) * Jvol;
        }

        /*!
         * Get the integration point coordinates and weight for a volume integral
         *
         * \param i: The integration point number
         * \param xi_begin: The starting iterator of the integration point in local coordinates
         * \param xi_end: The stopping iterator of the integration point in local coordinates
         * \param &weight: The weight of the integration point
         */
        template <class element_configuration>
        void QuadraticHex<element_configuration>::GetVolumeIntegrationPointData(
            const unsigned int i, typename element_configuration::local_point_out xi_begin,
            typename element_configuration::local_point_out                             xi_end,
            typename element_configuration::volume_integration_point_weight_value_type &weight) {
            TARDIGRADE_ERROR_TOOLS_CHECK(i < 8, "The integration point id " + std::to_string(i) +
                                                    " must be less than the number of integration points " +
                                                    std::to_string(8));

            std::copy(std::begin(volume_integration_points) + 3 * i,
                      std::begin(volume_integration_points) + 3 * (i + 1), xi_begin);

            weight = volume_integration_weights[i];
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
        void QuadraticHex<element_configuration>::GetSurfaceIntegrationPointData(
            const unsigned int s, const unsigned int i, typename element_configuration::local_point_out xi_begin,
            typename element_configuration::local_point_out                              xi_end,
            typename element_configuration::surface_integration_point_weight_value_type &weight) {
            std::copy(std::begin(surface_integration_points) + 4 * 3 * s + 3 * i,
                      std::begin(surface_integration_points) + 4 * 3 * s + 3 * (i + 1), xi_begin);

            weight = surface_integration_weights[4 * s + i];
        }

        /*!
         * Compute the global normal of a local normal
         * All of the vectors must be continuous
         *
         * \param &xi_begin: The pointer to the start of the local coordinates of the normal
         * \param &local_normal_begin: The pointer to the start of the local normal vector
         * \param &global_normal_begin: The pointer to the start of the global normal vector
         * \param configuration: Whether to compute the normal in the reference (0) or current (1) configuration
         */
        template <class element_configuration>
        void QuadraticHex<element_configuration>::GetGlobalNormal(
            const typename element_configuration::local_node_value_type *xi_begin,
            const typename element_configuration::local_node_value_type *local_normal_begin,
            typename element_configuration::local_node_value_type *global_normal_begin,
            const bool configuration){

            using dxdxi_type = typename std::iterator_traits<typename element_configuration::node_in>::value_type;

            std::array<dxdxi_type, 9> dxdxi, dxidx;

            if (configuration) {
                this->GetLocalQuantityGradient(xi_begin, xi_begin + element_configuration::local_dim, this->x_begin, this->x_end, std::begin(dxdxi),
                                               std::end(dxdxi));

            } else {
                this->GetLocalQuantityGradient(xi_begin, xi_begin + element_configuration::local_dim, this->X_begin, this->X_end, std::begin(dxdxi),
                                               std::end(dxdxi));
            }

            Eigen::Map<
                const Eigen::Matrix<dxdxi_type,
                                    3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            Eigen::Map<
                Eigen::Matrix<dxdxi_type,
                                    3, 3, Eigen::RowMajor> >
                _dxidx(dxidx.data());

            auto Jvol = _dxdxi.determinant();
            _dxidx = _dxdxi.inverse().eval();

            std::fill(global_normal_begin, global_normal_begin + element_configuration::dim, dxdxi_type());

            for ( unsigned int i = 0; i < element_configuration::dim; ++i ){

                for ( unsigned int I = 0; I < element_configuration::local_dim; ++I ){

                    *(global_normal_begin + i) += Jvol * dxidx[element_configuration::dim * I + i] * (*(local_normal_begin + I));

                }

            }

            auto norm = std::sqrt(std::inner_product(global_normal_begin, global_normal_begin + element_configuration::dim, global_normal_begin, dxdxi_type()));

            std::transform(global_normal_begin, global_normal_begin + element_configuration::dim, global_normal_begin, std::bind(std::divides<>(), std::placeholders::_1, norm));

        }

    }  // namespace finiteElement

}  // namespace tardigradeBalanceEquations
