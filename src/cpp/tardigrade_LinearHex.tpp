/**
 ******************************************************************************
 * \file tardigrade_LinearHex.tpp
 ******************************************************************************
 * The template file for the Linear Hexahedral element
 ******************************************************************************
 */

namespace tardigradeBalanceEquations {

    namespace finiteElement {

        /*!
         * Constructor for the linear 8 node hex finite element class
         *
         * \param &_x_begin: The starting iterator for the current node positions
         * \param &_x_end: The stopping iterator for the current node positions
         * \param &_X_begin: The starting iterator for the reference node positions
         * \param &_X_end: The stopping iterator for the reference node positions
         */
        template <class element_configuration, typename T, class node_in, class local_point_in>
        LinearHex<element_configuration, T, node_in, local_point_in>::LinearHex(
            const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin, const node_in &_X_end)
            : FiniteElementBase<LinearHexConfiguration, node_in, local_point_in>(_x_begin, _x_end, _X_begin, _X_end,
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
        template <class element_configuration, typename T, class node_in, class local_point_in>
        void LinearHex<element_configuration, T, node_in, local_point_in>::GetShapeFunctions(
            const local_point_in &xi_begin, const local_point_in &xi_end, typename element_configuration::shape_functions_out N_begin,
            typename element_configuration::shape_functions_out N_end) {

            TARDIGRADE_ERROR_TOOLS_CHECK((size_type)(N_end - N_begin) == 8,
                                         "The dimension of the shape-function iterator is " +
                                             std::to_string((size_type)(N_end - N_begin)));

            for (std::pair<size_type, typename element_configuration::shape_functions_out> i(0, N_begin); i.second != N_end; ++i.first, ++i.second) {
                *i.second = (1 + (*(this->local_node_xi_begin + 3 * i.first + 0)) * (*(xi_begin + 0))) *
                            (1 + (*(this->local_node_xi_begin + 3 * i.first + 1)) * (*(xi_begin + 1))) *
                            (1 + (*(this->local_node_xi_begin + 3 * i.first + 2)) * (*(xi_begin + 2))) / 8;
            }
        }

        /*!
         * Get the shape function gradients of the linear hexahedral element
         *
         * \param &xi_begin: The starting iterator of the local coordinates (must have dimension 3)
         * \param &xi_end: The stopping iterator of the local coordinates (must have dimension 3)
         * \param dNdxi_begin: The starting iterator of the local gradient of the shape functions (must have
         * dimension 24)
         * \param dNdxi_end: The stopping iterator of the local gradient of the shape functions (must have dimension
         * 24)
         */
        template <class element_configuration, typename T, class node_in, class local_point_in>
        void LinearHex<element_configuration, T, node_in, local_point_in>::GetLocalShapeFunctionGradients(const local_point_in    &xi_begin,
                                                                                 const local_point_in    &xi_end,
                                                                                 typename element_configuration::grad_shape_functions_out dNdxi_begin,
                                                                                 typename element_configuration::grad_shape_functions_out dNdxi_end) {

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
        template <class element_configuration, typename T, class node_in, class local_point_in>
        void LinearHex<element_configuration, T, node_in, local_point_in>::GetGlobalShapeFunctionGradients(const local_point_in &xi_begin,
                                                                                  const local_point_in &xi_end,
                                                                                  const node_in &node_positions_begin,
                                                                                  const node_in &node_positions_end,
                                                                                  typename element_configuration::grad_shape_functions_out value_begin,
                                                                                  typename element_configuration::grad_shape_functions_out value_end) {

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
        template <class element_configuration, typename T, class node_in, class local_point_in>
        void LinearHex<element_configuration, T, node_in, local_point_in>::
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
        template <class element_configuration, typename T, class node_in, class local_point_in>
        void LinearHex<element_configuration, T, node_in, local_point_in>
            ::GetSurfaceIntegralJacobianOfTransformation(
                const unsigned int s,
                const local_point_in &xi_begin, const local_point_in &xi_end,
                typename std::iterator_traits<node_in>::value_type &value, const bool configuration) {

            using dxdxi_type = typename std::iterator_traits<node_in>::value_type;

            std::array<dxdxi_type, 9> dxdxi;
            std::array<dxdxi_type, 9> A, Ainv;

            if (configuration) {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->x_begin, this->x_end, std::begin(dxdxi),
                                               std::end(dxdxi));

            } else {
                this->GetLocalQuantityGradient(xi_begin, xi_end, this->X_begin, this->X_end, std::begin(dxdxi),
                                               std::end(dxdxi));
            }

            Eigen::Map<const Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _dxdxi(dxdxi.data());

            auto Jvol = _dxdxi.determinant();

            std::fill(std::begin(A), std::end(A), dxdxi_type());
            for ( unsigned int i = 0; i < 3; ++i){
                for ( unsigned int I = 0; I < 3; ++I ){
                    for ( unsigned int J = 0; J < 3; ++J ){
                        A[3*I+J] += dxdxi[dim*i+I] * dxdxi[dim*i+J];
                    }
                }
            }

            Eigen::Map<const Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _A(A.data());
            Eigen::Map<Eigen::Matrix<typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor> >
                _Ainv(Ainv.data());

            _Ainv = _A.inverse().eval();

            value = dxdxi_type();

            for (unsigned int I = 0; I < 3; ++I){
                for (unsigned int J = 0; J < 3; ++J){
                    value += surface_normals[3*s+I] * Ainv[3*I+J] * surface_normals[3*s+J];
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
        template <class element_configuration, typename T, class node_in, class local_point_in>
        void LinearHex<element_configuration, T, node_in, local_point_in>
            ::GetVolumeIntegrationPointData(const unsigned int i, typename element_configuration::local_point_out xi_begin,
                                                       typename element_configuration::local_point_out xi_end, typename element_configuration::volume_integration_point_weight_value_type &weight){

            TARDIGRADE_ERROR_TOOLS_CHECK(i<8, "The integration point id " + std::to_string(i) + " must be less than the number of integration points " + std::to_string(8));

            std::copy(std::begin(volume_integration_points) + 3 * i, std::begin(volume_integration_points) + 3 * (i + 1), xi_begin);

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
        template <class element_configuration, typename T, class node_in, class local_point_in>
        void LinearHex<element_configuration, T, node_in, local_point_in>::
             GetSurfaceIntegrationPointData(const unsigned int s, const unsigned int i, typename element_configuration::local_point_out xi_begin,
                                                    typename element_configuration::local_point_out xi_end, typename element_configuration::surface_integration_point_weight_value_type &weight){

            std::copy(std::begin(surface_integration_points) + 4 * 3 * s + 3 * i, std::begin(surface_integration_points) + 4 * 3 * s + 3 * (i+1), xi_begin);

            weight = surface_integration_weights[4*s+i];

        }

    }

}
