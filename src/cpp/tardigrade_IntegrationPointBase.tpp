/*!
=====================================================================
|               tardigrade_IntegrationPointBase.tpp                 |
=====================================================================
| The template definition file for an integration point             |
=====================================================================
*/

#include <algorithm>
#include <functional>

namespace tardigradeBalanceEquations {

    namespace finiteElement {

        /*!
         * Assemble the response of an integration point
         */
        template <class configuration>
        void IntegrationPointBase<configuration>::assembleIntegrationPointResponse() {
            std::fill(std::begin(_response), std::end(_response), typename configuration::response_type());

            computeConstantPointResponse();

            for (_i = 0; _i < configuration::num_test_functions; ++_i) {
                std::fill(std::begin(_response_i), std::end(_response_i), typename configuration::response_type());

                computeVariablePointResponse();

                accumulateIntegrationPointResponse();
            }

            std::transform(std::begin(_response), std::end(_response), std::begin(_response),
                           std::bind(std::multiplies<>(), std::placeholders::_1, _Jxw));
        }

        /*!
         * Assemble the Jacobian of an integration point
         */
        template <class configuration>
        void IntegrationPointBase<configuration>::assembleIntegrationPointJacobian() {
            std::fill(std::begin(_jacobian), std::end(_jacobian), typename configuration::jacobian_type());

            computeConstantPointJacobian();

            for (_i = 0; _i < configuration::num_test_functions; ++_i) {
                for (_j = 0; _j < configuration::num_interpolation_functions; ++_j) {
                    std::fill(std::begin(_jacobian_ij), std::end(_jacobian_ij),
                              typename configuration::jacobian_type());

                    computeVariablePointJacobian();

                    accumulateIntegrationPointJacobian();
                }
            }

            std::transform(std::begin(_jacobian), std::end(_jacobian), std::begin(_jacobian),
                           std::bind(std::multiplies<>(), std::placeholders::_1, _Jxw));
        }

        /*!
         * Accumulate the response of the current test function into the net response
         */
        template <class configuration>
        void IntegrationPointBase<configuration>::accumulateIntegrationPointResponse() {
            std::transform(std::begin(_response_i), std::end(_response_i),
                           std::begin(_response) + _i * configuration::num_response_terms,
                           std::begin(_response) + _i * configuration::num_response_terms, std::plus<>());
        }

        /*!
         * Accumulate the Jacobian of the current test function into the net Jacobian
         */
        template <class configuration>
        void IntegrationPointBase<configuration>::accumulateIntegrationPointJacobian() {
            constexpr unsigned int stride = configuration::num_interpolation_functions * configuration::num_dof;

            const unsigned int row_offset = configuration::num_response_terms * stride * _i;

            const unsigned int column_offset = configuration::num_dof * _j;

            for (unsigned int row = 0; row < configuration::num_response_terms; ++row) {
                std::transform(std::begin(_jacobian_ij) + configuration::num_dof * row,
                               std::begin(_jacobian_ij) + configuration::num_dof * (row + 1),
                               std::begin(_jacobian) + row_offset + stride * row + column_offset,
                               std::begin(_jacobian) + row_offset + stride * row + column_offset, std::plus<>());
            }
        }

    }  // namespace finiteElement

}  // namespace tardigradeBalanceEquations
