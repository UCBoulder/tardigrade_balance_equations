/*!
=====================================================================
|                 tardigrade_IntegrationPointBase.h                 |
=====================================================================
| The header definition file for an integration point               |
=====================================================================
*/

#ifndef TARDIGRADE_INTEGRATIONPOINTBASE_H
#define TARDIGRADE_INTEGRATIONPOINTBASE_H

#include <array>

namespace tardigradeBalanceEquations {

    namespace finiteElement {

        namespace unit_test {

            class IntegrationPointBaseTester;

        }

        /*!
         * The configuration of the integration point
         */
        template <int _dim, int _num_response_terms, int _num_dof, int _num_test_functions,
                  int _num_interpolation_functions, typename _test_type = double, typename _interpolation_type = double,
                  typename _response_type = double, typename _jacobian_type = double, typename _dof_type = double,
                  typename _Jxw_type = double>
        class IntegrationPointConfiguration {
           public:
            //! The type of the test function
            using test_type = _test_type;

            //! The type of the interpolation function
            using interpolation_type = _interpolation_type;

            //! The type of the response
            using response_type = _response_type;

            //! The type of the Jacobian
            using jacobian_type = _jacobian_type;

            //! The type of the dof
            using dof_type = _dof_type;

            //! The type of the product of the Jacobian and the Integration point weight
            using Jxw_type = _Jxw_type;

            //! The spatial dimension
            static constexpr unsigned int dim = _dim;

            //! The number of response terms
            static constexpr unsigned int num_response_terms = _num_response_terms;

            //! The number of degrees of freedom
            static constexpr unsigned int num_dof = _num_dof;

            //! The number of test functions which are defined at the integration point
            static constexpr unsigned int num_test_functions = _num_test_functions;

            //! The number of interpolation functions which are defined at the integration point
            static constexpr unsigned int num_interpolation_functions = _num_interpolation_functions;
        };

        /*!
         * The base class for integration points
         */
        template <class configuration>
        class IntegrationPointBase {
           public:
            /*!
             * Default constructor
             */
            IntegrationPointBase()
                : _test({}),
                  _grad_test({}),
                  _interpolation({}),
                  _grad_interpolation({}),
                  _dof({}),
                  _grad_dof({}),
                  _previous_dof({}),
                  _previous_grad_dof({}),
                  _Jxw(typename configuration::Jxw_type()) {};

            /*!
             * Full constructor
             *
             * \param &test: The test function values
             * \param &grad_test: The spatial gradient of the test function values
             * \param &interpolation: The interpolation function values
             * \param &grad_interpolation: The spatial gradient of the interpolation function values
             * \param &dof: The dof values at the integration point
             * \param &grad_dof: The spatial gradient of the dof values at the integration point
             * \param &previous_dof: The previous dof values at the integration point
             * \param &previous_grad_dof: The previous spatial gradient of the dof values at the integration point
             * \param &Jxw: The product of the Jacobian of transformation and the integration point weight
             */
            IntegrationPointBase(
                const std::array<typename configuration::test_type, configuration::num_test_functions> &test,
                const std::array<typename configuration::test_type,
                                 configuration::num_test_functions * configuration::dim>               &grad_test,
                const std::array<typename configuration::test_type, configuration::num_interpolation_functions>
                                                                                                  &interpolation,
                const std::array<typename configuration::test_type,
                                 configuration::num_interpolation_functions * configuration::dim> &grad_interpolation,
                const std::array<typename configuration::test_type, configuration::num_dof>       &dof,
                const std::array<typename configuration::test_type, configuration::num_dof * configuration::dim>
                                                                                            &grad_dof,
                const std::array<typename configuration::test_type, configuration::num_dof> &previous_dof,
                const std::array<typename configuration::test_type, configuration::num_dof * configuration::dim>
                                                       &previous_grad_dof,
                const typename configuration::Jxw_type &Jxw)
                : _test(test),
                  _grad_test(grad_test),
                  _interpolation(interpolation),
                  _grad_interpolation(grad_interpolation),
                  _dof(dof),
                  _grad_dof(grad_dof),
                  _previous_dof(previous_dof),
                  _previous_grad_dof(previous_grad_dof),
                  _Jxw(Jxw) {};

            //! Get the net response
            const auto getResponse() { return &_response; }

            //! Get the net Jacobian
            const auto getJacobian() { return &_jacobian; }

           protected:
            /*!
             * A function that is called prior to looping through the test functions.
             */
            virtual void computeConstantPointResponse() {};

            /*!
             * A function that is called for every test function value
             *
             * Should define _response_i
             */
            virtual void computeVariablePointResponse() {};

            /*!
             * A function that is called prior to looping through the test and interpolation functions
             */
            virtual void computeConstantPointJacobian() {};

            /*!
             * A function that is called for every combindation of test and interpolation functions
             *
             * Should define _jacobian_ij
             */
            virtual void computeVariablePointJacobian() {};

            void assembleIntegrationPointResponse();

            void assembleIntegrationPointJacobian();

            void accumulateIntegrationPointResponse();

            void accumulateIntegrationPointJacobian();

            //! The current test function index
            unsigned int _i = 0;

            //! The current interpolation function index
            unsigned int _j = 0;

            //! The test functions
            const std::array<typename configuration::test_type, configuration::num_test_functions> _test;

            //! The gradient of the test functions
            const std::array<typename configuration::test_type, configuration::num_test_functions * configuration::dim>
                _grad_test;

            //! The interpolation functions
            const std::array<typename configuration::interpolation_type, configuration::num_interpolation_functions>
                _interpolation;

            //! The gradient of the interpolation functions
            const std::array<typename configuration::interpolation_type,
                             configuration::num_interpolation_functions * configuration::dim>
                _grad_interpolation;

            //! The DOF vector
            const std::array<typename configuration::dof_type, configuration::num_dof> _dof;

            //! The DOF gradient vector
            const std::array<typename configuration::dof_type, configuration::num_dof * configuration::dim> _grad_dof;

            //! The previous DOF vector
            const std::array<typename configuration::dof_type, configuration::num_dof> _previous_dof;

            //! The previous DOF gradient vector
            const std::array<typename configuration::dof_type, configuration::num_dof * configuration::dim>
                _previous_grad_dof;

            //! The product of the volumetric Jacobian and the integration point weight
            const typename configuration::Jxw_type _Jxw;

            //! The current response vector
            std::array<typename configuration::response_type, configuration::num_response_terms> _response_i;

            //! The current Jacobian row-major matrix
            std::array<typename configuration::response_type,
                       configuration::num_response_terms * configuration::num_dof>
                _jacobian_ij;

           private:
            friend class unit_test::IntegrationPointBaseTester;  //!< Friend class which allows modification of private
                                                                 //!< variables. ONLY TO BE USED FOR TESTING!

            //! The net response vector
            std::array<typename configuration::response_type,
                       configuration::num_response_terms * configuration::num_test_functions>
                _response;

            //! The net Jacobian row-major matrix
            std::array<typename configuration::jacobian_type,
                       configuration::num_response_terms * configuration::num_test_functions * configuration::num_dof *
                           configuration::num_interpolation_functions>
                _jacobian;
        };

    }  // namespace finiteElement

}  // namespace tardigradeBalanceEquations

#include "tardigrade_IntegrationPointBase.tpp"

#endif
