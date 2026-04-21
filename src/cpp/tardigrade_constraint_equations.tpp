/**
 ******************************************************************************
 * \file tardigrade_constraint_equations.tpp
 ******************************************************************************
 * The template file for the constraint equations associated with the balance
 * equations
 ******************************************************************************
 */

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>

#include "tardigrade_constraint_equations.h"

namespace tardigradeBalanceEquations {

    namespace constraintEquations {

        template <int predicted_internal_energy_index, typename internal_energy_type, class material_response_iter,
                  typename test_function_type, typename result_type>
        void computeInternalEnergyConstraint(const internal_energy_type   &internal_energy,
                                             const material_response_iter &material_response_begin,
                                             const material_response_iter &material_response_end,
                                             const test_function_type &test_function, result_type &result) {
            /*!
             * Compute the value of the constraint on the internal energy.
             *
             * Traditionally, the internal energy is related to the temperature via
             *
             * \f$ e = c \theta \f$
             *
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             *
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             *
             * \param &internal_energy: The internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &result: The resulting error between the internal energy DOF and the material response
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                predicted_internal_energy_index < (unsigned int)(material_response_end - material_response_begin),
                "The index for the predicted internal energy is out of range for the material response")

            result = (*(material_response_begin + predicted_internal_energy_index) - internal_energy) * test_function;
        }

        template <int predicted_internal_energy_index, typename internal_energy_type, typename density_type,
                  class material_response_iter, typename test_function_type, typename result_type>
        void computeInternalEnergyConstraint(const internal_energy_type &internal_energy, const density_type &density,
                                             const material_response_iter &material_response_begin,
                                             const material_response_iter &material_response_end,
                                             const test_function_type &test_function, result_type &result) {
            /*!
             * Compute the value of the constraint on the internal energy.
             *
             * Traditionally, the internal energy per unit mass is related to the temperature via
             *
             * \f$ e = c \theta \f$
             *
             * where \f$ e \f$ is the internal energy per unit mass, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             *
             * \f$ R = \rho e^{p} - \tilde{e} \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation and \f$ \tilde{e} \f$ is the internal energy
             * per unit volume we can force the connection and improve the numeric properties of the
             * solution by using the internal energy per unit volume as the degree of freedom.
             *
             * \param &internal_energy: The internal energy per unit volume degree of freedom
             * \param &density: The density of the material
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &result: The resulting error between the internal energy DOF and the material response
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                predicted_internal_energy_index < (unsigned int)(material_response_end - material_response_begin),
                "The index for the predicted internal energy is out of range for the material response")

            result = (density * (*(material_response_begin + predicted_internal_energy_index)) - internal_energy) *
                     test_function;
        }

        template <class configuration, int predicted_internal_energy_index,
                  typename internal_energy_type, class material_response_iter, class material_response_jacobian_iter,
                  typename test_function_type, typename interpolation_function_type,
                  class interpolation_function_gradient_iter, class full_material_response_dof_gradient_iter,
                  typename dUDotdU_type, typename result_type, class dRdRho_iter, class dRdU_iter, class dRdW_iter,
                  class dRdTheta_iter, class dRdE_iter, class dRdVF_iter, class dRdZ_iter, class dRdUMesh_iter>
        void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy, const material_response_iter &material_response_begin,
            const material_response_iter          &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function, const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter     &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter     &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU, const unsigned int &phase, result_type &result, dRdRho_iter dRdRho_begin,
            dRdRho_iter dRdRho_end, dRdU_iter dRdU_begin, dRdU_iter dRdU_end, dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end, dRdE_iter dRdE_begin, dRdE_iter dRdE_end,
            dRdVF_iter dRdVF_begin, dRdVF_iter dRdVF_end, dRdZ_iter dRdZ_begin, dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the value of the constraint on the internal energy.
             *
             * Traditionally, the internal energy is related to the temperature via
             *
             * \f$ e = c \theta \f$
             *
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             *
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             *
             * \param &internal_energy: The internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial
             * gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial
             * gradient of the interpolation function
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the
             * material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the
             * material response dof vector
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &phase: The current active phase
             * \param &result: The resulting error between the internal energy DOF and the material response
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase spatial DOF
             * associated with velocity
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase spatial DOF
             * associated with velocity
             * \param &dRdW_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement
             * dof
             * \param &dRdW_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement
             * dof
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the phase
             * temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the phase
             * temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the phase internal
             * energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the phase internal
             * energy
             * \param &dRdVF_begin: The starting iterator of the derivative of the residual w.r.t. the phase volume
             * fraction
             * \param &dRdVF_end: The stopping iterator of the derivative of the residual w.r.t. the phase volume
             * fraction
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the phase additional
             * dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh
             * displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh
             * displacement
             */

            const unsigned int     nphases            = (unsigned int)(dRdRho_end - dRdRho_begin);

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * configuration::material::dimension == (unsigned int)(dRdU_end - dRdU_begin),
                "dRdU must be a consistent size with the material response dimension and the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * configuration::material::dimension == (unsigned int)(dRdW_end - dRdW_begin),
                "dRdW must be a consistent size with the material response dimension and the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(dRdTheta_end - dRdTheta_begin),
                                         "dRdTheta must be the same size as the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(dRdE_end - dRdE_begin),
                                         "dRdE must be the same size as the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                configuration::material::dimension == (unsigned int)(dRdUMesh_end - dRdUMesh_begin),
                "dRdUMesh must be the same size as the spatial dimension of the material response")

            std::fill(dRdRho_begin, dRdRho_end, 0);
            std::fill(dRdU_begin, dRdU_end, 0);
            std::fill(dRdW_begin, dRdW_end, 0);
            std::fill(dRdTheta_begin, dRdTheta_end, 0);
            std::fill(dRdE_begin, dRdE_end, 0);
            std::fill(dRdVF_begin, dRdVF_end, 0);
            std::fill(dRdZ_begin, dRdZ_end, 0);
            std::fill(dRdUMesh_begin, dRdUMesh_end, 0);

            computeInternalEnergyConstraint<predicted_internal_energy_index>(internal_energy, material_response_begin,
                                                                             material_response_end, test_function,
                                                                             result);

            *(dRdE_begin + phase) -= test_function * interpolation_function;

            // PREDICTED INTERNAL ENERGY
            // density
            for (auto p = std::pair<unsigned int, dRdRho_iter>(0, dRdRho_begin); p.second != dRdRho_end;
                 ++p.first, ++p.second) {
                *p.second += test_function *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::density_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::density_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // velocity
            for (auto p = std::pair<unsigned int, dRdU_iter>(0, dRdU_begin); p.second != dRdU_end;
                 ++p.first, ++p.second) {
                *p.second += test_function *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::velocity_index + p.first)) *
                             interpolation_function * dUDotdU;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::velocity_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a)) * dUDotdU;
                }
            }

            // displacement
            for (auto p = std::pair<unsigned int, dRdW_iter>(0, dRdW_begin); p.second != dRdW_end;
                 ++p.first, ++p.second) {
                *p.second += test_function *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::displacement_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::displacement_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // temperature
            for (auto p = std::pair<unsigned int, dRdTheta_iter>(0, dRdTheta_begin); p.second != dRdTheta_end;
                 ++p.first, ++p.second) {
                *p.second += test_function *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::temperature_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::temperature_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // internal energy
            for (auto p = std::pair<unsigned int, dRdE_iter>(0, dRdE_begin); p.second != dRdE_end;
                 ++p.first, ++p.second) {
                *p.second += test_function *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::internal_energy_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::internal_energy_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // volume fraction
            for (auto p = std::pair<unsigned int, dRdVF_iter>(0, dRdVF_begin); p.second != dRdVF_end;
                 ++p.first, ++p.second) {
                *p.second += test_function *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::volume_fraction_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::volume_fraction_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // additional dof
            for (auto p = std::pair<unsigned int, dRdZ_iter>(0, dRdZ_begin); p.second != dRdZ_end;
                 ++p.first, ++p.second) {
                *p.second += test_function *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::additional_dof_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::additional_dof_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // mesh displacement
            for (unsigned int I = 0; I < (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof); ++I) {
                for (unsigned int k = 0; k < configuration::material::dimension; ++k) {
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *(dRdUMesh_begin + a) -=
                            test_function *
                            (*(material_response_jacobian_begin +
                               (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                   (predicted_internal_energy_index) +
                               (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) + configuration::material::dimension * I + k)) *
                            (*(full_material_response_dof_gradient_begin + configuration::material::dimension * I + a)) *
                            (*(interpolation_function_gradient_begin + k));
                    }
                }
            }

            for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                *(dRdUMesh_begin + a) += result * (*(interpolation_function_gradient_begin + a));
            }
        }

        template <class configuration, int predicted_internal_energy_index,
                  typename internal_energy_type, typename density_type, class material_response_iter,
                  class material_response_jacobian_iter, typename test_function_type,
                  typename interpolation_function_type, class interpolation_function_gradient_iter,
                  class full_material_response_dof_gradient_iter, typename dUDotdU_type, typename result_type,
                  class dRdRho_iter, class dRdU_iter, class dRdW_iter, class dRdTheta_iter, class dRdE_iter,
                  class dRdVF_iter, class dRdZ_iter, class dRdUMesh_iter>
        void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy, const density_type &density,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function, const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter     &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter     &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU, const unsigned int &phase, result_type &result, dRdRho_iter dRdRho_begin,
            dRdRho_iter dRdRho_end, dRdU_iter dRdU_begin, dRdU_iter dRdU_end, dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end, dRdE_iter dRdE_begin, dRdE_iter dRdE_end,
            dRdVF_iter dRdVF_begin, dRdVF_iter dRdVF_end, dRdZ_iter dRdZ_begin, dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the value of the constraint on the internal energy.
             *
             * Traditionally, the internal energy per unit mass is related to the temperature via
             *
             * \f$ e = c \theta \f$
             *
             * where \f$ e \f$ is the internal energy per unit mass, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             *
             * \f$ R = \rho e^{p} - \tilde{e} \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation and \f$ \tilde{e} \f$ is the internal energy
             * per unit volume we can force the connection and improve the numeric properties of the
             * solution by using the internal energy per unit volume as the degree of freedom.
             *
             * \param &internal_energy: The internal energy degree of freedom
             * \param &density: The material density
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial
             * gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial
             * gradient of the interpolation function
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the
             * material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the
             * material response dof vector
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &phase: The current active phase
             * \param &result: The resulting error between the internal energy DOF and the material response
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase spatial DOF
             * associated with velocity
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase spatial DOF
             * associated with velocity
             * \param &dRdW_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement
             * dof
             * \param &dRdW_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement
             * dof
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the phase
             * temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the phase
             * temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the phase internal
             * energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the phase internal
             * energy
             * \param &dRdVF_begin: The starting iterator of the derivative of the residual w.r.t. the phase volume
             * fraction
             * \param &dRdVF_end: The stopping iterator of the derivative of the residual w.r.t. the phase volume
             * fraction
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the phase additional
             * dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh
             * displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh
             * displacement
             */

            const unsigned int     nphases            = (unsigned int)(dRdRho_end - dRdRho_begin);

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * configuration::material::dimension == (unsigned int)(dRdU_end - dRdU_begin),
                "dRdU must be a consistent size with the material response dimension and the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * configuration::material::dimension == (unsigned int)(dRdW_end - dRdW_begin),
                "dRdW must be a consistent size with the material response dimension and the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(dRdTheta_end - dRdTheta_begin),
                                         "dRdTheta must be the same size as the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(dRdE_end - dRdE_begin),
                                         "dRdE must be the same size as the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                configuration::material::dimension == (unsigned int)(dRdUMesh_end - dRdUMesh_begin),
                "dRdUMesh must be the same size as the spatial dimension of the material response")

            std::fill(dRdRho_begin, dRdRho_end, 0);
            std::fill(dRdU_begin, dRdU_end, 0);
            std::fill(dRdW_begin, dRdW_end, 0);
            std::fill(dRdTheta_begin, dRdTheta_end, 0);
            std::fill(dRdE_begin, dRdE_end, 0);
            std::fill(dRdVF_begin, dRdVF_end, 0);
            std::fill(dRdZ_begin, dRdZ_end, 0);
            std::fill(dRdUMesh_begin, dRdUMesh_end, 0);

            computeInternalEnergyConstraint<predicted_internal_energy_index>(internal_energy, density,
                                                                             material_response_begin,
                                                                             material_response_end, test_function,
                                                                             result);

            *(dRdRho_begin + phase) +=
                test_function * (*(material_response_begin + predicted_internal_energy_index)) * interpolation_function;
            *(dRdE_begin + phase) -= test_function * interpolation_function;

            // PREDICTED INTERNAL ENERGY
            // density
            for (auto p = std::pair<unsigned int, dRdRho_iter>(0, dRdRho_begin); p.second != dRdRho_end;
                 ++p.first, ++p.second) {
                *p.second += test_function * density *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::density_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function * density *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::density_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // velocity
            for (auto p = std::pair<unsigned int, dRdU_iter>(0, dRdU_begin); p.second != dRdU_end;
                 ++p.first, ++p.second) {
                *p.second += test_function * density *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::velocity_index + p.first)) *
                             interpolation_function * dUDotdU;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function * density *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::velocity_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a)) * dUDotdU;
                }
            }

            // displacement
            for (auto p = std::pair<unsigned int, dRdW_iter>(0, dRdW_begin); p.second != dRdW_end;
                 ++p.first, ++p.second) {
                *p.second += test_function * density *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::displacement_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function * density *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::displacement_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // temperature
            for (auto p = std::pair<unsigned int, dRdTheta_iter>(0, dRdTheta_begin); p.second != dRdTheta_end;
                 ++p.first, ++p.second) {
                *p.second += test_function * density *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::temperature_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function * density *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::temperature_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // internal energy
            for (auto p = std::pair<unsigned int, dRdE_iter>(0, dRdE_begin); p.second != dRdE_end;
                 ++p.first, ++p.second) {
                *p.second += test_function * density *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::internal_energy_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function * density *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::internal_energy_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // internal energy
            for (auto p = std::pair<unsigned int, dRdVF_iter>(0, dRdVF_begin); p.second != dRdVF_end;
                 ++p.first, ++p.second) {
                *p.second += test_function * density *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::volume_fraction_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function * density *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::volume_fraction_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // additional dof
            for (auto p = std::pair<unsigned int, dRdZ_iter>(0, dRdZ_begin); p.second != dRdZ_end;
                 ++p.first, ++p.second) {
                *p.second += test_function * density *
                             (*(material_response_jacobian_begin +
                                (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                    (predicted_internal_energy_index) +
                                nphases * configuration::material::dof::additional_dof_index + p.first)) *
                             interpolation_function;

                for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                    *p.second += test_function * density *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (predicted_internal_energy_index) +
                                    (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) +
                                    configuration::material::dimension * (nphases * configuration::material::dof::additional_dof_index + p.first) + a)) *
                                 (*(interpolation_function_gradient_begin + a));
                }
            }

            // mesh displacement
            for (unsigned int I = 0; I < (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof); ++I) {
                for (unsigned int k = 0; k < configuration::material::dimension; ++k) {
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *(dRdUMesh_begin + a) -=
                            test_function * density *
                            (*(material_response_jacobian_begin +
                               (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) * (1 + configuration::material::dimension) *
                                   (predicted_internal_energy_index) +
                               (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) + configuration::material::dimension * I + k)) *
                            (*(full_material_response_dof_gradient_begin + configuration::material::dimension * I + a)) *
                            (*(interpolation_function_gradient_begin + k));
                    }
                }
            }

            for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                *(dRdUMesh_begin + a) += result * (*(interpolation_function_gradient_begin + a));
            }
        }

        template <int predicted_internal_energy_index, class internal_energy_iter, class material_response_iter,
                  typename test_function_type, class result_iter>
        void computeInternalEnergyConstraint(const internal_energy_iter   &internal_energy_begin,
                                             const internal_energy_iter   &internal_energy_end,
                                             const material_response_iter &material_response_begin,
                                             const material_response_iter &material_response_end,
                                             const test_function_type &test_function, result_iter result_begin,
                                             result_iter result_end) {
            /*!
             * Compute the value of the constraint on the internal energy for a multiphase problem.
             *
             * Traditionally, the internal energy is related to the temperature via
             *
             * \f$ e = c \theta \f$
             *
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             *
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             *
             * \param &internal_energy_begin: The starting iterator of the internal energy degree of freedom
             * \param &internal_energy_end: The stopping iterator of the internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &result_begin: The starting iterator of the resulting error between the internal energy DOF and
             * the material response
             * \param &result_end: The starting iterator of the resulting error between the internal energy DOF and the
             * material response
             */

            const unsigned int nphases = (unsigned int)(result_end - result_begin);
            const unsigned int material_response_size =
                (unsigned int)(material_response_end - material_response_begin) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(internal_energy_end - internal_energy_begin),
                                         "The number of internal energy values must be equal to the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == (unsigned int)(material_response_end - material_response_begin),
                "The material response vector must be a scalar multiple of the number of phases")

            for (auto v = std::pair<unsigned int, result_iter>(0, result_begin); v.second != result_end;
                 ++v.first, ++v.second) {
                computeInternalEnergyConstraint<predicted_internal_energy_index>(
                    *(internal_energy_begin + v.first), material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * (v.first + 1), test_function, *v.second);
            }
        }

        template <int predicted_internal_energy_index, class internal_energy_iter, class density_iter,
                  class material_response_iter, typename test_function_type, class result_iter>
        void computeInternalEnergyConstraint(const internal_energy_iter &internal_energy_begin,
                                             const internal_energy_iter &internal_energy_end,
                                             const density_iter &density_begin, const density_iter &density_end,
                                             const material_response_iter &material_response_begin,
                                             const material_response_iter &material_response_end,
                                             const test_function_type &test_function, result_iter result_begin,
                                             result_iter result_end) {
            /*!
             * Compute the value of the constraint on the internal energy.
             *
             * Traditionally, the internal energy per unit mass is related to the temperature via
             *
             * \f$ e = c \theta \f$
             *
             * where \f$ e \f$ is the internal energy per unit mass, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             *
             * \f$ R = \rho e^{p} - \tilde{e} \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation and \f$ \tilde{e} \f$ is the internal energy
             * per unit volume we can force the connection and improve the numeric properties of the
             * solution by using the internal energy per unit volume as the degree of freedom.
             *
             * \param &internal_energy_begin: The starting iterator of the internal energy degree of freedom
             * \param &internal_energy_end: The stopping iterator of the internal energy degree of freedom
             * \param &density_begin: The starting iterator of the density degree of freedom
             * \param &density_end: The stopping iterator of the density degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &result_begin: The starting iterator of the resulting error between the internal energy DOF and
             * the material response
             * \param &result_end: The starting iterator of the resulting error between the internal energy DOF and the
             * material response
             */

            const unsigned int nphases = (unsigned int)(result_end - result_begin);
            const unsigned int material_response_size =
                (unsigned int)(material_response_end - material_response_begin) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(internal_energy_end - internal_energy_begin),
                                         "The number of internal energy values must be equal to the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == (unsigned int)(material_response_end - material_response_begin),
                "The material response vector must be a scalar multiple of the number of phases")

            for (auto v = std::pair<unsigned int, result_iter>(0, result_begin); v.second != result_end;
                 ++v.first, ++v.second) {
                computeInternalEnergyConstraint<predicted_internal_energy_index>(
                    *(internal_energy_begin + v.first), *(density_begin + v.first),
                    material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * (v.first + 1), test_function, *v.second);
            }
        }

        template <class configuration, int predicted_internal_energy_index,
                  class internal_energy_iter, class material_response_iter, class material_response_jacobian_iter,
                  typename test_function_type, typename interpolation_function_type,
                  class interpolation_function_gradient_iter, class full_material_response_dof_gradient_iter,
                  typename dUDotdU_type, class result_iter, class dRdRho_iter, class dRdU_iter, class dRdW_iter,
                  class dRdTheta_iter, class dRdE_iter, class dRdVF_iter, class dRdZ_iter, class dRdUMesh_iter>
        void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin, const internal_energy_iter &internal_energy_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function, const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter     &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter     &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU, result_iter result_begin, result_iter result_end, dRdRho_iter dRdRho_begin,
            dRdRho_iter dRdRho_end, dRdU_iter dRdU_begin, dRdU_iter dRdU_end, dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end, dRdE_iter dRdE_begin, dRdE_iter dRdE_end,
            dRdVF_iter dRdVF_begin, dRdVF_iter dRdVF_end, dRdZ_iter dRdZ_begin, dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the value of the constraint on the internal energy.
             *
             * Traditionally, the internal energy is related to the temperature via
             *
             * \f$ e = c \theta \f$
             *
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             *
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             *
             * \param &internal_energy_begin: The starting iterator of the internal energy degree of freedom
             * \param &internal_energy_end: The stopping iterator internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial
             * gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial
             * gradient of the interpolation function
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the
             * material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the
             * material response dof vector
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &result_begin: The starting iterator of the resulting error between the internal energy DOF and
             * the material response
             * \param &result_end: The stopping starting iterator of the resulting error between the internal energy DOF
             * and the material response
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase spatial DOF
             * associated with velocity
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase spatial DOF
             * associated with velocity
             * \param &dRdW_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement
             * dof
             * \param &dRdW_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement
             * dof
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the phase
             * temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the phase
             * temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the phase internal
             * energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the phase internal
             * energy
             * \param &dRdVF_begin: The starting iterator of the derivative of the residual w.r.t. the phase volume
             * fraction
             * \param &dRdVF_end: The stopping iterator of the derivative of the residual w.r.t. the phase volume
             * fraction
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the phase additional
             * dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh
             * displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh
             * displacement
             */

            const unsigned int nphases = (unsigned int)(result_end - result_begin);
            const unsigned int material_response_size =
                (unsigned int)(material_response_end - material_response_begin) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(internal_energy_end - internal_energy_begin),
                                         "The number of internal energy values must be equal to the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == (unsigned int)(material_response_end - material_response_begin),
                "The material response vector must be a scalar multiple of the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size * (1 + configuration::material::dimension) *
                        (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) ==
                    (unsigned int)(material_response_jacobian_end - material_response_jacobian_begin),
                "The material response jacobian must have a consistent size with the material response vector and the "
                "configuration::material::num_dof\n  number of phases           : " +
                    std::to_string(nphases) +
                    "\n  material_response_size     : " + std::to_string(material_response_size) +
                    "\n  dof per phase              : " + std::to_string(configuration::material::num_phase_dof) +
                    "\n  additional dof             : " + std::to_string(configuration::material::num_additional_dof) +
                    "\n  material response dimension: " + std::to_string(configuration::material::dimension) +
                    "\n  expected jacobian size     : " +
                    std::to_string(nphases * material_response_size * (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) *
                                   (1 + configuration::material::dimension)) +
                    "\n  actual jacobian size       : " +
                    std::to_string((unsigned int)(material_response_jacobian_end - material_response_jacobian_begin)))

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::material::dimension * (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) ==
                                             (unsigned int)(full_material_response_dof_gradient_end -
                                                            full_material_response_dof_gradient_begin),
                                         "The full material response dof gradient have a size of the material response "
                                         "dimension times the number of dof in the material response")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * configuration::material::dimension == (unsigned int)(dRdU_end - dRdU_begin),
                "dRdU must be a consistent size with the material response dimension and the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * configuration::material::dimension == (unsigned int)(dRdW_end - dRdW_begin),
                "dRdW must be a consistent size with the material response dimension and the number of phases squared")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * nphases == (unsigned int)(dRdTheta_end - dRdTheta_begin),
                                         "dRdTheta must be the same size as the number of phases squared")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * nphases == (unsigned int)(dRdE_end - dRdE_begin),
                                         "dRdE must be the same size as the number of phases squared")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * configuration::material::dimension ==
                                             (unsigned int)(dRdUMesh_end - dRdUMesh_begin),
                                         "dRdUMesh must be the same size as the spatial dimension of the material "
                                         "response times the number of phases")

            for (auto v = std::pair<unsigned int, result_iter>(0, result_begin); v.second != result_end;
                 ++v.first, ++v.second) {
                computeInternalEnergyConstraint<
                    configuration, predicted_internal_energy_index>(
                    *(internal_energy_begin + v.first), material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * (v.first + 1),
                    material_response_jacobian_begin + material_response_size *
                                                           (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) *
                                                           (1 + configuration::material::dimension) * v.first,
                    material_response_jacobian_begin + material_response_size *
                                                           (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) *
                                                           (1 + configuration::material::dimension) * (v.first + 1),
                    test_function, interpolation_function, interpolation_function_gradient_begin,
                    interpolation_function_gradient_end, full_material_response_dof_gradient_begin,
                    full_material_response_dof_gradient_end, dUDotdU, v.first, *v.second,
                    dRdRho_begin + nphases * v.first, dRdRho_begin + nphases * (v.first + 1),
                    dRdU_begin + nphases * configuration::material::dimension * v.first,
                    dRdU_begin + nphases * configuration::material::dimension * (v.first + 1),
                    dRdW_begin + nphases * configuration::material::dimension * v.first,
                    dRdW_begin + nphases * configuration::material::dimension * (v.first + 1), dRdTheta_begin + nphases * v.first,
                    dRdTheta_begin + nphases * (v.first + 1), dRdE_begin + nphases * v.first,
                    dRdE_begin + nphases * (v.first + 1), dRdVF_begin + nphases * v.first,
                    dRdVF_begin + nphases * (v.first + 1), dRdZ_begin + configuration::material::num_additional_dof * v.first,
                    dRdZ_begin + configuration::material::num_additional_dof * (v.first + 1), dRdUMesh_begin + configuration::material::dimension * v.first,
                    dRdUMesh_begin + configuration::material::dimension * (v.first + 1));
            }
        }

        template <class configuration, int predicted_internal_energy_index,
                  class internal_energy_iter, class density_iter, class material_response_iter,
                  class material_response_jacobian_iter, typename test_function_type,
                  typename interpolation_function_type, class interpolation_function_gradient_iter,
                  class full_material_response_dof_gradient_iter, typename dUDotdU_type, class result_iter,
                  class dRdRho_iter, class dRdU_iter, class dRdW_iter, class dRdTheta_iter, class dRdE_iter,
                  class dRdVF_iter, class dRdZ_iter, class dRdUMesh_iter>
        void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin, const internal_energy_iter &internal_energy_end,
            const density_iter &density_begin, const density_iter &density_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function, const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter     &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter     &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU, result_iter result_begin, result_iter result_end, dRdRho_iter dRdRho_begin,
            dRdRho_iter dRdRho_end, dRdU_iter dRdU_begin, dRdU_iter dRdU_end, dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end, dRdE_iter dRdE_begin, dRdE_iter dRdE_end,
            dRdVF_iter dRdVF_begin, dRdVF_iter dRdVF_end, dRdZ_iter dRdZ_begin, dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the value of the constraint on the internal energy.
             *
             * Traditionally, the internal energy is related to the temperature via
             *
             * \f$ e = c \theta \f$
             *
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             *
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             *
             * \param &internal_energy_begin: The starting iterator of the internal energy degree of freedom
             * \param &internal_energy_end: The stopping iterator of the internal energy degree of freedom
             * \param &density_begin: The starting iterator of the material density
             * \param &density_end: The stopping iterator of the material density
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial
             * gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial
             * gradient of the interpolation function
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the
             * material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the
             * material response dof vector
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &result_begin: The starting iterator of the resulting error between the internal energy DOF and
             * the material response
             * \param &result_end: The stopping starting iterator of the resulting error between the internal energy DOF
             * and the material response
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase spatial DOF
             * associated with velocity
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase spatial DOF
             * associated with velocity
             * \param &dRdW_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement
             * dof
             * \param &dRdW_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement
             * dof
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the phase
             * temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the phase
             * temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the phase internal
             * energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the phase internal
             * energy
             * \param &dRdVF_begin: The starting iterator of the derivative of the residual w.r.t. the phase volume
             * fraction
             * \param &dRdVF_end: The stopping iterator of the derivative of the residual w.r.t. the phase volume
             * fraction
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the phase additional
             * dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh
             * displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh
             * displacement
             */

            const unsigned int nphases = (unsigned int)(result_end - result_begin);
            const unsigned int material_response_size =
                (unsigned int)(material_response_end - material_response_begin) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(internal_energy_end - internal_energy_begin),
                                         "The number of internal energy values must be equal to the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == (unsigned int)(material_response_end - material_response_begin),
                "The material response vector must be a scalar multiple of the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * material_response_size * (1 + configuration::material::dimension) *
                                                 (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) ==
                                             (unsigned int)(material_response_jacobian_end -
                                                            material_response_jacobian_begin),
                                         "The material response jacobian vector must be a scalar multiple of the "
                                         "number of phases, the material response size and 1 + the material response "
                                         "dimension times the number of dof in the material response")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::material::dimension * (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) ==
                                             (unsigned int)(full_material_response_dof_gradient_end -
                                                            full_material_response_dof_gradient_begin),
                                         "The full material response dof gradient have a size of the material response "
                                         "dimension times the number of dof in the material resionse")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * configuration::material::dimension == (unsigned int)(dRdU_end - dRdU_begin),
                "dRdU must be a consistent size with the material response dimension and the number of phases")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * configuration::material::dimension == (unsigned int)(dRdW_end - dRdW_begin),
                "dRdW must be a consistent size with the material response dimension and the number of phases squared")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * nphases == (unsigned int)(dRdTheta_end - dRdTheta_begin),
                                         "dRdTheta must be the same size as the number of phases squared")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * nphases == (unsigned int)(dRdE_end - dRdE_begin),
                                         "dRdE must be the same size as the number of phases squared")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * configuration::material::dimension ==
                                             (unsigned int)(dRdUMesh_end - dRdUMesh_begin),
                                         "dRdUMesh must be the same size as the spatial dimension of the material "
                                         "response times the number of phases")

            for (auto v = std::pair<unsigned int, result_iter>(0, result_begin); v.second != result_end;
                 ++v.first, ++v.second) {
                computeInternalEnergyConstraint<
                    configuration, predicted_internal_energy_index>(
                    *(internal_energy_begin + v.first), *(density_begin + v.first),
                                          material_response_begin + material_response_size * v.first,
                                          material_response_begin + material_response_size * (v.first + 1),
                                          material_response_jacobian_begin +
                                              material_response_size * (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) *
                                                  (1 + configuration::material::dimension) * v.first,
                                          material_response_jacobian_begin +
                                              material_response_size * (nphases * configuration::material::num_phase_dof + configuration::material::num_additional_dof) *
                                                  (1 + configuration::material::dimension) * (v.first + 1),
                                          test_function, interpolation_function, interpolation_function_gradient_begin,
                                          interpolation_function_gradient_end,
                                          full_material_response_dof_gradient_begin,
                                          full_material_response_dof_gradient_end, dUDotdU, v.first, *v.second,
                                          dRdRho_begin + nphases * v.first, dRdRho_begin + nphases * (v.first + 1),
                                          dRdU_begin + nphases * configuration::material::dimension * v.first,
                                          dRdU_begin + nphases * configuration::material::dimension * (v.first + 1),
                                          dRdW_begin + nphases * configuration::material::dimension * v.first,
                                          dRdW_begin + nphases * configuration::material::dimension * (v.first + 1),
                                          dRdTheta_begin + nphases * v.first, dRdTheta_begin + nphases * (v.first + 1),
                                          dRdE_begin + nphases * v.first, dRdE_begin + nphases * (v.first + 1),
                                          dRdVF_begin + nphases * v.first, dRdVF_begin + nphases * (v.first + 1),
                                          dRdZ_begin + configuration::material::num_additional_dof * v.first,
                                          dRdZ_begin + configuration::material::num_additional_dof * (v.first + 1),
                                          dRdUMesh_begin + configuration::material::dimension * v.first,
                                          dRdUMesh_begin + configuration::material::dimension * (v.first + 1));
            }
        }

        template <class displacement_dot_iter, class velocity_iter, typename test_function_type, class result_iter>
        void computeDisplacementConstraint(const displacement_dot_iter &displacement_dot_begin,
                                           const displacement_dot_iter &displacement_dot_end,
                                           const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
                                           const test_function_type &test_function, result_iter result_begin,
                                           result_iter result_end) {
            /*!
             * Compute the constraint for the displacement to the velocity DOF. If this is integrated over the volume
             * then the test function should be provided as normal. If this is directly applied at the nodes, then
             * setting the test function to 1 will achieve the desired result.
             *
             * \param displacement_dot_begin: The starting iterator for the time derivative of the displacement dof
             * \param displacement_dot_end: The stopping iterator for the time derivative of the displacement dof
             * \param velocity_begin: The starting iterator for the velocity dof
             * \param velocity_end: The stopping iterator for the velocity dof
             * \param test_function: The test function. If the constraint is being applied in a volume-integrated
             *     way then this should be set to the test function value. If it is being applied at a mesh's nodes
             *     then it should be set to 1.
             * \param result_begin: The starting iterator of the constraint violation result
             * \param result_end: The stopping iterator of the constraint violation result
             */

            using result_type = typename std::iterator_traits<result_iter>::value_type;

            TARDIGRADE_ERROR_TOOLS_EVAL(const unsigned int length =
                                            (unsigned int)(displacement_dot_end - displacement_dot_begin);)

            TARDIGRADE_ERROR_TOOLS_CHECK(length == (unsigned int)(velocity_end - velocity_begin),
                                         "The velocity and density dot vectors must be the same size")

            TARDIGRADE_ERROR_TOOLS_CHECK(length == (unsigned int)(result_end - result_begin),
                                         "The result and density dot vectors must be the same size")

            std::transform(displacement_dot_begin, displacement_dot_end, velocity_begin, result_begin,
                           std::minus<result_type>());

            std::transform(result_begin, result_end, result_begin,
                           std::bind(std::multiplies<result_type>(), std::placeholders::_1, test_function));
        }

        template <class configuration, class displacement_dot_iter, class velocity_iter, typename test_function_type,
                  typename interpolation_function_type, class interpolation_function_gradient_iter,
                  typename dDDotdD_type, class result_iter, class dRdD_iter, class dRdV_iter, class dRdUMesh_iter>
        void computeDisplacementConstraint(
            const displacement_dot_iter &displacement_dot_begin, const displacement_dot_iter &displacement_dot_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const test_function_type &test_function, const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const dDDotdD_type &dDDotdD, result_iter result_begin, result_iter result_end, dRdD_iter dRdD_begin,
            dRdD_iter dRdD_end, dRdV_iter dRdV_begin, dRdV_iter dRdV_end, dRdUMesh_iter dRdUMesh_begin,
            dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the constraint for the displacement to the velocity DOF. If this is integrated over the volume
             * then the test function, interpolation function and interpolation function gradient should be provided
             * as normal. If this is directly applied at the nodes, then setting the test function and interpolation
             * function to 1 and the interpolation function gradient to zero will achieve the desired result.
             *
             * \param displacement_dot_begin: The starting iterator for the time derivative of the displacement dof
             * \param displacement_dot_end: The stopping iterator for the time derivative of the displacement dof
             * \param velocity_begin: The starting iterator for the velocity dof
             * \param velocity_end: The stopping iterator for the velocity dof
             * \param test_function: The test function. If the constraint is being applied in a volume-integrated
             *     way then this should be set to the test function value. If it is being applied at a mesh's nodes
             *     then it should be set to 1.
             * \param interpolation_function: The interpolation function. If the constraint is being applied in a
             * volume-integrated way then this should be set to the interpolation function value. If it is being applied
             * at a mesh's nodes then it should be set to 1.
             * \param interpolation_function_gradient_begin: The starting iterator of the interpolation function. If the
             * constraint is being applied in a volume-integrated way then this should be set to the interpolation
             * function gradient value. If it is being applied at a mesh's nodes then it should be set to 0.
             * \param interpolation_function_gradient_end: The starting iterator of the interpolation function. If the
             * constraint is being applied in a volume-integrated way then this should be set to the interpolation
             * function gradient value. If it is being applied at a mesh's nodes then it should be set to 0.
             * \param dDDotdD: The derivative of the displacement dot w.r.t. the displacement
             * \param result_begin: The starting iterator of the constraint violation result
             * \param result_end: The stopping iterator of the constraint violation result
             * \param dRdD_begin: The starting iterator for the derivative of the residual w.r.t. the displacement (only
             * the diagonal is returned as everything else is zero)
             * \param dRdD_end:   The stopping iterator for the derivative of the residual w.r.t. the displacement (only
             * the diagonal is returned as everything else is zero)
             * \param dRdV_begin: The starting iterator for the derivative of the residual w.r.t. the velocity (only the
             * diagonal is returned as everything else is zero)
             * \param dRdV_end:   The stopping iterator for the derivative of the residual w.r.t. the velocity (only the
             * diagonal is returned as everything else is zero)
             * \param dRdUMesh_begin: The starting iterator for the derivative of the residual w.r.t. the mesh
             * displacement
             * \param dRdUMesh_end: The stopping iterator for the derivative of the residual w.r.t. the mesh
             * displacement
             */

            TARDIGRADE_ERROR_TOOLS_EVAL(const unsigned int length =
                                            (unsigned int)(displacement_dot_end - displacement_dot_begin);)

            TARDIGRADE_ERROR_TOOLS_CHECK(length == (unsigned int)(dRdD_end - dRdD_begin),
                                         "The dRdD and density dot vectors must be the same size")

            TARDIGRADE_ERROR_TOOLS_CHECK(length == (unsigned int)(dRdV_end - dRdV_begin),
                                         "The dRdV and density dot vectors must be the same size")

            TARDIGRADE_ERROR_TOOLS_CHECK(length * configuration::dimension == (unsigned int)(dRdUMesh_end - dRdUMesh_begin),
                                         "The dRdUMesh vector must be the length of the density dot vector time the "
                                         "length of the interpolation function gradient vector")

            computeDisplacementConstraint(displacement_dot_begin, displacement_dot_end, velocity_begin, velocity_end,
                                          test_function, result_begin, result_end);

            std::fill(dRdD_begin, dRdD_end, 0);
            std::fill(dRdV_begin, dRdV_end, 0);
            std::fill(dRdUMesh_begin, dRdUMesh_end, 0);

            for (auto v = std::pair<unsigned int, dRdD_iter>(0, dRdD_begin); v.second != dRdD_end;
                 ++v.first, ++v.second) {
                *v.second               = test_function * dDDotdD * interpolation_function;
                *(dRdV_begin + v.first) = -test_function * interpolation_function;

                for (unsigned int a = 0; a < configuration::dimension; ++a) {
                    *(dRdUMesh_begin + configuration::dimension * v.first + a) =
                        (*(result_begin + v.first)) * (*(interpolation_function_gradient_begin + a));
                }
            }
        }

        template <class configuration, int cauchy_stress_index, int predicted_internal_energy_index,
                  int body_force_index, int interphasic_force_index, int heat_flux_index,
                  int internal_heat_generation_index, int interphasic_heat_transfer_index,
                  int trace_mass_change_velocity_gradient_index, class density_iter, class volume_fraction_iter,
                  class material_response_iter, class material_response_jacobian_iter, class mixture_response_iter,
                  class mixture_jacobian_iter>
        inline void computeMixtureMaterialResponse(
            const density_iter &density_begin, const density_iter &density_end,
            const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const mixture_response_iter &mixture_response_begin, const mixture_response_iter &mixture_response_end,
            const mixture_jacobian_iter &mixture_jacobian_begin, const mixture_jacobian_iter &mixture_jacobian_end) {
            /*!
             * Assemble the mixture material response from the behavior of each phase
             *
             * We assume that the mixture model's DOF are included in the full Jacobian as the final "phase" but
             * the number of phases passed in (i.e., the number of densities) is the number of each of the true
             * phases. If there are \f$ n \f$ true phases, then there are \f$n+1\f$ phases defined in the DOF
             * vector. Similarly, there should be \f$n\f$ volume fractions, material responses, and material response
             * jacobians passed in.
             *
             * Some of the material responses are assembled using a volume fraction approach (e.g., the Cauchy stress),
             * others are assembled with a density-weighted approach (e.g., the body force), and others are summed
             * directly (e.g., the interphasic heat transfer). The difference is based on the balance equations.
             *
             * Direct Sum:
             * mass change, interphasic force, heat flux, interphasic heat transfer
             *
             * Volume Fraction Weighted:
             * Cauchy stress, trace mass change velocity gradient
             *
             * Density weighted
             * internal energy, body force, internal heat generation
             */

            using density_type = typename std::iterator_traits<density_iter>::value_type;

            auto material_response_size = (mixture_response_end - mixture_response_begin);

            auto num_dof = (mixture_jacobian_end - mixture_jacobian_begin) / material_response_size;

            auto num_phases = (material_response_end - material_response_begin) / material_response_size;

            TARDIGRADE_ERROR_TOOLS_CHECK(material_response_size * num_dof ==
                                             (unsigned int)(mixture_jacobian_end - mixture_jacobian_begin),
                                         "The mixture jacobian size must be an integer multiple of the material "
                                         "response size (i.e., the number of degrees of freedom)")

            TARDIGRADE_ERROR_TOOLS_CHECK(material_response_size * num_phases ==
                                             (unsigned int)(material_response_end - material_response_begin),
                                         "The material response size must be an integer multiple of the mixture "
                                         "response size (i.e., the number of true phases) )")

            TARDIGRADE_ERROR_TOOLS_CHECK(material_response_size * num_phases * num_dof ==
                                             (unsigned int)(material_response_jacobian_end -
                                                            material_response_jacobian_begin),
                                         "The material response jacobian size must be an integer multiple of the "
                                         "mixture response jacobian size (i.e., the number of true phases)")

            std::fill(mixture_response_begin, mixture_response_end, 0);

            std::fill(mixture_jacobian_begin, mixture_jacobian_end, 0);

            density_type density_sum = 0;
            for (auto v = density_begin; v != density_end; ++v) {
                density_sum += *v;
            };

            for (unsigned int phase = 0; phase < num_phases; ++phase) {
                //
                // Incorporate the direct summations
                //

                // mass change rate
                for (unsigned int i = 0; i < 1; ++i) {
                    *(mixture_response_begin + configuration::material::mass_change_index + i) +=
                        *(material_response_begin + phase * material_response_size + configuration::material::mass_change_index + i);

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + configuration::material::mass_change_index) + j) +=
                            *(material_response_jacobian_begin + material_response_size * num_dof * phase +
                              num_dof * (i + configuration::material::mass_change_index) + j);
                    }
                }

                // interphasic force
                for (unsigned int i = 0; i < configuration::dimension; ++i) {
                    *(mixture_response_begin + interphasic_force_index + i) +=
                        *(material_response_begin + phase * material_response_size + interphasic_force_index + i);

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + interphasic_force_index) + j) +=
                            *(material_response_jacobian_begin + material_response_size * num_dof * phase +
                              num_dof * (i + interphasic_force_index) + j);
                    }
                }

                // heat flux
                for (unsigned int i = 0; i < configuration::dimension; ++i) {
                    *(mixture_response_begin + heat_flux_index + i) +=
                        *(material_response_begin + phase * material_response_size + heat_flux_index + i);

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + heat_flux_index) + j) +=
                            *(material_response_jacobian_begin + material_response_size * num_dof * phase +
                              num_dof * (i + heat_flux_index) + j);
                    }
                }

                // interphasic heat transfer
                for (unsigned int i = 0; i < 1; ++i) {
                    *(mixture_response_begin + interphasic_heat_transfer_index + i) +=
                        *(material_response_begin + phase * material_response_size + interphasic_heat_transfer_index +
                          i);

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + interphasic_heat_transfer_index) + j) +=
                            *(material_response_jacobian_begin + material_response_size * num_dof * phase +
                              num_dof * (i + interphasic_heat_transfer_index) + j);
                    }
                }

                //
                // Incorporate volume-fraction summations
                //

                // Cauchy stress
                for (unsigned int i = 0; i < configuration::dimension * configuration::dimension; ++i) {
                    *(mixture_response_begin + cauchy_stress_index + i) +=
                        (*(volume_fraction_begin + phase)) *
                        (*(material_response_begin + phase * material_response_size + cauchy_stress_index + i));

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + cauchy_stress_index) + j) +=
                            (*(volume_fraction_begin + phase)) *
                            (*(material_response_jacobian_begin + material_response_size * num_dof * phase +
                               num_dof * (i + cauchy_stress_index) + j));
                    }

                    // Add contributions due to dependence on the volume fraction
                    *(mixture_jacobian_begin + num_dof * (i + cauchy_stress_index) +
                      (num_phases + 1) * configuration::material::dof::volume_fraction_index + phase) +=
                        (*(material_response_begin + material_response_size * phase + cauchy_stress_index + i));
                }

                // trace mass change velocity gradient
                for (unsigned int i = 0; i < 1; ++i) {
                    *(mixture_response_begin + trace_mass_change_velocity_gradient_index + i) +=
                        (*(volume_fraction_begin + phase)) *
                        (*(material_response_begin + phase * material_response_size +
                           trace_mass_change_velocity_gradient_index + i));

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + trace_mass_change_velocity_gradient_index) + j) +=
                            (*(volume_fraction_begin + phase)) *
                            (*(material_response_jacobian_begin + material_response_size * num_dof * phase +
                               num_dof * (i + trace_mass_change_velocity_gradient_index) + j));
                    }

                    // Add contributions due to dependence on the volume fraction
                    *(mixture_jacobian_begin + num_dof * (i + trace_mass_change_velocity_gradient_index) +
                      (num_phases + 1) * configuration::material::dof::volume_fraction_index + phase) +=
                        (*(material_response_begin + material_response_size * phase +
                           trace_mass_change_velocity_gradient_index + i));
                }

                //
                // Incorporate density-weighted summations
                //

                // body force
                for (unsigned int i = 0; i < configuration::dimension; ++i) {
                    *(mixture_response_begin + body_force_index + i) +=
                        (*(density_begin + phase)) *
                        (*(material_response_begin + phase * material_response_size + body_force_index + i)) /
                        density_sum;

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + body_force_index) + j) +=
                            (*(density_begin + phase)) *
                            (*(material_response_jacobian_begin + material_response_size * num_dof * phase +
                               num_dof * (i + body_force_index) + j)) /
                            density_sum;
                    }

                    // Add contributions due to dependence on the density
                    *(mixture_jacobian_begin + num_dof * (i + body_force_index) + (num_phases + 1) * configuration::material::dof::density_index +
                      phase) += (*(material_response_begin + material_response_size * phase + body_force_index + i)) /
                                density_sum;

                    for (unsigned int j = 0; j < num_phases; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + body_force_index) + (num_phases + 1) * configuration::material::dof::density_index +
                          j) += -(*(density_begin + phase)) *
                                (*(material_response_begin + material_response_size * phase + body_force_index + i)) /
                                (density_sum * density_sum);
                    }
                }

                // internal energy
                for (unsigned int i = 0; i < 1; ++i) {
                    *(mixture_response_begin + predicted_internal_energy_index + i) +=
                        (*(density_begin + phase)) *
                        (*(material_response_begin + phase * material_response_size + predicted_internal_energy_index + i)) /
                        density_sum;

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + predicted_internal_energy_index) + j) +=
                            (*(density_begin + phase)) *
                            (*(material_response_jacobian_begin + material_response_size * num_dof * phase +
                               num_dof * (i + predicted_internal_energy_index) + j)) /
                            density_sum;
                    }

                    // Add contributions due to dependence on the density
                    *(mixture_jacobian_begin + num_dof * (i + predicted_internal_energy_index) +
                      (num_phases + 1) * configuration::material::dof::density_index + phase) +=
                        (*(material_response_begin + material_response_size * phase + predicted_internal_energy_index + i)) /
                        density_sum;

                    for (unsigned int j = 0; j < num_phases; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + predicted_internal_energy_index) +
                          (num_phases + 1) * configuration::material::dof::density_index + j) +=
                            -(*(density_begin + phase)) *
                            (*(material_response_begin + material_response_size * phase + predicted_internal_energy_index + i)) /
                            (density_sum * density_sum);
                    }
                }

                // internal heat generation
                for (unsigned int i = 0; i < 1; ++i) {
                    *(mixture_response_begin + internal_heat_generation_index + i) +=
                        (*(density_begin + phase)) *
                        (*(material_response_begin + phase * material_response_size + internal_heat_generation_index +
                           i)) /
                        density_sum;

                    for (unsigned int j = 0; j < num_dof; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + internal_heat_generation_index) + j) +=
                            (*(density_begin + phase)) *
                            (*(material_response_jacobian_begin + phase * material_response_size * num_dof +
                               num_dof * (i + internal_heat_generation_index) + j)) /
                            density_sum;
                    }

                    // Add contributions due to dependence on the density
                    *(mixture_jacobian_begin + num_dof * (i + internal_heat_generation_index) +
                      (num_phases + 1) * configuration::material::dof::density_index + phase) +=
                        (*(material_response_begin + material_response_size * phase + internal_heat_generation_index +
                           i)) /
                        density_sum;

                    for (unsigned int j = 0; j < num_phases; ++j) {
                        *(mixture_jacobian_begin + num_dof * (i + internal_heat_generation_index) +
                          (num_phases + 1) * configuration::material::dof::density_index + j) +=
                            -(*(density_begin + phase)) *
                            (*(material_response_begin + material_response_size * phase +
                               internal_heat_generation_index + i)) /
                            (density_sum * density_sum);
                    }
                }
            }
        }

    }  // namespace constraintEquations

}  // namespace tardigradeBalanceEquations
