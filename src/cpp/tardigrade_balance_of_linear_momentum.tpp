/**
 ******************************************************************************
 * \file tardigrade_balance_of_linear_momentum.tpp
 ******************************************************************************
 * The template file for the equations associated with the balance of linear
 * momentum
 ******************************************************************************
 */

#include <algorithm>
#include <functional>
#include <numeric>

#include "tardigrade_balance_of_linear_momentum.h"

namespace tardigradeBalanceEquations {

    namespace balanceOfLinearMomentum {

        template <class configuration, typename density_type, typename density_dot_type, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter, class body_force_iter,
                  class result_iter>
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin, const body_force_iter &body_force_end, result_iter result_begin,
            result_iter result_end) {
            /*!
             * Compute the non-divergence part of the balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho
             * v_i \right) = -f_i \f$
             *
             * and the non-divergence part is
             *
             * \f$ -\left( \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) +
             * f_i = 0 \f$
             *
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. The expression
             * we implement is therefore
             *
             * \f$ -\rho_{,j} v_j v_i - \rho v_{i,j} v_j - \rho v_i v_{j,j} + \rho b_i - \frac{\partial \rho}{\partial
             * t} v_i - \rho \frac{\partial v_i}{\partial t} \f$
             *
             * \param &density: The mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot: The partial time derivative of the density \f$ \left( \frac{\partial}{\partial t}
             * \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right)
             * \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            // Compute the inner product of the velocity and the density gradient
            using density_gradient_type = typename std::iterator_traits<density_gradient_iter>::value_type;
            using velocity_type         = typename std::iterator_traits<velocity_iter>::value_type;
            using grad_rho_dot_type = decltype(std::declval<density_gradient_type &>() * std::declval<velocity_type>());
            using velocity_gradient_type = typename std::iterator_traits<velocity_gradient_iter>::value_type;

            grad_rho_dot_type grad_rho_dot_v =
                std::inner_product(density_gradient_begin, density_gradient_end, velocity_begin, 0.);

            // Compute the trace of the velocity gradient
            velocity_gradient_type trace_velocity_gradient = 0;

            for (unsigned int i = 0; i < configuration::dimension; ++i) {
                trace_velocity_gradient += *(velocity_gradient_begin + configuration::dimension * i + i);
            }

            // Compute the balance of linear momentum
            for (unsigned int i = 0; i < configuration::dimension; ++i) {
                *(result_begin + i) = density * ((*(body_force_begin + i)) - (*(velocity_dot_begin + i)) -
                                                 trace_velocity_gradient * (*(velocity_begin + i))) -
                                      grad_rho_dot_v * (*(velocity_begin + i)) - density_dot * (*(velocity_begin + i));

                for (unsigned int j = 0; j < configuration::dimension; ++j) {
                    *(result_begin + i) -= density * (*(velocity_gradient_begin + configuration::dimension * i + j)) *
                                           (*(velocity_begin + j));
                }
            }
        }

        template <class configuration, typename density_type, typename density_dot_type, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter, class body_force_iter,
                  class cauchy_stress_iter, typename volume_fraction_type, typename testFunction_type,
                  class testFunctionGradient_iter, class result_iter>
        void computeBalanceOfLinearMomentum(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin, const body_force_iter &body_force_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction, const testFunction_type &test_function,
            const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter &test_function_gradient_end, result_iter result_begin,
            result_iter result_end) {
            /*!
             * Compute the balance of linear momentum
             *
             * \param &density: The mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot: The partial time derivative of the density \f$ \left( \frac{\partial}{\partial t}
             * \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right)
             * \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the cauchy stress \f$ \left( \sigma_{ij} \right)
             * \f$
             * \param &cauchy_stress_end: The stopping iterator of the cauchy stress \f$ \left( \sigma_{ij} \right) \f$
             * \param &volume_fraction: The volume fraction of the phase. Only applied to the Cauchy stress because the
             * density is assumed to be the apparent density i.e., the mass of the phase per unit volume.
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$
             * \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            using result_type = typename std::iterator_traits<result_iter>::value_type;
            using final_type  = decltype(std::declval<result_type &>() * std::declval<testFunction_type &>());

            std::array<result_type, configuration::dimension> non_divergence_result;

            std::array<result_type, configuration::dimension> divergence_result;

            // Compute the non-divergence part of the balance of linear momentum
            computeBalanceOfLinearMomentumNonDivergence<configuration>(
                density, density_dot, density_gradient_begin, density_gradient_end, velocity_begin, velocity_end,
                velocity_dot_begin, velocity_dot_end, velocity_gradient_begin, velocity_gradient_end, body_force_begin,
                body_force_end, std::begin(non_divergence_result), std::end(non_divergence_result));

            std::transform(std::cbegin(non_divergence_result), std::cend(non_divergence_result),
                           std::begin(non_divergence_result),
                           std::bind(std::multiplies<final_type>(), std::placeholders::_1, test_function));

            // Compute the divergence part of the balance of linear momentum
            computeBalanceOfLinearMomentumDivergence<configuration>(test_function_gradient_begin,
                                                                    test_function_gradient_end, cauchy_stress_begin,
                                                                    cauchy_stress_end, volume_fraction,
                                                                    std::begin(divergence_result),
                                                                    std::end(divergence_result));

            // Assemble the full balance equation
            std::transform(std::cbegin(non_divergence_result), std::cend(non_divergence_result),
                           std::cbegin(divergence_result), result_begin, std::plus<final_type>());
        }

        template <class configuration, int body_force_index, int cauchy_stress_index,
                  int interphasic_force_index, typename density_type, typename density_dot_type,
                  class density_gradient_iter, class velocity_iter, class velocity_dot_iter,
                  class velocity_gradient_iter, class material_response_iter, typename volume_fraction_type,
                  typename testFunction_type, class testFunctionGradient_iter, class result_iter>
        void computeBalanceOfLinearMomentum(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const volume_fraction_type &volume_fraction, const testFunction_type &test_function,
            const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter &test_function_gradient_end, result_iter result_begin,
            result_iter result_end) {
            /*!
             * Compute the balance of linear momentum including the inter-phasic force
             *
             * body_force_index: The index of the material response vector where the body force is located
             * cauchy_stress_index: The index of the material response vector where the cauchy stress force is located
             * interphasic_force_index: The index of the material response vector where the net interphasic force is
             * located
             *
             * \param &density: The mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot: The partial time derivative of the density \f$ \left( \frac{\partial}{\partial t}
             * \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &volume_fraction: The volume fraction of the phase. Only applied to the Cauchy stress because the
             * density is assumed to be the apparent density i.e., the mass of the phase per unit volume.
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$
             * \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            computeBalanceOfLinearMomentum<configuration>(
                density, density_dot, density_gradient_begin, density_gradient_end, velocity_begin, velocity_end,
                velocity_dot_begin, velocity_dot_end, velocity_gradient_begin, velocity_gradient_end,
                material_response_begin + body_force_index,
                material_response_begin + body_force_index + configuration::material::dimension,
                material_response_begin + cauchy_stress_index,
                material_response_begin + cauchy_stress_index + configuration::material::dimension * configuration::material::dimension,
                volume_fraction, test_function, test_function_gradient_begin, test_function_gradient_end, result_begin,
                result_end);

            // Add the contribution from the interphasic force
            for (auto v = std::pair<unsigned int, result_iter>(0, result_begin); v.second != result_end;
                 ++v.first, ++v.second) {
                *v.second += test_function * (*(material_response_begin + interphasic_force_index + v.first));
            }
        }

        template <class configuration, int body_force_index, int cauchy_stress_index,
                  int interphasic_force_index, typename density_type,
                  typename density_dot_type, class density_gradient_iter, class velocity_iter, class velocity_dot_iter,
                  class velocity_gradient_iter, class material_response_iter, class material_response_jacobian_iter,
                  typename volume_fraction_type, typename testFunction_type, class testFunctionGradient_iter,
                  typename interpolationFunction_type, class interpolationFunctionGradient_iter,
                  class full_material_response_dof_gradient_iter, typename dDensityDotdDensity_type,
                  typename dUDotdU_type, typename dUDDotdU_type, class result_iter, class dRdRho_iter, class dRdU_iter,
                  class dRdW_iter, class dRdTheta_iter, class dRdE_iter, class dRdVolumeFraction_iter, class dRdZ_iter,
                  class dRdUMesh_iter>
        void computeBalanceOfLinearMomentum(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const volume_fraction_type &volume_fraction, const testFunction_type &test_function,
            const testFunctionGradient_iter                &test_function_gradient_begin,
            const testFunctionGradient_iter                &test_function_gradient_end,
            const interpolationFunction_type               &interpolation_function,
            const interpolationFunctionGradient_iter       &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter       &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dDensityDotdDensity_type &dDensityDotdDensity, const dUDotdU_type &dUDotdU,
            const dUDDotdU_type &dUDDotdU, const unsigned int phase, result_iter result_begin, result_iter result_end,
            dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end, dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin, dRdW_iter dRdW_end, dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end,
            dRdE_iter dRdE_begin, dRdE_iter dRdE_end, dRdVolumeFraction_iter dRdVolumeFraction_begin,
            dRdVolumeFraction_iter dRdVolumeFraction_end, dRdZ_iter dRdZ_begin, dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the balance of linear momentum including the inter-phasic force
             *
             * body_force_index: The index of the material response vector where the body force is located
             * cauchy_stress_index: The index of the material response vector where the cauchy stress force is located
             * interphasic_force_index: The index of the material response vector where the net interphasic force is
             * located
             *
             * \param &density: The mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot: The partial time derivative of the density \f$ \left( \frac{\partial}{\partial t}
             * \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response vector Jacobian
             * \param &material_response_jacobian_end: The stopping iterator of the material response vector Jacobian
             * \param &volume_fraction: The volume fraction of the phase. Only applied to the Cauchy stress because the
             * density is assumed to be the apparent density i.e., the mass of the phase per unit volume.
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$
             * \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &interpolation_function: The value of the interpolation function \f$ \left( \phi \right) \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the gradient of the interpolation
             * function \f$ \left( \phi_{,i} \right) \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the gradient of the interpolation
             * function \f$ \left( \phi_{,i} \right) \f$
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the
             * material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the
             * material response dof vector
             * \param &dDensityDotdDensity: The total derivative of the time derivative of the density w.r.t. the
             * density
             * \param &dUDotdU: The total derivative of the time derivative of the spatial dof w.r.t. the spatial dof (1
             * if the spatial DOF is the velocity)
             * \param &dUDDotdU: The total derivative of the second time derivative of the spatial dof w.r.t. the
             * spatial dof
             * \param &phase: The phase the balance equation applies to
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             * \param &dRdRho_begin: The starting iterator of the Jacobian of the result w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the Jacobian of the result w.r.t. the density
             * \param &dRdU_begin: The starting iterator of the Jacobian of the result w.r.t. the spatial dof
             * \param &dRdU_end: The stopping iterator of the Jacobian of the result w.r.t. the spatial dof
             * \param &dRdW_begin: The starting iterator of the Jacobian of the result w.r.t. the displacement
             * \param &dRdW_end: The stopping iterator of the Jacobian of the result w.r.t. the displacement
             * \param &dRdTheta_begin: The starting iterator of the Jacobian of the result w.r.t. the temperature
             * \param &dRdTheta_end: The stopping iterator of the Jacobian of the result w.r.t. the temperature
             * \param &dRdE_begin: The starting iterator of the Jacobian of the result w.r.t. the internal energy
             * \param &dRdE_end: The stopping iterator of the Jacobian of the result w.r.t. the internal energy
             * \param &dRdVolumeFraction_begin: The starting iterator of the Jacobian of the result w.r.t. the volume
             * fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the Jacobian of the result w.r.t. the volume
             * fraction
             * \param &dRdZ_begin: The starting iterator of the Jacobian of the result w.r.t. the additional dof
             * \param &dRdZ_end: The stopping iterator of the Jacobian of the result w.r.t. the additional dof
             * \param &dRdUMesh_begin: The starting iterator of the Jacobian of the result w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the Jacobian of the result w.r.t. the mesh displacement
             */

            using result_type = typename std::iterator_traits<result_iter>::value_type;

            std::array<result_type, configuration::dimension>                            dRdRho_phase;
            std::array<result_type, configuration::dimension * configuration::dimension> dRdU_phase;
            std::array<result_type, configuration::dimension * configuration::dimension> dRdB_phase;
            std::array<result_type, configuration::dimension * configuration::dimension * configuration::dimension>
                                                              dRdCauchy_phase;
            std::array<result_type, configuration::dimension> dRdVF_phase;

            computeBalanceOfLinearMomentum<configuration>(
                density, density_dot, density_gradient_begin, density_gradient_end, velocity_begin, velocity_end,
                velocity_dot_begin, velocity_dot_end, velocity_gradient_begin, velocity_gradient_end,
                material_response_begin + body_force_index,
                material_response_begin + body_force_index + configuration::material::dimension,
                material_response_begin + cauchy_stress_index,
                material_response_begin + cauchy_stress_index + configuration::material::dimension * configuration::material::dimension,
                volume_fraction, test_function, test_function_gradient_begin, test_function_gradient_end,
                interpolation_function, interpolation_function_gradient_begin, interpolation_function_gradient_end,
                dDensityDotdDensity, dUDotdU, dUDDotdU, result_begin, result_end, std::begin(dRdRho_phase),
                std::end(dRdRho_phase), std::begin(dRdU_phase), std::end(dRdU_phase), std::begin(dRdB_phase),
                std::end(dRdB_phase), std::begin(dRdCauchy_phase), std::end(dRdCauchy_phase), std::begin(dRdVF_phase),
                std::end(dRdVF_phase), dRdUMesh_begin, dRdUMesh_end);

            // Add the contribution from the interphasic force
            for (auto v = std::pair<unsigned int, result_iter>(0, result_begin); v.second != result_end;
                 ++v.first, ++v.second) {
                *v.second += test_function * (*(material_response_begin + interphasic_force_index + v.first));
            }

            // Set the number of phases
            const unsigned int     nphases       = (unsigned int)(dRdRho_end - dRdRho_begin) / configuration::dimension;

            std::fill(dRdRho_begin, dRdRho_end, 0);
            std::fill(dRdU_begin, dRdU_end, 0);
            std::fill(dRdW_begin, dRdW_end, 0);
            std::fill(dRdTheta_begin, dRdTheta_end, 0);
            std::fill(dRdE_begin, dRdE_end, 0);
            std::fill(dRdVolumeFraction_begin, dRdVolumeFraction_end, 0);
            std::fill(dRdZ_begin, dRdZ_end, 0);

            for (unsigned int i = 0; i < configuration::dimension; ++i) {
                // BODY FORCE CONTRIBUTIONS
                for (unsigned int j = 0; j < configuration::dimension; ++j) {
                    // density
                    for (auto p = std::pair<unsigned int, dRdRho_iter>(0, dRdRho_begin + nphases * i);
                         p.second != dRdRho_begin + nphases * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdB_phase[configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (body_force_index + j) +
                                        nphases * configuration::material::dof::density_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdB_phase[configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (body_force_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::density_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // velocity
                    for (auto p =
                             std::pair<unsigned int, dRdU_iter>(0, dRdU_begin + nphases * configuration::dimension * i);
                         p.second != dRdU_begin + nphases * configuration::dimension * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdB_phase[configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (body_force_index + j) +
                                        nphases * configuration::material::dof::velocity_index + p.first)) *
                                     interpolation_function * dUDotdU;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdB_phase[configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (body_force_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::velocity_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a)) * dUDotdU;
                        }
                    }

                    // displacement
                    for (auto p =
                             std::pair<unsigned int, dRdW_iter>(0, dRdW_begin + nphases * configuration::dimension * i);
                         p.second != dRdW_begin + nphases * configuration::dimension * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdB_phase[configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (body_force_index + j) +
                                        nphases * configuration::material::dof::displacement_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdB_phase[configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (body_force_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::displacement_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // temperature
                    for (auto p = std::pair<unsigned int, dRdTheta_iter>(0, dRdTheta_begin + nphases * i);
                         p.second != dRdTheta_begin + nphases * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdB_phase[configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (body_force_index + j) +
                                        nphases * configuration::material::dof::temperature_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdB_phase[configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (body_force_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::temperature_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // internal energy
                    for (auto p = std::pair<unsigned int, dRdE_iter>(0, dRdE_begin + nphases * i);
                         p.second != dRdE_begin + nphases * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdB_phase[configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (body_force_index + j) +
                                        nphases * configuration::material::dof::internal_energy_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdB_phase[configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (body_force_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::internal_energy_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // volume fraction
                    for (auto p =
                             std::pair<unsigned int, dRdVolumeFraction_iter>(0, dRdVolumeFraction_begin + nphases * i);
                         p.second != dRdVolumeFraction_begin + nphases * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdB_phase[configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (body_force_index + j) +
                                        nphases * configuration::material::dof::volume_fraction_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdB_phase[configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (body_force_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::volume_fraction_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // additional dof
                    for (auto p = std::pair<unsigned int, dRdZ_iter>(0, dRdZ_begin + configuration::material::dof::num_additional_dof * i);
                         p.second != dRdZ_begin + configuration::material::dof::num_additional_dof * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdB_phase[configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (body_force_index + j) +
                                        nphases * configuration::material::dof::additional_dof_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdB_phase[configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (body_force_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::additional_dof_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // mesh displacement
                    for (unsigned int I = 0; I < (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof); ++I) {
                        for (unsigned int k = 0; k < configuration::material::dimension; ++k) {
                            for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                                *(dRdUMesh_begin + configuration::dimension * i + a) -=
                                    dRdB_phase[configuration::dimension * i + j] *
                                    (*(material_response_jacobian_begin +
                                       (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                           (body_force_index + j) +
                                       (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) + configuration::material::dimension * I +
                                       k)) *
                                    (*(full_material_response_dof_gradient_begin + configuration::material::dimension * I + a)) *
                                    (*(interpolation_function_gradient_begin + k));
                            }
                        }
                    }
                }

                // CAUCHY STRESS CONTRIBUTIONS
                for (unsigned int j = 0; j < configuration::dimension * configuration::dimension; ++j) {
                    // density
                    for (auto p = std::pair<unsigned int, dRdRho_iter>(0, dRdRho_begin + nphases * i);
                         p.second != dRdRho_begin + nphases * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (cauchy_stress_index + j) +
                                        nphases * configuration::material::dof::density_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (cauchy_stress_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::density_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // velocity
                    for (auto p =
                             std::pair<unsigned int, dRdU_iter>(0, dRdU_begin + nphases * configuration::dimension * i);
                         p.second != dRdU_begin + nphases * configuration::dimension * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (cauchy_stress_index + j) +
                                        nphases * configuration::material::dof::velocity_index + p.first)) *
                                     interpolation_function * dUDotdU;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (cauchy_stress_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::velocity_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a)) * dUDotdU;
                        }
                    }

                    // displacement
                    for (auto p =
                             std::pair<unsigned int, dRdW_iter>(0, dRdW_begin + nphases * configuration::dimension * i);
                         p.second != dRdW_begin + nphases * configuration::dimension * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (cauchy_stress_index + j) +
                                        nphases * configuration::material::dof::displacement_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (cauchy_stress_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::displacement_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // temperature
                    for (auto p = std::pair<unsigned int, dRdTheta_iter>(0, dRdTheta_begin + nphases * i);
                         p.second != dRdTheta_begin + nphases * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (cauchy_stress_index + j) +
                                        nphases * configuration::material::dof::temperature_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (cauchy_stress_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::temperature_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // internal energy
                    for (auto p = std::pair<unsigned int, dRdE_iter>(0, dRdE_begin + nphases * i);
                         p.second != dRdE_begin + nphases * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (cauchy_stress_index + j) +
                                        nphases * configuration::material::dof::internal_energy_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (cauchy_stress_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::internal_energy_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // volume fraction
                    for (auto p =
                             std::pair<unsigned int, dRdVolumeFraction_iter>(0, dRdVolumeFraction_begin + nphases * i);
                         p.second != dRdVolumeFraction_begin + nphases * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (cauchy_stress_index + j) +
                                        nphases * configuration::material::dof::volume_fraction_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (cauchy_stress_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::volume_fraction_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // additional dof
                    for (auto p = std::pair<unsigned int, dRdZ_iter>(0, dRdZ_begin + configuration::material::dof::num_additional_dof * i);
                         p.second != dRdZ_begin + configuration::material::dof::num_additional_dof * (i + 1); ++p.first, ++p.second) {
                        *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (cauchy_stress_index + j) +
                                        nphases * configuration::material::dof::additional_dof_index + p.first)) *
                                     interpolation_function;

                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *p.second += dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                         (*(material_response_jacobian_begin +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                (1 + configuration::material::dimension) * (cauchy_stress_index + j) +
                                            (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                            configuration::material::dimension * (nphases * configuration::material::dof::additional_dof_index + p.first) + a)) *
                                         (*(interpolation_function_gradient_begin + a));
                        }
                    }

                    // mesh displacement
                    for (unsigned int I = 0; I < (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof); ++I) {
                        for (unsigned int k = 0; k < configuration::material::dimension; ++k) {
                            for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                                *(dRdUMesh_begin + configuration::dimension * i + a) -=
                                    dRdCauchy_phase[configuration::dimension * configuration::dimension * i + j] *
                                    (*(material_response_jacobian_begin +
                                       (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                           (cauchy_stress_index + j) +
                                       (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) + configuration::material::dimension * I +
                                       k)) *
                                    (*(full_material_response_dof_gradient_begin + configuration::material::dimension * I + a)) *
                                    (*(interpolation_function_gradient_begin + k));
                            }
                        }
                    }
                }

                // INTERPHASIC FORCE CONTRIBUTIONS
                // density
                for (auto p = std::pair<unsigned int, dRdRho_iter>(0, dRdRho_begin + nphases * i);
                     p.second != dRdRho_begin + nphases * (i + 1); ++p.first, ++p.second) {
                    // DOF value contributions
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (interphasic_force_index + i) +
                                    nphases * configuration::material::dof::density_index + p.first)) *
                                 interpolation_function;

                    // DOF spatial gradient contributions
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *p.second += test_function *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (interphasic_force_index + i) +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                        configuration::material::dimension * (nphases * configuration::material::dof::density_index + p.first) + a)) *
                                     (*(interpolation_function_gradient_begin + a));
                    }
                }

                *(dRdRho_begin + nphases * i + phase) += dRdRho_phase[i];

                // velocity
                for (auto p =
                         std::pair<unsigned int, dRdU_iter>(0, dRdU_begin + nphases * configuration::dimension * i);
                     p.second != dRdU_begin + nphases * configuration::dimension * (i + 1); ++p.first, ++p.second) {
                    // DOF value contributions
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (interphasic_force_index + i) +
                                    nphases * configuration::material::dof::velocity_index + p.first)) *
                                 interpolation_function * dUDotdU;

                    // DOF spatial gradient contributions
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *p.second += test_function *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (interphasic_force_index + i) +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                        configuration::material::dimension * (nphases * configuration::material::dof::velocity_index + p.first) + a)) *
                                     (*(interpolation_function_gradient_begin + a)) * dUDotdU;
                    }
                }

                std::transform(dRdU_begin + nphases * configuration::dimension * i + configuration::dimension * phase,
                               dRdU_begin + nphases * configuration::dimension * i +
                                   configuration::dimension * (phase + 1),
                               std::begin(dRdU_phase) + configuration::dimension * i,
                               dRdU_begin + nphases * configuration::dimension * i + configuration::dimension * phase,
                               std::plus<result_type>());

                // displacement
                for (auto p =
                         std::pair<unsigned int, dRdW_iter>(0, dRdW_begin + nphases * configuration::dimension * i);
                     p.second != dRdW_begin + nphases * configuration::dimension * (i + 1); ++p.first, ++p.second) {
                    // DOF value contributions
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (interphasic_force_index + i) +
                                    nphases * configuration::material::dof::displacement_index + p.first)) *
                                 interpolation_function;

                    // DOF spatial gradient contributions
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *p.second += test_function *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (interphasic_force_index + i) +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                        configuration::material::dimension * (nphases * configuration::material::dof::displacement_index + p.first) + a)) *
                                     (*(interpolation_function_gradient_begin + a));
                    }
                }

                // temperature
                for (auto p = std::pair<unsigned int, dRdTheta_iter>(0, dRdTheta_begin + nphases * i);
                     p.second != dRdTheta_begin + nphases * (i + 1); ++p.first, ++p.second) {
                    // DOF value contributions
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (interphasic_force_index + i) +
                                    nphases * configuration::material::dof::temperature_index + p.first)) *
                                 interpolation_function;

                    // DOF spatial gradient contributions
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *p.second += test_function *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (interphasic_force_index + i) +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                        configuration::material::dimension * (nphases * configuration::material::dof::temperature_index + p.first) + a)) *
                                     (*(interpolation_function_gradient_begin + a));
                    }
                }

                // internal energy
                for (auto p = std::pair<unsigned int, dRdE_iter>(0, dRdE_begin + nphases * i);
                     p.second != dRdE_begin + nphases * (i + 1); ++p.first, ++p.second) {
                    // DOF value contributions
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (interphasic_force_index + i) +
                                    nphases * configuration::material::dof::internal_energy_index + p.first)) *
                                 interpolation_function;

                    // DOF spatial gradient contributions
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *p.second += test_function *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (interphasic_force_index + i) +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                        configuration::material::dimension * (nphases * configuration::material::dof::internal_energy_index + p.first) + a)) *
                                     (*(interpolation_function_gradient_begin + a));
                    }
                }

                // volume fraction
                for (auto p = std::pair<unsigned int, dRdVolumeFraction_iter>(0, dRdVolumeFraction_begin + nphases * i);
                     p.second != dRdVolumeFraction_begin + nphases * (i + 1); ++p.first, ++p.second) {
                    // DOF value contributions
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (interphasic_force_index + i) +
                                    nphases * configuration::material::dof::volume_fraction_index + p.first)) *
                                 interpolation_function;

                    // DOF spatial gradient contributions
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *p.second += test_function *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (interphasic_force_index + i) +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                        configuration::material::dimension * (nphases * configuration::material::dof::volume_fraction_index + p.first) + a)) *
                                     (*(interpolation_function_gradient_begin + a));
                    }
                }
                *(dRdVolumeFraction_begin + nphases * i + phase) += dRdVF_phase[i] * interpolation_function;

                // additional dof
                for (auto p = std::pair<unsigned int, dRdZ_iter>(0, dRdZ_begin + configuration::material::dof::num_additional_dof * i);
                     p.second != dRdZ_begin + configuration::material::dof::num_additional_dof * (i + 1); ++p.first, ++p.second) {
                    // DOF value contributions
                    *p.second += test_function *
                                 (*(material_response_jacobian_begin +
                                    (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                        (interphasic_force_index + i) +
                                    nphases * configuration::material::dof::additional_dof_index + p.first)) *
                                 interpolation_function;

                    // DOF spatial gradient contributions
                    for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                        *p.second += test_function *
                                     (*(material_response_jacobian_begin +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                            (interphasic_force_index + i) +
                                        (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) +
                                        configuration::material::dimension * (nphases * configuration::material::dof::additional_dof_index + p.first) + a)) *
                                     (*(interpolation_function_gradient_begin + a));
                    }
                }

                // mesh displacement
                for (unsigned int I = 0; I < (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof); ++I) {
                    for (unsigned int j = 0; j < configuration::material::dimension; ++j) {
                        for (unsigned int a = 0; a < configuration::material::dimension; ++a) {
                            *(dRdUMesh_begin + configuration::dimension * i + a) -=
                                test_function *
                                (*(material_response_jacobian_begin +
                                   (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) * (1 + configuration::material::dimension) *
                                       (interphasic_force_index + i) +
                                   (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) + configuration::material::dimension * I + j)) *
                                (*(full_material_response_dof_gradient_begin + configuration::material::dimension * I + a)) *
                                (*(interpolation_function_gradient_begin + j));
                        }
                    }
                }

                for (unsigned int a = 0; a < configuration::dimension; ++a) {
                    *(dRdUMesh_begin + configuration::dimension * i + a) +=
                        test_function * (*(material_response_begin + interphasic_force_index + i)) *
                        (*(interpolation_function_gradient_begin + a));
                }
            }
        }

        template <class configuration, int body_force_index, int cauchy_stress_index,
                  int interphasic_force_index, class density_iter,
                  class density_dot_iter, class density_gradient_iter, class velocity_iter, class velocity_dot_iter,
                  class velocity_gradient_iter, class material_response_iter, class material_response_jacobian_iter,
                  class volume_fraction_iter, typename testFunction_type, class testFunctionGradient_iter,
                  typename interpolationFunction_type, class interpolationFunctionGradient_iter,
                  class full_material_response_dof_gradient_iter, typename dDensityDotdDensity_type,
                  typename dUDotdU_type, typename dUDDotdU_type, class result_iter, class dRdRho_iter, class dRdU_iter,
                  class dRdW_iter, class dRdTheta_iter, class dRdE_iter, class dRdVolumeFraction_iter, class dRdZ_iter,
                  class dRdUMesh_iter>
        void computeBalanceOfLinearMomentum(
            const density_iter &density_begin, const density_iter &density_end,
            const density_dot_iter &density_dot_begin, const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
            const testFunction_type &test_function, const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter                &test_function_gradient_end,
            const interpolationFunction_type               &interpolation_function,
            const interpolationFunctionGradient_iter       &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter       &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dDensityDotdDensity_type &dDensityDotdDensity, const dUDotdU_type &dUDotdU,
            const dUDDotdU_type &dUDDotdU, result_iter result_begin, result_iter result_end, dRdRho_iter dRdRho_begin,
            dRdRho_iter dRdRho_end, dRdU_iter dRdU_begin, dRdU_iter dRdU_end, dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end, dRdE_iter dRdE_begin, dRdE_iter dRdE_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdZ_iter dRdZ_begin, dRdZ_iter dRdZ_end, dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the balance of linear momentum including the inter-phasic force
             *
             * body_force_index: The index of the material response vector where the body force is located
             * cauchy_stress_index: The index of the material response vector where the cauchy stress force is located
             * interphasic_force_index: The index of the material response vector where the net interphasic force is
             * located
             *
             * \param &density_begin: The starting iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_end: The stopping iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response vector Jacobian
             * \param &material_response_jacobian_end: The stopping iterator of the material response vector Jacobian
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of the phase. Only applied to
             * the Cauchy stress because the density is assumed to be the apparent density i.e., the mass of the phase
             * per unit volume.
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of the phase. Only applied to
             * the Cauchy stress because the density is assumed to be the apparent density i.e., the mass of the phase
             * per unit volume.
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$
             * \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &interpolation_function: The value of the interpolation function \f$ \left( \phi \right) \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the gradient of the interpolation
             * function \f$ \left( \phi_{,i} \right) \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the gradient of the interpolation
             * function \f$ \left( \phi_{,i} \right) \f$
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the
             * material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the
             * material response dof vector
             * \param &dDensityDotdDensity: The total derivative of the time derivative of the density w.r.t. the
             * density
             * \param &dUDotdU: The total derivative of the time derivative of the spatial dof w.r.t. the spatial dof (1
             * if the spatial DOF is the velocity)
             * \param &dUDDotdU: The total derivative of the second time derivative of the spatial dof w.r.t. the
             * spatial dof
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             * \param &dRdRho_begin: The starting iterator of the Jacobian of the result w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the Jacobian of the result w.r.t. the density
             * \param &dRdU_begin: The starting iterator of the Jacobian of the result w.r.t. the spatial dof
             * \param &dRdU_end: The stopping iterator of the Jacobian of the result w.r.t. the spatial dof
             * \param &dRdW_begin: The starting iterator of the Jacobian of the result w.r.t. the displacement
             * \param &dRdW_end: The stopping iterator of the Jacobian of the result w.r.t. the displacement
             * \param &dRdTheta_begin: The starting iterator of the Jacobian of the result w.r.t. the temperature
             * \param &dRdTheta_end: The stopping iterator of the Jacobian of the result w.r.t. the temperature
             * \param &dRdE_begin: The starting iterator of the Jacobian of the result w.r.t. the internal energy
             * \param &dRdE_end: The stopping iterator of the Jacobian of the result w.r.t. the internal energy
             * \param &dRdVolumeFraction_begin: The starting iterator of the Jacobian of the result w.r.t. the volume
             * fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the Jacobian of the result w.r.t. the volume
             * fraction
             * \param &dRdZ_begin: The starting iterator of the Jacobian of the result w.r.t. the additional dof
             * \param &dRdZ_end: The stopping iterator of the Jacobian of the result w.r.t. the additional dof
             * \param &dRdUMesh_begin: The starting iterator of the Jacobian of the result w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the Jacobian of the result w.r.t. the mesh displacement
             */

            const unsigned int nphases = (unsigned int)(density_end - density_begin);

            const unsigned int material_response_size =
                (unsigned int)(material_response_end - material_response_begin) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(density_dot_end - density_dot_begin),
                                         "The length of density dot and density must be the same")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * configuration::dimension ==
                                             (unsigned int)(density_gradient_end - density_gradient_begin),
                                         "The length of the density gradient and the density must be consistent")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * configuration::dimension ==
                                             (unsigned int)(velocity_end - velocity_begin),
                                         "The length of the velocity and the density must be consistent")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * configuration::dimension ==
                                             (unsigned int)(velocity_dot_end - velocity_dot_begin),
                                         "The length of the velocity dot and the density must be consistent")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases * configuration::dimension * configuration::dimension ==
                                             (unsigned int)(velocity_gradient_end - velocity_gradient_begin),
                                         "The length of the velocity gradient and the density must be consistent")

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(volume_fraction_end - volume_fraction_begin),
                                         "The length of the volume fraction and the density must be consistent")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                (body_force_index + configuration::dimension) < material_response_size,
                "The material response vector must be larger than the body force index plus the dimension")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                (cauchy_stress_index + configuration::dimension * configuration::dimension) < material_response_size,
                "The material response vector must be larger than the Cauchy stress index plus the dimension squared")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                (interphasic_force_index + configuration::dimension) < material_response_size,
                "The material response vector must be larger than the interphasic force index plus the dimension")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size * (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                        (1 + configuration::material::dimension) ==
                    (unsigned int)(material_response_jacobian_end - material_response_jacobian_begin),

                "The material response jacobian must have a consistent size with the material response vector and the "
                "configuration::material::dof::num_dof\n  number of phases           : " +
                    std::to_string(nphases) +
                    "\n  material_response_size     : " + std::to_string(material_response_size) +
                    "\n  dof per phase              : " + std::to_string(configuration::material::dof::num_phase_dof) +
                    "\n  additional dof             : " + std::to_string(configuration::material::dof::num_additional_dof) +
                    "\n  material response dimension: " + std::to_string(configuration::material::dimension) +
                    "\n  expected jacobian size     : " +
                    std::to_string(nphases * material_response_size * (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                   (1 + configuration::material::dimension)) +
                    "\n  actual jacobian size       : " +
                    std::to_string((unsigned int)(material_response_jacobian_end - material_response_jacobian_begin)))

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension ==
                                             (unsigned int)(interpolation_function_gradient_end -
                                                            interpolation_function_gradient_begin),
                                         "The interpolation function gradient must have a size of dim")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases ==
                                             (unsigned int)(result_end - result_begin),
                                         "The result vector must be the same size as the density vector")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * nphases * 1 ==
                                             (unsigned int)(dRdRho_end - dRdRho_begin),
                                         "dRdRho must have a consistent size with the density vector")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * nphases * configuration::material::dimension ==
                                             (unsigned int)(dRdU_end - dRdU_begin),
                                         "dRdU must have a consistent size with the density vector")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * nphases * configuration::material::dimension ==
                                             (unsigned int)(dRdW_end - dRdW_begin),
                                         "dRdW must have a consistent size with the density vector")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * nphases * 1 ==
                                             (unsigned int)(dRdTheta_end - dRdTheta_begin),
                                         "dRdTheta must have a consistent size with the density vector")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * nphases * 1 ==
                                             (unsigned int)(dRdE_end - dRdE_begin),
                                         "dRdE must have a consistent size with the density vector")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * configuration::material::dof::num_additional_dof ==
                                             (unsigned int)(dRdZ_end - dRdZ_begin),
                                         "dRdZ must have a consistent size with the density vector")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * nphases * 1 ==
                                             (unsigned int)(dRdVolumeFraction_end - dRdVolumeFraction_begin),
                                         "dRdVF must have a consistent size with the density vector")

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * configuration::dimension ==
                                             (unsigned int)(dRdUMesh_end - dRdUMesh_begin),
                                         "dRdUMesh must have a consistent size with the density vector")

            for (auto v = std::pair<unsigned int, density_iter>(0, density_begin); v.second != density_end;
                 ++v.first, ++v.second) {
                computeBalanceOfLinearMomentum<
                    configuration, body_force_index, cauchy_stress_index,
                    interphasic_force_index>(
                    *(density_begin + v.first), *(density_dot_begin + v.first),
                    density_gradient_begin + configuration::dimension * v.first,
                    density_gradient_begin + configuration::dimension * (v.first + 1),
                    velocity_begin + configuration::dimension * v.first,
                    velocity_begin + configuration::dimension * (v.first + 1),
                    velocity_dot_begin + configuration::dimension * v.first,
                    velocity_dot_begin + configuration::dimension * (v.first + 1),
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * v.first,
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * (v.first + 1),
                    material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * (v.first + 1),
                    material_response_jacobian_begin + material_response_size *
                                                           (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                           (1 + configuration::material::dimension) * v.first,
                    material_response_jacobian_begin + material_response_size *
                                                           (nphases * configuration::material::dof::num_phase_dof + configuration::material::dof::num_additional_dof) *
                                                           (1 + configuration::material::dimension) * (v.first + 1),
                    *(volume_fraction_begin + v.first), test_function, test_function_gradient_begin,
                    test_function_gradient_end, interpolation_function, interpolation_function_gradient_begin,
                    interpolation_function_gradient_end, full_material_response_dof_gradient_begin,
                    full_material_response_dof_gradient_end, dDensityDotdDensity, dUDotdU, dUDDotdU, v.first,
                    result_begin + configuration::dimension * v.first,
                    result_begin + configuration::dimension * (v.first + 1),
                    dRdRho_begin + configuration::dimension * nphases * v.first,
                    dRdRho_begin + configuration::dimension * nphases * (v.first + 1),
                    dRdU_begin + configuration::dimension * configuration::dimension * nphases * v.first,
                    dRdU_begin + configuration::dimension * configuration::dimension * nphases * (v.first + 1),
                    dRdW_begin + configuration::dimension * configuration::dimension * nphases * v.first,
                    dRdW_begin + configuration::dimension * configuration::dimension * nphases * (v.first + 1),
                    dRdTheta_begin + configuration::dimension * nphases * v.first,
                    dRdTheta_begin + configuration::dimension * nphases * (v.first + 1),
                    dRdE_begin + configuration::dimension * nphases * v.first,
                    dRdE_begin + configuration::dimension * nphases * (v.first + 1),
                    dRdVolumeFraction_begin + configuration::dimension * nphases * v.first,
                    dRdVolumeFraction_begin + configuration::dimension * nphases * (v.first + 1),
                    dRdZ_begin + configuration::dimension * configuration::material::dof::num_additional_dof * v.first,
                    dRdZ_begin + configuration::dimension * configuration::material::dof::num_additional_dof * (v.first + 1),
                    dRdUMesh_begin + configuration::dimension * configuration::dimension * v.first,
                    dRdUMesh_begin + configuration::dimension * configuration::dimension * (v.first + 1));
            }
        }

        template <class configuration, int body_force_index, int cauchy_stress_index,
                  int interphasic_force_index, class density_iter, class density_dot_iter, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
                  class material_response_iter, class volume_fraction_iter, typename testFunction_type,
                  class testFunctionGradient_iter, class result_iter>
        void computeBalanceOfLinearMomentum(
            const density_iter &density_begin, const density_iter &density_end,
            const density_dot_iter &density_dot_begin, const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
            const testFunction_type &test_function, const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter &test_function_gradient_end, result_iter result_begin,
            result_iter result_end) {
            /*!
             * Compute the balance of linear momentum including the inter-phasic force
             *
             * body_force_index: The index of the material response vector where the body force is located
             * cauchy_stress_index: The index of the material response vector where the cauchy stress force is located
             * interphasic_force_index: The index of the material response vector where the net interphasic force is
             * located
             *
             * \param &density_begin: The starting iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_end: The stopping iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of the phase. Only applied to
             * the Cauchy stress because the density is assumed to be the apparent density i.e., the mass of the phase
             * per unit volume.
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of the phase. Only applied to
             * the Cauchy stress because the density is assumed to be the apparent density i.e., the mass of the phase
             * per unit volume.
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$
             * \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            const unsigned int nphases = (unsigned int)(density_end - density_begin);
            const unsigned int material_response_size =
                (unsigned int)(material_response_end - material_response_begin) / nphases;

            for (auto v = std::pair<unsigned int, density_iter>(0, density_begin); v.second != density_end;
                 ++v.first, ++v.second) {
                computeBalanceOfLinearMomentum<configuration, body_force_index,
                                               cauchy_stress_index, interphasic_force_index>(
                    *(density_begin + v.first), *(density_dot_begin + v.first),
                    density_gradient_begin + configuration::dimension * v.first,
                    density_gradient_begin + configuration::dimension * (v.first + 1),
                    velocity_begin + configuration::dimension * v.first,
                    velocity_begin + configuration::dimension * (v.first + 1),
                    velocity_dot_begin + configuration::dimension * v.first,
                    velocity_dot_begin + configuration::dimension * (v.first + 1),
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * v.first,
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * (v.first + 1),
                    material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * (v.first + 1),
                    *(volume_fraction_begin + v.first), test_function, test_function_gradient_begin,
                    test_function_gradient_end, result_begin + configuration::dimension * v.first,
                    result_begin + configuration::dimension * (v.first + 1));
            }
        }

        template <class configuration, typename density_type, typename density_dot_type, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter, class body_force_iter,
                  class cauchy_stress_iter, typename volume_fraction_type, typename testFunction_type,
                  class testFunctionGradient_iter, typename interpolationFunction_type,
                  class interpolationFunctionGradient_iter, class result_iter, class dRdRho_iter, class dRdU_iter,
                  class dRdB_iter, class dRdCauchy_iter, class dRdVolumeFraction_iter, class dRdUMesh_iter,
                  typename dRhoDotdRho_type, typename dUDotdU_type, typename dUDDotdU_type>
        void computeBalanceOfLinearMomentum(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin, const body_force_iter &body_force_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction, const testFunction_type &test_function,
            const testFunctionGradient_iter          &test_function_gradient_begin,
            const testFunctionGradient_iter          &test_function_gradient_end,
            const interpolationFunction_type         &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const dRhoDotdRho_type dRhoDotdRho, const dUDotdU_type dUDotdU, const dUDDotdU_type dUDDotdU,
            result_iter result_begin, result_iter result_end, dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end, dRdB_iter dRdB_begin, dRdB_iter dRdB_end,
            dRdCauchy_iter dRdCauchy_begin, dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the balance of linear momentum
             *
             * \param &density: The mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot: The partial time derivative of the density \f$ \left( \frac{\partial}{\partial t}
             * \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right)
             * \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the cauchy stress \f$ \left( \sigma_{ij} \right)
             * \f$
             * \param &cauchy_stress_end: The stopping iterator of the cauchy stress \f$ \left( \sigma_{ij} \right) \f$
             * \param &volume_fraction: The volume fraction of the phase. Only applied to the Cauchy stress because the
             * density is assumed to be the apparent density i.e., the mass of the phase per unit volume.
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$
             * \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &interpolation_function: The value of the inerpolation function \f$ \left( \phi \right) \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the gradient of the interpolation
             * function \f$ \left( \phi_{,i} \right) \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the gradient of the interpolation
             * function \f$ \left( \phi_{,i} \right) \f$
             * \param &dRhoDotdRho: The Jacobian of the time derivative of the density w.r.t. the density
             * \param &dUDotdU: The Jacobian of the time derivative of the displacement degree of freedom w.r.t. the
             * displacement degree of freedom
             * \param &dUDDotdU: The Jacobian of the second time derivative of the displacement degree of freedom w.r.t.
             * the displacement degree of freedom density
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             * \param &dRdRho_begin: The starting iterator of the Jacobian of the result w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the Jacobian of the result w.r.t. the density
             * \param &dRdU_begin: The starting iterator of the Jacobian of the result w.r.t. the displacement degree of
             * freedom
             * \param &dRdU_end: The stopping iterator of the Jacobian of the result w.r.t. the displacement degree of
             * freedom
             * \param &dRdB_begin: The starting iterator of the Jacobian of the result w.r.t. the body force
             * \param &dRdB_end: The stopping iterator of the Jacobian of the result w.r.t. the body force
             * \param &dRdCauchy_begin: The starting iterator of the Jacobian of the result w.r.t. the Cauchy stress
             * \param &dRdCauchy_end: The stopping iterator of the Jacobian of the result w.r.t. the Cauchy stress
             * \param &dRdVolumeFraction_begin: The starting iterator of the Jacobian of the result w.r.t. the volume
             * fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the Jacobian of the result w.r.t. the volume
             * fraction
             * \param &dRdUMesh_begin: The starting iterator of the Jacobian of the result w.r.t. the mesh displacement
             * (includes the terms due to the volume change)
             * \param &dRdUMesh_end: The stopping iterator of the Jacobian of the result w.r.t. the mesh displacement
             * (includes the terms due to the volume change)
             */

            using dRdRhoDot_type  = typename std::iterator_traits<dRdRho_iter>::value_type;
            using dRdGradRho_type = typename std::iterator_traits<dRdRho_iter>::value_type;

            using dRdVDot_type  = typename std::iterator_traits<dRdU_iter>::value_type;
            using dRdGradV_type = typename std::iterator_traits<dRdU_iter>::value_type;

            using dRdUMesh_type = typename std::iterator_traits<dRdUMesh_iter>::value_type;

            using result_type = typename std::iterator_traits<result_iter>::value_type;
            using final_type  = decltype(std::declval<result_type &>() * std::declval<testFunction_type &>());

            std::array<result_type, configuration::dimension> non_divergence_result;

            std::array<result_type, configuration::dimension> divergence_result;

            std::array<dRdRhoDot_type, configuration::dimension>                             dNonDivRdRhoDot;
            std::array<dRdGradRho_type, configuration::dimension * configuration::dimension> dNonDivRdGradRho;

            std::array<dRdVDot_type, configuration::dimension * configuration::dimension> dNonDivRdVDot;
            std::array<dRdGradV_type, configuration::dimension * configuration::dimension * configuration::dimension>
                dNonDivRdGradV;

            std::array<dRdUMesh_type, configuration::dimension * configuration::dimension> dDivRdTestFunctionGradient;

            // Compute the non-divergence part of the balance of linear momentum
            computeBalanceOfLinearMomentumNonDivergence<configuration>(
                density, density_dot, density_gradient_begin, density_gradient_end, velocity_begin, velocity_end,
                velocity_dot_begin, velocity_dot_end, velocity_gradient_begin, velocity_gradient_end, body_force_begin,
                body_force_end, std::begin(non_divergence_result), std::end(non_divergence_result), dRdRho_begin,
                dRdRho_end, std::begin(dNonDivRdRhoDot), std::end(dNonDivRdRhoDot), std::begin(dNonDivRdGradRho),
                std::end(dNonDivRdGradRho), dRdU_begin, dRdU_end, std::begin(dNonDivRdVDot), std::end(dNonDivRdVDot),
                std::begin(dNonDivRdGradV), std::end(dNonDivRdGradV), dRdB_begin, dRdB_end);

            std::transform(std::cbegin(non_divergence_result), std::cend(non_divergence_result),
                           std::begin(non_divergence_result),
                           std::bind(std::multiplies<final_type>(), std::placeholders::_1, test_function));

            // Compute the divergence part of the balance of linear momentum
            computeBalanceOfLinearMomentumDivergence<configuration>(
                test_function_gradient_begin, test_function_gradient_end, cauchy_stress_begin, cauchy_stress_end,
                volume_fraction, std::begin(divergence_result), std::end(divergence_result),
                std::begin(dDivRdTestFunctionGradient), std::end(dDivRdTestFunctionGradient), dRdCauchy_begin,
                dRdCauchy_end, dRdVolumeFraction_begin, dRdVolumeFraction_end);

            // Assemble the full balance equation
            std::transform(std::cbegin(non_divergence_result), std::cend(non_divergence_result),
                           std::cbegin(divergence_result), result_begin, std::plus<final_type>());

            // Assemble the Jacobians
            std::fill(dRdUMesh_begin, dRdUMesh_end, 0);
            for (unsigned int i = 0; i < configuration::dimension; ++i) {
                *(dRdRho_begin + i) =
                    test_function * ((*(dRdRho_begin + i)) + dNonDivRdRhoDot[i] * dRhoDotdRho) * interpolation_function;

                for (unsigned int j = 0; j < configuration::dimension; ++j) {
                    *(dRdRho_begin + i) += test_function * dNonDivRdGradRho[configuration::dimension * i + j] *
                                           (*(interpolation_function_gradient_begin + j));

                    *(dRdU_begin + configuration::dimension * i + j) =
                        test_function *
                        ((*(dRdU_begin + configuration::dimension * i + j)) * dUDotdU +
                         dNonDivRdVDot[configuration::dimension * i + j] * dUDDotdU) *
                        interpolation_function;

                    *(dRdB_begin + configuration::dimension * i + j) =
                        test_function * (*(dRdB_begin + configuration::dimension * i + j));

                    *(dRdUMesh_begin + configuration::dimension * i + j) +=
                        (*(result_begin + i)) * (*(interpolation_function_gradient_begin + j));

                    for (unsigned int k = 0; k < configuration::dimension; ++k) {
                        *(dRdU_begin + configuration::dimension * i + j) +=
                            test_function *
                            dNonDivRdGradV[configuration::dimension * configuration::dimension * i +
                                           configuration::dimension * j + k] *
                            dUDotdU * (*(interpolation_function_gradient_begin + k));

                        *(dRdUMesh_begin + configuration::dimension * i + k) -=
                            (test_function * dNonDivRdGradRho[configuration::dimension * i + j] *
                                 (*(density_gradient_begin + k)) +
                             dDivRdTestFunctionGradient[configuration::dimension * i + j] *
                                 (*(test_function_gradient_begin + k))) *
                            (*(interpolation_function_gradient_begin + j));

                        for (unsigned int a = 0; a < configuration::dimension; ++a) {
                            *(dRdUMesh_begin + configuration::dimension * i + a) -=
                                test_function *
                                dNonDivRdGradV[configuration::dimension * configuration::dimension * i +
                                               configuration::dimension * j + k] *
                                (*(velocity_gradient_begin + configuration::dimension * j + a)) *
                                (*(interpolation_function_gradient_begin + k));
                        }
                    }
                }
            }
        }

        template <class configuration, class density_iter, class density_dot_iter, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter, class body_force_iter,
                  class cauchy_stress_iter, class volume_fraction_iter, typename testFunction_type,
                  class testFunctionGradient_iter, typename interpolationFunction_type,
                  class interpolationFunctionGradient_iter, class result_iter, class dRdRho_iter, class dRdU_iter,
                  class dRdB_iter, class dRdCauchy_iter, class dRdVolumeFraction_iter, class dRdUMesh_iter,
                  typename dRhoDotdRho_type, typename dUDotdU_type, typename dUDDotdU_type>
        void computeBalanceOfLinearMomentum(
            const density_iter &density_begin, const density_iter &density_end,
            const density_dot_iter &density_dot_begin, const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin, const body_force_iter &body_force_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
            const testFunction_type &test_function, const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter          &test_function_gradient_end,
            const interpolationFunction_type         &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const dRhoDotdRho_type dRhoDotdRho, const dUDotdU_type dUDotdU, const dUDDotdU_type dUDDotdU,
            result_iter result_begin, result_iter result_end, dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end, dRdB_iter dRdB_begin, dRdB_iter dRdB_end,
            dRdCauchy_iter dRdCauchy_begin, dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end) {
            /*!
             * Compute the balance of linear momentum for the multiphase continuum
             *
             * \param &density_begin: The starting iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_end: The stopping iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right)
             * \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the cauchy stress \f$ \left( \sigma_{ij} \right)
             * \f$
             * \param &cauchy_stress_end: The stopping iterator of the cauchy stress \f$ \left( \sigma_{ij} \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of the phase. Only applied to
             * the Cauchy stress because the density is assumed to be the apparent density i.e., the mass of the phase
             * per unit volume.
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of the phase. Only applied to
             * the Cauchy stress because the density is assumed to be the apparent density i.e., the mass of the phase
             * per unit volume.
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$
             * \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &interpolation_function: The value of the inerpolation function \f$ \left( \phi \right) \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the gradient of the interpolation
             * function \f$ \left( \phi_{,i} \right) \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the gradient of the interpolation
             * function \f$ \left( \phi_{,i} \right) \f$
             * \param &dRhoDotdRho: The Jacobian of the time derivative of the density w.r.t. the density
             * \param &dUDotdU: The Jacobian of the time derivative of the displacement degree of freedom w.r.t. the
             * displacement degree of freedom
             * \param &dUDDotdU: The Jacobian of the second time derivative of the displacement degree of freedom w.r.t.
             * the displacement degree of freedom density
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             * \param &dRdRho_begin: The starting iterator of the Jacobian of the result w.r.t. the density
             *     The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's density
             * \param &dRdRho_end: The stopping iterator of the Jacobian of the result w.r.t. the density
             *     The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's density
             * \param &dRdU_begin: The starting iterator of the Jacobian of the result w.r.t. the displacement degree of
             * freedom The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's
             * displacement
             * \param &dRdU_end: The stopping iterator of the Jacobian of the result w.r.t. the displacement degree of
             * freedom The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's
             * displacement
             * \param &dRdB_begin: The starting iterator of the Jacobian of the result w.r.t. the body force
             *     The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's body force
             * \param &dRdB_end: The stopping iterator of the Jacobian of the result w.r.t. the body force
             *     The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's body force
             * \param &dRdCauchy_begin: The starting iterator of the Jacobian of the result w.r.t. the Cauchy stress
             *     The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's Cauchy
             * stress
             * \param &dRdCauchy_end: The stopping iterator of the Jacobian of the result w.r.t. the Cauchy stress
             *     The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's Cauchy
             * stress
             * \param &dRdVolumeFraction_begin: The starting iterator of the Jacobian of the result w.r.t. the volume
             * fraction The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's volume
             * fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the Jacobian of the result w.r.t. the volume
             * fraction The derivatives are stored as the derivative of each phase's residual w.r.t. the phase's volume
             * fraction
             * \param &dRdUMesh_begin: The starting iterator of the Jacobian of the result w.r.t. the mesh displacement
             * (includes the terms due to the volume change) The derivatives are stored as the derivative of each
             * phase's residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the Jacobian of the result w.r.t. the mesh displacement
             * (includes the terms due to the volume change) The derivatives are stored as the derivative of each
             * phase's residual w.r.t. the mesh displacement
             */

            TARDIGRADE_ERROR_TOOLS_EVAL(const unsigned int nphases = (unsigned int)(density_end - density_begin))

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(density_dot_end - density_dot_begin),
                                         "The density and the time derivative of the density must have the same size");

            TARDIGRADE_ERROR_TOOLS_CHECK(
                configuration::dimension * nphases == (unsigned int)(density_gradient_end - density_gradient_begin),
                "The density and the spatial gradient of the density must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases ==
                                             (unsigned int)(velocity_end - velocity_begin),
                                         "The density and the velocity must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(
                configuration::dimension * nphases == (unsigned int)(velocity_dot_end - velocity_dot_begin),
                "The density and the time derivative of the velocity must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension * nphases ==
                                             (unsigned int)(velocity_gradient_end - velocity_gradient_begin),
                                         "The density and the velocity gradient must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases ==
                                             (unsigned int)(body_force_end - body_force_begin),
                                         "The density and the body force must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension * nphases ==
                                             (unsigned int)(cauchy_stress_end - cauchy_stress_begin),
                                         "The density and the Cauchy stress must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(nphases == (unsigned int)(volume_fraction_end - volume_fraction_begin),
                                         "The density and the volume fraction must have the same size");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases ==
                                             (unsigned int)(dRdRho_end - dRdRho_begin),
                                         "The density and dRdRho must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * configuration::dimension ==
                                             (unsigned int)(dRdU_end - dRdU_begin),
                                         "The density and dRdU must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * configuration::dimension ==
                                             (unsigned int)(dRdB_end - dRdB_begin),
                                         "The density and dRdB must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * configuration::dimension *
                                                 configuration::dimension ==
                                             (unsigned int)(dRdCauchy_end - dRdCauchy_begin),
                                         "The density and dRdCauchy must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases ==
                                             (unsigned int)(dRdVolumeFraction_end - dRdVolumeFraction_begin),
                                         "The density and dRdVolumeFraction must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * nphases * configuration::dimension ==
                                             (unsigned int)(dRdUMesh_end - dRdUMesh_begin),
                                         "The density and dRdUMesh must have consistent sizes");

            for (auto rho = std::pair<unsigned int, density_iter>(0, density_begin); rho.second != density_end;
                 ++rho.first, ++rho.second) {
                computeBalanceOfLinearMomentum<configuration>(
                    *rho.second, *(density_dot_begin + rho.first),
                    density_gradient_begin + configuration::dimension * rho.first,
                    density_gradient_begin + configuration::dimension * (rho.first + 1),
                    velocity_begin + configuration::dimension * rho.first,
                    velocity_begin + configuration::dimension * (rho.first + 1),
                    velocity_dot_begin + configuration::dimension * rho.first,
                    velocity_dot_begin + configuration::dimension * (rho.first + 1),
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * rho.first,
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * (rho.first + 1),
                    body_force_begin + configuration::dimension * rho.first,
                    body_force_begin + configuration::dimension * (rho.first + 1),
                    cauchy_stress_begin + configuration::dimension * configuration::dimension * rho.first,
                    cauchy_stress_begin + configuration::dimension * configuration::dimension * (rho.first + 1),
                    *(volume_fraction_begin + rho.first), test_function, test_function_gradient_begin,
                    test_function_gradient_end, interpolation_function, interpolation_function_gradient_begin,
                    interpolation_function_gradient_end, dRhoDotdRho, dUDotdU, dUDDotdU,
                    result_begin + configuration::dimension * rho.first,
                    result_begin + configuration::dimension * (rho.first + 1),
                    dRdRho_begin + configuration::dimension * rho.first,
                    dRdRho_begin + configuration::dimension * (rho.first + 1),
                    dRdU_begin + configuration::dimension * configuration::dimension * rho.first,
                    dRdU_begin + configuration::dimension * configuration::dimension * (rho.first + 1),
                    dRdB_begin + configuration::dimension * configuration::dimension * rho.first,
                    dRdB_begin + configuration::dimension * configuration::dimension * (rho.first + 1),
                    dRdCauchy_begin +
                        configuration::dimension * configuration::dimension * configuration::dimension * rho.first,
                    dRdCauchy_begin + configuration::dimension * configuration::dimension * configuration::dimension *
                                          (rho.first + 1),
                    dRdVolumeFraction_begin + configuration::dimension * rho.first,
                    dRdVolumeFraction_begin + configuration::dimension * (rho.first + 1),
                    dRdUMesh_begin + configuration::dimension * configuration::dimension * rho.first,
                    dRdUMesh_begin + configuration::dimension * configuration::dimension * (rho.first + 1));
            }
        }

        template <class configuration, class density_iter, class density_dot_iter, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter, class body_force_iter,
                  class cauchy_stress_iter, class volume_fraction_iter, typename testFunction_type,
                  class testFunctionGradient_iter, class result_iter>
        void computeBalanceOfLinearMomentum(
            const density_iter &density_begin, const density_iter &density_end,
            const density_dot_iter &density_dot_begin, const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin, const body_force_iter &body_force_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
            const testFunction_type &test_function, const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter &test_function_gradient_end, result_iter result_begin,
            result_iter result_end) {
            /*!
             * Compute the balance of linear momentum for a multiphase material
             *
             * \param &density_begin: The starting iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_end: The stopping iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right)
             * \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the cauchy stress \f$ \left( \sigma_{ij} \right)
             * \f$
             * \param &cauchy_stress_end: The stopping iterator of the cauchy stress \f$ \left( \sigma_{ij} \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of the phase. Only applied to
             * the Cauchy stress because the density is assumed to be the apparent density
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of the phase. Only applied to
             * the Cauchy stress because the density is assumed to be the apparent density i.e., the mass of the phase
             * per unit volume.
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$
             * \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(density_dot_end - density_dot_begin),
                                         "The density and the time derivative of the density must have the same size");

            TARDIGRADE_ERROR_TOOLS_CHECK(
                configuration::dimension * (unsigned int)(density_end - density_begin) ==
                    (unsigned int)(density_gradient_end - density_gradient_begin),
                "The density and the spatial gradient of the density must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(velocity_end - velocity_begin),
                                         "The density and the velocity must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(
                configuration::dimension * (unsigned int)(density_end - density_begin) ==
                    (unsigned int)(velocity_dot_end - velocity_dot_begin),
                "The density and the time derivative of the velocity must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension *
                                                 (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(velocity_gradient_end - velocity_gradient_begin),
                                         "The density and the velocity gradient must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(body_force_end - body_force_begin),
                                         "The density and the body force must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension *
                                                 (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(cauchy_stress_end - cauchy_stress_begin),
                                         "The density and the Cauchy stress must have consistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(volume_fraction_end - volume_fraction_begin),
                                         "The density and the volume fraction must have the same size");

            for (auto rho = std::pair<unsigned int, density_iter>(0, density_begin); rho.second != density_end;
                 ++rho.first, ++rho.second) {
                computeBalanceOfLinearMomentum<configuration>(
                    *rho.second, *(density_dot_begin + rho.first),
                    density_gradient_begin + configuration::dimension * rho.first,
                    density_gradient_begin + configuration::dimension * (rho.first + 1),
                    velocity_begin + configuration::dimension * rho.first,
                    velocity_begin + configuration::dimension * (rho.first + 1),
                    velocity_dot_begin + configuration::dimension * rho.first,
                    velocity_dot_begin + configuration::dimension * (rho.first + 1),
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * rho.first,
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * (rho.first + 1),
                    body_force_begin + configuration::dimension * rho.first,
                    body_force_begin + configuration::dimension * (rho.first + 1),
                    cauchy_stress_begin + configuration::dimension * configuration::dimension * rho.first,
                    cauchy_stress_begin + configuration::dimension * configuration::dimension * (rho.first + 1),
                    *(volume_fraction_begin + rho.first), test_function, test_function_gradient_begin,
                    test_function_gradient_end, result_begin + configuration::dimension * rho.first,
                    result_begin + configuration::dimension * (rho.first + 1));
            }
        }

        template <class configuration, typename density_type, typename density_dot_type, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter, class body_force_iter,
                  class result_iter, class dRdRho_iter, class dRdRhoDot_iter, class dRdGradRho_iter, class dRdV_iter,
                  class dRdVDot_iter, class dRdGradV_iter, class dRdB_iter>
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin, const body_force_iter &body_force_end, result_iter result_begin,
            result_iter result_end, dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end, dRdRhoDot_iter dRdRhoDot_begin,
            dRdRhoDot_iter dRdRhoDot_end, dRdGradRho_iter dRdGradRho_begin, dRdGradRho_iter dRdGradRho_end,
            dRdV_iter dRdV_begin, dRdV_iter dRdV_end, dRdVDot_iter dRdVDot_begin, dRdVDot_iter dRdVDot_end,
            dRdGradV_iter dRdGradV_begin, dRdGradV_iter dRdGradV_end, dRdB_iter dRdB_begin, dRdB_iter dRdB_end) {
            /*!
             * Compute the non-divergence part of the balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho
             * v_i \right) = -f_i \f$
             *
             * and the non-divergence part is
             *
             * \f$ -\left( \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) +
             * f_i = 0 \f$
             *
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. The expression
             * we implement is therefore
             *
             * \f$ -\rho_{,j} v_j v_i - \rho v_{i,j} v_j - \rho v_i v_{j,j} + \rho b_i - \frac{\partial \rho}{\partial
             * t} v_i - \rho \frac{\partial v_i}{\partial t} \f$
             *
             * \param &density: The mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot: The partial time derivative of the density \f$ \left( \frac{\partial}{\partial t}
             * \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right)
             * \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             * \param &dRdRho_begin: The starting iterator of the derivative of the result w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the derivative of the result w.r.t. the density
             * \param &dRdRhoDot_begin: The starting iterator of the derivative of the result w.r.t. the partial
             * temporal derivative of the density
             * \param &dRdRhoDot_end: The stopping iterator of the derivative of the result w.r.t. the partial temporal
             * derivative of the density
             * \param &dRdGradRho_begin: The starting iterator of the derivative of the result w.r.t. the spatial
             * derivative of the density
             * \param &dRdGradRho_end: The stopping iterator of the derivative of the result w.r.t. the spatial
             * derivative of the density
             * \param &dRdV_begin: The starting iterator of the derivative of the result w.r.t. the velocity
             * \param &dRdV_end: The stopping iterator of the derivative of the result w.r.t. the velocity
             * \param &dRdVDot_begin: The starting iterator of the derivative of the result w.r.t. the partial temporal
             * derivative of the velocity
             * \param &dRdVDot_end: The stopping iterator of the derivative of the result w.r.t. the partial temporal
             * derivative of the velocity
             * \param &dRdGradV_begin: The starting iterator of the derivative of the result w.r.t. the spatial
             * derivative of the velocity
             * \param &dRdGradV_end: The stopping iterator of the derivative of the result w.r.t. the spatial derivative
             * of the velocity
             * \param &dRdB_begin: The starting iterator of the derivative of the result w.r.t. the body force
             * \param &dRdB_end: The stopping iterator of the derivative of the result w.r.t. the body force
             */

            // Compute the inner product of the velocity and the density gradient
            using density_gradient_type = typename std::iterator_traits<density_gradient_iter>::value_type;
            using velocity_type         = typename std::iterator_traits<velocity_iter>::value_type;
            using grad_rho_dot_type = decltype(std::declval<density_gradient_type &>() * std::declval<velocity_type>());
            using velocity_gradient_type = typename std::iterator_traits<velocity_gradient_iter>::value_type;

            grad_rho_dot_type grad_rho_dot_v =
                std::inner_product(density_gradient_begin, density_gradient_end, velocity_begin, 0.);

            // Compute the trace of the velocity gradient
            velocity_gradient_type trace_velocity_gradient = 0;

            for (unsigned int i = 0; i < configuration::dimension; ++i) {
                trace_velocity_gradient += *(velocity_gradient_begin + configuration::dimension * i + i);
            }

            // Initialize the gradients
            std::fill(dRdGradRho_begin, dRdGradRho_end, 0.);
            std::fill(dRdV_begin, dRdV_end, 0.);
            std::fill(dRdVDot_begin, dRdVDot_end, 0.);
            std::fill(dRdGradV_begin, dRdGradV_end, 0.);
            std::fill(dRdB_begin, dRdB_end, 0.);

            // Compute the balance of linear momentum

            for (unsigned int i = 0; i < configuration::dimension; ++i) {
                *(result_begin + i) = density * ((*(body_force_begin + i)) - (*(velocity_dot_begin + i)) -
                                                 trace_velocity_gradient * (*(velocity_begin + i))) -
                                      grad_rho_dot_v * (*(velocity_begin + i)) - density_dot * (*(velocity_begin + i));

                *(dRdRho_begin + i) = ((*(body_force_begin + i)) - (*(velocity_dot_begin + i)) -
                                       trace_velocity_gradient * (*(velocity_begin + i)));

                *(dRdRhoDot_begin + i) = -(*(velocity_begin + i));

                *(dRdV_begin + configuration::dimension * i + i) =
                    -density * trace_velocity_gradient - grad_rho_dot_v - density_dot;

                *(dRdVDot_begin + configuration::dimension * i + i) = -density;

                *(dRdB_begin + configuration::dimension * i + i) = density;

                for (unsigned int j = 0; j < configuration::dimension; ++j) {
                    *(result_begin + i) -= density * (*(velocity_gradient_begin + configuration::dimension * i + j)) *
                                           (*(velocity_begin + j));

                    *(dRdRho_begin + i) -=
                        (*(velocity_gradient_begin + configuration::dimension * i + j)) * (*(velocity_begin + j));

                    *(dRdGradRho_begin + configuration::dimension * i + j) =
                        -(*(velocity_begin + i)) * (*(velocity_begin + j));

                    *(dRdV_begin + configuration::dimension * i + j) -=
                        (*(density_gradient_begin + j)) * (*(velocity_begin + i)) +
                        density * (*(velocity_gradient_begin + configuration::dimension * i + j));

                    *(dRdGradV_begin + configuration::dimension * configuration::dimension * i +
                      configuration::dimension * j + j) -= density * (*(velocity_begin + i));

                    *(dRdGradV_begin + configuration::dimension * configuration::dimension * i +
                      configuration::dimension * i + j) -= density * (*(velocity_begin + j));
                }
            }
        }

        template <class configuration, class density_iter, class density_dot_iter, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter, class body_force_iter,
                  class result_iter>
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_iter &density_begin, const density_iter &density_end,
            const density_dot_iter &density_dot_begin, const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin, const body_force_iter &body_force_end, result_iter result_begin,
            result_iter result_end) {
            /*!
             * Compute the non-divergence part of the multiphase balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho
             * v_i \right) = -f_i \f$
             *
             * and the non-divergence part is
             *
             * \f$ -\left( \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) +
             * f_i = 0 \f$
             *
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. The expression
             * we implement is therefore
             *
             * \f$ -\rho_{,j} v_j v_i - \rho v_{i,j} v_j - \rho v_i v_{j,j} + \rho b_i - \frac{\partial \rho}{\partial
             * t} v_i - \rho \frac{\partial v_i}{\partial t} \f$
             *
             * \param &density_begin: The starting iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_end: The stopping iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right)
             * \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(density_dot_end - density_dot_begin),
                                         "The density and density dot vectors are not the same size");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(density_gradient_end - density_gradient_begin),
                                         "The density and density gradient vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(velocity_end - velocity_begin),
                                         "The density and velocity vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(velocity_dot_end - velocity_dot_begin),
                                         "The density and velocity dot vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension *
                                                 (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(velocity_gradient_end - velocity_gradient_begin),
                                         "The density and velocity gradient vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(body_force_end - body_force_begin),
                                         "The density and body force vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(result_end - result_begin),
                                         "The density and result vectors are of inconsistent sizes");

            for (auto rho = density_begin; rho != density_end; ++rho) {
                unsigned int phase = (unsigned int)(rho - density_begin);

                computeBalanceOfLinearMomentumNonDivergence<configuration>(
                    *(density_begin + phase), *(density_dot_begin + phase),
                    density_gradient_begin + configuration::dimension * phase,
                    density_gradient_begin + configuration::dimension * (phase + 1),
                    velocity_begin + configuration::dimension * phase,
                    velocity_begin + configuration::dimension * (phase + 1),
                    velocity_dot_begin + configuration::dimension * phase,
                    velocity_dot_begin + configuration::dimension * (phase + 1),
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * phase,
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    body_force_begin + configuration::dimension * phase,
                    body_force_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    result_begin + configuration::dimension * phase,
                    result_begin + configuration::dimension * (phase + 1));
            }
        }

        template <class configuration, class density_iter, class density_dot_iter, class density_gradient_iter,
                  class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter, class body_force_iter,
                  class result_iter, class dRdRho_iter, class dRdRhoDot_iter, class dRdGradRho_iter, class dRdV_iter,
                  class dRdVDot_iter, class dRdGradV_iter, class dRdB_iter>
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_iter &density_begin, const density_iter &density_end,
            const density_dot_iter &density_dot_begin, const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin, const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin, const body_force_iter &body_force_end, result_iter result_begin,
            result_iter result_end, dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end, dRdRhoDot_iter dRdRhoDot_begin,
            dRdRhoDot_iter dRdRhoDot_end, dRdGradRho_iter dRdGradRho_begin, dRdGradRho_iter dRdGradRho_end,
            dRdV_iter dRdV_begin, dRdV_iter dRdV_end, dRdVDot_iter dRdVDot_begin, dRdVDot_iter dRdVDot_end,
            dRdGradV_iter dRdGradV_begin, dRdGradV_iter dRdGradV_end, dRdB_iter dRdB_begin, dRdB_iter dRdB_end) {
            /*!
             * Compute the non-divergence part of the balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho
             * v_i \right) = -f_i \f$
             *
             * and the non-divergence part is
             *
             * \f$ -\left( \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) +
             * f_i = 0 \f$
             *
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. The expression
             * we implement is therefore
             *
             * \f$ -\rho_{,j} v_j v_i - \rho v_{i,j} v_j - \rho v_i v_{j,j} + \rho b_i - \frac{\partial \rho}{\partial
             * t} v_i - \rho \frac{\partial v_i}{\partial t} \f$
             *
             * \param &density_begin: The starting iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_end: The stopping iterator of the mass density per unit current volume \f$ \left( \rho
             * \right) \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \left(
             * \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left(
             * \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$
             * \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left(
             * \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left(
             * v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right)
             * \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             * \param &dRdRho_begin: The starting iterator of the derivative of the result w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the derivative of the result w.r.t. the density
             * \param &dRdRhoDot_begin: The starting iterator of the derivative of the result w.r.t. the partial
             * temporal derivative of the density
             * \param &dRdRhoDot_end: The stopping iterator of the derivative of the result w.r.t. the partial temporal
             * derivative of the density
             * \param &dRdGradRho_begin: The starting iterator of the derivative of the result w.r.t. the spatial
             * derivative of the density
             * \param &dRdGradRho_end: The stopping iterator of the derivative of the result w.r.t. the spatial
             * derivative of the density
             * \param &dRdV_begin: The starting iterator of the derivative of the result w.r.t. the velocity
             * \param &dRdV_end: The stopping iterator of the derivative of the result w.r.t. the velocity
             * \param &dRdVDot_begin: The starting iterator of the derivative of the result w.r.t. the partial temporal
             * derivative of the velocity
             * \param &dRdVDot_end: The stopping iterator of the derivative of the result w.r.t. the partial temporal
             * derivative of the velocity
             * \param &dRdGradV_begin: The starting iterator of the derivative of the result w.r.t. the spatial
             * derivative of the velocity
             * \param &dRdGradV_end: The stopping iterator of the derivative of the result w.r.t. the spatial derivative
             * of the velocity
             * \param &dRdB_begin: The starting iterator of the derivative of the result w.r.t. the body force
             * \param &dRdB_end: The stopping iterator of the derivative of the result w.r.t. the body force
             */

            TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(density_dot_end - density_dot_begin),
                                         "The density and density dot vectors are not the same size");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(density_gradient_end - density_gradient_begin),
                                         "The density and density gradient vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(velocity_end - velocity_begin),
                                         "The density and velocity vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(velocity_dot_end - velocity_dot_begin),
                                         "The density and velocity dot vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension *
                                                 (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(velocity_gradient_end - velocity_gradient_begin),
                                         "The density and velocity gradient vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(body_force_end - body_force_begin),
                                         "The density and body force vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(result_end - result_begin),
                                         "The density and result vectors are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(dRdRho_end - dRdRho_begin),
                                         "The density and dRdRho are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(dRdRhoDot_end - dRdRhoDot_begin),
                                         "The density and dRdRhoDot are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension *
                                                 (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(dRdGradRho_end - dRdGradRho_begin),
                                         "The density and dRdGradRho are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension *
                                                 (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(dRdV_end - dRdV_begin),
                                         "The density and dRdV are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension *
                                                 configuration::dimension *
                                                 (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(dRdGradV_end - dRdGradV_begin),
                                         "The density and dRdGradV are of inconsistent sizes");

            TARDIGRADE_ERROR_TOOLS_CHECK(configuration::dimension * configuration::dimension *
                                                 (unsigned int)(density_end - density_begin) ==
                                             (unsigned int)(dRdB_end - dRdB_begin),
                                         "The density and dRdB are of inconsistent sizes");

            for (auto rho = density_begin; rho != density_end; ++rho) {
                unsigned int phase = (unsigned int)(rho - density_begin);

                computeBalanceOfLinearMomentumNonDivergence<configuration>(
                    *(density_begin + phase), *(density_dot_begin + phase),
                    density_gradient_begin + configuration::dimension * phase,
                    density_gradient_begin + configuration::dimension * (phase + 1),
                    velocity_begin + configuration::dimension * phase,
                    velocity_begin + configuration::dimension * (phase + 1),
                    velocity_dot_begin + configuration::dimension * phase,
                    velocity_dot_begin + configuration::dimension * (phase + 1),
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * phase,
                    velocity_gradient_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    body_force_begin + configuration::dimension * phase,
                    body_force_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    result_begin + configuration::dimension * phase,
                    result_begin + configuration::dimension * (phase + 1),
                    dRdRho_begin + configuration::dimension * phase,
                    dRdRho_begin + configuration::dimension * (phase + 1),
                    dRdRhoDot_begin + configuration::dimension * phase,
                    dRdRhoDot_begin + configuration::dimension * (phase + 1),
                    dRdGradRho_begin + configuration::dimension * configuration::dimension * phase,
                    dRdGradRho_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    dRdV_begin + configuration::dimension * configuration::dimension * phase,
                    dRdV_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    dRdVDot_begin + configuration::dimension * configuration::dimension * phase,
                    dRdVDot_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    dRdGradV_begin +
                        configuration::dimension * configuration::dimension * configuration::dimension * phase,
                    dRdGradV_begin +
                        configuration::dimension * configuration::dimension * configuration::dimension * (phase + 1),
                    dRdB_begin + configuration::dimension * configuration::dimension * phase,
                    dRdB_begin + configuration::dimension * configuration::dimension * (phase + 1));
            }
        }

        template <class configuration, class testFunctionGradient_iter, class cauchyStress_iter,
                  typename volumeFraction_type, class result_iter>
        void computeBalanceOfLinearMomentumDivergence(const testFunctionGradient_iter &grad_test_function_begin,
                                                      const testFunctionGradient_iter &grad_test_function_end,
                                                      const cauchyStress_iter         &cauchy_stress_begin,
                                                      const cauchyStress_iter         &cauchy_stress_end,
                                                      const volumeFraction_type       &volume_fraction,
                                                      result_iter result_begin, result_iter result_end) {
            /*!
             * Compute the divergence part of the balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho
             * v_i \right) = -f_i \f$
             *
             * and the divergence part is
             *
             * \f$ \left( \sigma_{ji} \right)_{,j} \f$
             *
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. In the context
             * of a variationally-based solution technique we multiply this by a test function \f$ \psi \f$ and write
             *
             * \f$ \left( \psi \sigma_{ji} \right)_{,j} - \psi_{,j} \sigma_{ji} \f$. The divergence term is defined as
             * the second term. We also introduce a volume fraction
             * \f$ \phi \f$ to prepare for multi-phase continuum such that we will actually compute
             *
             * \f$ -\psi_{,j} \phi \bar{\sigma}_{ji} \f$
             *
             * where we made the definition that the effective Cauchy stress \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}}
             * \f$ where \f$\bar{\bf{\sigma}} \f$ is the true Cauchy stress of the phase.
             *
             * \param &grad_test_function_begin: The starting iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &grad_test_function_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}}
             * \right) \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}}
             * \right) \f$
             * \param &volume_fraction: The volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting iterator of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the divergence part of the balance of linear momentum
             */

            std::fill(result_begin, result_end, 0.);

            for (unsigned int i = 0; i < configuration::dimension; ++i) {
                for (unsigned int j = 0; j < configuration::dimension; ++j) {
                    *(result_begin + i) -= (*(grad_test_function_begin + j)) * volume_fraction *
                                           (*(cauchy_stress_begin + configuration::dimension * j + i));
                }
            }
        }

        template <class configuration, class testFunctionGradient_iter, class cauchyStress_iter,
                  typename volumeFraction_type, class result_iter, class dRdGradPsi_iter, class dRdCauchy_iter,
                  class dRdVolumeFraction_iter>
        void computeBalanceOfLinearMomentumDivergence(
            const testFunctionGradient_iter &grad_test_function_begin,
            const testFunctionGradient_iter &grad_test_function_end, const cauchyStress_iter &cauchy_stress_begin,
            const cauchyStress_iter &cauchy_stress_end, const volumeFraction_type &volume_fraction,
            result_iter result_begin, result_iter result_end, dRdGradPsi_iter dRdGradPsi_begin,
            dRdGradPsi_iter dRdGradPsi_end, dRdCauchy_iter dRdCauchy_begin, dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end) {
            /*!
             * Compute the divergence part of the balance of linear momentum and its Jacobians where the full equation
             * is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho
             * v_i \right) = -f_i \f$
             *
             * and the divergence part is
             *
             * \f$ \left( \sigma_{ji} \right)_{,j} \f$
             *
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. In the context
             * of a variationally-based solution technique we multiply this by a test function \f$ \psi \f$ and write
             *
             * \f$ \left( \psi \sigma_{ji} \right)_{,j} - \psi_{,j} \sigma_{ji} \f$. The divergence term is defined as
             * the second term. We also introduce a volume fraction
             * \f$ \phi \f$ to prepare for multi-phase continuum such that we will actually compute
             *
             * \f$ -\psi_{,j} \phi \bar{\sigma}_{ji} \f$
             *
             * where we made the definition that the effective Cauchy stress \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}}
             * \f$ where \f$\bar{\bf{\sigma}} \f$ is the true Cauchy stress of the phase.
             *
             * \param &grad_test_function_begin: The starting iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &grad_test_function_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}}
             * \right) \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}}
             * \right) \f$
             * \param &volume_fraction: The volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting iterator of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the divergence part of the balance of linear momentum
             * \param &dRdGradPsi_begin: The starting iterator of the derivative of the result w.r.t. the test function
             * spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdGradPsi_end: The stopping iterator of the derivative of the result w.r.t. the test function
             * spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdCauchy_begin: The starting iterator of the derivative of the result w.r.t. the Cauchy stress
             * \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdCauchy_end: The stopping iterator of the derivative of the result w.r.t. the Cauchy stress \f$
             * \left( \bf{\sigma} \right) \f$
             * \param &dRdVolumeFraction_begin: The starting iterator of the derivative of the result w.r.t. the volume
             * fraction \f$ \left( \bf{\phi} \right) \f$
             * \param &dRdVolumeFraction_end: The stopping iterator of the derivative of the result w.r.t. the volume
             * fraction \f$ \left( \bf{\phi} \right) \f$
             */

            std::fill(result_begin, result_end, 0.);
            std::fill(dRdGradPsi_begin, dRdGradPsi_end, 0.);
            std::fill(dRdCauchy_begin, dRdCauchy_end, 0.);
            std::fill(dRdVolumeFraction_begin, dRdVolumeFraction_end, 0.);

            for (unsigned int i = 0; i < configuration::dimension; ++i) {
                for (unsigned int j = 0; j < configuration::dimension; ++j) {
                    *(result_begin + i) -= (*(grad_test_function_begin + j)) * volume_fraction *
                                           (*(cauchy_stress_begin + configuration::dimension * j + i));

                    *(dRdGradPsi_begin + configuration::dimension * i + j) -=
                        volume_fraction * (*(cauchy_stress_begin + configuration::dimension * j + i));

                    *(dRdCauchy_begin + configuration::dimension * configuration::dimension * i +
                      configuration::dimension * j + i) -= (*(grad_test_function_begin + j)) * volume_fraction;

                    *(dRdVolumeFraction_begin + i) -=
                        (*(grad_test_function_begin + j)) * (*(cauchy_stress_begin + configuration::dimension * j + i));
                }
            }
        }

        template <class configuration, class testFunctionGradient_iter, class cauchyStress_iter,
                  class volumeFraction_iter, class result_iter>
        void computeBalanceOfLinearMomentumDivergence(const testFunctionGradient_iter &grad_test_function_begin,
                                                      const testFunctionGradient_iter &grad_test_function_end,
                                                      const cauchyStress_iter         &cauchy_stress_begin,
                                                      const cauchyStress_iter         &cauchy_stress_end,
                                                      const volumeFraction_iter       &volume_fraction_begin,
                                                      const volumeFraction_iter       &volume_fraction_end,
                                                      result_iter result_begin, result_iter result_end) {
            /*!
             * Compute the divergence part of the balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho
             * v_i \right) = -f_i \f$
             *
             * and the divergence part is
             *
             * \f$ \left( \sigma_{ji} \right)_{,j} \f$
             *
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. In the context
             * of a variationally-based solution technique we multiply this by a test function \f$ \psi \f$ and write
             *
             * \f$ \left( \psi \sigma_{ji} \right)_{,j} - \psi_{,j} \sigma_{ji} \f$. The divergence term is defined as
             * the second term. We also introduce a volume fraction
             * \f$ \phi \f$ to prepare for multi-phase continuum such that we will actually compute
             *
             * \f$ -\psi_{,j} \phi \bar{\sigma}_{ji} \f$
             *
             * where we made the definition that the effective Cauchy stress \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}}
             * \f$ where \f$\bar{\bf{\sigma}} \f$ is the true Cauchy stress of the phase.
             *
             * \param &grad_test_function_begin: The starting iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &grad_test_function_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}}
             * \right) \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}}
             * \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of the current phase \f$
             * \left( \phi \right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of the current phase \f$ \left(
             * \phi \right) \f$
             * \param &result_begin: The starting iterator of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the divergence part of the balance of linear momentum
             */

            for (auto volume_fraction = volume_fraction_begin; volume_fraction != volume_fraction_end;
                 ++volume_fraction) {
                unsigned int phase = (unsigned int)(volume_fraction - volume_fraction_begin);

                computeBalanceOfLinearMomentumDivergence<configuration>(
                    grad_test_function_begin, grad_test_function_end,
                    cauchy_stress_begin + configuration::dimension * configuration::dimension * phase,
                    cauchy_stress_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    *volume_fraction, result_begin + configuration::dimension * phase,
                    result_begin + configuration::dimension * (phase + 1));
            }
        }

        template <class configuration, class testFunctionGradient_iter, class cauchyStress_iter,
                  class volumeFraction_iter, class result_iter, class dRdGradPsi_iter, class dRdCauchy_iter,
                  class dRdVolumeFraction_iter>
        void computeBalanceOfLinearMomentumDivergence(
            const testFunctionGradient_iter &grad_test_function_begin,
            const testFunctionGradient_iter &grad_test_function_end, const cauchyStress_iter &cauchy_stress_begin,
            const cauchyStress_iter &cauchy_stress_end, const volumeFraction_iter &volume_fraction_begin,
            const volumeFraction_iter &volume_fraction_end, result_iter result_begin, result_iter result_end,
            dRdGradPsi_iter dRdGradPsi_begin, dRdGradPsi_iter dRdGradPsi_end, dRdCauchy_iter dRdCauchy_begin,
            dRdCauchy_iter dRdCauchy_end, dRdVolumeFraction_iter dRdVolumeFraction_begin,
            dRdVolumeFraction_iter dRdVolumeFraction_end) {
            /*!
             * Compute the divergence part of the balance of linear momentum and its Jacobians where the full equation
             * is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho
             * v_i \right) = -f_i \f$
             *
             * and the divergence part is
             *
             * \f$ \left( \sigma_{ji} \right)_{,j} \f$
             *
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. In the context
             * of a variationally-based solution technique we multiply this by a test function \f$ \psi \f$ and write
             *
             * \f$ \left( \psi \sigma_{ji} \right)_{,j} - \psi_{,j} \sigma_{ji} \f$. The divergence term is defined as
             * the second term. We also introduce a volume fraction
             * \f$ \phi \f$ to prepare for multi-phase continuum such that we will actually compute
             *
             * \f$ -\psi_{,j} \phi \bar{\sigma}_{ji} \f$
             *
             * where we made the definition that the effective Cauchy stress \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}}
             * \f$ where \f$\bar{\bf{\sigma}} \f$ is the true Cauchy stress of the phase.
             *
             * \param &grad_test_function_begin: The starting iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &grad_test_function_end: The stopping iterator of the gradient of the test function \f$ \left(
             * \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}}
             * \right) \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}}
             * \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of the current phase \f$
             * \left( \phi \right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of the current phase \f$ \left(
             * \phi \right) \f$
             * \param &result_begin: The starting iterator of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the divergence part of the balance of linear momentum
             * \param &dRdGradPsi_begin: The starting iterator of the derivative of the result w.r.t. the test function
             * spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdGradPsi_end: The stopping iterator of the derivative of the result w.r.t. the test function
             * spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdCauchy_begin: The starting iterator of the derivative of the result w.r.t. the Cauchy stress
             * \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdCauchy_end: The stopping iterator of the derivative of the result w.r.t. the Cauchy stress \f$
             * \left( \bf{\sigma} \right) \f$
             * \param &dRdVolumeFraction_begin: The starting iterator of the derivative of the result w.r.t. the volume
             * fraction \f$ \left( \bf{\phi} \right) \f$
             * \param &dRdVolumeFraction_end: The stopping iterator of the derivative of the result w.r.t. the volume
             * fraction \f$ \left( \bf{\phi} \right) \f$
             */

            for (auto volume_fraction = volume_fraction_begin; volume_fraction != volume_fraction_end;
                 ++volume_fraction) {
                unsigned int phase = (unsigned int)(volume_fraction - volume_fraction_begin);

                computeBalanceOfLinearMomentumDivergence<configuration>(
                    grad_test_function_begin, grad_test_function_end,
                    cauchy_stress_begin + configuration::dimension * configuration::dimension * phase,
                    cauchy_stress_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    *volume_fraction, result_begin + configuration::dimension * phase,
                    result_begin + configuration::dimension * (phase + 1),
                    dRdGradPsi_begin + configuration::dimension * configuration::dimension * phase,
                    dRdGradPsi_begin + configuration::dimension * configuration::dimension * (phase + 1),
                    dRdCauchy_begin +
                        configuration::dimension * configuration::dimension * configuration::dimension * phase,
                    dRdCauchy_begin +
                        configuration::dimension * configuration::dimension * configuration::dimension * (phase + 1),
                    dRdVolumeFraction_begin + configuration::dimension * phase,
                    dRdVolumeFraction_begin + configuration::dimension * (phase + 1));
            }
        }

    }  // namespace balanceOfLinearMomentum

}  // namespace tardigradeBalanceEquations
