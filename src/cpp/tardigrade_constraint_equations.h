/**
  ******************************************************************************
  * \file tardigrade_constraint_equations.h
  ******************************************************************************
  * The header file for the equations associated with the balance of mass
  ******************************************************************************
  */

#ifndef TARDIGRADE_CONSTRAINT_EQUATIONS_H
#define TARDIGRADE_CONSTRAINT_EQUATIONS_H

#define USE_EIGEN
#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace constraintEquations{

        template<
            int predicted_internal_energy_index,
            typename internal_energy_type,
            class material_response_iter,
            typename test_function_type,
            typename result_type
        >
        inline void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const test_function_type &test_function,
            result_type &result
        );

        template<
            int predicted_internal_energy_index,
            class internal_energy_iter,
            class material_response_iter,
            typename test_function_type,
            class result_iter
        >
        inline void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin,
            const internal_energy_iter &internal_energy_end,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const test_function_type &test_function,
            result_iter result_begin,
            result_iter result_end
        );

        template<
            int predicted_internal_energy_index,
            typename internal_energy_type,
            typename density_type,
            class material_response_iter,
            typename test_function_type,
            typename result_type
        >
        inline void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy,
            const density_type &density,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const test_function_type &test_function,
            result_type &result
        );

        template<
            int predicted_internal_energy_index,
            class internal_energy_iter,
            class density_iter,
            class material_response_iter,
            typename test_function_type,
            class result_iter
        >
        inline void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin,
            const internal_energy_iter &internal_energy_end,
            const density_iter &density_begin,
            const density_iter &density_end,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const test_function_type &test_function,
            result_iter result_begin,
            result_iter result_end
        );

        template<
            int material_response_dim, int predicted_internal_energy_index,
            int material_response_num_dof,
            typename internal_energy_type,
            class material_response_iter,
            class material_response_jacobian_iter,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dUDotdU_type,
            typename result_type,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter,
            class dRdTheta_iter, class dRdE_iter, class dRdVF_iter,
            class dRdZ_iter, class dRdUMesh_iter,
            int density_index         = 0,
            int displacement_index    = 1,
            int velocity_index        = 4,
            int temperature_index     = 7,
            int internal_energy_index = 8,
            int volume_fraction_index = 9,
            int additional_dof_index  = 10
        >
        inline void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU,
            const unsigned int &phase,
            result_type &result,
            dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end,
            dRdE_iter     dRdE_begin,     dRdE_iter     dRdE_end,
            dRdVF_iter    dRdVF_begin,    dRdVF_iter    dRdVF_end,
            dRdZ_iter     dRdZ_begin,     dRdZ_iter     dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        );

        template<
            int material_response_dim, int predicted_internal_energy_index,
            int material_response_num_dof,
            class internal_energy_iter,
            class material_response_iter,
            class material_response_jacobian_iter,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dUDotdU_type,
            class result_iter,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter,
            class dRdTheta_iter, class dRdE_iter, class dRdVF_iter,
            class dRdZ_iter, class dRdUMesh_iter,
            int density_index         = 0,
            int displacement_index    = 1,
            int velocity_index        = 4,
            int temperature_index     = 7,
            int internal_energy_index = 8,
            int volume_fraction_index = 9,
            int additional_dof_index  = 10
        >
        inline void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin,
            const internal_energy_iter &internal_energy_end,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU,
            result_iter result_begin,
            result_iter result_end,
            dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end,
            dRdE_iter     dRdE_begin,     dRdE_iter     dRdE_end,
            dRdVF_iter    dRdVF_begin,    dRdVF_iter    dRdVF_end,
            dRdZ_iter     dRdZ_begin,     dRdZ_iter     dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        );

        template<
            int material_response_dim, int predicted_internal_energy_index,
            int material_response_num_dof,
            typename internal_energy_type,
            typename density_type,
            class material_response_iter,
            class material_response_jacobian_iter,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dUDotdU_type,
            typename result_type,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter,
            class dRdTheta_iter, class dRdE_iter, class dRdVF_iter,
            class dRdZ_iter, class dRdUMesh_iter,
            int density_index         = 0,
            int displacement_index    = 1,
            int velocity_index        = 4,
            int temperature_index     = 7,
            int internal_energy_index = 8,
            int volume_fraction_index = 9,
            int additional_dof_index  = 10
        >
        inline void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy,
            const density_type &density,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU,
            const unsigned int &phase,
            result_type &result,
            dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end,
            dRdE_iter     dRdE_begin,     dRdE_iter     dRdE_end,
            dRdVF_iter    dRdVF_begin,    dRdVF_iter    dRdVF_end,
            dRdZ_iter     dRdZ_begin,     dRdZ_iter     dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        );

        template<
            int material_response_dim, int predicted_internal_energy_index,
            int material_response_num_dof,
            class internal_energy_iter,
            class density_iter,
            class material_response_iter,
            class material_response_jacobian_iter,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dUDotdU_type,
            class result_iter,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter,
            class dRdTheta_iter, class dRdE_iter, class dRdVF_iter,
            class dRdZ_iter, class dRdUMesh_iter,
            int density_index         = 0,
            int displacement_index    = 1,
            int velocity_index        = 4,
            int temperature_index     = 7,
            int internal_energy_index = 8,
            int volume_fraction_index = 9,
            int additional_dof_index  = 10
        >
        inline void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin,
            const internal_energy_iter &internal_energy_end,
            const density_iter &density_begin,
            const density_iter &density_end,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU,
            result_iter result_begin,
            result_iter result_end,
            dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end,
            dRdE_iter     dRdE_begin,     dRdE_iter     dRdE_end,
            dRdVF_iter    dRdVF_begin,    dRdVF_iter    dRdVF_end,
            dRdZ_iter     dRdZ_begin,     dRdZ_iter     dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        );

        template<
            int cauchy_stress_index,
            int internal_energy_index,
            int mass_change_index,
            int body_force_index,
            int interphasic_force_index,
            int heat_flux_index,
            int internal_heat_generation_index,
            int interphasic_heat_transfer_index,
            int trace_mass_change_velocity_gradient_index,
            class density_iter,
            class volume_fraction_iter,
            class material_response_iter,
            class material_response_jacobian_iter,
            class mixture_response_iter,
            class mixture_jacobian_iter
            int density_index                = 0,
            int volume_fraction_index        = 9
        >
        inline void assembleMixtureMaterialResponse(
            const density_iter &density_begin,
            const density_iter &density_end,
            const volume_fraction_iter &volume_fraction_begin,
            const volume_fraction_iter &volume_fraction_end,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const mixture_response_iter &mixture_response_begin,
            const mixture_response_iter &mixture_response_end,
            const mixture_jacobian_iter &mixture_jacobian_begin,
            const mixture_jacobian_iter &mixture_jacobian_end
        );

    }

}

#include "tardigrade_constraint_equations.cpp"

#endif
