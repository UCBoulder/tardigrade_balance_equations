/**
  ******************************************************************************
  * \file tardigrade_balance_surface_growth.h
  ******************************************************************************
  * The header file for the equations associated with calculating the surface
  * growth velocity
  ******************************************************************************
  */

#ifndef TARDIGRADE_SURFACE_GROWTH_BALANCE_H
#define TARDIGRADE_SURFACE_GROWTH_BALANCE_H

#include<array>

#define USE_EIGEN
#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace surfaceGrowth{

        template<
            class surfaceGrowthVelocity_iter,
            typename lagrangeMultiplier_type,
            class lagrangeMultiplierGradient_iter,
            typename testFunction_type,
            class testFunctionGradient_iter,
            class result_iter
        >
        void computeSurfaceGrowthBalance(
            const surfaceGrowthVelocity_iter &surfaceGrowthVelocity_begin, const surfaceGrowthVelocity_iter &surfaceGrowthVelocity_end,
            const lagrangeMultiplier_type &lagrangeMultiplier,
            const lagrangeMultiplierGradient_iter &lagrangeMultiplierGradient_begin, const lagrangeMultiplierGradient_iter &lagrangeMultiplierGradient_end,
            const testFunction_type &testFunction,
            const testFunctionGradient_iter &testFunctionGradient_begin, const testFunctionGradient_iter &testFunctionGradient_end,
            result_iter result_begin, result_iter result_end
        );

        template<
            class surfaceGrowthVelocity_iter,
            typename lagrangeMultiplier_type,
            class lagrangeMultiplierGradient_iter,
            typename testFunction_type,
            class testFunctionGradient_iter,
            class result_iter,
            typename interpolationFunction_type,
            class interpolationFunctionGradient_iter,
            class dRdV_iter, class dRdL_iter, class dRdUMesh_iter
        >
        void computeSurfaceGrowthBalance(
            const surfaceGrowthVelocity_iter &surfaceGrowthVelocity_begin, const surfaceGrowthVelocity_iter &surfaceGrowthVelocity_end,
            const lagrangeMultiplier_type &lagrangeMultiplier,
            const lagrangeMultiplierGradient_iter &lagrangeMultiplierGradient_begin, const lagrangeMultiplierGradient_iter &lagrangeMultiplierGradient_end,
            const testFunction_type &testFunction,
            const testFunctionGradient_iter &testFunctionGradient_begin, const testFunctionGradient_iter &testFunctionGradient_end,
            const interpolationFunction_type &interpolationFunction,
            const interpolationFunctionGradient_iter &interpolationFunctionGradient_begin, const interpolationFunctionGradient_iter &interpolationFunctionGradient_end,
            result_iter result_begin, result_iter result_end,
            dRdV_iter dRdV_begin, dRdV_iter dRdV_end,
            dRdL_iter dRdL_begin, dRdL_iter dRdL_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        );

        template<
            class surfaceGrowthVelocity_iter,
            class surfaceGrowthVelocityGradient_iter,
            typename testFunction_type,
            class testFunctionGradient_iter,
            typename result_type
        >
        void computeLagrangeMultiplierBalance(
            const surfaceGrowthVelocity_iter &surfaceGrowthVelocity_begin, const surfaceGrowthVelocity_iter &surfaceGrowthVelocity_end,
            const surfaceGrowthVelocityGradient_iter &surfaceGrowthVelocityGradient_begin, const surfaceGrowthVelocityGradient_iter &surfaceGrowthVelocityGradient_end,
            const testFunction_type &testFunction,
            const testFunctionGradient_iter &testFunctionGradient_begin, const testFunctionGradient_iter &testFunctionGradient_end,
            result_type &result
        );

        template<
            class surfaceGrowthVelocity_iter,
            class surfaceGrowthVelocityGradient_iter,
            typename testFunction_type,
            class testFunctionGradient_iter,
            typename result_type,
            typename interpolationFunction_type,
            class interpolationFunctionGradient_iter,
            class dRdV_iter, class dRdUMesh_iter
        >
        void computeLagrangeMultiplierBalance(
            const surfaceGrowthVelocity_iter &surfaceGrowthVelocity_begin, const surfaceGrowthVelocity_iter &surfaceGrowthVelocity_end,
            const surfaceGrowthVelocityGradient_iter &surfaceGrowthVelocityGradient_begin, const surfaceGrowthVelocityGradient_iter &surfaceGrowthVelocityGradient_end,
            const testFunction_type &testFunction,
            const testFunctionGradient_iter &testFunctionGradient_begin, const testFunctionGradient_iter &testFunctionGradient_end,
            const interpolationFunction_type &interpolationFunction,
            const interpolationFunctionGradient_iter &interpolationFunctionGradient_begin, const interpolationFunctionGradient_iter &interpolationFunctionGradient_end,
            result_type &result,
            dRdV_iter dRdV_begin, dRdV_iter dRdV_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        );

    }

}

#include "tardigrade_surface_growth_balance.cpp"

#endif
