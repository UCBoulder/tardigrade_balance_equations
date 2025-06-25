/**
  ******************************************************************************
  * \file tardigrade_surface_erosion_balance.h
  ******************************************************************************
  * The header file for the equations associated with calculating the surface
  * erosion velocity
  ******************************************************************************
  */

#ifndef TARDIGRADE_SURFACE_EROSION_BALANCE_H
#define TARDIGRADE_SURFACE_EROSION_BALANCE_H

#include<array>

#define USE_EIGEN
#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace surfaceErosion{

        template<
            class surfaceErosionVelocity_iter,
            typename lagrangeMultiplier_type,
            class lagrangeMultiplierGradient_iter,
            typename testFunction_type,
            class testFunctionGradient_iter,
            class result_iter
        >
        void computeSurfaceErosionBalance(
            const surfaceErosionVelocity_iter &surfaceErosionVelocity_begin, const surfaceErosionVelocity_iter &surfaceErosionVelocity_end,
            const lagrangeMultiplier_type &lagrangeMultiplier,
            const lagrangeMultiplierGradient_iter &lagrangeMultiplierGradient_begin, const lagrangeMultiplierGradient_iter &lagrangeMultiplierGradient_end,
            const testFunction_type &testFunction,
            const testFunctionGradient_iter &testFunctionGradient_begin, const testFunctionGradient_iter &testFunctionGradient_end,
            result_iter result_begin, result_iter result_end
        );

        template<
            class surfaceErosionVelocity_iter,
            class surfaceErosionVelocityGradient_iter,
            typename testFunction_type,
            class testFunctionGradient_iter,
            typename result_type
        >
        void computeLagrangeMultiplierBalance(
            const surfaceErosionVelocity_iter &surfaceErosionVelocity_begin, const surfaceErosionVelocity_iter &surfaceErosionVelocity_end,
            const surfaceErosionVelocityGradient_iter &surfaceErosionVelocityGradient_begin, const surfaceErosionVelocityGradient_iter &surfaceErosionVelocityGradient_end,
            const testFunction_type &testFunction,
            const testFunctionGradient_iter &testFunctionGradient_begin, const testFunctionGradient_iter &testFunctionGradient_end,
            result_type &result
        );

    }

}

#include "tardigrade_surface_erosion_balance.cpp"

#endif
