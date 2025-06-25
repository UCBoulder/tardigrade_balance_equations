/**
  ******************************************************************************
  * \file tardigrade_surface_erosion_balance.cpp
  ******************************************************************************
  * The header file for the equations associated with calculating the surface
  * erosion velocity
  ******************************************************************************
  */

#include "tardigrade_surface_erosion_balance.h"

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
        ){
            /*!
             * Compute the surface erosion balance equation
             * 
             * This equation arises from trying to minimize the magnitude of the surface erosion velocity subject
             * to the surface normal constraint on the velocity
             * 
             * \param &surfaceErosionVelocity_begin: The starting iterator for the surface erosion velocity
             * \param &surfaceErosionVelocity_end: The stopping iterator for the surface erosion velocity
             * \param &lagrangeMultiplier: The value of the lagrange multiplier
             * \param &lagrangeMultiplierGradient_begin: The starting iterator for the gradient of the lagrange multiplier
             * \param &lagrangeMultiplierGradient_end: The stopping iterator for the gradient of the lagrange multiplier
             * \param &testFunction: The value of the test function
             * \param &testFunctionGradient_begin: The starting iterator for the gradient of the test function
             * \param &testFunctionGradient_end: The stopping iterator for the gradient of the test function
             * \param &result_begin: The starting iterator for the result
             * \param &result_end: The stopping iterator for the result
             */

            // Definitions only used for error handling
            TARDIGRADE_ERROR_TOOLS_EVAL(
                unsigned int surface_erosion_velocity_size = ( unsigned int )( surfaceErosionVelocity_end - surfaceErosionVelocity_begin );
                unsigned int lagrange_multiplier_gradient_size = ( unsigned int )( lagrangeMultiplierGradient_end - lagrangeMultiplierGradient_begin );
                unsigned int test_function_gradient_size = ( unsigned int )( testFunctionGradient_end - testFunctionGradient_begin );
                unsigned int result_size = ( unsigned int )( result_end - result_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                surface_erosion_velocity_size == result_size,
                "The surface erosion velocity must be the same size as the result ( " + std::to_string( surface_erosion_velocity_size ) + " vs. " + std::to_string( result_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                lagrange_multiplier_gradient_size == result_size,
                "The lagrange multiplier gradient must be the same size as the result ( " + std::to_string( lagrange_multiplier_gradient_size ) + " vs. " + std::to_string( result_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                test_function_gradient_size == result_size,
                "The test function gradient must be the same size as the result ( " + std::to_string( test_function_gradient_size ) + " vs. " + std::to_string( result_size )
            )

            for (
                auto v = std::pair< unsigned int, result_iter >( 0, result_begin );
                v.second != result_end;
                ++v.first, ++v.second ){

                *v.second = testFunction * ( *( surfaceErosionVelocity_begin + v.first ) )
                          - lagrangeMultiplier * ( *( testFunctionGradient_begin + v.first ) )
                          - testFunction * ( *( lagrangeMultiplierGradient_begin + v.first ) );

            }

        }

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
        ){
            /*!
             * Compute the Lagrange Multiplier balance equation
             *
             * This equation arises from trying to minimize the magnitude of the surface erosion velocity subject
             * to the surface normal constraint on the velocity
             * 
             * \param &surfaceErosionVelocity_begin: The starting iterator for the surface erosion velocity
             * \param &surfaceErosionVelocity_end: The stopping iterator for the surface erosion velocity
             * \param &surfaceErosionVelocityGradient_begin: The starting iterator for the surface erosion velocity gradient
             * \param &surfaceErosionVelocityGradient_end: The stopping iterator for the surface erosion velocity gradient
             * \param &testFunction: The value of the test function
             * \param &testFunctionGradient_begin: The starting iterator for the gradient of the test function
             * \param &testFunctionGradient_end: The stopping iterator for the gradient of the test function
             * \param &result
             */

            const unsigned int surface_erosion_velocity_size = ( unsigned int )( surfaceErosionVelocity_end - surfaceErosionVelocity_begin );

            // Definitions only used for error handling
            TARDIGRADE_ERROR_TOOLS_EVAL(
                const unsigned int surface_erosion_velocity_gradient_size = ( unsigned int )( surfaceErosionVelocityGradient_end - surfaceErosionVelocityGradient_begin );
                unsigned int test_function_gradient_size = ( unsigned int )( testFunctionGradient_end - testFunctionGradient_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                surface_erosion_velocity_gradient_size == surface_erosion_velocity_size * surface_erosion_velocity_size,
                "The surface erosion velocity gradient must be the square of the size of the surface erosion velocity ( " + std::to_string( surface_erosion_velocity_gradient_size ) + " vs. " + std::to_string( surface_erosion_velocity_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                test_function_gradient_size == surface_erosion_velocity_size,
                "The test function gradient must be the same size as the surface erosion velocity ( " + std::to_string( test_function_gradient_size ) + " vs. " + std::to_string( surface_erosion_velocity_size )
            )

            result = 0;
            for (
                auto v = std::pair< unsigned int, surfaceErosionVelocity_iter >( 0, surfaceErosionVelocity_begin );
                v.second != surfaceErosionVelocity_end;
                ++v.first, ++v.second
            ){

                result += testFunction * ( *( surfaceErosionVelocityGradient_begin + surface_erosion_velocity_size * v.first + v.first ) )
                        + ( *( testFunctionGradient_begin + v.first ) ) * ( *v.second );

            }

        }

    }

}
