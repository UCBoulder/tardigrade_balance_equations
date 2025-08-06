/**
  ******************************************************************************
  * \file tardigrade_balance_of_surface_growth.cpp
  ******************************************************************************
  * The header file for the equations associated with calculating the surface
  * growth velocity
  ******************************************************************************
  */

#include "tardigrade_balance_of_surface_growth.h"

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
        ){
            /*!
             * Compute the surface growth balance equation
             * 
             * This equation arises from trying to minimize the magnitude of the surface growth velocity subject
             * to the surface normal constraint on the velocity
             * 
             * \param &surfaceGrowthVelocity_begin: The starting iterator for the surface growth velocity
             * \param &surfaceGrowthVelocity_end: The stopping iterator for the surface growth velocity
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
                const unsigned int surface_growth_velocity_size = ( unsigned int )( surfaceGrowthVelocity_end - surfaceGrowthVelocity_begin );
                const unsigned int lagrange_multiplier_gradient_size = ( unsigned int )( lagrangeMultiplierGradient_end - lagrangeMultiplierGradient_begin );
                const unsigned int test_function_gradient_size = ( unsigned int )( testFunctionGradient_end - testFunctionGradient_begin );
                const unsigned int result_size = ( unsigned int )( result_end - result_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                surface_growth_velocity_size == result_size,
                "The surface growth velocity must be the same size as the result ( " + std::to_string( surface_growth_velocity_size ) + " vs. " + std::to_string( result_size )
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

                *v.second = testFunction * ( *( surfaceGrowthVelocity_begin + v.first ) )
                          - lagrangeMultiplier * ( *( testFunctionGradient_begin + v.first ) )
                          - testFunction * ( *( lagrangeMultiplierGradient_begin + v.first ) );

            }

        }

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
        ){
            /*!
             * Compute the surface growth balance equation and Jacobians
             * 
             * This equation arises from trying to minimize the magnitude of the surface growth velocity subject
             * to the surface normal constraint on the velocity
             * 
             * \param &surfaceGrowthVelocity_begin: The starting iterator for the surface growth velocity
             * \param &surfaceGrowthVelocity_end: The stopping iterator for the surface growth velocity
             * \param &lagrangeMultiplier: The value of the lagrange multiplier
             * \param &lagrangeMultiplierGradient_begin: The starting iterator for the gradient of the lagrange multiplier
             * \param &lagrangeMultiplierGradient_end: The stopping iterator for the gradient of the lagrange multiplier
             * \param &testFunction: The value of the test function
             * \param &testFunctionGradient_begin: The starting iterator for the gradient of the test function
             * \param &testFunctionGradient_end: The stopping iterator for the gradient of the test function
             * \param &interpolationFunction: The value of the interpolation function
             * \param &interpolationFunctionGradient_begin: The starting iterator for the gradient of the interpolation function
             * \param &interpolationFunctionGradient_end: The stopping iterator for the gradient of the interpolation function
             * \param &result_begin: The starting iterator for the result
             * \param &result_end: The stopping iterator for the result
             * \param &dRdV_begin: The starting iterator of the derivative of the result w.r.t. the surface growth velocity
             * \param &dRdV_end: The stopping iterator of the derivative of the result w.r.t. the surface growth velocity
             * \param &dRdL_begin: The starting iterator of the derivative of the result w.r.t. the Lagrange multiplier
             * \param &dRdL_end: The stopping iterator of the derivative of the result w.r.t. the Lagrange multiplier
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the result w.r.t. the mesh displacement.
             *     Includes the derivative of the current differential volume
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the result w.r.t. the mesh displacement
             *     Includes the derivative of the current differential volume
             */

            const unsigned int surface_growth_velocity_size = ( unsigned int )( surfaceGrowthVelocity_end - surfaceGrowthVelocity_begin );
            const unsigned int interpolation_function_gradient_size = ( unsigned int )( interpolationFunctionGradient_end - interpolationFunctionGradient_begin );

            // Definitions only used for error handling
            TARDIGRADE_ERROR_TOOLS_EVAL(
                const unsigned int test_function_gradient_size = ( unsigned int )( testFunctionGradient_end - testFunctionGradient_begin );
                const unsigned int result_size = ( unsigned int )( result_end - result_begin );
                const unsigned int dRdV_size = ( unsigned int )( dRdV_end - dRdV_begin );
                const unsigned int dRdL_size = ( unsigned int )( dRdL_end - dRdL_begin );
                const unsigned int dRdUMesh_size = ( unsigned int )( dRdUMesh_end - dRdUMesh_begin );
                
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                interpolation_function_gradient_size == test_function_gradient_size,
                "dRdV must have a size of " + std::to_string( result_size * surface_growth_velocity_size ) + " rather than " + std::to_string( dRdV_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dRdV_size == result_size * surface_growth_velocity_size,
                "dRdV must have a size of " + std::to_string( result_size * surface_growth_velocity_size ) + " rather than " + std::to_string( dRdV_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dRdL_size == result_size * 1,
                "dRdL must have a size of " + std::to_string( result_size * 1 ) + " rather than " + std::to_string( dRdL_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dRdUMesh_size == result_size * interpolation_function_gradient_size,
                "dRdUMesh must have a size of " + std::to_string( result_size * interpolation_function_gradient_size ) + " rather than " + std::to_string( dRdUMesh_size )
            )

            computeSurfaceGrowthBalance(
                surfaceGrowthVelocity_begin, surfaceGrowthVelocity_end,
                lagrangeMultiplier, lagrangeMultiplierGradient_begin, lagrangeMultiplierGradient_end,
                testFunction, testFunctionGradient_begin, testFunctionGradient_end,
                result_begin, result_end
            );

            std::fill(
                dRdV_begin, dRdV_end, 0
            );

            std::fill(
                dRdL_begin, dRdL_end, 0
            );

            std::fill(
                dRdUMesh_begin, dRdUMesh_end, 0
            );

            for (
                auto v = std::pair< unsigned int, result_iter >( 0, result_begin );
                v.second != result_end;
                ++v.first, ++v.second
            ){

                *( dRdV_begin + surface_growth_velocity_size * v.first + v.first ) += testFunction * interpolationFunction; 

                *( dRdL_begin + v.first ) -= testFunction * ( *( interpolationFunctionGradient_begin + v.first ) )
                                           + ( *( testFunctionGradient_begin + v.first ) ) * interpolationFunction;

                for (
                    auto w = std::pair< unsigned int, interpolationFunctionGradient_iter >( 0, interpolationFunctionGradient_begin );
                    w.second != interpolationFunctionGradient_end;
                    ++w.first, ++w.second
                ){

                    *( dRdUMesh_begin + interpolation_function_gradient_size * v.first + w.first )
                    += (
                           ( *( testFunctionGradient_begin + w.first ) ) * lagrangeMultiplier + ( *( lagrangeMultiplierGradient_begin + w.first ) ) * testFunction
                       ) * ( *( interpolationFunctionGradient_begin + v.first ) )
                     + ( *( result_begin + v.first ) ) * ( *w.second );

                }

            }

        }

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
        ){
            /*!
             * Compute the Lagrange Multiplier balance equation
             *
             * This equation arises from trying to minimize the magnitude of the surface growth velocity subject
             * to the surface normal constraint on the velocity
             * 
             * \param &surfaceGrowthVelocity_begin: The starting iterator for the surface growth velocity
             * \param &surfaceGrowthVelocity_end: The stopping iterator for the surface growth velocity
             * \param &surfaceGrowthVelocityGradient_begin: The starting iterator for the surface growth velocity gradient
             * \param &surfaceGrowthVelocityGradient_end: The stopping iterator for the surface growth velocity gradient
             * \param &testFunction: The value of the test function
             * \param &testFunctionGradient_begin: The starting iterator for the gradient of the test function
             * \param &testFunctionGradient_end: The stopping iterator for the gradient of the test function
             * \param &result: The result
             */

            const unsigned int surface_growth_velocity_size = ( unsigned int )( surfaceGrowthVelocity_end - surfaceGrowthVelocity_begin );

            // Definitions only used for error handling
            TARDIGRADE_ERROR_TOOLS_EVAL(
                const unsigned int surface_growth_velocity_gradient_size = ( unsigned int )( surfaceGrowthVelocityGradient_end - surfaceGrowthVelocityGradient_begin );
                const unsigned int test_function_gradient_size = ( unsigned int )( testFunctionGradient_end - testFunctionGradient_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                surface_growth_velocity_gradient_size == surface_growth_velocity_size * surface_growth_velocity_size,
                "The surface growth velocity gradient must be the square of the size of the surface growth velocity ( " + std::to_string( surface_growth_velocity_gradient_size ) + " vs. " + std::to_string( surface_growth_velocity_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                test_function_gradient_size == surface_growth_velocity_size,
                "The test function gradient must be the same size as the surface growth velocity ( " + std::to_string( test_function_gradient_size ) + " vs. " + std::to_string( surface_growth_velocity_size )
            )

            result = 0;
            for (
                auto v = std::pair< unsigned int, surfaceGrowthVelocity_iter >( 0, surfaceGrowthVelocity_begin );
                v.second != surfaceGrowthVelocity_end;
                ++v.first, ++v.second
            ){

                result += testFunction * ( *( surfaceGrowthVelocityGradient_begin + surface_growth_velocity_size * v.first + v.first ) )
                        + ( *( testFunctionGradient_begin + v.first ) ) * ( *v.second );

            }

        }

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
        ){
            /*!
             * Compute the Lagrange Multiplier balance equation
             *
             * This equation arises from trying to minimize the magnitude of the surface growth velocity subject
             * to the surface normal constraint on the velocity
             * 
             * \param &surfaceGrowthVelocity_begin: The starting iterator for the surface growth velocity
             * \param &surfaceGrowthVelocity_end: The stopping iterator for the surface growth velocity
             * \param &surfaceGrowthVelocityGradient_begin: The starting iterator for the surface growth velocity gradient
             * \param &surfaceGrowthVelocityGradient_end: The stopping iterator for the surface growth velocity gradient
             * \param &testFunction: The value of the test function
             * \param &testFunctionGradient_begin: The starting iterator for the gradient of the test function
             * \param &testFunctionGradient_end: The stopping iterator for the gradient of the test function
             * \param &interpolationFunction: The value of the interpolation function
             * \param &interpolationFunctionGradient_begin: The starting iterator for the gradient of the interpolation function
             * \param &interpolationFunctionGradient_end: The stopping iterator for the gradient of the interpolation function
             * \param &result: The result
             * \param &dRdV_begin: The starting iterator of the derivative of the result w.r.t. the surface growth velocity
             * \param &dRdV_end: The stopping iterator of the derivative of the result w.r.t. the surface growth velocity
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the result w.r.t. the mesh displacement.
             *     Includes the derivative of the current differential volume
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the result w.r.t. the mesh displacement
             *     Includes the derivative of the current differential volume
             */

            const unsigned int surface_growth_velocity_size = ( unsigned int )( surfaceGrowthVelocity_end - surfaceGrowthVelocity_begin );
            const unsigned int interpolation_function_gradient_size = ( unsigned int )( interpolationFunctionGradient_end - interpolationFunctionGradient_begin );

            // Definitions only used for error handling
            TARDIGRADE_ERROR_TOOLS_EVAL(
                const unsigned int test_function_gradient_size = ( unsigned int )( testFunctionGradient_end - testFunctionGradient_begin );
                const unsigned int result_size = 1;
                const unsigned int dRdV_size = ( unsigned int )( dRdV_end - dRdV_begin );
                const unsigned int dRdUMesh_size = ( unsigned int )( dRdUMesh_end - dRdUMesh_begin );
                
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                interpolation_function_gradient_size == test_function_gradient_size,
                "dRdV must have a size of " + std::to_string( result_size * surface_growth_velocity_size ) + " rather than " + std::to_string( dRdV_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dRdV_size == result_size * surface_growth_velocity_size,
                "dRdV must have a size of " + std::to_string( result_size * surface_growth_velocity_size ) + " rather than " + std::to_string( dRdV_size )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dRdUMesh_size == result_size * interpolation_function_gradient_size,
                "dRdUMesh must have a size of " + std::to_string( result_size * interpolation_function_gradient_size ) + " rather than " + std::to_string( dRdUMesh_size )
            )

            computeLagrangeMultiplierBalance(
                surfaceGrowthVelocity_begin,         surfaceGrowthVelocity_end,
                surfaceGrowthVelocityGradient_begin, surfaceGrowthVelocityGradient_end,
                testFunction, testFunctionGradient_begin, testFunctionGradient_end,
                result
            );

            std::fill( dRdV_begin, dRdV_end, 0 );

            std::fill( dRdUMesh_begin, dRdUMesh_end, 0 );

            for (
                auto v = std::pair< unsigned int, surfaceGrowthVelocity_iter >( 0, surfaceGrowthVelocity_begin );
                v.second != surfaceGrowthVelocity_end;
                ++v.first, ++v.second
            ){

                *( dRdV_begin + v.first ) += ( *( testFunctionGradient_begin + v.first ) ) * interpolationFunction
                                           + testFunction * ( *( interpolationFunctionGradient_begin + v.first ) );

                *( dRdUMesh_begin + v.first ) += result * ( *( interpolationFunctionGradient_begin + v.first ) );

                for (
                    auto w = std::pair< unsigned int, interpolationFunctionGradient_iter >( 0, interpolationFunctionGradient_begin );
                    w.second != interpolationFunctionGradient_end;
                    ++w.first, ++w.second
                ){

                    *( dRdUMesh_begin + v.first ) -= (
                                                         ( *( testFunctionGradient_begin + v.first ) ) * ( *( surfaceGrowthVelocity_begin + w.first ) )
                                                   +     testFunction * ( *( surfaceGrowthVelocityGradient_begin + test_function_gradient_size * w.first + v.first ) )
                                                     ) * ( *( interpolationFunctionGradient_begin + w.first ) );


                }

            }

        }

    }

}
