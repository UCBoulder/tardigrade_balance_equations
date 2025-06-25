/**
  * \file test_tardigrade_surface_erosion_balance.cpp
  *
  * Tests for tardigrade_surface_erosion_balance.cpp
  */

#include<tardigrade_surface_erosion_balance.h>
#include<tardigrade_finite_element_utilities.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_hydraLinearTestMaterial.h>
#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_tardigrade_balance_equations_balance_of_energy
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer )
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer )
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

typedef double floatType; //!< Define the float type

typedef std::array< floatType, 3 > floatVector; //!< Define the float vector type

typedef std::array< floatType, 9 > secondOrderTensor; //!< Define the second order tensor type

BOOST_AUTO_TEST_CASE( test_computeSurfaceErosionBalance, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the computation of the surface erosion balance
     */

    std::array< floatType, 3 > v = { 1, 2, 3 };

    floatType lambda = 2.45;

    std::array< floatType, 3 > grad_lambda = { 0.11, 0.22, 0.33 };

    floatType test = 0.125;

    std::array< floatType, 3 > grad_test = { 0.1, -0.2, 0.3 };

    std::array< floatType, 3 > answer = { -0.13375,  0.7125 , -0.40125 };

    std::array< floatType, 3 > result;

    tardigradeBalanceEquations::surfaceErosion::computeSurfaceErosionBalance(
        std::begin( v ), std::end( v ),
        lambda, std::begin( grad_lambda ), std::end( grad_lambda ),
        test,   std::begin( grad_test ),   std::end( grad_test ),
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_computeLagrangeMultiplierBalance, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the computation of the Lagrange multiplier
     */

    std::array< floatType, 3 > v = { 1, 2, 3 };

    std::array< floatType, 9 > grad_v = { 0.11, -0.22, 0.33, -0.44, 0.55, -0.66, 0.77, -0.88, 0.99 };

    floatType test = 0.125;

    std::array< floatType, 3 > grad_test = { 0.1, -0.2, 0.3 };

    floatType answer = 0.80625;

    floatType result;

    tardigradeBalanceEquations::surfaceErosion::computeLagrangeMultiplierBalance(
        std::begin( v ), std::end( v ), std::begin( grad_v ), std::end( grad_v ),
        test,   std::begin( grad_test ),   std::end( grad_test ),
        result
    );

    BOOST_TEST( answer == result );

}
