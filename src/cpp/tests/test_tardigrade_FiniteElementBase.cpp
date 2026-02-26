/**
 * \file test_tardigrade_FiniteElementBase.cpp
 *
 * Tests for tardigrade_FiniteElementBase
 */

#include <tardigrade_FiniteElementBase.h>

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE test_tardigrade_FiniteElementBase
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element()

struct cout_redirect {
    cout_redirect(std::streambuf *new_buffer) : old(std::cout.rdbuf(new_buffer)) {}

    ~cout_redirect() { std::cout.rdbuf(old); }

   private:
    std::streambuf *old;
};

struct cerr_redirect {
    cerr_redirect(std::streambuf *new_buffer) : old(std::cerr.rdbuf(new_buffer)) {}

    ~cerr_redirect() { std::cerr.rdbuf(old); }

   private:
    std::streambuf *old;
};

typedef tardigradeBalanceEquations::finiteElement::floatType
    floatType;  //!< Define the float type to be the same as in the balance of mass

typedef tardigradeBalanceEquations::finiteElement::floatVector
    floatVector;  //!< Define the float vector type to be the same as in the balance of mass

typedef tardigradeBalanceEquations::finiteElement::secondOrderTensor
    secondOrderTensor;  //!< Define the second order tensor type to be the same as in the balance of mass

BOOST_AUTO_TEST_CASE(test_placeholder, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {}
