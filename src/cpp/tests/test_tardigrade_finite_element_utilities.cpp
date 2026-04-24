/**
 * \file test_tardigrade_finite_element_utilities.cpp
 *
 * Tests for tardigrade_finite_element_utilities
 */

#include <tardigrade_finite_element_utilities.h>

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE test_tardigrade_finite_element_utilities
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

BOOST_AUTO_TEST_CASE(test_computeGradientSpatialJacobian, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the jacobian of the gradient of a quantity w.r.t. the spatial displacement
     */

    std::array<floatType, 12> grad_a = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    floatVector grad_test = {0.1, 0.2, 0.3};

    std::array<floatType, 12> answer = {-0.2, -0.4, -0.6, -0.5, -1., -1.5, -0.8, -1.6, -2.4, -1.1, -2.2, -3.3};

    std::array<floatType, 12> result;
    std::fill(std::begin(result), std::end(result), 0);

    tardigradeBalanceEquations::finiteElement::computeGradientSpatialJacobian(std::begin(grad_a), 12, grad_test, 2,
                                                                              std::begin(result));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_compute_current_rate_of_change, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    double dt = 1.34;

    std::vector<double> v_t = {.11, .22, .33};

    std::vector<double> v_tp1 = {.44, .55, .66};

    std::vector<double> vDot_t = {-0.123, -0.14, 0.5};

    std::vector<double> vDot_tp1(3, 0);

    double alpha = 0.67;

    double dVDotdV;

    tardigradeBalanceEquations::finiteElement::compute_current_rate_of_change(dt, std::begin(v_t), std::end(v_t),
                                                                              std::begin(v_tp1), std::end(v_tp1),
                                                                              std::begin(vDot_t), std::end(vDot_t),
                                                                              alpha, std::begin(vDot_tp1),
                                                                              std::end(vDot_tp1), dVDotdV);

    std::vector<double> result(3, 0);
    for (unsigned int i = 0; i < 3; ++i) {
        result[i] = v_t[i] + dt * ((1 - alpha) * vDot_t[i] + alpha * vDot_tp1[i]);
    }

    BOOST_TEST(v_tp1 == result, CHECK_PER_ELEMENT);

    {
        double                 eps     = 1e-6;
        constexpr unsigned int NUM_VAR = 3;
        constexpr unsigned int NUM_OUT = 3;
        std::vector<double>    x       = v_tp1;
        std::vector<double>    jacobian(NUM_VAR * NUM_OUT, 0);

        for (unsigned int i = 0; i < NUM_VAR; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            double tempp, tempm;

            std::vector<double> rp(NUM_VAR, 0);
            std::vector<double> rm(NUM_VAR, 0);

            tardigradeBalanceEquations::finiteElement::compute_current_rate_of_change(
                dt, std::begin(v_t), std::end(v_t), std::begin(xp), std::end(xp), std::begin(vDot_t), std::end(vDot_t),
                alpha, std::begin(rp), std::end(rp), tempp);

            tardigradeBalanceEquations::finiteElement::compute_current_rate_of_change(
                dt, std::begin(v_t), std::end(v_t), std::begin(xm), std::end(xm), std::begin(vDot_t), std::end(vDot_t),
                alpha, std::begin(rm), std::end(rm), tempm);

            for (unsigned int j = 0; j < NUM_OUT; ++j) {
                jacobian[NUM_VAR * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }

        std::vector<double> jacobian_answer = {dVDotdV, 0, 0, 0, dVDotdV, 0, 0, 0, dVDotdV};

        BOOST_TEST(jacobian_answer == jacobian, CHECK_PER_ELEMENT);
    }
}

BOOST_AUTO_TEST_CASE(test_compute_current_acceleration, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    double dt = 1.34;

    std::vector<double> v_t = {.11, .22, .33};

    std::vector<double> v_tp1 = {.44, .55, .66};

    std::vector<double> vDot_t = {-0.123, -0.14, 0.50};

    std::vector<double> vDDot_t = {0.203, 1.014, 0.75};

    std::vector<double> vDDot_tp1(3, 0);

    double alpha = 0.67;

    double beta = 0.46;

    double dVDDotdV;

    tardigradeBalanceEquations::finiteElement::compute_current_acceleration(
        dt, std::begin(v_t), std::end(v_t), std::begin(v_tp1), std::end(v_tp1), std::begin(vDot_t), std::end(vDot_t),
        std::begin(vDDot_t), std::end(vDDot_t), alpha, beta, std::begin(vDDot_tp1), std::end(vDDot_tp1), dVDDotdV);

    std::vector<double> vDot_tp1(3, 0);
    std::vector<double> result(3, 0);
    for (unsigned int i = 0; i < 3; ++i) {
        vDot_tp1[i] = vDot_t[i] + dt * ((1 - beta) * vDDot_t[i] + beta * vDDot_tp1[i]);
        result[i]   = v_t[i] + dt * ((1 - alpha) * vDot_t[i] + alpha * vDot_tp1[i]);
    }

    BOOST_TEST(v_tp1 == result, CHECK_PER_ELEMENT);

    {
        double                 eps     = 1e-6;
        constexpr unsigned int NUM_VAR = 3;
        constexpr unsigned int NUM_OUT = 3;
        std::vector<double>    x       = v_tp1;
        std::vector<double>    jacobian(NUM_VAR * NUM_OUT, 0);

        for (unsigned int i = 0; i < NUM_VAR; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            double tempp, tempm;

            std::vector<double> rp(NUM_VAR, 0);
            std::vector<double> rm(NUM_VAR, 0);

            tardigradeBalanceEquations::finiteElement::compute_current_acceleration(
                dt, std::begin(v_t), std::end(v_t), std::begin(xp), std::end(xp), std::begin(vDot_t), std::end(vDot_t),
                std::begin(vDDot_t), std::end(vDDot_t), alpha, beta, std::begin(rp), std::end(rp), tempp);

            tardigradeBalanceEquations::finiteElement::compute_current_acceleration(
                dt, std::begin(v_t), std::end(v_t), std::begin(xm), std::end(xm), std::begin(vDot_t), std::end(vDot_t),
                std::begin(vDDot_t), std::end(vDDot_t), alpha, beta, std::begin(rm), std::end(rm), tempm);

            for (unsigned int j = 0; j < NUM_OUT; ++j) {
                jacobian[NUM_VAR * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }

        std::vector<double> jacobian_answer = {dVDDotdV, 0, 0, 0, dVDDotdV, 0, 0, 0, dVDDotdV};

        BOOST_TEST(jacobian_answer == jacobian, CHECK_PER_ELEMENT);
    }
}
