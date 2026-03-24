/**
 * \file test_tardigrade_QuadraticHex.cpp
 *
 * Tests for tardigrade_QuadraticHex
 */

#include <tardigrade_QuadraticHex.h>

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE test_tardigrade_QuadraticHex
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

BOOST_AUTO_TEST_CASE(test_QuadraticHex, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> xi = {
        -1.000000000e+00, -1.000000000e+00, -1.000000000e+00, +1.000000000e+00, -1.000000000e+00, -1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, -1.000000000e+00, -1.000000000e+00, +1.000000000e+00, -1.000000000e+00,
        -1.000000000e+00, -1.000000000e+00, +1.000000000e+00, +1.000000000e+00, -1.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, -1.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +0.000000000e+00, -1.000000000e+00, -1.000000000e+00, +1.000000000e+00, +0.000000000e+00, -1.000000000e+00,
        +0.000000000e+00, +1.000000000e+00, -1.000000000e+00, -1.000000000e+00, +0.000000000e+00, -1.000000000e+00,
        +0.000000000e+00, -1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, -1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        -1.000000000e+00, -1.000000000e+00, +0.000000000e+00, +1.000000000e+00, -1.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, -1.000000000e+00, +1.000000000e+00, +0.000000000e+00};

    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    std::array<floatType, 20> result;

    std::array<floatType, 20> answer;

    for (unsigned int i = 0; i < 8; i++) {
        e.GetShapeFunctions(std::cbegin(xi) + 3 * i, std::cbegin(xi) + 3 * (i + 1), std::begin(result),
                            std::end(result));

        std::fill(std::begin(answer), std::end(answer), 0);

        answer[i] = 1;

        BOOST_TEST(answer == result, CHECK_PER_ELEMENT);
    }

    std::array<floatType, 3> point = {0, 0, 0};

    e.GetShapeFunctions(std::cbegin(point), std::cend(point), std::begin(result), std::end(result));

    answer = {-0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, 0.25, 0.25,
              0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25, 0.25};

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    std::array<floatType, 20 * 3> dNdxi;

    e.GetLocalShapeFunctionGradients(std::cbegin(point), std::cend(point), std::begin(dNdxi), std::end(dNdxi));

    floatType eps = 1e-6;

    std::array<floatType, 20> vp, vm;

    for (unsigned int i = 0; i < 3; ++i) {
        floatType delta = eps * std::fabs(point[i]) + eps;

        std::array<floatType, 3> pp = point;
        std::array<floatType, 3> pm = point;

        pp[i] += delta;
        pm[i] -= delta;

        e.GetShapeFunctions(std::cbegin(pp), std::cend(pp), std::begin(vp), std::end(vp));
        e.GetShapeFunctions(std::cbegin(pm), std::cend(pm), std::begin(vm), std::end(vm));

        for (unsigned int j = 0; j < 20; j++) {
            BOOST_TEST(((vp[j] - vm[j]) / (2 * delta)) == dNdxi[3 * j + i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_QuadraticHex2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    std::array<floatType, 3> point = {+3.929383712e-01, -4.277213301e-01, -5.462970929e-01};

    std::array<floatType, 20 * 4> quantities = {
        -2.581465141e-01, +3.297684532e-01, -6.683990931e-03, +2.041580954e-01, -5.116301113e-01, -1.050387725e-02,
        +4.210808295e-01, -9.721523782e-02, -9.059246807e-01, +6.195543359e-02, +1.061158088e-01, -7.711300713e-01,
        -5.090515143e-01, +6.734037985e-01, +8.115901223e-01, -4.968685792e-01, -5.946891962e-01, +2.666253030e-01,
        -9.519217106e-01, +4.322266506e-01, +7.960701698e-01, -6.494674834e-01, +4.290239639e-02, -5.600843628e-01,
        +7.191786725e-01, +1.747178591e-01, -2.948180939e-01, -6.041881881e-02, -4.472549030e-01, -7.042590916e-01,
        -3.255481872e-01, -5.727464800e-01, +6.439562259e-01, -6.610541328e-01, +8.786398339e-01, +7.562994040e-01,
        -6.350302971e-01, -5.609524198e-01, +3.188754883e-01, +9.038765763e-01, -9.806126697e-01, +5.825500203e-01,
        +7.640208052e-01, -6.252148470e-01, +6.846030636e-01, -2.126818044e-02, -3.881186564e-01, -7.119225455e-01,
        -7.366425047e-01, +5.293732966e-01, -1.080907759e-01, -3.442208103e-01, +6.564608055e-01, -7.686106625e-01,
        -9.191327509e-01, +8.221293722e-01, -1.650269651e-01, -4.654258360e-01, +3.929470640e-01, +7.849815680e-01,
        -8.129219521e-01, -2.200162818e-01, +1.317284240e-01, +3.523961174e-01, -8.367614598e-01, -9.211575885e-01,
        +2.818603460e-01, +5.221061869e-01, +9.918942317e-01, +2.797250987e-01, -4.055630010e-01, -1.573836867e-01,
        +2.238366118e-01, +2.762281524e-01, +3.711713723e-01, +1.956957552e-01, -1.196331340e-01, +7.296372214e-01,
        -2.629868301e-01, +5.585957159e-01};

    std::array<floatType, 4> answer = {+4.764313111e-01, -5.475129130e-01, +3.676472387e-01, +1.033506397e+00};

    std::array<floatType, 4> result;

    e.InterpolateQuantity(std::cbegin(point), std::cend(point), std::cbegin(quantities), std::cend(quantities),
                          std::begin(result), std::end(result));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    std::array<floatType, 12> dQuantitydxi;

    e.GetLocalQuantityGradient(std::cbegin(point), std::cend(point), std::cbegin(quantities), std::cend(quantities),
                               std::begin(dQuantitydxi), std::end(dQuantitydxi));

    floatType eps = 1e-6;

    std::array<floatType, 4> vp, vm;

    for (unsigned int i = 0; i < 3; ++i) {
        floatType delta = eps * std::fabs(point[i]) + eps;

        std::array<floatType, 3> pp = point;
        std::array<floatType, 3> pm = point;

        pp[i] += delta;
        pm[i] -= delta;

        e.InterpolateQuantity(std::cbegin(pp), std::cend(pp), std::begin(quantities), std::end(quantities),
                              std::begin(vp), std::end(vp));
        e.InterpolateQuantity(std::cbegin(pm), std::cend(pm), std::begin(quantities), std::end(quantities),
                              std::begin(vm), std::end(vm));

        for (unsigned int j = 0; j < 4; j++) {
            BOOST_TEST(((vp[j] - vm[j]) / (2 * delta)) == dQuantitydxi[3 * j + i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_QuadraticHex3, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    std::array<floatType, 3> point = {+3.929383712e-01, -4.277213301e-01, -5.462970929e-01};

    std::array<floatType, 20 * 3> answer = {
        +8.458348050e-01, +3.632219699e-01, -4.688684357e-01, +1.463823677e+00, +1.198490629e-01, -1.130343261e+00,
        +4.681267498e-01, -4.856208634e-01, -1.129042247e-01, +1.051851273e+00, -1.211559191e-01, -3.620896831e-01,
        +7.440023149e-01, +5.116796729e-01, -8.929739763e-01, +2.067033117e-01, +2.821793576e-01, -8.086240235e-01,
        +4.169490940e-01, -1.784155895e-01, -5.857161174e-01, +6.681412321e-01, +6.802268622e-02, -5.665370610e-01,
        -6.020028291e-01, -1.044935908e+00, -1.542915989e-01, +5.512509503e-01, +1.634329385e+00, -1.171880270e+00,
        -1.170860567e+00, +4.611636116e-01, +1.559517186e-01, -1.436379797e+00, -2.308104438e-03, +3.603737428e-01,
        -6.137805961e-01, -6.515482615e-01, +1.030780936e+00, -2.503542839e-01, +1.543464278e-01, +6.705528762e-01,
        -5.187677104e-01, -2.956775695e-03, +4.770759820e-01, -6.010498433e-01, -1.423978525e-01, +5.478263648e-01,
        -9.747325517e-01, -8.083679309e-01, +1.082907821e+00, +8.134092215e-01, -5.549512902e-01, +9.001620459e-01,
        -3.686867186e-01, +4.352539484e-01, +5.235912803e-01, -6.934777324e-01, -3.738762767e-02, +5.050058842e-01};

    std::array<floatType, 20 * 3> result;

    e.GetGlobalShapeFunctionGradients(std::cbegin(point), std::cend(point), std::cbegin(x), std::cend(x),
                                      std::begin(result), std::end(result));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    answer = {
        +4.480811161e-01, -2.065459081e-03, -2.759991526e-02, +4.194015946e-01, -4.279136183e-01, -4.540533539e-01,
        -2.113780134e-02, -4.933538806e-01, -1.152103924e-02, +3.688540127e-01, -3.994355406e-01, +6.323402287e-02,
        +3.084076662e-01, +7.462607781e-02, -4.458838801e-01, -5.387738899e-02, +4.706926217e-02, -6.323828749e-01,
        -7.712361202e-02, -3.173806690e-01, -4.239590974e-01, +1.791478975e-01, -1.924314928e-01, -2.530222336e-01,
        -8.674827107e-01, -6.537739752e-01, -6.036401761e-01, +6.317044710e-01, +9.212674989e-01, -5.690532570e-01,
        -3.477162114e-01, +6.537739752e-01, -2.419592604e-01, -6.317044710e-01, +4.015009996e-01, -2.480012068e-01,
        -2.545302772e-01, -1.918254613e-01, +6.036401761e-01, +1.853499927e-01, +2.703114068e-01, +5.690532570e-01,
        -1.020242854e-01, +1.918254613e-01, +2.419592604e-01, -1.853499927e-01, +1.178054150e-01, +2.480012068e-01,
        -5.008157215e-01, -2.129449222e-01, +4.734837953e-01, +5.008157215e-01, -4.886145641e-01, +1.086436229e+00,
        +2.007437648e-01, +4.886145641e-01, +4.354801367e-01, -2.007437648e-01, +2.129449222e-01, +1.897882107e-01};

    e.GetGlobalShapeFunctionGradients(std::cbegin(point), std::cend(point), std::cbegin(X), std::cend(X),
                                      std::begin(result), std::end(result));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_QuadraticHex4, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    std::array<floatType, 3> point = {+3.929383712e-01, -4.277213301e-01, -5.462970929e-01};

    std::array<floatType, 15> answer = {1, 2, 3, 4, 5, 6, 7, 8, 9, -6, -5, -4, -1, -2, -3};

    std::array<floatType, 5> b2 = {-0.1, -0.2, -0.3, -0.4, -0.5};

    std::array<floatType, 100> quantity;
    std::fill(std::begin(quantity), std::end(quantity), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 5; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                quantity[5 * i + j] += answer[3 * j + k] * x[3 * i + k];
            }
            quantity[5 * i + j] += b2[j];
        }
    }

    std::array<floatType, 15> result;

    e.GetGlobalQuantityGradient(std::cbegin(point), std::cend(point), std::cbegin(quantity), std::cend(quantity),
                                std::begin(result), std::end(result));

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    std::fill(std::begin(quantity), std::end(quantity), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 5; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                quantity[5 * i + j] += answer[3 * j + k] * X[3 * i + k];
            }
            quantity[5 * i + j] += b2[j];
        }
    }

    e.GetGlobalQuantityGradient(std::cbegin(point), std::cend(point), std::cbegin(quantity), std::cend(quantity),
                                std::begin(result), std::end(result), 0);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    floatType eps = 1e-6;

    std::array<floatType, 60> dNdx;
    e.GetGlobalShapeFunctionGradients(std::cbegin(point), std::cend(point), std::cbegin(x), std::cend(x),
                                      std::begin(dNdx), std::end(dNdx));

    for (unsigned int i = 0; i < 100; ++i) {
        floatType delta = eps * std::fabs(quantity[i]) + eps;

        std::array<floatType, 100> qp, qm;

        qp = quantity;
        qm = quantity;

        qp[i] += delta;
        qm[i] -= delta;

        std::array<floatType, 15> vp, vm;

        e.GetGlobalQuantityGradient(std::cbegin(point), std::cend(point), std::cbegin(qp), std::cend(qp),
                                    std::begin(vp), std::end(vp));

        e.GetGlobalQuantityGradient(std::cbegin(point), std::cend(point), std::cbegin(qm), std::cend(qm),
                                    std::begin(vm), std::end(vm));

        unsigned int node      = i / 5;
        unsigned int component = i - node * 5;

        for (unsigned int j = 0; j < 15; ++j) {
            unsigned int row = j / 3;
            unsigned int col = j - 3 * row;

            if (component == row) {
                BOOST_TEST(dNdx[3 * node + col] == ((vp[j] - vm[j]) / (2 * delta)));

            } else {
                BOOST_TEST(((vp[j] - vm[j]) / (2 * delta)) == 0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_QuadraticHex5, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    std::array<floatType, 3> point = {+3.929383712e-01, -4.277213301e-01, -5.462970929e-01};

    floatType answer = 0.125;

    floatType result;

    e.GetVolumeIntegralJacobianOfTransformation(std::cbegin(point), std::cend(point), result, 0);

    BOOST_TEST(result == answer);

    answer = 0.125 * 0.3714657360390555;

    e.GetVolumeIntegralJacobianOfTransformation(std::cbegin(point), std::cend(point), result);

    BOOST_TEST(result == answer);
}

BOOST_AUTO_TEST_CASE(test_QuadraticHex6, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    std::vector<std::vector<floatType>> integration_point_answers = {
        {-0.57735027, -0.57735027, -0.57735027},
        {0.57735027,  -0.57735027, -0.57735027},
        {0.57735027,  0.57735027,  -0.57735027},
        {-0.57735027, 0.57735027,  -0.57735027},
        {-0.57735027, -0.57735027, 0.57735027 },
        {0.57735027,  -0.57735027, 0.57735027 },
        {0.57735027,  0.57735027,  0.57735027 },
        {-0.57735027, 0.57735027,  0.57735027 }
    };

    std::array<floatType, 8> weight_answers = {1, 1, 1, 1, 1, 1, 1, 1};

    for (unsigned int i = 0; i < 8; ++i) {
        std::array<floatType, 3> integration_point_result;
        floatType                weight_result;

        e.GetVolumeIntegrationPointData(i, std::begin(integration_point_result), std::end(integration_point_result),
                                        weight_result);

        BOOST_TEST(integration_point_result == integration_point_answers[i], CHECK_PER_ELEMENT);

        BOOST_TEST(weight_result == weight_answers[i]);
    }
}

BOOST_AUTO_TEST_CASE(test_QuadraticHex7, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    floatType answer = 1.0;
    floatType result = 0.0;

    for (unsigned int i = 0; i < 8; ++i) {
        std::array<floatType, 3> xi;
        floatType                weight;
        floatType                J;
        e.GetVolumeIntegrationPointData(i, std::begin(xi), std::end(xi), weight);
        e.GetVolumeIntegralJacobianOfTransformation(std::begin(xi), std::end(xi), J, 0);
        result += weight * J;
    }

    BOOST_TEST(result == answer);
}

BOOST_AUTO_TEST_CASE(test_QuadraticHex8, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    std::vector<std::vector<std::vector<floatType>>> integration_point_answers = {
        {{{-1, -0.57735027, -0.57735027},
          {-1, 0.57735027, -0.57735027},
          {-1, 0.57735027, 0.57735027},
          {-1, -0.57735027, 0.57735027}},
         {{1, -0.57735027, -0.57735027},
          {1, 0.57735027, -0.57735027},
          {1, 0.57735027, 0.57735027},
          {1, -0.57735027, 0.57735027}},
         {{-0.57735027, -1, -0.57735027},
          {0.57735027, -1, -0.57735027},
          {0.57735027, -1, 0.57735027},
          {-0.57735027, -1, 0.57735027}},
         {{-0.57735027, 1, -0.57735027},
          {0.57735027, 1, -0.57735027},
          {0.57735027, 1, 0.57735027},
          {-0.57735027, 1, 0.57735027}},
         {{-0.57735027, -0.57735027, -1},
          {0.57735027, -0.57735027, -1},
          {0.57735027, 0.57735027, -1},
          {-0.57735027, 0.57735027, -1}},
         {{-0.57735027, -0.57735027, 1},
          {0.57735027, -0.57735027, 1},
          {0.57735027, 0.57735027, 1},
          {-0.57735027, 0.57735027, 1}}}
    };

    std::vector<std::vector<floatType>> weight_answers = {
        {1, 1, 1, 1},
        {1, 1, 1, 1},
        {1, 1, 1, 1},
        {1, 1, 1, 1},
        {1, 1, 1, 1},
        {1, 1, 1, 1}
    };

    for (unsigned int s = 0; s < 6; ++s) {
        for (unsigned int i = 0; i < 4; ++i) {
            std::array<floatType, 3> xi;
            floatType                weight;
            e.GetSurfaceIntegrationPointData(s, i, std::begin(xi), std::end(xi), weight);
            BOOST_TEST(xi == integration_point_answers[s][i], CHECK_PER_ELEMENT);
            BOOST_TEST(weight == weight_answers[s][i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_QuadraticHex9, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    std::array<floatType, 60> X = {
        +0.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +0.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00,
        +1.000000000e+00, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +0.000000000e+00,
        +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +5.000000000e-01, +1.000000000e+00, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00,
        +0.000000000e+00, +0.000000000e+00, +5.000000000e-01, +1.000000000e+00, +0.000000000e+00, +5.000000000e-01,
        +1.000000000e+00, +1.000000000e+00, +5.000000000e-01, +0.000000000e+00, +1.000000000e+00, +5.000000000e-01};

    std::array<floatType, 9> A = {+5.309673435e-01, -2.457179268e-01, +4.124004435e-01,
                                  +4.701741978e-01, +7.466197759e-01, +1.988074143e-01,
                                  +3.664293192e-01, +1.395218471e-01, +9.568436484e-01};

    std::array<floatType, 3> b = {+4.987947519e-01, +7.204750069e-01, +2.712531897e-01};

    std::array<floatType, 60> x;

    std::fill(std::begin(x), std::end(x), 0);

    for (unsigned int i = 0; i < 20; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                x[3 * i + j] += X[3 * i + k] * A[3 * j + k];
            }
            x[3 * i + j] += b[j];
        }
    }

    using element_configuration = tardigradeBalanceEquations::finiteElement::QuadraticHexConfiguration;
    tardigradeBalanceEquations::finiteElement::QuadraticHex<element_configuration>
        e(std::cbegin(x), std::cend(x), std::cbegin(X), std::cend(X));

    std::array<floatType, 6> answers = {1, 1, 1, 1, 1, 1};
    std::array<floatType, 6> results;

    for (unsigned int s = 0; s < 6; ++s) {
        results[s] = 0.;
        for (unsigned int i = 0; i < 4; ++i) {
            std::array<floatType, 3> xi;
            floatType                weight;
            floatType                J;
            e.GetSurfaceIntegrationPointData(s, i, std::begin(xi), std::end(xi), weight);
            e.GetSurfaceIntegralJacobianOfTransformation(s, std::begin(xi), std::end(xi), J, 0);

            results[s] += J * weight;
        }
    }
    BOOST_TEST(results == answers, CHECK_PER_ELEMENT);
}
