/**
 * \file test_tardigrade_balance_of_surface_growth.cpp
 *
 * Tests for tardigrade_balance_of_surface_growth.cpp
 */

#include <tardigrade_balance_of_surface_growth.h>
#include <tardigrade_constitutive_tools.h>
#include <tardigrade_LinearHex.h>
#include <tardigrade_hydraLinearTestMaterial.h>
#define USE_EIGEN
#include <tardigrade_vector_tools.h>

#include <fstream>
#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE test_tardigrade_balance_equations_balance_of_surface_growth
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

typedef double floatType;  //!< Define the float type

typedef std::array<floatType, 3> floatVector;  //!< Define the float vector type

typedef std::array<floatType, 9> secondOrderTensor;  //!< Define the second order tensor type

template <int dim, int num_nodes, class xi_in, class v_tp1_in, class l_tp1_in, class umesh_tp1_in, class X_in,
          class value_out>
void evaluate_at_nodes(const xi_in &xi_begin, const xi_in &xi_end, const v_tp1_in &v_tp1_begin,
                       const v_tp1_in &v_tp1_end, const l_tp1_in &l_tp1_begin, const l_tp1_in &l_tp1_end,
                       const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end, const X_in &X_begin,
                       const X_in &X_end, value_out value_begin, value_out value_end) {
    // Update the mesh nodes
    std::array<typename std::iterator_traits<umesh_tp1_in>::value_type, dim * num_nodes> x_tp1;

    std::transform(X_begin, X_end, umesh_tp1_begin, std::begin(x_tp1),
                   std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>());

    // Instantiate the element
    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<
        floatType, typename std::array<floatType, 24>::const_iterator,
        typename std::array<floatType, 3>::const_iterator, typename std::array<floatType, 8>::iterator,
        typename std::array<floatType, 24>::iterator>
        e(std::cbegin(x_tp1), std::cend(x_tp1), X_begin, X_end);

    // Interpolate the quantities to the local point
    std::array<typename std::iterator_traits<v_tp1_in>::value_type, dim> v_tp1_p;

    std::array<typename std::iterator_traits<l_tp1_in>::value_type, 1> l_tp1_p;

    std::array<typename std::iterator_traits<l_tp1_in>::value_type, 1 * dim> grad_l_tp1_p;

    e.InterpolateQuantity(xi_begin, xi_end, v_tp1_begin, v_tp1_end, std::begin(v_tp1_p), std::end(v_tp1_p));

    e.InterpolateQuantity(xi_begin, xi_end, l_tp1_begin, l_tp1_end, std::begin(l_tp1_p), std::end(l_tp1_p));

    e.GetGlobalQuantityGradient(xi_begin, xi_end, l_tp1_begin, l_tp1_end, std::begin(grad_l_tp1_p),
                                std::end(grad_l_tp1_p));

    // Get the Jacobian of transformation
    std::array<floatType, dim * dim> dxdxi;
    e.GetLocalQuantityGradient(xi_begin, xi_end, std::cbegin(x_tp1), std::cend(x_tp1), std::begin(dxdxi),
                               std::end(dxdxi));

    floatType J =
        tardigradeVectorTools::determinant<typename std::array<floatType, dim * dim>::const_iterator, floatType, 3, 3>(
            std::cbegin(dxdxi), std::cend(dxdxi), 3, 3);

    std::array<floatType, num_nodes> Ns;
    e.GetShapeFunctions(xi_begin, xi_end, std::begin(Ns), std::end(Ns));

    std::array<floatType, num_nodes * dim> dNdx;
    e.GetGlobalShapeFunctionGradients(xi_begin, xi_end, std::cbegin(x_tp1), std::cend(x_tp1), std::begin(dNdx),
                                      std::end(dNdx));

    for (unsigned int i = 0; i < num_nodes; ++i) {
        tardigradeBalanceEquations::surfaceGrowth::computeSurfaceGrowthBalance(
            std::cbegin(v_tp1_p), std::cend(v_tp1_p), l_tp1_p[0], std::cbegin(grad_l_tp1_p), std::cend(grad_l_tp1_p),
            Ns[i], std::cbegin(dNdx) + dim * i, std::cbegin(dNdx) + dim * (i + 1), value_begin + dim * i,
            value_begin + dim * (i + 1));

        std::transform(value_begin + dim * i, value_begin + dim * (i + 1), value_begin + dim * i,
                       std::bind(std::multiplies<>(), std::placeholders::_1, J));
    }
}

template <int dim, int num_nodes, class xi_in, class v_tp1_in, class l_tp1_in, class umesh_tp1_in, class X_in,
          class value_out, class dRdV_iter, class dRdL_iter, class dRdUMesh_iter>
void evaluate_at_nodes(const xi_in &xi_begin, const xi_in &xi_end, const v_tp1_in &v_tp1_begin,
                       const v_tp1_in &v_tp1_end, const l_tp1_in &l_tp1_begin, const l_tp1_in &l_tp1_end,
                       const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end, const X_in &X_begin,
                       const X_in &X_end, value_out value_begin, value_out value_end, dRdV_iter dRdV_begin,
                       dRdV_iter dRdV_end, dRdL_iter dRdL_begin, dRdL_iter dRdL_end, dRdUMesh_iter dRdUMesh_begin,
                       dRdUMesh_iter dRdUMesh_end) {
    // Update the mesh nodes
    std::array<typename std::iterator_traits<umesh_tp1_in>::value_type, dim * num_nodes> x_tp1;

    std::transform(X_begin, X_end, umesh_tp1_begin, std::begin(x_tp1),
                   std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>());

    // Instantiate the element
    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<
        floatType, typename std::array<floatType, 24>::const_iterator,
        typename std::array<floatType, 3>::const_iterator, typename std::array<floatType, 8>::iterator,
        typename std::array<floatType, 24>::iterator>
        e(std::cbegin(x_tp1), std::cend(x_tp1), X_begin, X_end);

    // Interpolate the quantities to the local point
    std::array<typename std::iterator_traits<v_tp1_in>::value_type, dim> v_tp1_p;

    std::array<typename std::iterator_traits<l_tp1_in>::value_type, 1> l_tp1_p;

    std::array<typename std::iterator_traits<l_tp1_in>::value_type, 1 * dim> grad_l_tp1_p;

    e.InterpolateQuantity(xi_begin, xi_end, v_tp1_begin, v_tp1_end, std::begin(v_tp1_p), std::end(v_tp1_p));

    e.InterpolateQuantity(xi_begin, xi_end, l_tp1_begin, l_tp1_end, std::begin(l_tp1_p), std::end(l_tp1_p));

    e.GetGlobalQuantityGradient(xi_begin, xi_end, l_tp1_begin, l_tp1_end, std::begin(grad_l_tp1_p),
                                std::end(grad_l_tp1_p));

    std::array<floatType, dim> value_n;

    std::array<floatType, dim * dim> dRdV_n;

    std::array<floatType, dim * 1> dRdL_n;

    std::array<floatType, dim * dim> dRdUMesh_n;

    std::fill(dRdV_begin, dRdV_end, 0);

    std::fill(dRdL_begin, dRdL_end, 0);

    std::fill(dRdUMesh_begin, dRdUMesh_end, 0);

    // Get the Jacobian of transformation
    std::array<floatType, dim * dim> dxdxi;
    e.GetLocalQuantityGradient(xi_begin, xi_end, std::cbegin(x_tp1), std::cend(x_tp1), std::begin(dxdxi),
                               std::end(dxdxi));

    floatType J =
        tardigradeVectorTools::determinant<typename std::array<floatType, dim * dim>::const_iterator, floatType, 3, 3>(
            std::cbegin(dxdxi), std::cend(dxdxi), 3, 3);

    std::array<floatType, num_nodes> Ns;
    e.GetShapeFunctions(xi_begin, xi_end, std::begin(Ns), std::end(Ns));

    std::array<floatType, num_nodes * dim> dNdx;
    e.GetGlobalShapeFunctionGradients(xi_begin, xi_end, std::cbegin(x_tp1), std::cend(x_tp1), std::begin(dNdx),
                                      std::end(dNdx));

    for (unsigned int i = 0; i < num_nodes; ++i) {
        tardigradeBalanceEquations::surfaceGrowth::computeSurfaceGrowthBalance(
            std::cbegin(v_tp1_p), std::cend(v_tp1_p), l_tp1_p[0], std::cbegin(grad_l_tp1_p), std::cend(grad_l_tp1_p),
            Ns[i], std::cbegin(dNdx) + dim * i, std::cbegin(dNdx) + dim * (i + 1), value_begin + dim * i,
            value_begin + dim * (i + 1));

        std::transform(value_begin + dim * i, value_begin + dim * (i + 1), value_begin + dim * i,
                       std::bind(std::multiplies<>(), std::placeholders::_1, J));
    }

    for (unsigned int i = 0; i < num_nodes; ++i) {
        for (unsigned int j = 0; j < num_nodes; ++j) {
            tardigradeBalanceEquations::surfaceGrowth::computeSurfaceGrowthBalance(
                std::cbegin(v_tp1_p), std::cend(v_tp1_p), l_tp1_p[0], std::cbegin(grad_l_tp1_p),
                std::cend(grad_l_tp1_p), Ns[i], std::cbegin(dNdx) + dim * i, std::cbegin(dNdx) + dim * (i + 1), Ns[j],
                std::cbegin(dNdx) + dim * j, std::cbegin(dNdx) + dim * (j + 1), std::begin(value_n), std::end(value_n),
                std::begin(dRdV_n), std::end(dRdV_n), std::begin(dRdL_n), std::end(dRdL_n), std::begin(dRdUMesh_n),
                std::end(dRdUMesh_n));

            std::transform(std::begin(value_n), std::end(value_n), std::begin(value_n),
                           std::bind(std::multiplies<>(), std::placeholders::_1, J));

            for (unsigned int k = 0; k < dim; ++k) {
                BOOST_TEST(value_n[k] == *(value_begin + dim * i + k));

                *(dRdL_begin + dim * num_nodes * 1 * i + num_nodes * 1 * k + 1 * j + 0) += dRdL_n[k] * J;

                for (unsigned int l = 0; l < dim; ++l) {
                    *(dRdV_begin + dim * num_nodes * dim * i + num_nodes * dim * k + dim * j + l) +=
                        dRdV_n[dim * k + l] * J;

                    *(dRdUMesh_begin + dim * num_nodes * dim * i + num_nodes * dim * k + dim * j + l) +=
                        dRdUMesh_n[dim * k + l] * J;
                }
            }
        }
    }
}

template <int dim, int num_nodes, class xi_in, class v_tp1_in, class umesh_tp1_in, class X_in, class value_out>
void evaluate_at_nodes(const xi_in &xi_begin, const xi_in &xi_end, const v_tp1_in &v_tp1_begin,
                       const v_tp1_in &v_tp1_end, const umesh_tp1_in &umesh_tp1_begin,
                       const umesh_tp1_in &umesh_tp1_end, const X_in &X_begin, const X_in &X_end, value_out value_begin,
                       value_out value_end) {
    // Update the mesh nodes
    std::array<typename std::iterator_traits<umesh_tp1_in>::value_type, dim * num_nodes> x_tp1;

    std::transform(X_begin, X_end, umesh_tp1_begin, std::begin(x_tp1),
                   std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>());

    // Instantiate the element
    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<
        floatType, typename std::array<floatType, 24>::const_iterator,
        typename std::array<floatType, 3>::const_iterator, typename std::array<floatType, 8>::iterator,
        typename std::array<floatType, 24>::iterator>
        e(std::cbegin(x_tp1), std::cend(x_tp1), X_begin, X_end);

    // Interpolate the quantities to the local point
    std::array<typename std::iterator_traits<v_tp1_in>::value_type, dim> v_tp1_p;

    std::array<typename std::iterator_traits<v_tp1_in>::value_type, dim * dim> grad_v_tp1_p;

    e.InterpolateQuantity(xi_begin, xi_end, v_tp1_begin, v_tp1_end, std::begin(v_tp1_p), std::end(v_tp1_p));

    e.GetGlobalQuantityGradient(xi_begin, xi_end, v_tp1_begin, v_tp1_end, std::begin(grad_v_tp1_p),
                                std::end(grad_v_tp1_p));

    // Get the Jacobian of transformation
    std::array<floatType, dim * dim> dxdxi;
    e.GetLocalQuantityGradient(xi_begin, xi_end, std::cbegin(x_tp1), std::cend(x_tp1), std::begin(dxdxi),
                               std::end(dxdxi));

    floatType J =
        tardigradeVectorTools::determinant<typename std::array<floatType, dim * dim>::const_iterator, floatType, 3, 3>(
            std::cbegin(dxdxi), std::cend(dxdxi), 3, 3);

    std::array<floatType, num_nodes> Ns;
    e.GetShapeFunctions(xi_begin, xi_end, std::begin(Ns), std::end(Ns));

    std::array<floatType, num_nodes * dim> dNdx;
    e.GetGlobalShapeFunctionGradients(xi_begin, xi_end, std::cbegin(x_tp1), std::cend(x_tp1), std::begin(dNdx),
                                      std::end(dNdx));

    for (unsigned int i = 0; i < num_nodes; ++i) {
        tardigradeBalanceEquations::surfaceGrowth::computeLagrangeMultiplierBalance(
            std::cbegin(v_tp1_p), std::cend(v_tp1_p), std::cbegin(grad_v_tp1_p), std::cend(grad_v_tp1_p), Ns[i],
            std::cbegin(dNdx) + dim * i, std::cbegin(dNdx) + dim * (i + 1), *(value_begin + 1 * i));

        std::transform(value_begin + 1 * i, value_begin + 1 * (i + 1), value_begin + 1 * i,
                       std::bind(std::multiplies<>(), std::placeholders::_1, J));
    }
}

template <int dim, int num_nodes, class xi_in, class v_tp1_in, class umesh_tp1_in, class X_in, class value_out,
          class dRdV_iter, class dRdUMesh_iter>
void evaluate_at_nodes(const xi_in &xi_begin, const xi_in &xi_end, const v_tp1_in &v_tp1_begin,
                       const v_tp1_in &v_tp1_end, const umesh_tp1_in &umesh_tp1_begin,
                       const umesh_tp1_in &umesh_tp1_end, const X_in &X_begin, const X_in &X_end, value_out value_begin,
                       value_out value_end, dRdV_iter dRdV_begin, dRdV_iter dRdV_end, dRdUMesh_iter dRdUMesh_begin,
                       dRdUMesh_iter dRdUMesh_end) {
    // Update the mesh nodes
    std::array<typename std::iterator_traits<umesh_tp1_in>::value_type, dim * num_nodes> x_tp1;

    std::transform(X_begin, X_end, umesh_tp1_begin, std::begin(x_tp1),
                   std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>());

    // Instantiate the element
    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<
        floatType, typename std::array<floatType, 24>::const_iterator,
        typename std::array<floatType, 3>::const_iterator, typename std::array<floatType, 8>::iterator,
        typename std::array<floatType, 24>::iterator>
        e(std::cbegin(x_tp1), std::cend(x_tp1), X_begin, X_end);

    // Interpolate the quantities to the local point
    std::array<typename std::iterator_traits<v_tp1_in>::value_type, dim> v_tp1_p;

    std::array<typename std::iterator_traits<v_tp1_in>::value_type, dim * dim> grad_v_tp1_p;

    e.InterpolateQuantity(xi_begin, xi_end, v_tp1_begin, v_tp1_end, std::begin(v_tp1_p), std::end(v_tp1_p));

    e.GetGlobalQuantityGradient(xi_begin, xi_end, v_tp1_begin, v_tp1_end, std::begin(grad_v_tp1_p),
                                std::end(grad_v_tp1_p));

    floatType value_n;

    std::array<floatType, 1 * dim> dRdV_n;

    std::array<floatType, 1 * dim> dRdUMesh_n;

    std::fill(dRdV_begin, dRdV_end, 0);

    std::fill(dRdUMesh_begin, dRdUMesh_end, 0);

    // Get the Jacobian of transformation
    std::array<floatType, dim * dim> dxdxi;
    e.GetLocalQuantityGradient(xi_begin, xi_end, std::cbegin(x_tp1), std::cend(x_tp1), std::begin(dxdxi),
                               std::end(dxdxi));

    floatType J =
        tardigradeVectorTools::determinant<typename std::array<floatType, dim * dim>::const_iterator, floatType, 3, 3>(
            std::cbegin(dxdxi), std::cend(dxdxi), 3, 3);

    std::array<floatType, num_nodes> Ns;
    e.GetShapeFunctions(xi_begin, xi_end, std::begin(Ns), std::end(Ns));

    std::array<floatType, num_nodes * dim> dNdx;
    e.GetGlobalShapeFunctionGradients(xi_begin, xi_end, std::cbegin(x_tp1), std::cend(x_tp1), std::begin(dNdx),
                                      std::end(dNdx));

    for (unsigned int i = 0; i < num_nodes; ++i) {
        tardigradeBalanceEquations::surfaceGrowth::computeLagrangeMultiplierBalance(
            std::cbegin(v_tp1_p), std::cend(v_tp1_p), std::cbegin(grad_v_tp1_p), std::cend(grad_v_tp1_p), Ns[i],
            std::cbegin(dNdx) + dim * i, std::cbegin(dNdx) + dim * (i + 1), *(value_begin + 1 * i));

        std::transform(value_begin + 1 * i, value_begin + 1 * (i + 1), value_begin + 1 * i,
                       std::bind(std::multiplies<>(), std::placeholders::_1, J));
    }

    for (unsigned int i = 0; i < num_nodes; ++i) {
        for (unsigned int j = 0; j < num_nodes; ++j) {
            tardigradeBalanceEquations::surfaceGrowth::computeLagrangeMultiplierBalance(
                std::cbegin(v_tp1_p), std::cend(v_tp1_p), std::cbegin(grad_v_tp1_p), std::cend(grad_v_tp1_p), Ns[i],
                std::cbegin(dNdx) + dim * i, std::cbegin(dNdx) + dim * (i + 1), Ns[j], std::cbegin(dNdx) + dim * j,
                std::cbegin(dNdx) + dim * (j + 1), value_n, std::begin(dRdV_n), std::end(dRdV_n),
                std::begin(dRdUMesh_n), std::end(dRdUMesh_n));

            value_n *= J;

            for (unsigned int k = 0; k < 1; ++k) {
                BOOST_TEST(value_n == *(value_begin + 1 * i + k));

                for (unsigned int l = 0; l < dim; ++l) {
                    *(dRdV_begin + 1 * num_nodes * dim * i + num_nodes * dim * k + dim * j + l) +=
                        dRdV_n[dim * k + l] * J;

                    *(dRdUMesh_begin + 1 * num_nodes * dim * i + num_nodes * dim * k + dim * j + l) +=
                        dRdUMesh_n[dim * k + l] * J;
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_computeSurfaceGrowthBalance, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the surface growth balance
     */

    std::array<floatType, 3> v = {1, 2, 3};

    floatType lambda = 2.45;

    std::array<floatType, 3> grad_lambda = {0.11, 0.22, 0.33};

    floatType test = 0.125;

    std::array<floatType, 3> grad_test = {0.1, -0.2, 0.3};

    std::array<floatType, 3> answer = {-0.13375, 0.7125, -0.40125};

    std::array<floatType, 3> result;

    tardigradeBalanceEquations::surfaceGrowth::computeSurfaceGrowthBalance(std::begin(v), std::end(v), lambda,
                                                                           std::begin(grad_lambda),
                                                                           std::end(grad_lambda), test,
                                                                           std::begin(grad_test), std::end(grad_test),
                                                                           std::begin(result), std::end(result));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_computeSurfaceGrowthBalance2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the surface growth balance
     */

    constexpr unsigned int dim       = 3;
    constexpr unsigned int num_nodes = 8;

    std::array<floatType, dim * num_nodes> v = {
        0.2479059,   -0.76876321, -0.36542904, -0.17034758, 0.73261832, -0.49908927, -0.03393147, 0.97111957,
        0.03897024,  0.22578905,  -0.75874267, 0.6526816,   0.20612026, 0.09013601,  -0.31447233, -0.39175842,
        -0.16595558, 0.36260153,  0.75091368,  0.02084467,  0.33862757, 0.17187311,  0.249807,    0.3493781};

    std::array<floatType, 1 * num_nodes> l = {0.68468488,  -0.83361002, 0.52736568,  -0.51266725,
                                              -0.61155408, 0.14491391,  -0.80857497, 0.77065365};

    std::array<floatType, dim * num_nodes> u_mesh = {
        0.34317802, 0.72904971, 0.43857224, 0.32179195, 0.7509966,  0.45705522, 0.2944771,  0.74330725,
        0.45514841, 0.31586316, 0.72136035, 0.43666543, 0.34830949, 0.77712613, 0.427784,   0.32692343,
        0.79907302, 0.44626697, 0.29960857, 0.79138367, 0.44436016, 0.32099464, 0.76943677, 0.42587719};

    std::array<floatType, dim * num_nodes> X = {0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0.,
                                                0., 0., 1., 1., 0., 1., 1., 1., 1., 0., 1., 1.};

    std::array<floatType, dim> local_point = {0.78677833, 0.88832004, 0.00367335};

    std::array<floatType, dim * num_nodes> answer = {
        -1.978679682e-04, -7.313817536e-04, +3.717514247e-04, +2.322889029e-03, -6.019379618e-03, +3.089268667e-03,
        +3.669697337e-02, +1.812398816e-02, +4.642040395e-02, -3.653393777e-03, +1.941782906e-03, +5.591897041e-03,
        -2.027411640e-04, -7.365220033e-04, +5.549611719e-04, +2.311407182e-03, -6.061646770e-03, +4.624358234e-03,
        +3.648381646e-02, +1.829345074e-02, +7.233331169e-02, -3.738060810e-03, +1.960375726e-03, +8.684549327e-03};

    std::array<floatType, dim * num_nodes> result;

    evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v), std::cend(v),
                                      std::cbegin(l), std::cend(l), std::cbegin(u_mesh), std::cend(u_mesh),
                                      std::cbegin(X), std::cend(X), std::begin(result), std::end(result));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    std::array<floatType, num_nodes * dim * num_nodes * dim> dRdV;

    std::array<floatType, num_nodes * dim * num_nodes * 1> dRdL;

    std::array<floatType, num_nodes * dim * num_nodes * dim> dRdUMesh;

    std::fill(std::begin(result), std::end(result), 0);

    evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v), std::cend(v),
                                      std::cbegin(l), std::cend(l), std::cbegin(u_mesh), std::cend(u_mesh),
                                      std::cbegin(X), std::cend(X), std::begin(result), std::end(result),
                                      std::begin(dRdV), std::end(dRdV), std::begin(dRdL), std::end(dRdL),
                                      std::begin(dRdUMesh), std::end(dRdUMesh));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    {
        floatType                       eps      = 1e-6;
        constexpr unsigned int          VAR_SIZE = dim * num_nodes;
        constexpr unsigned int          OUT_SIZE = dim * num_nodes;
        std::array<floatType, VAR_SIZE> x        = v;

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            floatType delta = eps * std::fabs(x[i]) + eps;

            std::array<floatType, VAR_SIZE> xp = x;
            std::array<floatType, VAR_SIZE> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            std::array<floatType, OUT_SIZE> rp, rm;

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(xp),
                                              std::cend(xp), std::cbegin(l), std::cend(l), std::cbegin(u_mesh),
                                              std::cend(u_mesh), std::cbegin(X), std::cend(X), std::begin(rp),
                                              std::end(rp));

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(xm),
                                              std::cend(xm), std::cbegin(l), std::cend(l), std::cbegin(u_mesh),
                                              std::cend(u_mesh), std::cbegin(X), std::cend(X), std::begin(rm),
                                              std::end(rm));

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                BOOST_TEST(dRdV[VAR_SIZE * j + i] == (rp[j] - rm[j]) / (2 * delta));
            }
        }
    }

    {
        floatType                       eps      = 1e-6;
        constexpr unsigned int          VAR_SIZE = 1 * num_nodes;
        constexpr unsigned int          OUT_SIZE = dim * num_nodes;
        std::array<floatType, VAR_SIZE> x        = l;

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            floatType delta = eps * std::fabs(x[i]) + eps;

            std::array<floatType, VAR_SIZE> xp = x;
            std::array<floatType, VAR_SIZE> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            std::array<floatType, OUT_SIZE> rp, rm;

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v),
                                              std::cend(v), std::cbegin(xp), std::cend(xp), std::cbegin(u_mesh),
                                              std::cend(u_mesh), std::cbegin(X), std::cend(X), std::begin(rp),
                                              std::end(rp));

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v),
                                              std::cend(v), std::cbegin(xm), std::cend(xm), std::cbegin(u_mesh),
                                              std::cend(u_mesh), std::cbegin(X), std::cend(X), std::begin(rm),
                                              std::end(rm));

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                BOOST_TEST(dRdL[VAR_SIZE * j + i] == (rp[j] - rm[j]) / (2 * delta));
            }
        }
    }

    {
        floatType                       eps      = 1e-6;
        constexpr unsigned int          VAR_SIZE = dim * num_nodes;
        constexpr unsigned int          OUT_SIZE = dim * num_nodes;
        std::array<floatType, VAR_SIZE> x        = u_mesh;

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            floatType delta = eps * std::fabs(x[i]) + eps;

            std::array<floatType, VAR_SIZE> xp = x;
            std::array<floatType, VAR_SIZE> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            std::array<floatType, OUT_SIZE> rp, rm;

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v),
                                              std::cend(v), std::cbegin(l), std::cend(l), std::cbegin(xp),
                                              std::cend(xp), std::cbegin(X), std::cend(X), std::begin(rp),
                                              std::end(rp));

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v),
                                              std::cend(v), std::cbegin(l), std::cend(l), std::cbegin(xm),
                                              std::cend(xm), std::cbegin(X), std::cend(X), std::begin(rm),
                                              std::end(rm));

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                BOOST_TEST(dRdUMesh[VAR_SIZE * j + i] == (rp[j] - rm[j]) / (2 * delta));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_computeLagrangeMultiplierBalance, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the Lagrange multiplier
     */

    std::array<floatType, 3> v = {1, 2, 3};

    std::array<floatType, 9> grad_v = {0.11, -0.22, 0.33, -0.44, 0.55, -0.66, 0.77, -0.88, 0.99};

    floatType test = 0.125;

    std::array<floatType, 3> grad_test = {0.1, -0.2, 0.3};

    floatType answer = 0.80625;

    floatType result;

    tardigradeBalanceEquations::surfaceGrowth::computeLagrangeMultiplierBalance(std::begin(v), std::end(v),
                                                                                std::begin(grad_v), std::end(grad_v),
                                                                                test, std::begin(grad_test),
                                                                                std::end(grad_test), result);

    BOOST_TEST(answer == result);
}

BOOST_AUTO_TEST_CASE(test_computeLagrangeMultiplierBalance2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the surface growth balance
     */

    constexpr unsigned int dim       = 3;
    constexpr unsigned int num_nodes = 8;

    std::array<floatType, dim * num_nodes> v = {
        0.2479059,   -0.76876321, -0.36542904, -0.17034758, 0.73261832, -0.49908927, -0.03393147, 0.97111957,
        0.03897024,  0.22578905,  -0.75874267, 0.6526816,   0.20612026, 0.09013601,  -0.31447233, -0.39175842,
        -0.16595558, 0.36260153,  0.75091368,  0.02084467,  0.33862757, 0.17187311,  0.249807,    0.3493781};

    std::array<floatType, dim * num_nodes> u_mesh = {
        0.34317802, 0.72904971, 0.43857224, 0.32179195, 0.7509966,  0.45705522, 0.2944771,  0.74330725,
        0.45514841, 0.31586316, 0.72136035, 0.43666543, 0.34830949, 0.77712613, 0.427784,   0.32692343,
        0.79907302, 0.44626697, 0.29960857, 0.79138367, 0.44436016, 0.32099464, 0.76943677, 0.42587719};

    std::array<floatType, dim * num_nodes> X = {0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0.,
                                                0., 0., 1., 1., 0., 1., 1., 1., 1., 0., 1., 1.};

    std::array<floatType, dim> local_point = {0.78677833, 0.88832004, 0.00367335};

    std::array<floatType, 1 * num_nodes> answer = {-3.519554977e-03, -1.922568427e-02, +4.802186299e-02,
                                                   -1.498709318e-02, -3.256041881e-03, -1.694175527e-02,
                                                   +8.939038467e-02, -1.020322998e-02};

    std::array<floatType, 1 * num_nodes> result;

    evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v), std::cend(v),
                                      std::cbegin(u_mesh), std::cend(u_mesh), std::cbegin(X), std::cend(X),
                                      std::begin(result), std::end(result));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    std::array<floatType, num_nodes * 1 * num_nodes * dim> dRdV;

    std::array<floatType, num_nodes * 1 * num_nodes * dim> dRdUMesh;

    std::fill(std::begin(result), std::end(result), 0);

    evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v), std::cend(v),
                                      std::cbegin(u_mesh), std::cend(u_mesh), std::cbegin(X), std::cend(X),
                                      std::begin(result), std::end(result), std::begin(dRdV), std::end(dRdV),
                                      std::begin(dRdUMesh), std::end(dRdUMesh));

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    {
        floatType                       eps      = 1e-6;
        constexpr unsigned int          VAR_SIZE = dim * num_nodes;
        constexpr unsigned int          OUT_SIZE = 1 * num_nodes;
        std::array<floatType, VAR_SIZE> x        = v;

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            floatType delta = eps * std::fabs(x[i]) + eps;

            std::array<floatType, VAR_SIZE> xp = x;
            std::array<floatType, VAR_SIZE> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            std::array<floatType, OUT_SIZE> rp, rm;

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(xp),
                                              std::cend(xp), std::cbegin(u_mesh), std::cend(u_mesh), std::cbegin(X),
                                              std::cend(X), std::begin(rp), std::end(rp));

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(xm),
                                              std::cend(xm), std::cbegin(u_mesh), std::cend(u_mesh), std::cbegin(X),
                                              std::cend(X), std::begin(rm), std::end(rm));

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                BOOST_TEST(dRdV[VAR_SIZE * j + i] == (rp[j] - rm[j]) / (2 * delta));
            }
        }
    }

    {
        floatType                       eps      = 1e-6;
        constexpr unsigned int          VAR_SIZE = dim * num_nodes;
        constexpr unsigned int          OUT_SIZE = 1 * num_nodes;
        std::array<floatType, VAR_SIZE> x        = u_mesh;

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            floatType delta = eps * std::fabs(x[i]) + eps;

            std::array<floatType, VAR_SIZE> xp = x;
            std::array<floatType, VAR_SIZE> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            std::array<floatType, OUT_SIZE> rp, rm;

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v),
                                              std::cend(v), std::cbegin(xp), std::cend(xp), std::cbegin(X),
                                              std::cend(X), std::begin(rp), std::end(rp));

            evaluate_at_nodes<dim, num_nodes>(std::cbegin(local_point), std::cend(local_point), std::cbegin(v),
                                              std::cend(v), std::cbegin(xm), std::cend(xm), std::cbegin(X),
                                              std::cend(X), std::begin(rm), std::end(rm));

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                BOOST_TEST(dRdUMesh[VAR_SIZE * j + i] == (rp[j] - rm[j]) / (2 * delta));
            }
        }
    }
}
