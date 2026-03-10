/**
 * \file test_tardigrade_IntegrationPointBase.cpp
 *
 * Tests for tardigrade_IntegrationPointBase
 */

#include <tardigrade_IntegrationPointBase.h>

#include <fstream>
#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE test_tardigrade_IntegrationPointBase
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

BOOST_AUTO_TEST_CASE(tardigrade_IntegrationPointBase_assembleIntegrationPointResponse,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class IntegrationPointConfigurationTest
        : public tardigradeBalanceEquations::finiteElement::IntegrationPointConfiguration<3, 2, 3, 4, 5> {};

    class IntegrationPointBaseMock
        : public tardigradeBalanceEquations::finiteElement::IntegrationPointBase<IntegrationPointConfigurationTest> {
       public:
        using tardigradeBalanceEquations::finiteElement::IntegrationPointBase<
            IntegrationPointConfigurationTest>::IntegrationPointBase;

        unsigned int num_computeConstantPointResponse = 0;

        unsigned int num_computeVariablePointResponse = 0;

        unsigned int num_computeConstantPointJacobian = 0;

        unsigned int num_computeVariablePointJacobian = 0;

        void public_assembleIntegrationPointResponse() { assembleIntegrationPointResponse(); }

        void public_assembleIntegrationPointJacobian() { assembleIntegrationPointJacobian(); }

        double constant_response = 0;

        double constant_jacobian = 0;

       protected:
        void computeConstantPointResponse() override {
            num_computeConstantPointResponse++;
            constant_response = 0.123;
        }

        void computeVariablePointResponse() override {
            num_computeVariablePointResponse++;
            _response_i = {constant_response * (_i + 1), constant_response * (_i + 2)};
        }

        void computeConstantPointJacobian() override {
            num_computeConstantPointJacobian++;
            constant_jacobian = 0.234;
        }

        void computeVariablePointJacobian() override {
            num_computeVariablePointJacobian++;
            _jacobian_ij = {constant_jacobian * (.1 * _i + .2 * _j + 1), constant_jacobian * (.1 * _i + .2 * _j + 2),
                            constant_jacobian * (.1 * _i + .2 * _j + 3), constant_jacobian * (.1 * _i + .2 * _j + 4),
                            constant_jacobian * (.1 * _i + .2 * _j + 5), constant_jacobian * (.1 * _i + .2 * _j + 6)};
        }
    };

    std::array<double, 4> test = {1, 2, 3, 4};

    std::array<double, 4 * 3> grad_test = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    std::array<double, 5> interpolation = {1, 2, 3, 4, 5};

    std::array<double, 5 * 3> grad_interpolation = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    std::array<double, 3> dof = {1, 2, 3};

    std::array<double, 3 * 3> grad_dof = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    double Jxw = .345;

    std::array<double, 8> response_answer = {0.042435, 0.08487, 0.08487, 0.127305,
                                             0.127305, 0.16974, 0.16974, 0.212175};

    std::array<double, 120> jacobian_answer = {
        0.08073,  0.16146,  0.24219,  0.096876, 0.177606, 0.258336, 0.113022, 0.193752, 0.274482, 0.129168, 0.209898,
        0.290628, 0.145314, 0.226044, 0.306774, 0.32292,  0.40365,  0.48438,  0.339066, 0.419796, 0.500526, 0.355212,
        0.435942, 0.516672, 0.371358, 0.452088, 0.532818, 0.387504, 0.468234, 0.548964, 0.088803, 0.169533, 0.250263,
        0.104949, 0.185679, 0.266409, 0.121095, 0.201825, 0.282555, 0.137241, 0.217971, 0.298701, 0.153387, 0.234117,
        0.314847, 0.330993, 0.411723, 0.492453, 0.347139, 0.427869, 0.508599, 0.363285, 0.444015, 0.524745, 0.379431,
        0.460161, 0.540891, 0.395577, 0.476307, 0.557037, 0.096876, 0.177606, 0.258336, 0.113022, 0.193752, 0.274482,
        0.129168, 0.209898, 0.290628, 0.145314, 0.226044, 0.306774, 0.16146,  0.24219,  0.32292,  0.339066, 0.419796,
        0.500526, 0.355212, 0.435942, 0.516672, 0.371358, 0.452088, 0.532818, 0.387504, 0.468234, 0.548964, 0.40365,
        0.48438,  0.56511,  0.104949, 0.185679, 0.266409, 0.121095, 0.201825, 0.282555, 0.137241, 0.217971, 0.298701,
        0.153387, 0.234117, 0.314847, 0.169533, 0.250263, 0.330993, 0.347139, 0.427869, 0.508599, 0.363285, 0.444015,
        0.524745, 0.379431, 0.460161, 0.540891, 0.395577, 0.476307, 0.557037, 0.411723, 0.492453, 0.573183};

    IntegrationPointBaseMock point(test, grad_test, interpolation, grad_interpolation, dof, grad_dof, Jxw);

    point.public_assembleIntegrationPointResponse();

    BOOST_TEST(point.num_computeConstantPointResponse == 1);

    BOOST_TEST(point.num_computeVariablePointResponse == 4);

    BOOST_TEST(point.num_computeConstantPointJacobian == 0);

    BOOST_TEST(point.num_computeVariablePointJacobian == 0);

    BOOST_TEST(response_answer == *point.getResponse(), CHECK_PER_ELEMENT);

    point.public_assembleIntegrationPointJacobian();

    BOOST_TEST(point.num_computeConstantPointResponse == 1);

    BOOST_TEST(point.num_computeVariablePointResponse == 4);

    BOOST_TEST(point.num_computeConstantPointJacobian == 1);

    BOOST_TEST(point.num_computeVariablePointJacobian == 20);

    BOOST_TEST(jacobian_answer == *point.getJacobian(), CHECK_PER_ELEMENT);
}
