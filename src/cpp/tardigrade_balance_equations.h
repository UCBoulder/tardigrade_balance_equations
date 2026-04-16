/**
 ******************************************************************************
 * \file tardigrade_balance_equations.h
 ******************************************************************************
 * The header file which defines common values for the balance equations
 ******************************************************************************
 */

#ifndef TARDIGRADE_BALANCE_EQUATIONS_H
#define TARDIGRADE_BALANCE_EQUATIONS_H

namespace tardigradeBalanceEquations {

    template <unsigned int _dimension = 3, unsigned int _material_response_dimension = _dimension>
    class BalanceEquationConfigurationBase {
       public:
        //! The spatial dimension of the balance equations
        static constexpr unsigned int dimension = _dimension;

        //! The spatial dimension of the material
        static constexpr unsigned int material_response_dimension = _material_response_dimension;
    };

}  // namespace tardigradeBalanceEquations

#include "tardigrade_balance_equations.tpp"

#endif
