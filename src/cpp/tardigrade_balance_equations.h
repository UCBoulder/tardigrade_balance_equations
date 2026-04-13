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

    template<unsigned int _dimension = 3>
    class BalanceEquationConfigurationBase{

        public:

            //! The spatial dimension
            static constexpr unsigned int dimension = _dimension;

    };

}

#include "tardigrade_balance_equations.tpp"

#endif
