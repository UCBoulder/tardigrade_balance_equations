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

    /*!
     * A class which defines the configuration of the material response
     */
    template<int _dimension = 3, int _mass_change_index = 10>
    class MaterialResponseConfigurationBase {
        public:
            //! The spatial dimension of the material model
            static constexpr int dimension = _dimension;

            //! The index of the mass-change rate in the material response vector
            static constexpr int mass_change_index = _mass_change_index;

    };

    /*!
     * A class which defines the configuration of the the balance equations
     */
    template <unsigned int _dimension = 3, class _material = MaterialResponseConfigurationBase<_dimension>>
    class BalanceEquationConfigurationBase {
       public:
        //! The spatial dimension of the balance equations
        static constexpr unsigned int dimension = _dimension;

        //! The configuration of the material response
        using material = _material;

    };

}  // namespace tardigradeBalanceEquations

#include "tardigrade_balance_equations.tpp"

#endif
