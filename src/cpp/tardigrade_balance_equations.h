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
    template<int _num_phase_dof, int _num_additional_dof, int _dimension = 3, int _mass_change_index = 10>
    class MaterialResponseConfigurationBase {
        public:
            //! The spatial dimension of the material model
            static constexpr int dimension = _dimension;

            //! The index of the mass-change rate in the material response vector
            static constexpr int mass_change_index = _mass_change_index;

            //! The number of degrees of freedom associated with a specific phase of material
            static constexpr int num_phase_dof = _num_phase_dof;

            //! The number of additional degrees of freedom common to all phases
            static constexpr int num_additional_dof = _num_additional_dof;

            //! The number of degrees of freedom required in total for each material response
            static constexpr int num_dof = _num_phase_dof + _num_additional_dof;

    };

    /*!
     * A class which defines the configuration of the the balance equations
     */
    template <class _material, unsigned int _dimension = 3, bool _energy_is_per_unit_volume=true>
    class BalanceEquationConfigurationBase {
       public:
        //! The spatial dimension of the balance equations
        static constexpr unsigned int dimension = _dimension;

        //! Whether the energy is per unit volume or not
        static constexpr bool energy_is_per_unit_volume = _energy_is_per_unit_volume;

        //! The configuration of the material response
        using material = _material;

    };

}  // namespace tardigradeBalanceEquations

#include "tardigrade_balance_equations.tpp"

#endif
