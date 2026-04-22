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
     * A class which defines the configuration of the material response output
     */
    template <int _cauchy_stress_index = 0, int _predicted_internal_energy_index = 9, int _mass_change_rate_index = 10,
              int _body_force_index = 11, int _interphasic_force_index = 14, int _heat_flux_index = 17,
              int _internal_heat_generation_index = 20, int _interphasic_heat_transfer_index = 21,
              int _trace_mass_change_velocity_gradient_index = 22, int _mass_diffusion_index = 23>
    class MaterialResponseOutputConfigurationBase {
       public:
        // MATERIAL RESPONSE VECTOR STRUCTURE
        //! The index of the Cauchy stress in the material response vector
        static constexpr int cauchy_stress_index = _cauchy_stress_index;

        //! The index of the predicted internal energy in the material response vector
        static constexpr int predicted_internal_energy_index = _predicted_internal_energy_index;

        //! The index of the mass-change rate in the material response vector
        static constexpr int mass_change_rate_index = _mass_change_rate_index;

        //! The index of the body force in the material response vector
        static constexpr int body_force_index = _body_force_index;

        //! The index of the net interphasic force in the material response vector
        static constexpr int interphasic_force_index = _interphasic_force_index;

        //! The index of the heat flux in the material response vector
        static constexpr int heat_flux_index = _heat_flux_index;

        //! The index of the internal heat generation in the material response vector
        static constexpr int internal_heat_generation_index = _internal_heat_generation_index;

        //! The index of the interphasic heat transfer in the material response vector
        static constexpr int interphasic_heat_transfer_index = _interphasic_heat_transfer_index;

        //! The index of the trace of the mass-change velocity gradient in the material response vector
        static constexpr int trace_mass_change_velocity_gradient_index = _trace_mass_change_velocity_gradient_index;

        //! The index of the mass-diffusion in the material response vector
        static constexpr int mass_diffusion_index = _mass_diffusion_index;
    };

    /*!
     * A class which defines the DOF configuration for the material response
     */
    template <int _num_phase_dof, int _num_additional_dof, int _density_index = 0, int _displacement_index = 1,
              int _velocity_index = 4, int _temperature_index = 7, int _internal_energy_index = 8,
              int _volume_fraction_index = 9, int _additional_dof_index = 10>
    class MaterialResponseDOFConfigurationBase {
       public:
        //! The number of degrees of freedom associated with a specific phase of material
        static constexpr int num_phase_dof = _num_phase_dof;

        //! The number of additional degrees of freedom common to all phases
        static constexpr int num_additional_dof = _num_additional_dof;

        //! The number of degrees of freedom required in total for each material response
        static constexpr int num_dof = _num_phase_dof + _num_additional_dof;

        // MATERIAL DOF VECTOR STRUCTURE
        //! The index of the density in the material dof vector
        static constexpr int density_index = _density_index;

        //! The index of the displacement in the material dof vector
        static constexpr int displacement_index = _displacement_index;

        //! The index of the velocity in the material dof vector
        static constexpr int velocity_index = _velocity_index;

        //! The index of the temperature in the material dof vector
        static constexpr int temperature_index = _temperature_index;

        //! The index of the internal energy in the material dof vector
        static constexpr int internal_energy_index = _internal_energy_index;

        //! The index of the volume fraction in the material dof vector
        static constexpr int volume_fraction_index = _volume_fraction_index;

        //! The index of the additional degrees of freedom in the material dof vector
        static constexpr int additional_dof_index = _additional_dof_index;
    };

    /*!
     * A class which defines the configuration of the material response
     */
    template <int _num_phase_dof, int _num_additional_dof, int _dimension = 3,
              class _output = MaterialResponseOutputConfigurationBase<>,
              class _dof    = MaterialResponseDOFConfigurationBase<_num_phase_dof, _num_additional_dof>>
    class MaterialResponseConfigurationBase {
       public:
        //! The spatial dimension of the material model
        static constexpr int dimension = _dimension;

        //! The configuration of the output for the material
        using output = _output;

        //! The configuration of the DOF for the material
        using dof = _dof;
    };

    /*!
     * A class which defines the configuration of the the balance equations
     */
    template <class _material, unsigned int _dimension = 3, bool _energy_is_per_unit_volume = true>
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
