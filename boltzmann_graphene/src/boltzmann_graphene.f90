!===============================================================================
! boltzmann_graphene.x
!
!    Copyright (C) 2025 Andrea Tomadin, Camillo Tassi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!===============================================================================
! Dimensional units:
! Energy:   eV
! Length:   nm
! Time:     fs
!===============================================================================

! Conditional compilation.
#define WITHCOULOMB     1
#define WITHPHONONS     0
#define WITHSUBSTRATE   0

#define STOPAFTER60S    0
#define EXTRACHECKS     0
#define IEEECHECKS      0

! Replacements in the code.
#define LOGME           print *,
#define VALUNIT         42
#define VALFILE         "out-values.txt"
#define MYOMPNUMTHREADS 24

! Source files with modules.
#include "utils.f90"
#include "electrons.f90"
#include "phonons.f90"
#include "optics.f90"
#include "time_evolution.f90"


program boltzmann_graphene
    use omp_lib
    use color_msg
    use band_structure
    use differential_transmission
    use electron_thermodynamics
    use electron_distribution
    use phonon_modes
    use phonon_distribution
    use electron_scattering
    use electron_phonon_scattering
    use time_evolution

    implicit none
!$    character(len=50) :: msg
    ! Switches to perform a task alternative to the time-evolution.
    ! Used for checks or auxialiary calculations.
    ! Check electron thermodynamic quantities.
    logical :: just_check_electron_thermo
    ! Check the equilibrium optical conductivity for various parameters.
    logical :: just_check_eq_opt_cond
    ! Check the spectral density for electron-electron scattering.
    logical :: just_check_spectral_density
    ! Check energy conservation of the electron distribution.
    logical :: just_check_energy_few_steps
    namelist /boltzmann_graphene_params/ just_check_electron_thermo, &
        & just_check_eq_opt_cond, just_check_spectral_density, &
        & just_check_energy_few_steps

    ! Setup the OMP library, if enabled at compilation.
!$    call omp_set_num_threads(MYOMPNUMTHREADS)
!$    write(msg, '(A,I3)') "boltzmann_graphene: Max number of OMP threads: ", omp_get_max_threads()
!$    LOGME msg

    ! Print alerts on the hard-coded parameters
    if (WITHCOULOMB .eq. 0) then
        LOGME trim(msg_red("boltzmann_graphene: Neglecting Coulomb interactions."))
    end if
    if (WITHPHONONS .eq. 0) then
        LOGME trim(msg_red("boltzmann_graphene: Neglecting electron-phonon scattering."))
    end if
    if (WITHSUBSTRATE .eq. 0) then
        LOGME trim(msg_red("boltzmann_graphene: Neglecting phonon relaxation."))
    end if

    ! Initialize a file where data is appended.
    open (unit=VALUNIT, file=VALFILE, status="unknown")
    write(VALUNIT, '(A)') "! Calculated parameters."
    close(unit=VALUNIT)

    ! Read parameters from input file.
    open (unit=37, file="input.txt", status="old")
    read (37, nml=boltzmann_graphene_params)
    close (37)

    LOGME "boltzmann_graphene: Start program."
    
    ! Copy input file so that it is will be saved in the results folder.
    call copy_file("input.txt", "out-input.txt")

    ! Initialize modules.

    call init_band_structure
    LOGME "boltzmann_graphene: Band structure initialized."

    call init_electron_distribution
    LOGME "boltzmann_graphene: Electron distribution initialized."

    call init_phonon_modes
    LOGME "boltzmann_graphene: Phonon modes initialized."

    call init_phonon_distribution
    LOGME "boltzmann_graphene: Phonon distribution initialized."

    call init_electron_scattering
    LOGME "boltzmann_graphene: Electron-Electron scattering initialized."

    call init_electron_phonon_scattering
    LOGME "boltzmann_graphene: Electron-Phonon scattering initialized."

    call init_time_evolution
    LOGME "boltzmann_graphene: Time-evolution initialized."

    ! Run specific tasks.
    
    if (just_check_electron_thermo) then
        call electron_thermodynamics_test_1
        LOGME "boltzmann_graphene: Electron thermodynamics tested."
        stop
    end if

    if (just_check_eq_opt_cond) then
        call equilibrium_optical_conductivity_test
        LOGME "boltzmann_graphene: Equilibrium optical conductivity tested."
        stop
    end if

    if (just_check_spectral_density) then
        call analyze_spectral_density_support
        LOGME "boltzmann_graphene: Spectral density tested."
        stop
    end if

    if (just_check_energy_few_steps) then
        call analyze_energy_conservation
        LOGME "boltzmann_graphene: Electron energy conservation tested."
        stop
    end if

    ! Run the time-evolution corresponding to the relaxation after the pump.
    
    call calculate_time_evolution

contains

    subroutine copy_file(src_name, dest_name)
        implicit none
        character(len=*), intent(in) :: src_name, dest_name
        integer :: in_unit, out_unit, ios
        character(len=1024) :: line
      
        in_unit  = 10
        out_unit = 20
      
        open(unit=in_unit, file=src_name, status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Error opening source file.'
      
        open(unit=out_unit, file=dest_name, status='replace', action='write', iostat=ios)
        if (ios /= 0) stop 'Error opening destination file.'
      
        do
          read(in_unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          write(out_unit, '(A)') trim(line)
        end do
      
        close(in_unit)
        close(out_unit)
    end subroutine copy_file

end program boltzmann_graphene
