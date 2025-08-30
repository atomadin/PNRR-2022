module differential_transmission
    use math
    use phys_const
    use band_structure
    use electron_distribution
    use electron_scattering
    implicit none

    ! Allocate a maximum number of probe frequencies.
    integer, parameter :: wp_n_max = 30
    ! Probe frequencies [eV].
    double precision, dimension(:), allocatable :: wp
    ! Actual number of probe frequencies.
    integer :: wp_n
    ! Broadening [eV].
    double precision :: eta, eta_sq
    ! Coefficient for the differential transmission [fs/nm].
    double precision :: coeff_dt

    ! Optical conductivity [nm/fs].
    double precision, dimension(:), allocatable :: oc_wp, oc_eq_wp

    ! Whether to use the dielectric function due to electron-electron
    ! interaction.  This introduce an additional screening factor.
    logical :: use_ee_screening

    ! Dielectric function, which depends on frequency.
    ! Only meaningful when considering ee screening.
    double precision, dimension(:), allocatable :: eps_wp

    private

    public :: wp_n
    public :: init_differential_transmission, calculate_diff_transm

    public :: equilibrium_optical_conductivity_test

contains

    subroutine init_differential_transmission(mu_eq, tempK_eq)
        implicit none
        double precision, intent(in) :: mu_eq, tempK_eq
        ! character(len=50) :: fmtStr
        logical :: use_band_defined_wavelengths
        integer :: l_sp, l_ba, l_pc, l_wp
        double precision :: temp_eq
        double precision, dimension(:,:,:), allocatable :: f_pc, h_pc
        double precision, dimension(wp_n_max) :: probe_wavelengths_nm
        double precision :: broadening_eV
        namelist /diff_transm_params/ probe_wavelengths_nm, &
            & use_band_defined_wavelengths, broadening_eV, use_ee_screening

        probe_wavelengths_nm = -1.0

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=diff_transm_params)
        close (37)

        eta = broadening_eV
        eta_sq = eta*eta
        temp_eq = k_boltz * tempK_eq

        ! Use probe frequencies defined in the band dispersion.
        ! This might be useful to study specific transitions or to optimize the
        ! calculation to the discretization of the Brillouin Zone.
        if (use_band_defined_wavelengths) then
            wp_n = ot_n
            allocate( wp(wp_n) )
            wp(1:wp_n) = en_ot(1:ot_n)
        else
            ! Count the effective number of given probe frequencies.
            ! Assume there is at least one.
            do l_wp = 1, wp_n_max
                if (probe_wavelengths_nm(l_wp) .gt. 0.0) then
                    ! Update the number of probe frequencies.
                    wp_n = l_wp
                else
                    exit
                end if
            end do
            ! Allocate the array of frequency of the correct size and
            ! populate it.
            allocate( wp(wp_n) )
            ! Convert from nm to eV.
            wp(:) = hc1240 / probe_wavelengths_nm(:)
        end if

        ! Allocate the array for the quantities that have to be evaluated on
        ! the probe frequency mesh.
        allocate( eps_wp(wp_n) )
        allocate( oc_wp(wp_n) )
        allocate( oc_eq_wp(wp_n) )

        ! Coefficient between the differential transmission and the
        ! variation of the optical conductivity.
        coeff_dt = - (4.0 * pi / c_light ) * 2.0 / (refr_ind_top + refr_ind_bottom)

        ! Initialize equilibrium distributions.
        allocate( f_pc(sp_n, ba_n, pc_n) )
        allocate( h_pc(sp_n, ba_n, pc_n) )
        
        do l_sp = 1, sp_n
            do l_ba = 1, ba_n
                do l_pc = 1, pc_n
                    f_pc(l_sp, l_ba, l_pc) = fermi_dirac(en_bz(l_ba, l_pc), mu_eq, temp_eq)
                end do
            end do
        end do
        h_pc = 1.0 - f_pc

        ! Calculate the optical conductivity at equilibrium, for all frequencies.
        call calculate_optical_cond_re(f_pc)
        oc_eq_wp(:) = oc_wp(:)

        deallocate( h_pc )
        deallocate( f_pc )

        ! Save the probe frequencies.
        open (unit=30, file="out-opt-diff-transm-freqs.csv", status="unknown")
        do l_wp = 1, wp_n
            write (30, '(F15.8)') wp(l_wp)
        end do
        close (unit=30)

        ! Save the optical conductivity at equilibrium.
        open (unit=30, file="out-opt-cond-eq.csv", status="unknown")
        do l_wp = 1, wp_n
            write (30, '(F15.8)') oc_eq_wp(l_wp)
        end do
        close (unit=30)

    end subroutine init_differential_transmission

    function optical_cond_unscr_re(hv, f_pc) result( oc )
        ! Calculate the optical conductivity *for graphene*, using the
        ! tight-binding expression in
        ! S. Stauber et al., Phys. Rev. B 78, 085432 (2008).
        ! A general calculation should use the eigenvectors of the tight-binding
        ! model to evaluate the matrix elements of the current operator,
        ! and the result may not be compact, see e.g.
        ! J. Esteve-Paredes et al., SciPost Phys. Core 6, 002 (2023).
        implicit none
        ! Optical conductivity [nm/fs].
        double precision :: oc
        ! Frequency [eV].
        double precision :: hv
        ! Electron distribution.
        double precision, dimension(:,:,:), intent(in) :: f_pc
        ! Wave vector and spin indices.
        integer :: l_pc, l_sp
        ! Band indices.
        integer, parameter :: l_v_ba = 1  ! valence
        integer, parameter :: l_c_ba = 2  ! conduction
        complex(kind=8) :: z
        double precision :: z2, cs

        oc = 0.0
        ! Iterate over the spin and the wave vectors in the Brillouin Zone.
        do l_sp = 1, sp_n
            do l_pc = 1, pc_n
                z = fconn_stauber_bz(l_pc)
                z2 = abs(z)**2
                cs = 18.0 - 4.0 * z2 + 18.0 * (real(z)**2-imag(z)**2) / z2
                oc = oc + cs * &
                    & (f_pc(l_sp, l_v_ba, l_pc) - f_pc(l_sp, l_c_ba, l_pc)) * ( &
                    & delta_light( hv - (en_bz(l_c_ba,l_pc) - en_bz(l_v_ba,l_pc)) ) - &
                    & delta_light( hv - (en_bz(l_v_ba,l_pc) - en_bz(l_c_ba,l_pc)) ) &
                    & )
            end do
        end do
        
        ! Multiply by an overall factor.
        oc = (oc / pc_n) * (e_sq /4.0/hbar) * (pi/6.0/s3) * ((t_hopp**2)/hv)
    end function optical_cond_unscr_re

    function delta_light(en) result (d)
        ! Delta function to set the energy of the transition to that of the
        ! incoming photon.
        implicit none
        ! Result [1/eV].
        double precision :: d
        ! Argument [eV].
        double precision, intent(in) :: en
        double precision :: en_sq

        en_sq = en*en
        ! Avoid very large values which produce numerical exceptions.
        if (en_sq .lt. 100.0*eta_sq) then
            d = exp(-en_sq/eta_sq)/pi_sqrt/eta
        else
            d = 0.0
        endif
    end function delta_light

    subroutine calculate_optical_cond_re(f_pc)
        ! Run the function optical_conductivity_re over an array of
        ! frequencies.
        implicit none
        ! Electron distribution.
        double precision, dimension(:,:,:), intent(in) :: f_pc
        integer :: l_wp

        ! Define the electron-electron screening factor.
        eps_wp = 1.0
        if (use_ee_screening) then
            call optical_eps(eps_wp, wp, wp_n, f_pc)
        end if

        ! Calculate the optical conductivity for all probe frequencies
        ! and divide by the screening factor.
        do l_wp = 1, wp_n
            oc_wp(l_wp) = optical_cond_unscr_re(wp(l_wp), f_pc) / eps_wp(l_wp)
        end do
    end subroutine calculate_optical_cond_re

    subroutine calculate_diff_transm(dt_wp, f_pc)
        ! Calculate the differential transmission according to the formula
        ! DT/T (t) = -4pi/c * 2/(n1+n2) * Re[Ds(t)]
        ! where n1 and n2 are the refractive indices and the Ds is the
        ! optical photoconductivity, i.e. Ds(t) = s(t) - s(equilibrium).
        ! (The conductivity is in Gaussian CGS, e.g. nm/fs, in this formula.)
        implicit none
        double precision, dimension(:), intent(out) :: dt_wp
        double precision, dimension(:,:,:), intent(in) :: f_pc

        ! Calculate the optical conductivity for all probe frequencies.
        call calculate_optical_cond_re(f_pc)

        ! Calculate the photoconductivity and multiply by a constant.
        dt_wp = coeff_dt * (oc_wp - oc_eq_wp)
    end subroutine calculate_diff_transm

    subroutine equilibrium_optical_conductivity_test
        ! Calculate the optical conductivity for low temperature and several
        ! frequencies and broadenings.
        implicit none
        integer :: l_pc, l_ba, l_sp, i_ww, i_dd
        double precision :: mu_eq, tempK_eq, temp_eq
        ! Extremes of the frequency mesh.
        double precision :: wmin, wmax
        ! Extremes of the broadenind mesh.
        double precision :: dmin, dmax
        integer, parameter :: ww_n = 301
        integer, parameter :: dd_n = 7
        double precision, dimension(ww_n) :: ww  ! frequency mesh
        double precision, dimension(dd_n) :: dd  ! broadening mesh
        double precision, dimension(:,:,:), allocatable :: f_pc
        ! Conductivity [nm/fs].
        double precision, dimension(ww_n, dd_n) :: cond_ww
        ! Conductivity divided by the universal value.
        double precision, dimension(ww_n, dd_n) :: cond_uni_ww
        ! Dielectric function.
        double precision, dimension(ww_n) :: eps_ww
        character(len=50) :: fmtStr

        tempK_eq = 300.0
        temp_eq = k_boltz * tempK_eq

        allocate( f_pc(sp_n, ba_n, pc_n) )

        do l_sp = 1, sp_n
            do l_ba = 1, ba_n
                do l_pc = 1, pc_n
                    f_pc(l_sp, l_ba, l_pc) = fermi_dirac(en_bz(l_ba, l_pc), mu_eq, temp_eq)
                end do
            end do
        end do

        ! Calculate the frequency mesh.
        wmin = 0.500
        wmax = 5.000

        do i_ww = 1, ww_n
            ww(i_ww) = wmin + (wmax - wmin) * (i_ww - 1.0) / (ww_n - 1.0)
        end do

        ! Calculate the broadening mesh.
        dmin = 0.050
        dmax = 0.200

        do i_dd = 1, dd_n
            dd(i_dd) = dmin + (dmax - dmin) * (i_dd - 1.0) / (dd_n - 1.0)
        end do

        ! Calculate the unscreened conductivity.
        do i_dd = 1, dd_n
            do i_ww = 1, ww_n
                eta = dd(i_dd)
                eta_sq = eta * eta
                cond_ww(i_ww, i_dd) = optical_cond_unscr_re(ww(i_ww), f_pc)
            end do
        end do

        ! Calculate the dielectric function.
        call print_lindhard(f_pc, "out-test-cond-lindhard")
        eps_ww = 1.0
        if (use_ee_screening) then
            call optical_eps(eps_ww, ww, ww_n, f_pc)
        end if

        ! Factor in screening and make dimensionless.
        do i_dd = 1, dd_n
            cond_uni_ww(:, i_dd) = cond_ww(:, i_dd) / eps_ww(:) / (e_sq / 4.0 / hbar)
        end do

        ! Save results.
        open (unit=30, file="out-test-cond-freq.csv", status="unknown")
        do i_ww = 1, ww_n
            write (30, '(F15.8)') ww(i_ww)
        end do
        close (unit=30)

        open (unit=30, file="out-test-cond-eta.csv", status="unknown")
        do i_dd = 1, dd_n
            write (30, '(F15.8)') dd(i_dd)
        end do
        close (unit=30)

        write (fmtStr, '(A,I2,A)') "(", dd_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-test-cond-unscr.csv", status="unknown")
        do i_ww = 1, ww_n
            write (30, fmtStr) (cond_ww(i_ww, i_dd) / (e_sq / 4.0 / hbar), i_dd = 1, dd_n)
        end do
        close (unit=30)

        write (fmtStr, '(A,I2,A)') "(", dd_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-test-cond.csv", status="unknown")
        do i_ww = 1, ww_n
            write (30, fmtStr) (cond_uni_ww(i_ww, i_dd), i_dd = 1, dd_n)
        end do
        close (unit=30)

        if (use_ee_screening) then
            open (unit=30, file="out-test-cond-eps.csv", status="unknown")
            do i_ww = 1, ww_n
                write (30, '(F15.8)') eps_ww(i_ww)
            end do
            close (unit=30)
        end if

        call print_particle_distro(f_pc, 'out-test-cond-distro.csv')

    end subroutine equilibrium_optical_conductivity_test

end module differential_transmission




