module PHONON_MODES
    use phys_const
    use color_msg
    use band_structure
    use electron_distribution
    implicit none

    ! Phonon modes (optical):
    !   - longitudinal at Gamma (LOG)
    !   - transverse at Gamma (TOG)
    !   - A' at K
    !   - A' at K1
#define L_LOG 1
#define L_TOG 2
#define L_A1K 3
#define L_AK1 4
    ! Number of phonon modes.
    integer, parameter :: ph_n = 4
    ! Support of each mode in the PC (or BZ).
    ! For each wave vector, mark if the wave vector belongs to a mode.
    ! The supports of the modes on the BZ do not need to be disjoint.
    logical, dimension(:,:), allocatable :: is_ph_pc
    ! Index of the high-symmetry point of each mode.
    integer, dimension(:), allocatable :: l_x_ph
    ! Phonon energies in the BZ mesh [eV].
    double precision, dimension(:, :), allocatable :: en_ph_pc

    ! IPs measured with respect to the appropriate high-symmetry point for
    ! each mode.
    integer, dimension(:, :, :), allocatable :: ip_ph_pc
    ! Wave vectors measured for each mode [nm^-1].
    double precision, dimension(:, :, :), allocatable :: wv_ph_pc

    ! Maximum norm [nm^-1] of the wave vectors measured with respect to the
    ! high-symmetry points.  At most it is half the distance between Gamma-K
    ! (or equivalently K-K1).
    double precision :: d_cutoff

    ! Constants that define the phonon energies [eV].
    double precision, parameter :: ph_en_G = 0.196
    double precision, parameter :: ph_en_K = 0.161

contains

    subroutine init_phonon_modes
        implicit none
        logical :: phonon_wavevector_cutoff_use
        double precision :: phonon_wavevector_cutoff_ratio, r, d
        namelist /phonon_modes_params/ phonon_wavevector_cutoff_use, &
            phonon_wavevector_cutoff_ratio

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=phonon_modes_params)
        close (37)

        ! Calculate the wave vector cutoff around the high-symmetry points.
        d = sqrt(wv_pc(1,l_pc_kappa(1))**2 + wv_pc(2,l_pc_kappa(1))**2)
        if (phonon_wavevector_cutoff_use) then
            r = phonon_wavevector_cutoff_ratio
            if ((r .gt. 1.0) .or. (r .lt. 0.0)) then
                LOGME trim(msg_red("init_phonon_modes: Cutoff ratio not between 0 and 1."))
                stop
            else
                ! Take a fraction of the distance between the K and the Gamma point.
                d_cutoff = r * 0.499 * d
            end if
        else
            ! Set the cutoff to a large value that does not affect the partition of the Brillouin zone.
            d_cutoff = d
        end if

        ! Initialize the support of the modes in the BZ.
        allocate (l_x_ph(ph_n))
        allocate (is_ph_pc(ph_n,pc_n))
        call init_modes_support

        ! Integer projections from high-symmetry points.
        allocate (ip_ph_pc(ndim,ph_n,pc_n))
        allocate (wv_ph_pc(ndim,ph_n,pc_n))
        call init_reduced_wave_vectors

        ! Initialize the phonon dispersion in the BZ.
        allocate (en_ph_pc(ph_n,pc_n))
        call init_modes_dispersion

    end subroutine init_phonon_modes

    subroutine init_modes_support
        ! Initialize the array that states if a wave vector belongs to a mode.
        ! NOTE: Several partitions of the Brillouin Zone can be used here.
        ! In particular, not all wave vectors need to belong to a mode and/or
        ! one vector can belong to more than one mode.
        implicit none
        character(len=50) :: fmt_string
        integer :: l_pc
        double precision :: d_gamma, d_kappa, d_kappa1

        ! High-symmetry point for each mode.
        !           LOG         TOG         A1,K           A1,K1
        l_x_ph = (/ l_pc_gamma, l_pc_gamma, l_pc_kappa(1), l_pc_kappa(2) /)

        ! Iterate over all wave vectors and assign it to a mode if it is
        ! closer than a given quantity to the proper high-symmetry point.
        is_ph_pc = .false.
        do l_pc = 1, pc_n
            ! Distance from Gamma.
            d_gamma = sqrt(wv_bz(1,l_pc)**2 + wv_bz(2,l_pc)**2)
            ! Distance from Kappa.
            d_kappa = sqrt(wv_va_pc(1,1,l_pc)**2 + wv_va_pc(2,1,l_pc)**2)
            ! Distance from Kappa1.
            d_kappa1 = sqrt(wv_va_pc(1,2,l_pc)**2 + wv_va_pc(2,2,l_pc)**2)
            ! Assign the wave vector.
            if (d_gamma .lt. d_kappa .and. d_gamma .lt. d_kappa1) then
                ! Gamma
                if (d_gamma .lt. d_cutoff) then
                    !                         LOG      TOG     A1,K    A1,K1
                    is_ph_pc(:,l_pc) = (/  .true.,  .true., .false., .false. /)
                end if
            else if (d_kappa .lt. d_kappa1) then
                ! Kappa
                if (d_kappa .lt. d_cutoff) then
                    is_ph_pc(:,l_pc) = (/ .false., .false.,  .true., .false. /)
                end if
            else
                ! Kappa1
                if (d_kappa1 .lt. d_cutoff) then
                    is_ph_pc(:,l_pc) = (/ .false., .false., .false.,  .true. /)
                end if
            endif
        end do

        ! Modes support.  Print logicals as 1 or 0.
        ! Print this array immediately to allow debugging if something goes
        ! wrong later on in the code.
        write (fmt_string, "(A,I1,A)") "(", ph_n - 1, "(I1,', '),I1)"
        open (unit=30, file="out-geo-phonon-support.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) merge(1, 0, is_ph_pc(:,l_pc))
        end do
        close (unit=30)
    end subroutine init_modes_support

    subroutine init_reduced_wave_vectors
        ! Calculate the "reduced" phonon wave vectors, i.e. wave vectors
        ! relative to the high-symmetry point of each mode.
        implicit none
        character(len=200) :: msgStr
        integer :: l_ph, l_pc
        integer, dimension(ndim) :: ii_diff, ii_opp

        ! Assign by default large values to wave vectors that do not belong
        ! to any mode, to help noticing if these values enter calculations
        ! by mistake.
        ip_ph_pc = -1000000
        do l_ph = 1, ph_n
            do l_pc = 1, pc_n
                if (is_ph_pc(l_ph,l_pc)) then
                    ! Around the Gamma point we take values in the BZ.
                    if (l_ph .eq. L_LOG .or. l_ph .eq. L_TOG) then
                        ip_ph_pc(:,l_ph,l_pc) = ip_bz(:,l_pc)
                    ! Around the K and K1 point the values are the same as for
                    ! the electron valleys.
                    else if (l_ph .eq. L_A1K) then
                        ip_ph_pc(:,l_ph,l_pc) = ip_va_pc(:,1,l_pc)
                    else
                        ip_ph_pc(:,l_ph,l_pc) = ip_va_pc(:,2,l_pc)
                    end if
                    ! Check that the opposite wave vector is in the support,
                    ! otherwise we have a problem with the support.
                    ! To calculate the opposite vector we must subtract the
                    ! reduced wave vector from the appropriate high-symmetry
                    ! point.
                    ii_diff = ip_pc(:,l_x_ph(l_ph)) - ip_ph_pc(:,l_ph,l_pc)
                    call ip_fold(ii_opp(1), ii_opp(2), ii_diff)
                    if (.not. is_ph_pc(l_ph,l_from_ip(ii_opp(1),ii_opp(2),"init_reduced_wave_vectors"))) then
                        write(msgStr, '(A,I1)') "init_reduced_wave_vectors: Opposite vector not found in support for mode ", l_ph
                        LOGME trim(msg_red(msgStr))
                        LOGME ii_opp
                        stop
                    end if
                end if
            end do
        end do

        ! Calculate the wave vector in physical units [nm^-1].
        do l_ph = 1, ph_n
            do l_pc = 1, ph_n
                wv_ph_pc(:,l_ph,l_pc) = &
                    (ip_ph_pc(1,l_ph,l_pc)/(pri1_n*1.0)) * wv_pri1 &
                    + (ip_ph_pc(2,l_ph,l_pc)/(pri1_n*1.0)) * wv_pri2
            end do
        end do
    end subroutine init_reduced_wave_vectors

    subroutine init_modes_dispersion
        ! Initialize the array with the mode dispersion.
        implicit none
        integer :: l_pc, l_ph

        ! Mode energy is constant.
        ! NOTE: A simple improvement here could be to use linear or quadratic
        ! approximations to the dispersion.  Write the assignment to all modes
        ! explicitly, to make modifications easier.
        ! Assign a large value outside the modes, to hopefully notice if
        ! these values enter calculations by mistake.
        en_ph_pc = -1000000.0
        do l_pc = 1, pc_n
            ! LOG
            l_ph = L_LOG
            if (is_ph_pc(l_ph,l_pc)) then
                en_ph_pc(l_ph,l_pc) = ph_en_G
            end if
            ! TOG
            l_ph = L_TOG
            if (is_ph_pc(l_ph,l_pc)) then
                en_ph_pc(l_ph,l_pc) = ph_en_G
            end if
            ! A' at K
            l_ph = L_A1K
            if (is_ph_pc(l_ph,l_pc)) then
                en_ph_pc(l_ph,l_pc) = ph_en_K
            end if
            ! A' at K1
            l_ph = L_AK1
            if (is_ph_pc(l_ph,l_pc)) then
                en_ph_pc(l_ph,l_pc) = ph_en_K
            end if
        end do
    end subroutine init_modes_dispersion

end module PHONON_MODES

module PHONON_DISTRIBUTION
    use color_msg
    use phonon_modes
    implicit none

    ! Phonon distribution functions.  The first index corresponds to the phonon
    ! mode and the second to the wave vector.

    ! Phonon distribution at equilibrium.
    double precision, dimension(:, :), allocatable :: u_eq_pc

    ! Phonon distribution at the initial time.
    double precision, dimension(:, :), allocatable :: u_i_pc

    ! Timescale of relaxation [fs] to the equilibrium temperature.
    ! This is a phenomenological description of the nonlinear interaction
    ! between optical and acoustic phonons in the substrate.
    ! In principle it is a scattering process, but we can view it here as a
    ! property of the phonon state, hence we include it in this module.
    double precision :: sbs_tau

    private

    public :: u_i_pc

    public :: init_phonon_distribution, calculate_sbs_rate, &
        calculate_phonon_energy_density, &
        calculate_phonon_avg_occup, print_phonon_distro

contains

    subroutine init_phonon_distribution
        implicit none

        integer :: l_pc, l_ph
        double precision :: ph_mu, ph_temp
        ! User-friendly parameters for the input file.
        double precision :: temperature_ph_K, substrate_relaxation_fs
        namelist /phonon_distribution_params/ temperature_ph_K, substrate_relaxation_fs

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=phonon_distribution_params)
        close (37)

        ! Set parameters.
        sbs_tau = substrate_relaxation_fs
        ph_temp = k_boltz*temperature_ph_K
        ph_mu = 0.0  ! phonon chemical potential is always zero

        ! Allocate the arrays.
        allocate (u_i_pc(ph_n, pc_n))
        allocate (u_eq_pc(ph_n, pc_n))

        ! Initialize the equilibrium distribution.
        ! Initialize to a non-physical negative value the distribution on
        ! wave vectors which do not belong to the support of the modes.
        u_eq_pc = -1.0
        do l_ph = 1, ph_n
            do l_pc = 1, pc_n
                if (is_ph_pc(l_ph,l_pc)) then
                    u_eq_pc(l_ph, l_pc) = bose_einstein(en_ph_pc(l_ph, l_pc), ph_mu, ph_temp)
                end if
            end do
        end do

        ! The phonons are initially at equilibrium.
        u_i_pc = u_eq_pc
    end subroutine init_phonon_distribution

    subroutine calculate_sbs_rate(u1_pc, u_pc)
        ! Calculate phenomenological rate of relaxation to equilibrium.
        implicit none
        ! Time-derivative of the phonon distribution.
        double precision, dimension(:, :), intent(out) :: u1_pc
        ! Phonon distribution.
        double precision, dimension(:, :), intent(in) :: u_pc
        integer :: l_ph, l_pc

        u1_pc = 0.0
        do l_ph = 1, ph_n
            do l_pc = 1, pc_n
                if (is_ph_pc(l_ph, l_pc)) then
                    u1_pc(l_ph,l_pc) = -(u_pc(l_ph,l_pc) - u_eq_pc(l_ph,l_pc)) / sbs_tau
                end if
            end do 
        end do
    end subroutine calculate_sbs_rate

    subroutine calculate_phonon_avg_occup(ua, u_pc)
        ! Calculate mode-resolved phonon average occupation [dimensionless].
        implicit none
        ! Particle density [nm^-2].
        double precision, dimension(:), intent(out) :: ua
        ! Phonon distribution.
        double precision, dimension(:, :), intent(in) :: u_pc
        integer :: l_ph, l_pc, n

        ua = 0.0
        do l_ph = 1, ph_n
            n = 0  ! count the number of wave vectors in the mode support
            do l_pc = 1, pc_n
                if (is_ph_pc(l_ph,l_pc)) then
                    ua(l_ph) = ua(l_ph) + u_pc(l_ph,l_pc)
                    n = n + 1
                end if
            end do
            ua(l_ph) = ua(l_ph) / n  ! average
        end do
    end subroutine calculate_phonon_avg_occup

    subroutine calculate_phonon_energy_density(ue, u_pc)
        ! Calculate mode-resolved phonon energy density [eV / nm^2].
        implicit none
        ! Phonon energy density, mode-resolved [nm^-2].
        double precision, dimension(:), intent(out) :: ue
        ! Phonon distribution.
        double precision, dimension(:, :), intent(in) :: u_pc
        integer :: l_ph, l_pc

        ue = 0.0
        do l_ph = 1, ph_n
            do l_pc = 1, pc_n
                if (is_ph_pc(l_ph,l_pc)) then
                    ue(l_ph) = ue(l_ph) + u_pc(l_ph,l_pc) * en_ph_pc(l_ph, l_pc)
                end if
            end do
            ue(l_ph) = ue(l_ph) * length_m_2
        end do
    end subroutine calculate_phonon_energy_density

    subroutine print_phonon_distro(u_pc, fileName)
        ! Print the phonon distribution.
        implicit none
        character(len=*), intent(in) :: fileName
        double precision, dimension(:, :), intent(in) :: u_pc
        character(len=50) :: fmtString
        integer :: l_pc, l_ph, col_n

        col_n = ph_n
        write (fmtString, '(A,I2,A)') "(", col_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file=fileName, status="unknown")
        ! Notice that not all values are meaningful, i.e. not all wave vectors
        ! are in the support of all modes.
        do l_pc = 1, pc_n
            write (30, fmtString) (u_pc(l_ph, l_pc), l_ph=1, ph_n)
        end do
        close (unit=30)
    end subroutine print_phonon_distro

end module PHONON_DISTRIBUTION

module ELECTRON_PHONON_SCATTERING
#if IEEECHECKS
    use, intrinsic :: ieee_arithmetic
#endif
    use math
    use color_msg
    use band_structure
    use phonon_modes
    use phonon_distribution
    implicit none

    ! Constants that define the electron-phonon coupling [eV^2].
    ! See S. Piscanec et al., Phys. Rev. Lett. 93, 185503 (2004).
    double precision, parameter :: epc_g_gamma_sqr_eV_2 = 0.0405
    double precision, parameter :: epc_g_kappa_sqr_eV_2 = 0.0994

    ! Support of electrons that take part in the scattering by phonons.
    ! In principle, all states.  However, it might be efficient to only
    ! consider states around the K and K1 points.
    logical, dimension(:), allocatable :: is_el_pc

    ! Absorption and emission processes.
    integer, parameter :: ae_n = 2
    ! Sign in +/- for absorption/emission of a phonon with wave vector +q.
    integer, dimension(ae_n), parameter :: s_ae = (/ +1, -1 /)

    ! Sum or difference of phonon and electron wave vectors, constrained by
    ! the request that the wave vectors are in the respective supports.
    logical, dimension(:, :, :, :), allocatable :: is_sod_ph_pc
    integer, dimension(:, :, :, :), allocatable :: l_sod_ph_pc

    ! Constant electron-phonon coefficients for all scattering processes,
    ! multiplied by the approximated delta function of energy conservation.
    double precision, dimension(:, :, :, :, :, :), allocatable :: g2delta_mat

    ! Scattering rates.
    ! Electron-by-phonon, depend on the electron and phonon distribution.
    double precision, dimension(:, :, :), allocatable :: gamma_out, gamma_in
    ! Phonon-by-electrons, depend on the electron distribution.
    double precision, dimension(:, :), allocatable :: gamma_em, gamma_abs

    ! Which broadening model to use for energy conservation:
    ! 1 = uniform, constant
    ! 3 = process-dependent, constant
    integer :: broadening_ph_model
#define BPM_UCO 1
#define BPM_PDC 3
    ! Target value for the broadening [eV], used in the uniform constant model.
    double precision :: broadening_ph_eV

    ! Parameter to use an average value of the EPC around the Dirac point.
    ! Needed when the mesh is so sparse that no other states are available
    ! for transitions except the K and K1 points.
    logical :: dirac_point_average

    ! Fine-tuning parameter for the transition rate.
    double precision :: epc_speedup

    ! Fine-tuning parameter for the transition rate around the K point.
    double precision :: dirac_point_speedup

    ! Whether to suppress inter- or intra-band processese.
    ! Might be useful for debugging purposes, or to asses the relative
    ! importance of different relaxation channels.
    logical :: suppress_interband, suppress_intraband

    private

    public :: init_electron_phonon_scattering, calculate_ebp_rate, calculate_pbe_rate
    public :: print_ebp_gammas

contains

    subroutine init_electron_phonon_scattering
        implicit none
        character(len=50) :: fmt_string
        integer :: l_pc
        ! Whether to use a reduced support with respect to the whole PC.
        logical :: electron_wavevec_cutoff_use
        ! Fraction of the maximum cutoff.
        double precision :: electron_wavevec_cutoff_ratio, r
        ! Maximum norm [nm^-2] of wave vectors measured from K or K1.
        double precision :: d_el_cutoff
        double precision :: d, d_kappa, d_kappa1

        namelist /electron_phonon_scatt_params/ &
            suppress_interband, suppress_intraband, &
            electron_wavevec_cutoff_use, electron_wavevec_cutoff_ratio, &
            epc_speedup, dirac_point_average, dirac_point_speedup, &
            broadening_ph_model, broadening_ph_eV
            
        open (unit=37, file="input.txt", status="old")
        read (37, nml=electron_phonon_scatt_params)
        close (37)

        ! Support of electrons that take part in the scattering by phonons.
        allocate (is_el_pc(pc_n))
        ! By default, all wave vectors are included.
        is_el_pc = .true.
        if (electron_wavevec_cutoff_use) then
            ! Calculate the wave vector cutoff around the K and K1 points.
            d = sqrt(wv_pc(1,l_pc_kappa(1))**2 + wv_pc(2,l_pc_kappa(1))**2)
            r = electron_wavevec_cutoff_ratio
            if ((r .gt. 1.0) .or. (r .lt. 0.0)) then
                LOGME trim(msg_red("init_electron_phonon_scattering: Cutoff ratio not between 0 and 1."))
                stop
            else
                ! Take a fraction of the distance between the K and the Gamma
                ! point.
                d_el_cutoff = r * d
            end if
            do l_pc = 1, pc_n
                ! Distance from Kappa.
                d_kappa = sqrt(wv_va_pc(1,1,l_pc)**2 + wv_va_pc(2,1,l_pc)**2)
                ! Distance from Kappa1.
                d_kappa1 = sqrt(wv_va_pc(1,2,l_pc)**2 + wv_va_pc(2,2,l_pc)**2)
                ! Exclude wave vectors which are away from both K and K1.
                if ((d_kappa .gt. d_el_cutoff) .and. (d_kappa1 .gt. d_el_cutoff)) then
                    is_el_pc(l_pc) = .false.
                end if
            end do
        end if

        ! Constant arrays, need to be calculated only once.
        ! Lookup arrays for k +/- q, i.e. whether the result is in the support
        ! for the electron wave vectors and the corresponding index.
        ! Indices:
        !   1. absorption, k+q (1); emission, k-q (2)
        !   2. phonon mode index
        !   3. phonon wave vector index
        !   4. incoming electron wave vector index
        allocate (is_sod_ph_pc(ae_n, ph_n, pc_n, pc_n))
        allocate (l_sod_ph_pc(ae_n, ph_n, pc_n, pc_n))
        ! Coefficient in the Boltzmann summation.
        ! Indices:
        !   1. absorption, k+q (1); emission, k-q (2)
        !   2. phonon mode index
        !   3. phonon wave vector index
        !   4. outgoing electron band index
        !   5. incoming electron band index
        !   6. incoming electron wave vector index
        allocate (g2delta_mat(ae_n, ph_n, pc_n, ba_n, ba_n, pc_n))
        call calculate_EPC_matrix_elements

        ! Arrays with the scattering rates.
        ! Same size as the phonon distribution.
        allocate (gamma_abs(ph_n, pc_n))
        allocate (gamma_em(ph_n, pc_n))
        ! Same size as the electron distribution.
        allocate (gamma_in(sp_n, ba_n, pc_n))
        allocate (gamma_out(sp_n, ba_n, pc_n))

        ! Set-up the calculation of the broadening in the delta functions
        ! of energy conservation.
        select case (broadening_ph_model)
            case (BPM_UCO)
                LOGME trim(msg_aqua("init_electron_phonon_scattering: Using uniform constant broadening."))
            case (BPM_PDC)
                LOGME trim(msg_aqua("init_electron_phonon_scattering: Using process-dependent broadening."))
            case default
                LOGME trim(msg_red("init_electron_phonon_scattering: Undefined broadening model."))
                stop
        end select

        ! Save the arrays.
        ! Electrons support.  Print logicals as 1 or 0.
        fmt_string = "(I1)"
        open (unit=30, file="out-geo-phonon-scatt-el-supp.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) merge(1, 0, is_el_pc(l_pc))
        end do
        close (unit=30)
    end subroutine init_electron_phonon_scattering

    function angle_betw_vec(ii_1, ii_2) result(theta)
        ! Calculate the angle between two wave vectors, given their integer
        ! projections.
        ! The angle is between 0 and pi.
        implicit none
#if IEEECHECKS
        character(len=200) :: msg
#endif
        integer :: s1, s2
        integer, dimension(ndim), intent(in) :: ii_1, ii_2
        double precision :: theta, costheta

        ! Check that the direction of the vectors is well defined.
#if EXTRACHECKS
        if (((ii_1(1) .eq. 0) .and. (ii_1(2) .eq. 0)) .or. &
            ((ii_2(1) .eq. 0) .and. (ii_2(2) .eq. 0))) then
            LOGME trim(msg_red("angle_betw_vec: Zero vector."))
            stop 
        end if
#endif
        ! Avoid boundary problems when vectors are parallel or antiparallel.
        ! The condition for being parallel or antiparallel is that the vector
        ! product is zero, or equivalently that the projection of 1 onto the
        ! orthogonal of 2 is zero.
        if ((ii_1(1)*ii_2(2) - ii_1(2)*ii_2(1)) .eq. 0) then
            ! Get the sign of a (nonzero) component.
            if (ii_1(1) .ne. 0) then
                s1 = sign(1, ii_1(1))
                s2 = sign(1, ii_2(1))
            else
                s1 = sign(1, ii_1(2))
                s2 = sign(1, ii_2(2))
            end if
            if (s1 .eq. s2) then
                ! Vectors are parallel.
                costheta = 1.0
                theta = 0.0
            else
                ! Vectors are antiparallel.
                costheta = -1.0
                theta = pi
            end if
        else
            costheta = float(dot_product(ii_1,ii_2)) &
                / sqrt(float(dot_product(ii_1,ii_1))) &
                / sqrt(float(dot_product(ii_2,ii_2)))
            theta = acos(costheta)
        end if
#if IEEECHECKS
        if (ieee_is_nan(theta)) then
            write(msg, '(A,F15.8,4I4)') "angle_betw_vec: Theta NaN for ", costheta, ii_1, ii_2
            LOGME trim(msg_red(msg))
            stop
        end if
#endif
    end function angle_betw_vec

    function g2_func(l_ph, ii_ph_q, l_o_ba, l_i_ba, ii_va_i) result(g2)
        ! Calculate the EPC coefficient g^2 given the indices of the states
        ! which define the scattering *coefficient*:
        !
        ! l_ph      phonon mode
        ! ii_ph_q   absorbed phonon wavevector q, measured from the
        !           high-symmetry point of the mode
        ! l_out_ba  incoming electron band
        ! l_in_ba   outgoingelectron band
        ! ii_va_i   incoming electron wave vector, measured from the center of
        !           the valley
        !
        ! Notice that the wave vectors of the scattering *event* might be
        ! different because the coefficient of several events can be
        ! simplified in the form g_{\pm q}^{k_{1}} and we take advantage of this
        ! fact to simplify the algorithm.
        !
        implicit none
        double precision :: g2
        integer, intent(in) :: l_ph, l_o_ba, l_i_ba
        integer, dimension(:), intent(in) :: ii_ph_q, ii_va_i
        logical :: is_intraband, is_interband
        ! Reduced integer projections of the outgoing electron.
        integer, dimension(ndim) :: ii_va_o
        ! Angles between wave vectors.
        double precision :: theta_iq, theta_oq
        double precision :: cos_p, cos_m
        ! EPC coefficients multiplied by a factor.
        double precision :: epc2_G, epc2_K

        ! Categorize the electron scattering.
        if (l_o_ba .eq. l_i_ba) then
            is_intraband = .true.
        else
            is_intraband = .false.
        end if
        is_interband = .not. is_intraband

        ! Suppress specific transition families.
        if (is_intraband .and. suppress_intraband) then
            g2 = 0.0
            ! Exit the function here.
            return
        else if (is_intraband .and. suppress_interband) then
            g2 = 0.0
            ! Exit the function here.
            return
        end if

        ! Calculate the integer projections of the "out" wave vector.
        ! This is the sum of the "reduced" wave vectors, measured with respect
        ! to the high-symmetry point of the phonon supports, or from the center
        ! of the electron valleys.
        ! The actual final electron wave vector could be in another valley.
        ! Here we do not fold the result into the BZ or the PC, because we want
        ! the coupling coefficient to be continuous.
        ii_va_o = ii_va_i + ii_ph_q

        ! Multiply the EPC coefficients by a factor.
        ! Used everywhere except for Dirac-point averaging, where a specific
        ! coefficient is used.
        epc2_G = epc_g_gamma_sqr_eV_2 * epc_speedup
        epc2_K = epc_g_kappa_sqr_eV_2 * epc_speedup

        ! Pathological cases due to vanishing wave vector(s).
        if ( (ii_ph_q(1) .eq. 0) .and. (ii_ph_q(2) .eq. 0) .and. &
            (ii_va_i(1) .eq. 0) .and. (ii_va_i(2) .eq. 0) ) then
            ! q=0 and k=0

            ! PROBLEM WITH A MESH WHICH IS NOT FINE ENOUGH:
            ! States in conduction band and states in valence band are
            ! separated by an energy which is larger than the phonon energy.
            ! So no phonon, with any wave vector, leads to recombination.
            ! The only exception are the Dirac points, where recombination
            ! could take place with q=0 phonons.
            ! However, if we set the coefficient to zero, we inhibit
            ! recombination entirely.  Hence some average value must be used.
            ! In other words, k=0 wave vectors should represent the whole patch
            ! of wave vectors around that point of the mesh, and using the
            ! (pathological) value appropriate to the Dirac point leads to a
            ! qualitatively different global dynamics.
            if (dirac_point_average) then
                if ((l_ph .eq. L_LOG) .or. (l_ph .eq. L_TOG)) then
                    g2 = epc_g_gamma_sqr_eV_2 * dirac_point_speedup
                else  ! L_A1K or L_AK1
                    g2 = epc_g_kappa_sqr_eV_2 * 2.0 * dirac_point_speedup
                end if
            else
                g2 = 0.0
            end if
            ! Exit the function here.
            return
        else if ( (ii_ph_q(1) .eq. 0) .and. (ii_ph_q(2) .eq. 0) ) then
            ! q=0, k!=0
            if (is_intraband) then
                g2 = 0.0
            else
                if ((l_ph .eq. L_LOG) .or. (l_ph .eq. L_TOG)) then
                    g2 = epc2_G
                else  ! L_A1K or L_AK1
                    g2 = epc2_K * 2.0
                end if
            endif
            ! Exit the function here.
            return
        else if ( ((ii_va_i(1) .eq. 0) .and. (ii_va_i(2) .eq. 0)) & 
            ! k=0, q!=0
            .or. ((ii_va_o(1) .eq. 0) .and. (ii_va_o(2) .eq. 0)) ) then
            ! k+q=0, q!=0
            ! Substitute the cosine with its vanishing average.
            cos_p = 0.0
            cos_m = 0.0
        else
            ! The angles between the (reduced) wave vectors are well-defined.
            ! print *, ii_va_i, " | ", ii_ph_q, " | ", ii_va_o  ! Debug.
            theta_iq = angle_betw_vec(ii_va_i,ii_ph_q)
            theta_oq = angle_betw_vec(ii_va_o,ii_ph_q)
            cos_p = cos(theta_oq + theta_iq)
            cos_m = cos(theta_oq - theta_iq)
        end if

        ! Formulas in Piscanec et al., PRL 2004.
        if (l_ph .eq. L_LOG) then
            if (is_interband) then
                g2 = epc2_G * (1.0 + cos_p)
            else
                g2 = epc2_G * (1.0 - cos_p)
            end if
        else if (l_ph .eq. L_TOG) then
            if (is_interband) then
                g2 = epc2_G * (1.0 - cos_p)
            else
                g2 = epc2_G * (1.0 + cos_p)
            end if
        else if (l_ph .eq. L_A1K) then
            if (is_interband) then
                g2 = epc2_K * (1.0 + cos_m)
            else
                g2 = epc2_K * (1.0 - cos_m)
            end if
        else  ! (l_ph .eq. L_AK1)
            if (is_interband) then
                g2 = epc2_K * (1.0 + cos_m)
            else
                g2 = epc2_K * (1.0 - cos_m)
            end if
        end if
    end function g2_func

    function delta_eph(l_o_ba, l_o_pc, l_i_ba, l_i_pc, en_d) result(d)
        ! Gaussian approximation to the Dirac delta function of energy
        ! conservation.  The approximation has to be Gaussian to work:
        ! See C. Illg et al., J Theor Appl Phys 10, 1 (2016)
        ! DOI: 10.1007/s40094-015-0193-5
        implicit none
        ! Result [1/eV].
        double precision :: d
        ! Indices of the scattering electron's states.
        integer, intent(in) :: l_o_ba, l_o_pc, l_i_ba, l_i_pc
        ! Positive/negative energy of the absorbed/emitted phonon.
        double precision, intent(in) :: en_d
        double precision :: en, en_sq, eta_loc_ph, eta_loc_ph_sq
        ! Minimum allowed broadening [eV].
        double precision, parameter :: eta_min = 0.005

        ! Broadening.
        select case (broadening_ph_model)
            case (BPM_UCO)
                ! Uniform value.
                eta_loc_ph = broadening_ph_eV
            case (BPM_PDC)
                ! Process-dependent value.
                ! Normalize the joint level spacing using the target value.
                eta_loc_ph = broadening_ph_eV * jde_bz(l_o_ba,l_o_pc,l_i_ba,l_i_pc) / jde_avg
                ! Correct for very small broadening, to be on the safe side.
                eta_loc_ph = max(eta_loc_ph, eta_min)
        end select
        eta_loc_ph_sq = eta_loc_ph * eta_loc_ph
        ! Energy combination that nullifies at the delta peak.
        en = en_bz(l_o_ba, l_o_pc) - en_bz(l_i_ba, l_i_pc) - en_d
        en_sq = en * en
        ! Avoid very large values which produce numerical exceptions.
        if (en_sq .lt. 25.0*eta_loc_ph_sq) then
            d = exp(-en_sq/eta_loc_ph_sq)/pi_sqrt/eta_loc_ph
        else
            d = 0.0
        endif
    end function delta_eph


    subroutine calculate_EPC_matrix_elements
        implicit none
#if IEEECHECKS
        character(len=300) :: msg
#endif
        ! Indices.
        integer :: l_ph, l_q_pc, l_i_pc, l_o_pc, l_i_ba, l_o_ba, l_ae
        ! Integer projections of reduced wave vectors.
        integer, dimension(ndim) :: ii_ph_q, ii_va_i
        ! Matrix element squared.
        double precision :: g2
        ! Dirac delta function of energy conservation.
        double precision :: a0
        ! IPs of wave vectors.
        integer, dimension(ndim) :: ii_o_unfolded, ii_o

        ! Determine whether the sum of a phonon and electron wave vectors
        ! is in the electron support.  Of course, the initial wave vectors
        ! must be in their respective supports as well.
        is_sod_ph_pc = .false.
        l_sod_ph_pc = -1  ! invalid index assigned to forbidden processes
        ! Iterate over the phonon mode.
        do l_ph = 1, ph_n
        ! Iterate over the phonon wave vector.
        do l_q_pc = 1, pc_n
        if (is_ph_pc(l_ph,l_q_pc)) then
            ! Iterate over the electron wave vectors.
            do l_i_pc = 1, pc_n
            if (is_el_pc(l_i_pc)) then
                ! Iterate over absorption/emission processes.
                do l_ae = 1, ae_n
                ii_o_unfolded = ip_pc(:,l_i_pc) + s_ae(l_ae) * ip_pc(:,l_q_pc)
                ! Here we fold into the PC, but do not need reduced wave
                ! vectors, because we need the indices in the arrays of the
                ! distributions.  The reduced wave vectors are used only in the
                ! calculation of the EPC coefficient.
                call ip_fold(ii_o(1), ii_o(2), ii_o_unfolded)
                ! Index of the outgoing electron wave vector.
                l_o_pc = l_from_ip(ii_o(1), ii_o(2), "calculate_EPC_matrix_elements")
                ! Check that the wave vector is in the support.
                if (is_el_pc(l_o_pc)) then
                    is_sod_ph_pc(l_ae,l_ph,l_q_pc,l_i_pc) = .true.
                    l_sod_ph_pc(l_ae,l_ph,l_q_pc,l_i_pc) = l_o_pc
                end if
                end do  ! absorption/emission processes
            end if
            end do  ! electron wave vector
        end if
        end do  ! phonon wave vector
        end do  ! phonon mode

        ! Calculate the matrix elements.
        !-----------------------------------------------------------------------
        ! NOTE: Since this matrix is not small and the calculation of each
        ! element is fast, it is not clear whether it is more efficient to
        ! calculate the elements here once and for all or to recalculate them
        ! at each step when needed.  A test should be done.
        !-----------------------------------------------------------------------
        ! Phonon mode.
        do l_ph = 1, ph_n
        ! Incoming band.
        do l_i_ba = 1, ba_n
        ! Outgoing band.
        do l_o_ba = 1, ba_n
        ! Phonon wave vector.
        do l_q_pc = 1, pc_n
        if (is_ph_pc(l_ph,l_q_pc)) then
            ! Incoming electron wave vector.
            do l_i_pc = 1, pc_n
            if (is_el_pc(l_i_pc)) then
                ! Iterate over absorption/emission processes.
                do l_ae = 1, ae_n
                ! Check that the outgoing wave vector is in the support.
                if (is_sod_ph_pc(l_ae,l_ph,l_q_pc,l_i_pc)) then
                    ! The reduced phonon wave vector multiplied by a +/- sign.
                    ii_ph_q = s_ae(l_ae) * ip_ph_pc(:,l_ph,l_q_pc)
                    ! The reduced electron wave vector depends on its valley.
                    ii_va_i = ip_va_pc(:, va_pc(l_i_pc),l_i_pc)
                    ! The scattering coefficient needs the "reduced" phonon
                    ! and electron wave vectors.
                    g2 = g2_func(l_ph, ii_ph_q, l_o_ba, l_i_ba, ii_va_i)
                    ! The delta function of the energy conservation.
                    ! The sign of the phonon energy is +hw for +q
                    ! and -hw for -q, corresponding to
                    ! absorption: k -> k+q, e_(k+q) = e_k + hw
                    ! emission: k -> k-q, e_(k-q) = e_k - hw
                    l_o_pc = l_sod_ph_pc(l_ae,l_ph,l_q_pc,l_i_pc)
                    ! delta( e_(k +/- q) - e_k -/+ hw )
                    a0 = delta_eph(l_o_ba, l_o_pc, l_i_ba, l_i_pc, s_ae(l_ae)*en_ph_pc(l_ph,l_q_pc))
                    !-----------------------------------------------------------
                    g2delta_mat(l_ae, l_ph, l_q_pc, l_o_ba, l_i_ba, l_i_pc) = g2 * a0
#if IEEECHECKS
                    ! Check that values are not singular.
                    if (ieee_is_nan(g2)) then
                        write(msg, '(A,6I4)') "calculate_EPC_matrix_elements: EPC NaN", &
                            l_ae, l_ph, l_q_pc, l_o_ba, l_i_ba, l_i_pc
                        LOGME trim(msg_red(msg))
                        stop
                    end if
                    if (ieee_is_nan(a0)) then
                        write(msg, '(A,6I4)') "calculate_EPC_matrix_elements: Delta NaN", &
                            l_ae, l_ph, l_q_pc, l_o_ba, l_i_ba, l_i_pc
                        LOGME trim(msg_red(msg))
                        stop
                    end if
#endif
                end if
                end do  ! absorpion/emission processes
            end if
            end do  ! incoming electron wave vector
        end if
        end do  ! phonon wave vector
        end do  ! outgoing band
        end do  ! incoming band
        end do  ! phonon mode
    end subroutine calculate_EPC_matrix_elements

    subroutine calculate_pbe_gammas(f_pc, h_pc)
        implicit none
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        ! Indices of the phonon state we calculate the rates for.
        integer :: l_ph, l_q_pc
        ! Indices for the summations.
        integer :: l_i_pc, l_i_ba, l_o_ba, l_sp, l_ae
        ! Coefficients.
        double precision :: cx
        ! Variable where the sum is accumulated.
        double precision, dimension(ae_n) :: total_ae

        cx = 2.0*pi/hbar/pc_n
        gamma_abs = 0.0
        gamma_em = 0.0

!$OMP PARALLEL DO PRIVATE(l_ph, l_q_pc, total_ae, l_ae, l_i_ba, l_o_ba, l_sp, l_i_pc)
        ! Iterate over phonon mode.
        do l_ph = 1, ph_n
        ! Iterate over phonon wave vector.
        do l_q_pc = 1, pc_n
        if (is_ph_pc(l_ph,l_q_pc)) then
            ! Initialize the result of the summation.
            total_ae = 0.0
            ! Iterate over absorption/emission.
            do l_ae = 1, ae_n
            !-------------------------------------------------------
            ! Sum over incoming band.
            do l_i_ba = 1, ba_n
            ! Sum over outgoing band.
            do l_o_ba = 1, ba_n
            ! Sum over spin, conserved.
            do l_sp = 1, sp_n
            ! Sum over incoming electron wave vector.
            do l_i_pc = 1, pc_n
            if (is_el_pc(l_i_pc)) then
                ! Check that the outgoing wave vector is in the support.
                if (is_sod_ph_pc(l_ae,l_ph,l_q_pc,l_i_pc)) then
                    total_ae(l_ae) = total_ae(l_ae) &
                        + g2delta_mat(l_ae, l_ph, l_q_pc, l_o_ba, l_i_ba, l_i_pc) &
                        * f_pc(l_sp, l_i_ba, l_i_pc) &
                        * h_pc(l_sp, l_o_ba, l_sod_ph_pc(l_ae, l_ph, l_q_pc, l_i_pc))
                end if
            end if
            end do  ! incoming electron wave vector
            end do  ! spin
            end do  ! outgoing band
            end do  ! incoming band
            !-------------------------------------------------------
            end do  ! absorption/emission
            gamma_abs(l_ph, l_q_pc) = cx * total_ae(1)
            gamma_em(l_ph, l_q_pc) = cx * total_ae(2)
        end if
        end do  ! phonon wave vector
        end do  ! phonon mode
!$OMP END PARALLEL DO
    end subroutine calculate_pbe_gammas

    subroutine calculate_ebp_gammas(u_pc, f_pc, h_pc)
        implicit none
        ! Phonon distribution.
        double precision, dimension(:, :), intent(in) :: u_pc
        ! Electron and hole distribution.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        ! Indices of the electron state (k) we calculate the rates for.
        integer :: l_k_pc, l_k_ba, l_sp
        ! Indices of the phonon (q) taking part in the scattering, summed over.
        integer :: l_ph, l_q_pc
        ! The other electron state corresponds to 4 possible processes:
        ! - outgoing after absorbing a phonon (k+q)
        ! - outgoing after emitting a phonon (k-q)
        ! - incoming before absorbing a phonon (k-q)
        ! - incoming before emitting a phonon (k+q)
        ! Indices for the summations.
        integer :: l_o_pc, l_o_ba
        ! Notice that due to the exploitation of time-reversal symmetry and
        ! other properties of the electron-phonon coupling coefficients,
        ! the coefficients of the absorption and emission terms are exchanged
        ! in the "in" term.
        ! Coefficients.
        double precision :: cx
        ! Variables to perform the summation.
        double precision :: total_out, total_in

        cx = 2.0*pi/hbar/pc_n
        gamma_out = 0.0
        gamma_in = 0.0

!$OMP PARALLEL DO PRIVATE(l_k_ba, l_sp, l_k_pc, total_out, total_in, l_o_ba, l_ph, l_q_pc, l_o_pc)
        ! Iterate over incoming electron wave vector.
        do l_k_pc = 1, pc_n
        if (is_el_pc(l_k_pc)) then
            ! Iterate over electron band.
            do l_k_ba = 1, ba_n
            ! Iterate over spin, conserved.
            do l_sp = 1, sp_n
            ! Initialize the result of the summation.
            total_out = 0.0
            total_in = 0.0
            !-------------------------------------------------------------------
            ! Sum over other electron band.
            do l_o_ba = 1, ba_n
            ! Sum over phonon mode.
            do l_ph = 1, ph_n
            ! Sum over phonon wave vector.
            do l_q_pc = 1, pc_n
            if (is_ph_pc(l_ph,l_q_pc)) then
                ! Outgoing rate.
                ! absorption k -> k+q
                if (is_sod_ph_pc(1,l_ph,l_q_pc,l_k_pc)) then
                    l_o_pc = l_sod_ph_pc(1,l_ph,l_q_pc,l_k_pc)  ! k+q
                    total_out = total_out + &  ! g_(q) delta(e_(k+q) - e_k - hw)
                        g2delta_mat(1, l_ph, l_q_pc, l_o_ba, l_k_ba, l_k_pc) &
                        * h_pc(l_sp, l_o_ba, l_o_pc) &  ! 1-f(k+q)
                        * u_pc(l_ph, l_q_pc)  ! n
                end if
                ! emission k -> k-q
                if (is_sod_ph_pc(2,l_ph,l_q_pc,l_k_pc)) then
                    l_o_pc = l_sod_ph_pc(2,l_ph,l_q_pc,l_k_pc)  ! k-q
                    total_out = total_out + &  ! g_(-q) delta(e_(k-q) - e_k + hw)
                        g2delta_mat(2, l_ph, l_q_pc, l_o_ba, l_k_ba, l_k_pc) &
                        * h_pc(l_sp, l_o_ba, l_o_pc) &  ! 1-f(k-q)
                        * (u_pc(l_ph, l_q_pc) + 1)  ! n+1
                end if

                ! Incoming rate.
                ! absorption k-q -> k
                if (is_sod_ph_pc(2,l_ph,l_q_pc,l_k_pc)) then
                    l_o_pc = l_sod_ph_pc(2,l_ph,l_q_pc,l_k_pc)  ! k-q
                    total_in = total_in + &  ! g_(-q) delta(e_(k-q) - e_k + hw)
                        g2delta_mat(2, l_ph, l_q_pc, l_o_ba, l_k_ba, l_k_pc) &
                        * f_pc(l_sp, l_o_ba, l_o_pc) &  ! f(k-q)
                        * u_pc(l_ph, l_q_pc)  ! n
                end if
                ! emission k+q -> k
                if (is_sod_ph_pc(1,l_ph,l_q_pc,l_k_pc)) then
                    l_o_pc = l_sod_ph_pc(1,l_ph,l_q_pc,l_k_pc)  ! k+q
                    total_in = total_in + &  ! g_(q) delta(e_(k+q) - e_k - hw)
                        g2delta_mat(1, l_ph, l_q_pc, l_o_ba, l_k_ba, l_k_pc) &
                        * f_pc(l_sp, l_o_ba, l_o_pc) &  ! f(k+q)
                        * (u_pc(l_ph, l_q_pc) + 1)  ! n+1
                end if
            end if
            end do  ! phonon wave vector
            end do  ! phonon mode
            end do  ! outgoing mode
            !-------------------------------------------------------------------
            gamma_out(l_sp, l_k_ba, l_k_pc) = cx * total_out
            gamma_in(l_sp, l_k_ba, l_k_pc) = cx * total_in
            end do  ! spin
            end do  ! incoming band
        end if
        end do  ! incoming electron wave vector
!$OMP END PARALLEL DO

    end subroutine calculate_ebp_gammas

    subroutine calculate_pbe_rate(u1_pc, u_pc, f_pc, h_pc)
        ! Calculate the phonon scattering rate due to el-ph interactions.
        implicit none
        ! Time-derivative of the phonon distribution.
        double precision, dimension(:, :), intent(out) :: u1_pc
        ! Electron, hole and phonon distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        double precision, dimension(:, :), intent(in) :: u_pc

        ! Calculate distribution-dependent quantities.
        call calculate_pbe_gammas(f_pc, h_pc)

        ! Total scattering rate at given wave vector.
        ! Notice that not all the entries in u_pc, u1_pc are meaningful.
        ! We should only calculate the rate where the modes have support.
        ! However, using broadcasting is more efficient and it yields the same
        ! result as long as the gamma matrices are set to zero for those values
        ! which are outside the support.
        u1_pc = -u_pc*gamma_abs + (u_pc + 1)*gamma_em
    end subroutine calculate_pbe_rate

    subroutine calculate_ebp_rate(f1_pc, u_pc, f_pc, h_pc)
        ! Calculate the electron scattering rate due to el-ph interactions.
        implicit none
        ! Time-derivative of the electron distribution.
        double precision, dimension(:, :, :), intent(out) :: f1_pc
        ! Electron, hole and phonon distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        double precision, dimension(:, :), intent(in) :: u_pc

        ! Calculate distribution-dependent quantities.
        call calculate_ebp_gammas(u_pc, f_pc, h_pc)

        ! Total scattering rate at given wave vector.
        f1_pc = -f_pc*gamma_out + h_pc*gamma_in
    end subroutine calculate_ebp_rate

    subroutine print_ebp_gammas(baseName)
        implicit none
        character(len=*), intent(in) :: baseName
        integer :: l_sp, l_ba, l_pc, col_n
        character(len=50) :: fmtString, fileName

        col_n = ba_n*sp_n
        write (fmtString, '(A,I2,A)') "(", col_n - 1, "(F15.8, ', '),F15.8)"

        write(fileName, '(2A)') baseName, "-out.csv"
        open (unit=30, file=fileName, status="unknown")
        do l_pc = 1, pc_n
            write (30, fmtString) ((gamma_out(l_sp,l_ba,l_pc), l_ba=1,ba_n), l_sp=1,sp_n)
        end do
        close (unit=30)

        write(fileName, '(2A)') baseName, "-in.csv"
        open (unit=30, file=fileName, status="unknown")
        do l_pc = 1, pc_n
            write (30, fmtString) ((gamma_in(l_sp,l_ba,l_pc), l_ba=1,ba_n), l_sp=1,sp_n)
        end do
        close (unit=30)
    end subroutine print_ebp_gammas

end module ELECTRON_PHONON_SCATTERING
