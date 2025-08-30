module BAND_STRUCTURE
    use color_msg
    use phys_const
    use math
    implicit none
    ! Number of dimensions.
    integer, parameter :: ndim = 2
    ! Number of spin states.
    integer, parameter :: sp_n = 2
    ! Number of bands.
    ! Band indices increase from bottom to top: valence = 1, conduction = 2.
    integer, parameter :: ba_n = 2
    ! Number of valleys, useful to make contact to the continuous model.
    integer, parameter :: va_n = 2
    ! Carbon-carbon distance [nm].
    double precision, parameter :: a_CC = 0.142
    ! Connecting vectors in real space [nm].
    double precision, dimension(ndim), parameter :: delta1 = (/0.5*a_CC, 0.5*a_CC*s3/)
    double precision, dimension(ndim), parameter :: delta2 = (/0.5*a_CC, -0.5*a_CC*s3/)
    double precision, dimension(ndim), parameter :: delta3 = (/-1.0*a_CC, 0.0d0/)
    ! Reciprocal primitive vectors [nm^-1].
    double precision, dimension(ndim), parameter :: wv_pri1 = (/(2.0*pi)/(3.0*a_CC), (2.0*pi)/(s3*a_CC)/)
    double precision, dimension(ndim), parameter :: wv_pri2 = (/(2.0*pi)/(3.0*a_CC), -(2.0*pi)/(s3*a_CC)/)
    ! Mesh size along each reciprocal primitive vector.
    integer :: pri1_n
    ! Mesh step along each reciprocal primitive vector [nm^-1].
    double precision :: wv_dk1, wv_dk2
    ! Mesh size in the reciprocal primitive cell (PC).
    integer :: pc_n
    ! Each vector in the PC is of the form
    ! kk = (i_1 / pri1_n) * wv_pri1 + (i_2 / pri1_n) * wv_pri2
    ! Let us call i_1 and i_2 the integers projections (IPs) of kk.
    ! IPs in the PC mesh.
    integer, dimension(:, :), allocatable :: ip_pc
    ! Wave vectors in the PC mesh [nm^-1].
    ! The two components are x and y.
    double precision, dimension(:, :), allocatable :: wv_pc
    ! IPs in the Brillouin Zone (BZ) mesh.
    integer, dimension(:, :), allocatable :: ip_bz
    ! Wave vectors in the BZ mesh [nm^-1].
    ! The two components are x and y.
    double precision, dimension(:, :), allocatable :: wv_bz
    ! Length of the wave vectors in the BZ mesh [nm^-1].
    double precision, dimension(:), allocatable :: norm_bz
    ! IPs in the "symmetrized" BZ mesh, where wave vectors are summed to their
    ! equivalent replicas which are not in the BZ, but have the same distance
    ! from the Gamma point.
    integer, dimension(:, :), allocatable :: ip_sym_bz
    ! Wave vectors in the "symmetrized" BZ mesh [nm^-1].
    ! The two components are x and y.
    double precision, dimension(:, :), allocatable :: wv_sym_bz
    ! True if the wave vector is a symmetrized one.
    logical, dimension(:), allocatable :: is_sym_bz
    ! Band velocity.
    double precision, dimension(:, :, :), allocatable :: vel_bz
    ! Connection function f(k) = \sum_{\alpha} e^{i k \cdot \delta_\alpha}
    ! calculated on the BZ.
    complex(kind=8), dimension(:), allocatable :: fconn_bz
    ! Same thing but with a different phase, to make contact with
    ! T. Stauber et al., Phys. Rev. B 78, 085432 (2008).
    ! NOTE: This seems redundant, is it possible to simplify the formulas?
    complex(kind=8), dimension(:), allocatable :: fconn_stauber_bz
    ! Phase f(k)/|f(k)|.
    complex(kind=8), dimension(:), allocatable :: fconn_ph_bz
    ! Hopping amplitude [eV].
    double precision :: t_hopp
    ! Non-zero band energy at the Dirac point, for numerical reasons [eV].
    double precision :: en_dirac
    ! Phase function used to calculate the dispersion.
    complex(kind=8), dimension(:), allocatable :: ph_nn_bz
    ! Band dispersion on the BZ mesh [eV].
    double precision, dimension(:, :), allocatable :: en_bz
    ! Minimum and maximum energy [eV].
    double precision :: en_bz_min, en_bz_max
    ! Average level distance [eV].
    double precision :: de_avg
    ! Maximum level distance [eV].
    double precision :: de_max
    ! Gradient of the band dispersion on the BZ mesh [eV/nm].
    double precision, dimension(:,:,:), allocatable :: ge_bz
    ! Effective "joint" level spacing at a scattering event [eV].
    double precision, dimension(:,:,:,:), allocatable :: jde_bz
    double precision :: jde_avg
    ! Frequencies corresponding to optical transitions of interest [eV].
    integer :: ot_n
    double precision, dimension(:), allocatable :: en_ot

    ! Index of high-symmetry points in the PC (or BZ) mesh.
    integer :: l_pc_gamma
    integer, dimension(va_n) :: l_pc_kappa  ! K, K1
    integer :: l_pc_em  ! M point between K and K1
    ! Path in the BZ.
    integer :: l_pc_path_n    
    integer, dimension(:), allocatable :: l_pc_path

    ! IPs measured with respect to K or K' in the two valleys, respectively.
    integer, dimension(:, :, :), allocatable :: ip_va_pc
    ! Wave vectors measured in the two valleys [nm^-1].
    ! The two components are x and y.
    double precision, dimension(:, :, :), allocatable :: wv_va_pc
    ! Valley labels of the wave vectors.
    ! All states are labelled, even those away from the K and K1 points.
    ! Useful to make contact with the continuous model.
    integer, dimension(:), allocatable :: va_pc

#if 0
    ! Integer projections of high-symmetry points in the reciprocal space.
    ! Notice that here we include equivalent points which are not all included
    ! in the BZ.
    integer, dimension(ndim) :: ii_gamma
    integer, dimension(ndim,3,va_n):: ii_kappa_eqv
#endif

    ! Whether to define a BZ mesh which is symmetric with respect to the Gamma
    ! point.
    logical :: define_symmetric_bz

    ! Area of the BZ [nm^-2].
    double precision :: area_bz
    ! Inverse of the area of the system [nm^-2].
    double precision :: length_m_2

    ! Energy bins to calculate the density of states.
    ! Minimum number of levels per bin.
    integer :: min_level_num_per_bin
    ! Number of energy bins.
    integer :: eb_n
    ! Parameters for the energy bins.
    double precision :: eb_en_min, eb_en_max, eb_en_w
    ! Extremes of the energy bins.  The number is eb_n+1.
    double precision, dimension(:), allocatable :: en_eb
    ! Number of wave vectors in each energy bin.
    integer, dimension(:), allocatable :: lk_n_eb
    ! Indices of the wave vectors in each energy bin.
    integer, dimension(:, :, :), allocatable :: lk_pc_eb
    ! Index of the energy bin for each wave vector.
    integer, dimension(:, :), allocatable :: l_eb_pc
    ! Density of states.
    double precision, dimension(:), allocatable :: dos_eb

    ! For each level, the distance in energy to the closest other level.
    ! Useful to calculate the average level spacing.
    double precision, dimension(:,:), allocatable :: de_pc_ba

    private

    public :: ndim, ba_n, va_n, sp_n, pri1_n, pc_n, ot_n, &
        ip_pc, ip_bz, l_pc_gamma, l_pc_kappa, l_pc_em, l_pc_path_n, l_pc_path, &
#if 0
        ii_gamma, ii_kappa_eqv, &
#endif
        ip_va_pc, wv_va_pc, va_pc, &
        en_bz, t_hopp, en_bz_max, en_bz_min, de_avg, de_max, &
        fconn_stauber_bz, fconn_ph_bz, &
        vel_bz, ge_bz, jde_bz, jde_avg, en_ot, &
        wv_dk1, wv_dk2, &
        wv_pc, wv_bz, wv_sym_bz, wv_pri1, wv_pri2, area_bz, length_m_2, &
        eb_n, lk_n_eb, lk_pc_eb, eb_en_w, dos_eb, &
        de_pc_ba

    public :: l_from_ip, ip_fold, ip_fold_f, dsq_from_ip, wv_from_ip
    public :: init_band_structure

contains

    subroutine init_band_structure
        implicit none

        character(len=50) :: fmt_string
        integer :: i1, i2, l_pc, l_ba, l_eb, l_va, l_i, l_ot
        double precision, dimension(va_n) :: dsq_x_va
        ! User-friendly parameters for the input file.
        integer :: factor_wavevector_mesh
        double precision :: absorption_M_point_eV, pseudo_gap_eV
        namelist /band_structure_params/ factor_wavevector_mesh, &
            & define_symmetric_bz, min_level_num_per_bin, absorption_M_point_eV, &
            & pseudo_gap_eV
        namelist /band_structure_values/ pc_n, area_bz, length_m_2, en_bz_min, &
            & en_bz_max, de_avg, eb_n, t_hopp, wv_dk1, wv_dk2

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=band_structure_params)
        close (37)

        ! Set some parameters.
        pri1_n = 6*factor_wavevector_mesh
        pc_n = pri1_n*pri1_n
        wv_dk1 = sqrt(wv_pri1(1)**2 + wv_pri1(2)**2) / pri1_n
        wv_dk2 = sqrt(wv_pri2(1)**2 + wv_pri2(2)**2) / pri1_n
        en_dirac = pseudo_gap_eV / 2.0

        ! Derived quantities.
        area_bz = abs(wv_pri1(1)*wv_pri2(2) - wv_pri1(2)*wv_pri2(1))
        length_m_2 = (area_bz/pc_n)/((2*pi)*(2*pi))

        ! Initialize the integer projections of the wave vectors in the PC mesh.
        allocate (ip_pc(ndim, pc_n))
        do i1 = 0, pri1_n - 1
            do i2 = 0, pri1_n - 1
                ip_pc(:, l_from_ip(i1, i2)) = (/i1, i2/)
            end do
        end do

        ! Save the indices of a high-symmetry path in the BZ.
        l_pc_path_n = 1 + pri1_n/3 + pri1_n/6 + pri1_n/2
        allocate (l_pc_path(l_pc_path_n))
        l_i = 1
        ! Gamma - K
        do i1 = 0, pri1_n/3
            i2 = 2 * i1
            l_pc_path(l_i) = l_from_ip(i1,i2)
            l_i = l_i + 1
        end do
        ! K - M
        do i1 = pri1_n/3 + 1, pri1_n/2
            i2 = pri1_n - i1
            l_pc_path(l_i) = l_from_ip(i1,i2)
            l_i = l_i + 1
        end do
        ! M - Gamma
        do i1 = pri1_n/2 - 1, 0, -1
            i2 = i1
            l_pc_path(l_i) = l_from_ip(i1,i2)
            l_i = l_i + 1
        end do

        ! Initialize the wave vectors in the PC mesh.
        allocate (wv_pc(ndim, pc_n))
        do i1 = 0, pri1_n - 1
            do i2 = 0, pri1_n - 1
                ! wv_pc(:, l_from_ip(i1, i2)) = (i1/(pri1_n*1.0))*wv_pri1 + (i2/(pri1_n*1.0))*wv_pri2
                wv_pc(:, l_from_ip(i1, i2)) = wv_from_ip(i1,i2)
            end do
        end do

        ! High-simmetry point in the wave vector mesh.
        ! Notice that these are not orthogonal projections, but the components
        ! along the primitive vectors.
        l_pc_gamma = l_from_ip(0, 0)
        l_pc_kappa(1) = l_from_ip(2*pri1_n/3, pri1_n/3)
        l_pc_kappa(2) = l_from_ip(pri1_n/3, 2*pri1_n/3)
        l_pc_em = l_from_ip(pri1_n/2, pri1_n/2)

#if 0
        ! Integer projections of high-symmetry points, possibly not included
        ! in the wave vector mesh.
        ! Notice that these are not orthogonal projections, but the components
        ! along the primitive vectors.
        ii_gamma = (/ 0, 0 /)
        ii_kappa_eqv(:,1,1) = (/ 2*pri1_n/3, pri1_n/3 /)
        ii_kappa_eqv(:,1,2) = (/ pri1_n/3, 2*pri1_n/3 /)
        ii_kappa_eqv(:,2,1) = -ii_kappa_eqv(:,1,2)
        ii_kappa_eqv(:,2,2) = ii_kappa_eqv(:,1,1) + ii_kappa_eqv(:,2,1)
        ii_kappa_eqv(:,3,1) = -ii_kappa_eqv(:,2,2)
        ii_kappa_eqv(:,3,2) = -ii_kappa_eqv(:,1,1)
#endif

        ! Integer projections measured with respect to the K and K1 points.
        ! The difference has to be taken in the PC, not the BZ.
        allocate (ip_va_pc(ndim,va_n,pc_n))
        do l_va = 1, va_n
            do l_pc = 1, pc_n
                ip_va_pc(:,l_va,l_pc) = ip_pc(:,l_pc) - ip_pc(:,l_pc_kappa(l_va))
            end do
        end do

        ! Wave vectors measured with respect to the K or K1 point.
        allocate (wv_va_pc(ndim, va_n, pc_n))
        do l_va = 1, va_n
            do l_pc = 1, pc_n
                i1 = ip_va_pc(1,l_va,l_pc)
                i2 = ip_va_pc(2,l_va,l_pc)
                wv_va_pc(:,l_va,l_pc) = wv_from_ip(i1,i2)
            end do
        end do

        ! Assign a valley label to each wave vector, based on the closer
        ! high-symmetry point.
        allocate (va_pc(pc_n))
        do l_pc = 1, pc_n
            ! Quantity proportional to the norm squared of the integer
            ! projections measured with respect to the K and K1 points.
            do l_va = 1, va_n
                dsq_x_va(l_va) = dsq_from_ip(ip_va_pc(1,l_va,l_pc), ip_va_pc(2,l_va,l_pc))
            end do
            ! The valley index corresponds to the index of the minimum.
            va_pc(l_pc) = minloc(dsq_x_va,1)
        end do

        ! Initialize the BZ mesh.
        call init_brillouin_zone

        ! Function f(k).
        call init_connection_function

        ! Initialize the hopping amplitude of the tigh-binding model based on
        ! the value of the absorption at the M point.
        t_hopp = 1.0
        t_hopp = (absorption_M_point_eV / 2.0) / band_energy(2, wv_pc(:, l_pc_em))

        ! Initialize the band dispersion.
        allocate (en_bz(ba_n, pc_n))
        do l_pc = 1, pc_n
            do l_ba = 1, ba_n
                en_bz(l_ba, l_pc) = band_energy(l_ba, wv_pc(:, l_pc))
                ! Numerical correction due to the paucity of wave vectors around
                ! the Dirac point.
                en_bz(l_ba, l_pc) = sign(sqrt(en_bz(l_ba, l_pc)**2 + en_dirac**2), en_bz(l_ba, l_pc))
            end do
        end do

        ! Band extremes.
        en_bz_min = minval(reshape(en_bz, (/ba_n*pc_n/)))
        en_bz_max = maxval(reshape(en_bz, (/ba_n*pc_n/)))

        ! Calculate the average level distance.  This energy scale is used to
        ! define energy meshes in the rest of the program.
        allocate(de_pc_ba(pc_n, ba_n))
        call calculate_avg_level_distance

        ! Calculate the gradient of the band dispersion.
        ! This is used to calculate the local level spacing when evaluating
        ! scattering integrals.
        allocate(ge_bz(ndim, ba_n, pc_n))
        call calculate_band_gradient

        ! Effective "joint" level spacing at a scattering event [eV].
        ! Calculated as the different of the gradient time the mesh step.
        ! Used to estimate the broadening of delta functions in the scattering terms
        ! of the Boltzmann equation.
        allocate(jde_bz(ba_n, pc_n, ba_n, pc_n))
        call calculate_joint_level_spacing

        ! Bins in the energy space to calculate energy-resolved quantities.
        call init_energy_bins

        ! Density of states.
        call init_dos

        ! Initialize the band velocity.
        allocate (vel_bz(ndim, ba_n, pc_n))
        do l_pc = 1, pc_n
            ! Singular K and K1 points.
            if ((l_pc .eq. l_pc_kappa(1)) .or. (l_pc .eq. l_pc_kappa(2))) then
                vel_bz(:, :, l_pc) = 0.0
            else
                do l_ba = 1, ba_n
                    call calculate_band_velocity(vel_bz(:, l_ba, l_pc), l_ba, wv_pc(:, l_pc))
                end do
            end if
        end do

        ! Define a set of frequencies corresponding to optical transitions of
        ! interest.  This list of frequency might be used by a routine that
        ! calculates the optical conductivity.
        ! Graphene: interband transitions along the Kappa-M segment.
        ot_n = pri1_n/6 - 1
        allocate( en_ot(ot_n) )
        l_ot = 1
        do i1 = pri1_n/3 + 1, pri1_n/2 - 1
            i2 = pri1_n - i1
            en_ot(l_ot) = 2.0 * en_bz(2,l_from_ip(i1,i2))
            l_ot = l_ot + 1
        end do

        ! Save derived values.
        open (unit=VALUNIT, file=VALFILE, status="old", position="append")
        write (VALUNIT, nml=band_structure_values)
        close (VALUNIT)

        ! Print the wave vectors of a path in the BZ.
        call print_path_bz

        ! Save the arrays.
        ! Reciprocal primitive vectors.
        open (unit=30, file="out-geo-rec-lat-basis.csv", status="unknown")
        write (30, "(F15.8,', ',F15.8)") wv_pri1(:)
        write (30, "(F15.8,', ',F15.8)") wv_pri2(:)
        close (unit=30)
        ! PC mesh.
        write (fmt_string, "(A,I1,A)") "(", ndim - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-primitive-cell.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) wv_pc(:, l_pc)
        end do
        close (unit=30)
        ! Path between high-symmetry points.
        open (unit=30, file="out-geo-path.csv", status="unknown")
        write (30, "(I5)") l_pc_path
        close (unit=30)
        ! Valley labels.
        open (unit=30, file="out-geo-valley-labels.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, '(I1)') va_pc(l_pc)
        end do
        close (unit=30)
        ! Wave vectors of the K and K1 points.
        write (fmt_string, "(A,I1,A)") "(", ndim - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-wavevectors-K-K1.csv", status="unknown")
        write (30, fmt_string) wv_pc(:, l_pc_kappa(1))
        write (30, fmt_string) wv_pc(:, l_pc_kappa(2))
        close (unit=30)
        ! Wave vectors measured from the K or K1 point.
        write (fmt_string, "(A,I1,A)") "(", ndim - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-valley-K.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) wv_va_pc(:,1,l_pc)
        end do
        close (unit=30)
        open (unit=30, file="out-geo-valley-K1.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) wv_va_pc(:,2,l_pc)
        end do
        close (unit=30)
        ! BZ mesh.
        write (fmt_string, "(A,I1,A)") "(", ndim - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-brillouin-zone.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) wv_bz(:, l_pc)
        end do
        close (unit=30)
        ! Averaged wave vectors in the BZ mesh.
        write (fmt_string, "(A,I1,A)") "(", ndim - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-brillouin-zone-sym.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) wv_sym_bz(:, l_pc)
        end do
        close (unit=30)
        ! List of symmetrized vectors.
        open (unit=30, file="out-geo-brillouin-zone-is-sym.csv", status="unknown")
        do l_pc = 1, pc_n
            if (is_sym_bz(l_pc)) then
                write (30, "(I1)") 1
            else
                write (30, "(I1)") 0
            end if
        end do
        close (unit=30)
        ! Function f.  Modulus and phase.
        open (unit=30, file="out-geo-overlap.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, "(F15.8,', ', F15.8)") abs(fconn_bz(l_pc)), &
                & atan2(aimag(fconn_bz(l_pc)), real(fconn_bz(l_pc)))
        end do
        close (unit=30)
        ! Band dispersion.
        write (fmt_string, "(A,I1,A)") "(", ba_n - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-bands.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) en_bz(:, l_pc)
        end do
        close (unit=30)
        ! Gradient of the band dispersion.
        write (fmt_string, "(A,I1,A)") "(", ndim * ba_n - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-gradient.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) (ge_bz(:, l_ba, l_pc), l_ba=1, ba_n)
        end do
        close (unit=30)
        ! Band velocity.
        write (fmt_string, "(A,I1,A)") "(", ba_n - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-band_vel_x.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) (vel_bz(1, l_ba, l_pc), l_ba=1, ba_n)
        end do
        close (unit=30)
        write (fmt_string, "(A,I1,A)") "(", ba_n - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-band_vel_y.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) (vel_bz(2, l_ba, l_pc), l_ba=1, ba_n)
        end do
        close (unit=30)
        ! Energy bins.
        open (unit=30, file="out-geo-energy-bins.csv", status="unknown")
        do l_eb = 1, eb_n + 1
            write (30, "(F15.8)") en_eb(l_eb)
        end do
        close (unit=30)
        ! Energy bin for each state.
        write (fmt_string, "(A,I02,A)") "(", ba_n - 1, "(I03,', '),I03)"
        open (unit=30, file="out-geo-binning.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmt_string) l_eb_pc(:, l_pc)
        end do
        close (unit=30)
        ! Density of states.
        open (unit=30, file="out-geo-dos.csv", status="unknown")
        do l_eb = 1, eb_n
            write (30, "(F15.8)") dos_eb(l_eb)
        end do
        close (unit=30)

    end subroutine init_band_structure

    subroutine init_energy_bins
        ! Initialize the bins to calculate the energy density.
        implicit none
        logical :: found
        integer :: l_pc, l_eb, ib

        ! Set the minimum and maximum energy so that the whole band dispersion
        ! is included, with a little extra padding.
        eb_en_min = en_bz_min - 0.010
        eb_en_max = en_bz_max + 0.010

        ! Iterative procedure.  Start from a large number of bins and decrease
        ! it until there are at least (min_level_num_per_bin) values per bin.

        ! Start from a large numer of bins.  The number be a power of two plus
        ! one.
        eb_n = 2**8 + 1
        found = .false.

        do while (.not. found)
            ! The array has one element more, to make it easier to iterate over
            ! the bins, instead than over the bin extremes.
            allocate (en_eb(eb_n + 1))
            do l_eb = 1, eb_n + 1
                en_eb(l_eb) = eb_en_min + (eb_en_max - eb_en_min)*(l_eb - 1.0)/(1.0*eb_n)
            end do
            eb_en_w = en_eb(2) - en_eb(1)

            ! Initialize the list of wave vectors in each energy bin.
            allocate (lk_n_eb(eb_n))
            allocate (lk_pc_eb(2, pc_n, eb_n))
            allocate (l_eb_pc(ba_n, pc_n))
            ! For each wave vector and band, find the bin where it falls, add
            ! its index to the correponding column, and increment the number of
            ! elements in that column.
            lk_n_eb = 0
            lk_pc_eb = 0
            do ib = 1, ba_n
                do l_pc = 1, pc_n
                    do l_eb = 1, eb_n
                        if ((en_bz(ib, l_pc) .ge. en_eb(l_eb)) .and. (en_bz(ib, l_pc) .lt. en_eb(l_eb + 1))) then
                            lk_n_eb(l_eb) = lk_n_eb(l_eb) + 1
                            lk_pc_eb(:, lk_n_eb(l_eb), l_eb) = (/l_pc, ib/)
                            l_eb_pc(ib, l_pc) = l_eb
                            ! Each wave vector can be only in one bin, so skip
                            ! the others.
                            exit
                        end if
                    end do
                end do
            end do

            if (minval(lk_n_eb) .ge. min_level_num_per_bin) then
                ! The condition on the minimum number of levels is fulfilled.
                found = .true.
            else
                ! Prepare a new iteration with a different number of bins.
                deallocate (en_eb)
                deallocate (lk_n_eb)
                deallocate (lk_pc_eb)
                deallocate (l_eb_pc)
                eb_n = ((eb_n - 1)/2) + 1
                ! Stop if the number of bins is too small.
                if (eb_n .le. 5) then
                    LOGME trim(msg_red("init_energy_bins: Number of bins too small."))
                    stop
                end if
            end if
        end do

    end subroutine init_energy_bins

    subroutine init_dos
        implicit none
        ! The DOS is just lk_n_eb times a factor, but let us calculate it
        ! explicitly as a model for integrations in the wave vector space.
        integer :: l_eb, l_lk
        allocate (dos_eb(eb_n))
        dos_eb = 0.0
        ! Calculate density for each energy bin.
        do l_eb = 1, eb_n
            ! Integrate over the counter-image of the bin.
            do l_lk = 1, lk_n_eb(l_eb)
                dos_eb(l_eb) = dos_eb(l_eb) + 1.0
            end do
        end do
        dos_eb = dos_eb*length_m_2/eb_en_w
    end subroutine init_dos

    function dsq_from_ip(i1, i2) result (dsq)
        ! From the components to an integer proportional to the modulus squared.
        ! The basis vectors are proportional to (1, \pm \sqrt(3)).
        implicit none
        integer, intent(in) :: i1, i2
        integer :: dsq
        dsq = (i1 + i2)*(i1 + i2) + 3*(i1 - i2)*(i1 - i2)
    end function dsq_from_ip

    function wv_from_ip(i1, i2) result (wv)
        ! Return the x and y components [nm^-1] of the wave vector, given its
        ! integer projections, which are taken in the non-orthogonal basis.
        implicit none
        integer, intent(in) :: i1, i2
        double precision, dimension(ndim) :: wv
        wv = (i1/(pri1_n*1.0))*wv_pri1 + (i2/(pri1_n*1.0))*wv_pri2
    end function wv_from_ip

    subroutine init_brillouin_zone
        ! Initialize the BZ, by folding as appropriate the wave vectors of the
        ! PC.
        implicit none
        integer :: i1, i2, l0_cr, l1_cr, l_pc, l_sy
        ! Number of corners of the BZ.
        integer, parameter :: cr_n = 4
        ! Corners of the BZ.  The order decides which one "takes" a wave vector
        ! that is equally distant from two corners.
        integer, dimension(ndim, cr_n) :: ii_cr
        ! Distance vector.
        integer, dimension(ndim) :: ii_d
        ! Integer proportional to the square modulus of the relative position.
        integer, dimension(cr_n) :: dsq_cr
        ! Save which corner "takes" which point.
        integer, dimension(:), allocatable :: l_cr_pc
        ! Initialize the components of the corners of the BZ.
        ii_cr(:, 1) = (/0, 0/)
        ii_cr(:, 2) = (/pri1_n, 0/)
        ii_cr(:, 3) = (/0, pri1_n/)
        ii_cr(:, 4) = (/pri1_n, pri1_n/)
        ! Allocate the wave vectors in the WS cell.
        allocate (ip_bz(ndim, pc_n))
        allocate (wv_bz(ndim, pc_n))
        allocate (norm_bz(pc_n))
        allocate (l_cr_pc(pc_n))
        allocate (ip_sym_bz(ndim, pc_n))
        allocate (wv_sym_bz(ndim, pc_n))
        allocate (is_sym_bz(pc_n))
        ! Most of wave vectors need no symmetrization.
        is_sym_bz = .false.
        ! Calculate the components in the WS cell.
        ! Iterate over the wave vectors in the BZ mesh.
        do i1 = 0, pri1_n - 1
            do i2 = 0, pri1_n - 1
                ! Calculate the distance from the corners of the BZ.
                do l1_cr = 1, 4
                    ! Calculate the components of the distance vector.
                    ii_d = (/i1, i2/) - ii_cr(:, l1_cr)
                    ! From the components to an integer proportional to the
                    ! modulus.  The basis vectors are proportional to
                    ! (1, \pm \sqrt(3)).
                    dsq_cr(l1_cr) = dsq_from_ip(ii_d(1), ii_d(2))
                end do
                ! Find the corner with the minumum distance.
                ! In the first instance, resolve same distance with the order
                ! of the corners given above.
                l0_cr = cr_n
                do l1_cr = cr_n - 1, 1, -1
                    if (dsq_cr(l1_cr) .le. dsq_cr(l0_cr)) then
                        l0_cr = l1_cr
                    end if
                end do
                l_cr_pc(l_from_ip(i1, i2)) = l0_cr
                ! Symmetrize with respect to the Gamma point.
                if (define_symmetric_bz) then
                    ! Second instance, correct for points shared between corners
                    ! 1 and 4.
                    if ((l0_cr .eq. 1) .and. (dsq_cr(1) .eq. dsq_cr(4))) then
                        l0_cr = 4
                    end if
                    ! Finally, correct for the K point.
                    if ((dsq_cr(1) .eq. dsq_cr(4)) .and. (dsq_cr(1) .eq. dsq_cr(3))) then
                        l0_cr = 1
                    end if
                end if
                ! Fold the position by subtracting the components of the
                ! selected corner.
                l_pc = l_from_ip(i1, i2)
                ip_bz(:, l_pc) = ip_pc(:, l_pc) - ii_cr(:, l0_cr)
                ! Average of the possible foldings, i.e. corners at
                ! the same distance.
                ip_sym_bz(:, l_pc) = 0
                l_sy = 0
                do l1_cr = 1, cr_n
                    if (dsq_cr(l1_cr) .eq. dsq_cr(l0_cr)) then
                        ip_sym_bz(:, l_pc) = ip_sym_bz(:, l_pc) + ip_pc(:, l_pc) - ii_cr(:, l1_cr)
                        l_sy = l_sy + 1
                    end if
                end do
                ! If more than one corner is at the same distance, the wave
                ! vector has been symmetrized.
                if (l_sy .gt. 1) then
                    is_sym_bz(l_pc) = .true.
                end if
            end do
        end do
        ! Calculate the wave vectors in the WS cell.
        do i1 = 0, pri1_n - 1
            do i2 = 0, pri1_n - 1
                l_pc = l_from_ip(i1, i2)
                wv_bz(:, l_pc) = (ip_bz(1, l_pc)/(pri1_n*1.0))*wv_pri1 + (ip_bz(2, l_pc)/(pri1_n*1.0))*wv_pri2
            end do
        end do
        ! Calculate the wave vectors in the WS cell with symmetry.
        do i1 = 0, pri1_n - 1
            do i2 = 0, pri1_n - 1
                l_pc = l_from_ip(i1, i2)
                wv_sym_bz(:, l_pc) = (ip_sym_bz(1, l_pc)/(pri1_n*1.0))*wv_pri1 + (ip_sym_bz(2, l_pc)/(pri1_n*1.0))*wv_pri2
            end do
        end do
        ! Calculate the norm of the wave vectors.
        do l_pc = 1, pc_n
            norm_bz(l_pc) = sqrt(wv_bz(1, l_pc)**2 + wv_bz(2, l_pc)**2)
        end do
        ! Print the corners with respect to which each wave vector is folded.
        open (unit=30, file="out-geo-icc.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, "(I5)") l_cr_pc(l_pc)
        end do
        close (unit=30)
    end subroutine init_brillouin_zone

    subroutine init_connection_function
        implicit none
        integer :: l_pc
        double precision :: en
        double precision, dimension(ndim) :: kk, kk_a
        allocate (fconn_bz(pc_n))
        allocate (fconn_stauber_bz(pc_n))
        allocate (fconn_ph_bz(pc_n))
        do l_pc = 1, pc_n
            fconn_bz(l_pc) = band_phase_nn(wv_bz(:,l_pc))
            fconn_stauber_bz(l_pc) = fconn_bz(l_pc) * &
                & exp(cmplx(0.0, -real(dot_product(wv_bz(:,l_pc), delta3), kind=4)))
            if (abs(fconn_bz(l_pc)) .lt. 1.0e-9) then
                fconn_ph_bz(l_pc) = 0.0
            else
                fconn_ph_bz(l_pc) = fconn_bz(l_pc) / abs(fconn_bz(l_pc))
            end if
            ! Check with the analytical result.
            kk = wv_pc(:, l_pc)
            kk_a = kk*a_CC
            en = sqrt(3.0 + 2.0*cos(s3*kk_a(2)) + 4.0*cos(s3*kk_a(2)/2.0)*cos(3.0*kk_a(1)/2.0))
            if (abs(en - abs(fconn_bz(l_pc))) .gt. 1.0d-3) then
                LOGME trim(msg_red("init_connection_function: Discrepancy with analytical result."))
                stop
            end if
        end do
    end subroutine init_connection_function

    subroutine ip_fold(ii_folded_1, ii_folded_2, ii)
        ! Fold the integer projections of a wave vector into the PC.
        implicit none
        ! integer, dimension(ndim), intent(out) :: ii_folded
        integer, intent(out) :: ii_folded_1, ii_folded_2
        integer, dimension(ndim), intent(in) :: ii
        ii_folded_1 = modulo(ii(1), pri1_n)
        ii_folded_2 = modulo(ii(2), pri1_n)
    end subroutine ip_fold

    function ip_fold_f(ii) result( ii_folded )
        ! Fold the integer projections of a wave vector into the PC.
        implicit none
        integer, dimension(ndim) :: ii_folded
        integer, dimension(ndim), intent(in) :: ii
        ii_folded(1) = modulo(ii(1), pri1_n)
        ii_folded(2) = modulo(ii(2), pri1_n)
    end function ip_fold_f

    subroutine l_from_ip_log_error_and_stop(i1, i2, msg)
        ! Separate the logging routine from the main function, so these lines
        ! are called only in the case of error, because the routine l_from_ip
        ! is called many many times.
        implicit none
        integer, intent(in) :: i1, i2
        character(len=*), intent(in), optional :: msg
        character(len=100) :: errString
        print *, i1, pri1_n, i2, pri1_n
        if (present(msg)) then
            write(errString, '(A,2I3,2A)') "l_from_ip: Wrong indices:", i1, i2, ".  Message: ", msg
        else
            write(errString, '(A,2I3,A)') "l_from_ip: Wrong indices:", i1, i2, "."
        end if
        LOGME trim(msg_red(errString))
        stop
    end subroutine

    function l_from_ip(i1, i2, msg)
        ! Calculate the index of a wave vector in the PC mesh, given its integer
        ! projections.
        implicit none
        integer :: l_from_ip
        ! Integer components, going from 0 to the mesh size minus one.
        integer, intent(in) :: i1, i2
        character(len=*), intent(in), optional :: msg

        if (present(msg)) then
#if EXTRACHECKS
            if ((i1 .lt. 0) .or. (i1 .ge. pri1_n) .or. (i2 .lt. 0) .or. (i2 .ge. pri1_n)) then
                print *, i1, pri1_n, i2, pri1_n
                call l_from_ip_log_error_and_stop(i1, i2, msg)
            end if
#endif
        end if
        l_from_ip = (i2 + 1) + i1*pri1_n
    end function l_from_ip

    function ip_in_bz(l1_pc, ii1, l0_pc)
        ! Check if a doublet of integer projections is part of the BZ mesh.
        ! The function returns true/false, and if true the first variable
        ! returns the index of the vector in the BZ.
        ! Brute force search: this may take a long time.
        implicit none
        logical :: ip_in_bz
        integer, intent(out) :: l1_pc
        ! Start the search from this index.
        integer, intent(in) :: l0_pc
        integer, dimension(:), intent(in) :: ii1
        integer :: l2_pc
        integer, dimension(ndim) :: ii2

        l1_pc = -1
        ip_in_bz = .false.
        ! Iterate over all the remaining wave vectors.
        do l2_pc = l0_pc, pc_n
            ii2 = ip_bz(:, l2_pc)
            if ((ii2(1) .eq. ii1(1)) .and. (ii2(2) .eq. ii1(2))) then
                ! Wave vector found.
                l1_pc = l2_pc
                ip_in_bz = .true.
                ! No need to continue searching.
                exit
            end if
        end do
    end function ip_in_bz

    subroutine calculate_avg_level_distance
        ! Calculate the average level distance.
        implicit none
        integer :: l_ba, l1_pc, l2_pc, l2_min
        double precision :: de, de_min, de_tot

        de_tot = 0.0
        de_max = 0.0
        ! Analyze levels separately in the two bands.
        do l_ba = 1, ba_n
            ! For each level, find the closest level.  Exclude degeneracies.
            do l1_pc = 1, pc_n
                de_min = 100.0
                do l2_pc = 1, pc_n
                    de = abs(en_bz(l_ba, l1_pc) - en_bz(l_ba, l2_pc))
                    if ((de .lt. de_min) .and. (de .gt. 0.001)) then
                        de_min = de
                        l2_min = l2_pc
                    end if
                end do
                ! Add the level distance to the total, save the closest level.
                de_tot = de_tot + de_min
                de_pc_ba(l1_pc, l_ba) = de_min
                ! Update the value of the maximum level distance.
                ! Notice that this is the maximum value that the minimum level
                ! distance (except degeneracies) achieves in the band.
                ! The naming of variables might lead to confusion in the
                ! following line, but should be clear in the rest of the
                ! program.
                if (de_min .gt. de_max) then
                    de_max = de_min
                end if
            end do
        end do
        ! The average level distance is the total divided by the number of
        ! levels.
        de_avg = de_tot / (ba_n * pc_n)
    end subroutine calculate_avg_level_distance

    function i_prev(i) result(j)
        implicit none
        integer, intent(in) :: i
        integer :: j
        if (i .eq. 0) then
            j = pri1_n - 1
        else
            j = i - 1
        end if
    end function

    function i_next(i) result(j)
        implicit none
        integer, intent(in) :: i
        integer :: j
        if (i .eq. pri1_n - 1) then
            j = 0
        else
            j = i + 1
        end if
    end function


#if 1
    ! Calculate the gradient using the discretized energies on the mesh, or
    ! the analytical formula.

    subroutine calculate_band_gradient
        ! Calculate the gradient of the band dispersion.
        ! The two components of the gradient vector correspond to the x and y
        ! direction.
        implicit none
        integer :: i1c, i2c  ! center, previous, next
        integer :: ib, lcc_pc, lpc_pc, lnc_pc, lcp_pc, lcn_pc
        double precision :: e_pc, e_nc, e_cp, e_cn, de1, de2
        double precision :: b1x, b1y, b2x, b2y, J

        ! Versors of the reciprocal primitive vectors.
        b1x = wv_pri1(1) / sqrt(wv_pri1(1)**2 + wv_pri1(2)**2)
        b1y = wv_pri1(2) / sqrt(wv_pri1(1)**2 + wv_pri1(2)**2)
        b2x = wv_pri2(1) / sqrt(wv_pri2(1)**2 + wv_pri2(2)**2)
        b2y = wv_pri2(2) / sqrt(wv_pri2(1)**2 + wv_pri2(2)**2)
        ! Determinant of the Jacobian of the coordinate transformation.
        J = b1x * b2y - b1y * b2x

        do ib = 1, ba_n
            ! Iterate over all points in the mesh.
            do i1c = 0, pri1_n - 1
                do i2c = 0, pri1_n - 1
                    ! Indices of the points to calculate the derivatives.
                    lcc_pc = l_from_ip(i1c, i2c, "calculate_band_gradient")
                    lpc_pc = l_from_ip(i_prev(i1c), i2c, "calculate_band_gradient")
                    lnc_pc = l_from_ip(i_next(i1c), i2c, "calculate_band_gradient")
                    lcp_pc = l_from_ip(i1c, i_prev(i2c), "calculate_band_gradient")
                    lcn_pc = l_from_ip(i1c, i_next(i2c), "calculate_band_gradient")
                    ! Band energy derivatives along the reciprocal vectors.
                    e_pc = en_bz(ib, lpc_pc)
                    e_nc = en_bz(ib, lnc_pc)
                    e_cp = en_bz(ib, lcp_pc)
                    e_cn = en_bz(ib, lcn_pc)
                    de1 = (e_nc - e_pc) / 2.0 / wv_dk1
                    de2 = (e_cn - e_cp) / 2.0 / wv_dk2
                    ! Components of the gradient vector.
                    ge_bz(1, ib, lcc_pc) = (b2y / J) * de1 - (b1y / J) * de2
                    ge_bz(2, ib, lcc_pc) = -(b2x / J) * de1 + (b1x / J) * de2
                end do
            end do
        end do
    end subroutine calculate_band_gradient

#else

    subroutine calculate_band_gradient
        implicit none
        integer :: l_pc, l_ba
        double precision :: e_xp, e_xm, e_yp, e_ym, dex, dey
        double precision, dimension(ndim) :: k, kxp, kxm, kyp, kym
        ! Wave vector increment to calculate the derivative [nm^-1].
        double precision, parameter :: dk = 0.1

        do l_ba = 1, ba_n
            do l_pc = 1, pc_n
                k = wv_bz(:, l_pc)
                kxp = k + (/ dk , 0.0d0 /)
                kxm = k - (/ dk , 0.0d0 /)
                kyp = k + (/ 0.0d0 , dk /)
                kym = k - (/ 0.0d0 , dk /)
                e_xp = band_energy(l_ba, kxp)
                e_xm = band_energy(l_ba, kxm)
                e_yp = band_energy(l_ba, kyp)
                e_ym = band_energy(l_ba, kym)
                dex = (e_xp - e_xm) / 2.0 / dk
                dey = (e_yp - e_ym) / 2.0 / dk
                ge_bz(1, l_ba, l_pc) = dex
                ge_bz(2, l_ba, l_pc) = dey
            end do
        end do
    end subroutine calculate_band_gradient

#endif

    subroutine calculate_joint_level_spacing
        ! Used to calculate the broadening of a delta function with the
        ! Gaussian approximation and a process-dependent broadening.
        ! See J.R. Yates et al., Phys. Rev. B 75, 195121 (2007)
        ! DOI: 10.1103/PhysRevB.75.195121
        implicit none
        integer :: l1_ba, l1_pc, l2_ba, l2_pc, jde_num
        double precision :: wv_dk_avg
        double precision :: jde_min, jde_max, x
        namelist /joint_density_of_states/ jde_min, jde_max, jde_avg

        ! Variable for the histogram of the values.
        integer, parameter :: n_bins = 201
        integer :: i_bin
        integer, dimension(n_bins) :: jde_hist
        double precision :: bin_w
        double precision, dimension(n_bins) :: bins_l

        bin_w = 1.5 / (n_bins - 1.0)
        bins_l(1) = 0.0
        do i_bin = 2, n_bins
            bins_l(i_bin) = bins_l(i_bin - 1) + bin_w
        end do
        jde_hist = 0

        jde_min = 100.0
        jde_max = 0.0
        jde_avg = 0.0
        jde_num = 0
        wv_dk_avg = sqrt(wv_dk1**2 + wv_dk2**2)
        do l1_ba = 1, ba_n
            do l1_pc = 1, pc_n
                do l2_ba = 1, ba_n
                    do l2_pc = 1, pc_n
                        x = sqrt((ge_bz(1,l1_ba,l1_pc) - ge_bz(1,l2_ba,l2_pc))**2 &
                            + (ge_bz(2,l1_ba,l1_pc) - ge_bz(2,l2_ba,l2_pc))**2) * wv_dk_avg
                        jde_bz(l1_ba, l1_pc, l2_ba, l2_pc) = x
                        ! For the statistics, ignore the zeros when the states
                        ! are the same.
                        if ((l1_ba .ne. l2_ba) .or. (l1_pc .ne. l2_pc)) then
                            ! Statistics of the level spacing.
                            jde_max = max(jde_max, x)
                            jde_min = min(jde_min, x)
                            jde_num = jde_num + 1
                            jde_avg = jde_avg + (x - jde_avg) / jde_num
                            ! Bin the value in the correct bin.
                            ! Start from the rightmost bin and proceed to the left.
                            do i_bin = n_bins, 0, -1
                                if (i_bin .eq. 0) then
                                    ! If we reach this index, it means that
                                    ! we have not found a bin, that is,
                                    ! the value is negative, which is wrong.
                                    LOGME trim(msg_red("calculate_joint_level_spacing: Negative eta."))
                                    print *, x
                                    stop
                                else if (x .gt. bins_l(i_bin)) then
                                    jde_hist(i_bin) = jde_hist(i_bin) + 1
                                    exit
                                end if
                            end do
                        end if
                    end do
                end do
            end do
        end do
        ! Print statistical figures.
        open (unit=VALUNIT, file=VALFILE, status="old", position="append")
        write (VALUNIT, nml=joint_density_of_states)
        close (VALUNIT)
        ! Print histogram.
        open (unit=37, file="out-geo-jde-hist.csv", status="unknown")
        do i_bin = 1, n_bins
            write (37, '(F15.8, ", ", I10)') bins_l(i_bin), jde_hist(i_bin)
        end do
        close (37)
    end subroutine calculate_joint_level_spacing

    function band_phase_nn(kk) result (f)
        implicit none
        complex(kind=8) :: f
        double precision, dimension(:), intent(in) ::kk
        f = exp(cmplx(0.0, real(dot_product(kk, delta1), kind=4))) &
            & + exp(cmplx(0.0, real(dot_product(kk, delta2), kind=4))) &
            & + exp(cmplx(0.0, real(dot_product(kk, delta3), kind=4)))
    end function band_phase_nn

    function band_energy(l_ba, kk) result (en)
        ! Expression for the band dispersion.
        implicit none
        double precision :: en
        integer, intent(in) :: l_ba
        double precision, dimension(:), intent(in) :: kk
        if (l_ba .eq. 1) then
            ! Valence band.
            en = -t_hopp*abs(band_phase_nn(kk))
        else if (l_ba .eq. 2) then
            ! Conduction band.
            en = t_hopp*abs(band_phase_nn(kk))
        end if
    end function band_energy

    subroutine calculate_band_velocity(vv, l_ba, kk)
        ! Calculate the band velocity.
        implicit none
        double precision, dimension(:), intent(out) :: vv
        integer, intent(in) :: l_ba
        double precision, dimension(:), intent(in) :: kk
        double precision, dimension(ndim) :: kk_a
        double precision :: r
        kk_a = kk*a_CC
        r = sqrt(3.0 + 2.0*cos(s3*kk_a(2)) + 4.0*cos(s3*kk_a(2)/2.0)*cos(3.0*kk_a(1)/2.0))
        vv(1) = 1.0/hbar*3.0*t_hopp*a_CC*cos(s3*kk_a(2)/2.0)*sin(3.0*kk_a(1)/2.0)
        vv(2) = 1.0/hbar*s3*t_hopp*a_CC*(sin(s3*kk_a(2)/2.0)*cos(3.0*kk_a(1)/2.0) + sin(s3*kk_a(2)))
        if (l_ba .eq. 1) then
            ! Valence band.
            vv = vv/r
        else if (l_ba .eq. 2) then
            ! Conduction band.
            vv = -vv/r
        end if
    end subroutine calculate_band_velocity

    subroutine print_path_bz
        ! Print the indices of a path in the Wigner-Seitz cell.
        ! The path is: Gamma - M.
        integer :: i1, i2

        open (unit=30, file="out-geo-path.csv", status="unknown")
        do i1 = 0, pri1_n/2
            i2 = 0
            write (30, "(I5)") l_from_ip(i1, i2, "print_path_bz")
        end do
        close (unit=30)
    end subroutine print_path_bz

end module BAND_STRUCTURE

module ELECTRON_THERMODYNAMICS
    use band_structure
    use phys_const
    use color_msg
    implicit none

    ! Extremes of the chemical potential used in the bisection algorithms.
    double precision, parameter :: mu_min = -4.0
    double precision, parameter :: mu_max = 4.0
    ! Required precision in the calculation of the chemical potential [eV].
    double precision, parameter :: mu_prec = 0.001
    ! Maximum number of bisections to calculate the chemical potential.
    integer, parameter :: bis_max = 100

    ! States contribute positive carrier density (i.e. electrons) if the
    ! band energy is positive and negative carrier density (i.e. holes) if the
    ! band energy is negative.

    ! Note that distributions are assumed to be SPIN-BALANCED and the
    ! band dispersion is spin-independent.
    
    ! Note that the bisection algorithm to calculate the chemical potential
    ! assumes that the chemical potential INCREASES with the carrier density,
    ! d mu / d n_C > 0.

    private

    public :: carrier_density_one_band, chemical_potential_one_band, &
        & carrier_density_all_bands, chemical_potential_all_bands

    public :: electron_thermodynamics_test_1

contains

    function carrier_density_one_band(mu, tempK, l_ba) result( nC )
        implicit none
        ! Calculate the carrier density in one band, given the
        ! electron chemical potentials and temperature.
        ! The smaller one in absolute value is used.
        ! Carrier density [nm^-2].
        double precision :: nC
        ! Electron chemical potential [eV].
        double precision, intent(in) :: mu
        ! Electron temperature [K].
        double precision, intent(in) :: tempK
        ! Band index.
        integer, intent(in) :: l_ba

        integer :: l_pc
        double precision :: temp
        double precision, dimension(:), allocatable :: f_pc

        allocate( f_pc(pc_n) )
        temp = k_boltz * tempK

        ! Define the electron distribution function.
        do l_pc = 1, pc_n
            f_pc(l_pc) = fermi_dirac(en_bz(l_ba,l_pc), mu, temp)
        end do
        ! Sum electron density and subtract hole density.
        nC = 0.0
        do l_pc = 1, pc_n
            if (en_bz(l_ba,l_pc) .gt. 0.0) then
                nC = nC + f_pc(l_pc)
            else
                nC = nC - (1.0 - f_pc(l_pc))
            end if
        end do
        ! Multiply by the integration element and the spin degeneracy.
        nC = nC * length_m_2 * sp_n
        
        deallocate(f_pc)
    end function carrier_density_one_band

    function chemical_potential_one_band(nC, tempK, l_ba, msg) result( mu )
        implicit none
        ! Calculate the electron chemical potential in one band, given the
        ! carrier density and the temperature.
        ! Electron chemical potential [eV].
        double precision :: mu
        ! Carrier density [nm^-2].
        double precision, intent(in) :: nC
        ! Electron temperature [K].
        double precision, intent(in) :: tempK
        ! Band index.
        integer, intent(in) :: l_ba
        ! Message from the calling routine.
        character(len=*), intent(in) :: msg
        ! Work quantities for the bisection algorithm.
        integer :: bis_count
        double precision :: mu_left, mu_right, mu_mid
        double precision :: n_left, n_right, n_mid

        ! Apply the bisection until convergence is achieved or the maximum
        ! number of iterations is reached.
        n_left = carrier_density_one_band(mu_min, tempK, l_ba)
        n_right = carrier_density_one_band(mu_max, tempK, l_ba)
        if ((n_left .gt. nC) .or. (n_right .lt. nC)) then
            LOGME trim(msg_red("chemical_potential_one_band: Value outside boundaries."))
            LOGME trim(msg_red(msg))
            write(*, '(3F8.3)') nC, n_left, n_right
            stop
        end if
        mu_left = mu_min
        mu_right = mu_max
        bis_count = 0
        !!! write(*, '(A,F8.3)') "DEBUG: chemical_potential_one_band for density", nC
        do while (((mu_right - mu_left) .gt. mu_prec) .and. (bis_count .lt. bis_max))            
            mu_mid = 0.5 * (mu_left + mu_right)
            n_mid = carrier_density_one_band(mu_mid, tempK, l_ba)
            !!! write(*, '(I3,3F8.3)') bis_count, mu_left, mu_right, n_mid  ! DEBUG
            if (n_mid .lt. nC) then
                mu_left = mu_mid
            else
                mu_right = mu_mid
            end if
            bis_count = bis_count + 1
        end do
        if ((mu_right - mu_left) .gt. mu_prec) then
            LOGME trim(msg_red("chemical_potential_one_band: Bisection did not converge."))
            LOGME trim(msg_red(msg))
            stop
        else
            ! We take the mid point of the final bracket.
            mu = 0.5 * (mu_left + mu_right)
        end if
    end function chemical_potential_one_band

    function carrier_density_all_bands(mu, tempK) result( nC )
        implicit none
        ! Calculate the total carrier density, given the electron chemical
        ! potential and temperature.
        ! Carrier density [nm^-2].
        double precision :: nC
        ! Electron chemical potential [eV].
        double precision, intent(in) :: mu
        ! Electron temperature [K].
        double precision, intent(in) :: tempK

        integer :: l_ba
        double precision :: nC_ba

        ! Sum over all the bands.
        nC = 0.0
        do l_ba = 1, ba_n
            nC_ba = carrier_density_one_band(mu, tempK, l_ba)
            nC = nC + nC_ba
        end do
    end function carrier_density_all_bands

    function chemical_potential_all_bands(nC, tempK) result( mu )
        ! Calculate the electron chemical potential common to all bands,
        ! given the carrier density and the temperature.
        ! Electron chemical potential [eV].
        double precision :: mu
        ! Carrier density [nm^-2].
        double precision, intent(in) :: nC
        ! Electron temperature [K].
        double precision, intent(in) :: tempK

        ! Work quantities for the bisection algorithm.
        integer :: bis_count
        double precision :: mu_left, mu_right, mu_mid
        double precision :: n_left, n_right, n_mid

        ! Apply the bisection until convergence is achieved or maximum number
        ! of iterations is reached.
        n_left = carrier_density_all_bands(mu_min, tempK)
        n_right = carrier_density_all_bands(mu_max, tempK)
        if ((n_left .gt. nC) .or. (n_right .lt. nC)) then
            LOGME trim(msg_red("chemical_potential_all_bands: Value outside boundaries."))
            stop
        end if
        mu_left = mu_min
        mu_right = mu_max
        bis_count = 0
        !!! write(*, '(A,F8.3)') "DEBUG: chemical_potential_all_bands for density", nC
        do while (((mu_right - mu_left) .gt. mu_prec) .and. (bis_count .lt. bis_max))            
            mu_mid = 0.5 * (mu_left + mu_right)
            n_mid = carrier_density_all_bands(mu_mid, tempK)
            !!! write(*, '(I3,3F8.3)') bis_count, mu_left, mu_right, n_mid  ! DEBUG
            if (n_mid .lt. nC) then
                mu_left = mu_mid
            else
                mu_right = mu_mid
            end if
            bis_count = bis_count + 1
        end do
        if ((mu_right - mu_left) .gt. mu_prec) then
            LOGME trim(msg_red("chemical_potential_all_bands: Bisection did not converge."))
            stop
        else
            ! We take the mid point of the final bracket.
            mu = 0.5 * (mu_left + mu_right)
        end if
    end function chemical_potential_all_bands

    subroutine electron_thermodynamics_test_1
        implicit none

        integer :: i_eF
        integer, parameter :: eF_n = 31
        double precision :: eF, eFmin, eFmax, deF, tK_0, tK_r, tK_h, nC, mu_r, mu_h
        double precision, dimension(:,:), allocatable :: res
        
        tK_0 = 77.0  ! "almost" zero temperature
        tK_r = 300.0  ! room temperature
        tK_h = 2000.0  ! hot temperature
        eFmin = 0.100
        eFmax = 3.0
        allocate(res(eF_n, 4))

        ! Calculate carrier density and chemical potential on a mesh of
        ! Fermi energies.
        deF = (log10(eFmax) - log10(eFmin))/(eF_n - 1.0)
        do i_eF = 1, eF_n
            eF = 10.0**(log10(eFmin) + (i_eF - 1.0) * deF)
            nC = carrier_density_all_bands(eF, tK_0)
            mu_r = chemical_potential_all_bands(nC, tK_r)
            mu_h = chemical_potential_all_bands(nC, tK_h)
            res(i_eF,:) = (/ eF, nC, mu_r, mu_h /)
        end do

        open (unit=30, file="out-thrm-el-1.csv", status="unknown")
        do i_eF = 1, eF_n
            write (30, '(3(F15.8, ", "),F15.8)') res(i_eF, :)
        end do
        close (unit=30)

        deallocate(res)
    end subroutine electron_thermodynamics_test_1
end module ELECTRON_THERMODYNAMICS

module ELECTRON_DISTRIBUTION
    use band_structure
    use phys_const
    use electron_thermodynamics
    use color_msg
    implicit none

    ! Distribution functions.  The first index correspond to the quantum labels
    ! (e.g. spin), the second the band index, and the third to the wave vector.
    !                    spin : band : wave vector
    double precision, dimension(:, :, :), allocatable :: f_i_pc, h_i_pc

    ! "Almost" zero temperature to convert between Fermi energy and
    ! carrier density.
    double precision, parameter :: tempK_0 = 77.0

    ! Parameters before photoexcitation.
    ! These are needed to calculate e.g. the differential transmission.
    ! Carrier density [nm^-2].
    double precision :: nC_eq
    ! Temperature [K].
    double precision :: tempK_eq, temp_eq
    ! Chemical potential at equilibrium [eV].
    double precision :: mu_eq
    ! Electron and hole distributions at equilibrium.
    double precision, dimension(:, :, :), allocatable :: f_eq_pc, h_eq_pc

    private :: tempK_0, nC_eq

contains

    subroutine init_electron_distribution
        implicit none

        integer :: l_pc, l_ba, l_sp
        ! User-friendly parameters for the input file.
        integer :: initial_state
        ! Values at equilibrium before photoexcitation.
        double precision :: carrier_density_eq_cm
        double precision :: nC_eq
        double precision :: fermi_energy
        double precision :: temperature_eq_K
        namelist /electron_distribution_params/ initial_state, &
            & carrier_density_eq_cm, temperature_eq_K
        namelist /initial_distribution_values/ fermi_energy, mu_eq

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=electron_distribution_params)
        close (37)

        tempK_eq = temperature_eq_K
        temp_eq = k_boltz * tempK_eq
        ! Convert carrier density from cm^-2 to nm^-2.
        nC_eq = carrier_density_eq_cm * 1.0e-12 * 0.01
        ! Chemical potential at the equilibrium temperature.
        mu_eq = chemical_potential_all_bands(nC_eq, tempK_eq)
        ! Calculate the Fermi energy as the chemical potential at "almost"
        ! zero temperature.
        fermi_energy = chemical_potential_all_bands(nC_eq, tempK_0)

        ! Distribution at equilibrium, before photoexcitation.
        allocate (f_eq_pc(sp_n, ba_n, pc_n))
        allocate (h_eq_pc(sp_n, ba_n, pc_n))
        
        ! Calculate the equilibrium distributions.
        do l_pc = 1, pc_n
            do l_ba = 1, ba_n
                do l_sp = 1, sp_n
                    f_eq_pc(l_sp, l_ba, l_pc) = fermi_dirac(en_bz(l_ba, l_pc), mu_eq, temp_eq)
                end do
            end do
        end do
        h_eq_pc = 1.0 - f_eq_pc

        ! Distributions immediately after photoexcitation.
        allocate (f_i_pc(sp_n, ba_n, pc_n))
        allocate (h_i_pc(sp_n, ba_n, pc_n))

        ! Initialize the distributions.
        if (initial_state .eq. 1) then
            call init_fermi_dirac
        else if (initial_state .eq. 2) then
            call init_hot_electrons
        else if (initial_state .eq. 3) then
            call init_photoexcited
        else
            LOGME trim(msg_red("init_electron_distribution: Initial state unknown."))
            stop
        end if

        ! Save derived values.
        open (unit=VALUNIT, file=VALFILE, status="old", position="append")
        write (VALUNIT, nml=initial_distribution_values)
        close (VALUNIT)

    end subroutine init_electron_distribution

    subroutine init_fermi_dirac
        implicit none

        integer :: l_pc, l_ba, l_sp
        double precision :: mu, temp
        ! User-friendly parameters for the input file.
        double precision :: temperature_K
        namelist /distro_fermi_dirac_params/ temperature_K

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=distro_fermi_dirac_params)
        close (37)

        LOGME "init_fermi_dirac: Electron state initialized as Fermi-Dirac."

        ! Set parameters.
        ! Assume no spin polarization.
        mu = chemical_potential_all_bands(nC_eq, temperature_K)
        temp = k_boltz * temperature_K

        do l_pc = 1, pc_n
            do l_ba = 1, ba_n
                do l_sp = 1, sp_n
                    f_i_pc(l_sp, l_ba, l_pc) = fermi_dirac(en_bz(l_ba, l_pc), mu, temp)
                end do
            end do
        end do
        h_i_pc = 1.0 - f_i_pc
    end subroutine init_fermi_dirac

    subroutine init_hot_electrons
        implicit none

        integer :: l_pc, l_ba, l_sp
        integer, dimension(ba_n) :: ns_ba  ! sign of the carriers
        double precision :: nE_exc, temp
        double precision, dimension(ba_n) :: nC_ba, mu_ba
        character(len=50) :: msg
        ! User-friendly parameters for the input file.
        double precision :: photoexc_density_cm
        double precision :: temperature_K
        namelist /distro_hot_electrons_params/ photoexc_density_cm, &
            & temperature_K

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=distro_hot_electrons_params)
        close (37)

        LOGME "init_hot_electrons: Electron state initialized as a hot Fermi-Dirac."

        ! Convert photoexcited density from 10^12cm^-2 to nm^-2.
        nE_exc = 0.01 * photoexc_density_cm
        ! Photoexcitation goes from valence to conduction.
        ! This could be generalized to an input parameter if needed.
        ns_ba = (/ -1, 1 /)
        temp = k_boltz * temperature_K
        
        ! Set parameters.
        ! Carrier density in each band, which is the sum of the equilibrium
        ! density (at equilibrium temperature) and the photoexcited density.
        ! Assume no spin polarization.
        do l_ba = 1, ba_n
            nC_ba(l_ba) = ns_ba(l_ba) * nE_exc + carrier_density_one_band(mu_eq, tempK_eq, l_ba)
        end do
        ! Chemical potential in each band.
        ! In each band, the density is now kept constant but the temperature has
        ! its photoexcited value.
        do l_ba = 1, ba_n
            write(msg, '(A,I3)') "init_hot_electrons, band ", l_ba
            mu_ba(l_ba) = chemical_potential_one_band(nC_ba(l_ba), temperature_K, l_ba, msg)
        end do

        ! Initialize a Fermi-Dirac in each band with appropriate chemical
        ! potential and common temperature.
        do l_pc = 1, pc_n
            do l_sp = 1, sp_n
                do l_ba = 1, ba_n
                    f_i_pc(l_sp, l_ba, l_pc) = fermi_dirac(en_bz(l_ba, l_pc), mu_ba(l_ba), temp)
                end do
            end do
        end do

        h_i_pc = 1.0 - f_i_pc
    end subroutine init_hot_electrons

    subroutine init_photoexcited
        implicit none

        integer :: l_pc, l_ba, l_sp
        double precision :: temp, en0, de, nE_exc, nE_gauss, fl_exc, fl_gauss
        double precision, dimension(ba_n, pc_n) :: df_pm
        ! User-friendly parameters for the input file.
        logical :: use_fluence
        double precision :: photoexc_density_cm
        double precision :: pump_wavelength_nm
        double precision :: pump_broadening_eV
        double precision :: abs_fluence_uJ_cm
        namelist /distro_photoexcited_params/ pump_wavelength_nm, pump_broadening_eV, &
            & photoexc_density_cm, use_fluence, abs_fluence_uJ_cm

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=distro_photoexcited_params)
        close (37)

        LOGME "init_photoexcited: Electron state initialized as a photoexcited state."

        ! Set parameters.
        temp = k_boltz * tempK_eq
        en0 = 0.5 * hc1240 / pump_wavelength_nm
        de = 0.5 * pump_broadening_eV
        nE_exc = 0.01 * photoexc_density_cm  ! convert 10^12 cm^-2 to nm^-2
        fl_exc =  1.0d-14 * Joule_eV * 1.0d-6 * abs_fluence_uJ_cm  ! convert (uJ cm^-2) to (eV nm^-2)

        ! Produce the functional form of two gaussian distributions,
        ! in valence and conduction band, and then normalize them
        ! to have the desired integral.
        
        ! Gaussian in conduction band.
        l_ba = 2
        do l_pc = 1, pc_n
            df_pm(l_ba, l_pc) = exp(-((en_bz(l_ba, l_pc) - en0)/de)**2)
        end do
        ! Gaussian in valence band.
        l_ba = 1
        do l_pc = 1, pc_n
            df_pm(l_ba, l_pc) = -exp(-((en_bz(l_ba, l_pc) + en0)/de)**2)
        end do

        if (use_fluence) then
            ! Calculate the absorbed fluence corresponding to one gaussian.
            fl_gauss = 0.0
            l_ba = 2
            do l_pc = 1, pc_n
                ! Transitions are vertical, fluence is twice the energy density
                ! of the photoexcited electrons.
                fl_gauss = fl_gauss + df_pm(l_ba, l_pc) * en_bz(l_ba, l_pc) * 2.0
            end do
            fl_gauss = fl_gauss * length_m_2
            ! Normalize the gaussians to have the desired fluence.
            ! Assume no spin polarization.
            df_pm = df_pm  * (fl_exc / sp_n / fl_gauss)
        else
            ! Calculate the density corresponding to one gaussian.
            nE_gauss = 0.0
            l_ba = 2
            do l_pc = 1, pc_n
                nE_gauss = nE_gauss + df_pm(l_ba, l_pc)
            end do
            nE_gauss = nE_gauss * length_m_2
            ! Normalize the gaussians to have the desired density.
            ! Assume no spin polarization.
            df_pm = df_pm  * (nE_exc / sp_n / nE_gauss) 
        end if

        ! Equilibrium distribution.
        do l_pc = 1, pc_n
            do l_ba = 1, ba_n
                do l_sp = 1, sp_n
                    f_i_pc(l_sp, l_ba, l_pc) = fermi_dirac(en_bz(l_ba, l_pc), mu_eq, temp)
                end do
            end do
        end do

        ! Add the gaussians to the equilibrium distribution.
        do l_pc = 1, pc_n
            do l_ba = 1, ba_n
                do l_sp = 1, sp_n
                    f_i_pc(l_sp, l_ba, l_pc) = f_i_pc(l_sp, l_ba, l_pc) + df_pm(l_ba, l_pc)
                    ! Check that the distribution remains smaller than unity and positive.
                    if (f_i_pc(l_sp, l_ba, l_pc) .gt. 1.0) then
                        LOGME trim(msg_red("init_photoexcited: Distribution larger than one."))
                        stop
                    end if
                    if (f_i_pc(l_sp, l_ba, l_pc) .lt. 0.0) then
                        LOGME trim(msg_red("init_photoexcited: Distribution negative."))
                        stop
                    end if
                end do
            end do
        end do

        h_i_pc = 1.0 - f_i_pc
    end subroutine init_photoexcited

    subroutine calculate_electron_density(rf, rh, f_pc, h_pc)
        ! Calculate spin- and band-resolved electron and hole densities.
        implicit none
        ! Particle density [nm^-2].
        double precision, dimension(:, :), intent(out) :: rf, rh
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc

        ! Sum over the index in the BZ zone, resolve spin and band indices.
        rf = sum(f_pc, 3)*length_m_2
        rh = sum(h_pc, 3)*length_m_2
    end subroutine calculate_electron_density

    subroutine calculate_energy_density(re, f_pc, h_pc)
        ! Calculate energy density.
        implicit none
        ! Energy density [eV nm^-2].
        double precision, intent(out) :: re
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        integer :: l_sp, l_ba, l_pc

        re = 0.0
        do l_sp = 1, sp_n
            do l_pc = 1, pc_n
                ! Hole energy in valence band.
                l_ba = 1
                re = re + h_pc(l_sp, l_ba, l_pc)*(-en_bz(l_ba, l_pc))
                ! Electron energy in conduction band.
                l_ba = 2
                re = re + f_pc(l_sp, l_ba, l_pc)*en_bz(l_ba, l_pc)
            end do
        end do
        re = re*length_m_2
    end subroutine calculate_energy_density

    subroutine calculate_current_density(rj, f_pc)
        ! Calculate spin- and band-resolved particle current density.
        implicit none
        ! Particle current density [fs^-1 nm^-1].
        double precision, dimension(:, :, :), intent(out) :: rj
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc
        integer :: l_sp, l_ba, l_pc

        rj = 0.0
        do l_sp = 1, sp_n
            do l_ba = 1, ba_n
                do l_pc = 1, pc_n
                    if (l_pc .ne. l_pc_gamma) then
                        rj(:, l_sp, l_ba) = rj(:, l_sp, l_ba) + f_pc(l_sp, l_ba, l_pc)*vel_bz(:, l_ba, l_pc)
                    end if
                end do
            end do
        end do
        rj = rj*length_m_2
    end subroutine calculate_current_density

    subroutine calculate_momentum_density(rp, f_pc)
        ! Calculate spin- and band-resolved momentum density.
        implicit none
        ! Momentum density [eV fs nm^-1].
        double precision, dimension(:, :, :), intent(out) :: rp
        ! Electron distribution.
        double precision, dimension(:, :, :), intent(in) :: f_pc
        integer :: l_sp, l_ba, l_pc

        rp = 0.0
        do l_sp = 1, sp_n
            do l_ba = 1, ba_n
                do l_pc = 1, pc_n
                    if (l_pc .ne. l_pc_gamma) then
                        rp(:, l_sp, l_ba) = rp(:, l_sp, l_ba) + f_pc(l_sp, l_ba, l_pc)*wv_sym_bz(:, l_pc)*hbar
                    end if
                end do
            end do
        end do
        rp = rp*length_m_2
    end subroutine calculate_momentum_density

    subroutine calculate_el_distro_energy(rf_eb, rh_eb, f_pc, h_pc)
        ! Calculate spin-resolved particle distribution vs energy.
        implicit none
        ! Particle energy density [nm^-2 eV^-1].
        double precision, dimension(:, :), intent(out) :: rf_eb, rh_eb
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        integer :: l_sp, l_eb, l_lk

        rf_eb = 0.0
        rh_eb = 0.0
        do l_sp = 1, sp_n
            do l_eb = 1, eb_n
                do l_lk = 1, lk_n_eb(l_eb)
                    rf_eb(l_sp, l_eb) = rf_eb(l_sp, l_eb) + f_pc(l_sp, lk_pc_eb(2, l_lk, l_eb), lk_pc_eb(1, l_lk, l_eb))
                    rh_eb(l_sp, l_eb) = rh_eb(l_sp, l_eb) + h_pc(l_sp, lk_pc_eb(2, l_lk, l_eb), lk_pc_eb(1, l_lk, l_eb))
                end do
            end do
        end do
        rf_eb = rf_eb*length_m_2/eb_en_w
        rh_eb = rh_eb*length_m_2/eb_en_w
    end subroutine calculate_el_distro_energy

    subroutine calculate_distro_average(f_eb, h_eb, f_pc, h_pc)
        ! Calculate the spin-resolved particle distribution averaged over
        ! states with the same energy.
        implicit none
        ! Particle distribution [dimensionless].
        double precision, dimension(:, :), intent(out) :: f_eb, h_eb
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        integer :: l_sp, l_eb, l_lk

        f_eb = 0.0
        h_eb = 0.0
        do l_sp = 1, sp_n
            do l_eb = 1, eb_n
                do l_lk = 1, lk_n_eb(l_eb)
                    f_eb(l_sp, l_eb) = f_eb(l_sp, l_eb) + f_pc(l_sp, lk_pc_eb(2, l_lk, l_eb), lk_pc_eb(1, l_lk, l_eb))
                    h_eb(l_sp, l_eb) = h_eb(l_sp, l_eb) + h_pc(l_sp, lk_pc_eb(2, l_lk, l_eb), lk_pc_eb(1, l_lk, l_eb))
                end do
                if (lk_n_eb(l_eb) .ne. 0) then
                    f_eb(l_sp, l_eb) = f_eb(l_sp, l_eb)/lk_n_eb(l_eb)
                    h_eb(l_sp, l_eb) = h_eb(l_sp, l_eb)/lk_n_eb(l_eb)
                end if
            end do
        end do
    end subroutine calculate_distro_average

    subroutine print_particle_distro(f_pc, fileName)
        ! Print the electron distribution.
        implicit none
        double precision, dimension(:, :, :), intent(in) :: f_pc
        character(len=*), intent(in) :: fileName
        character(len=50) :: fmtString
        integer :: l_pc, l_ba, l_sp, col_n

        col_n = ba_n*sp_n
        write (fmtString, '(A,I2,A)') "(", col_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file=fileName, status="unknown")
        do l_pc = 1, pc_n
            write (30, fmtString) ((f_pc(l_sp, l_ba, l_pc), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)
    end subroutine print_particle_distro

    subroutine print_particle_energy_distro(rf_eb, fileName)
        ! Print the spin-resolved particle energy density.
        implicit none
        double precision, dimension(:, :), intent(in) :: rf_eb
        character(len=*), intent(in) :: fileName
        character(len=50) :: fmtString
        integer :: l_sp, l_eb

        write (fmtString, '(A,I2,A)') "(", sp_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file=fileName, status="unknown")
        do l_eb = 1, eb_n
            write (30, fmtString) (rf_eb(l_sp, l_eb), l_sp=1, sp_n)
        end do
        close (unit=30)
    end subroutine print_particle_energy_distro

    subroutine print_averaged_particle_distro(f_eb, fileName)
        ! Print the spin-resolved average particle distribution.
        implicit none
        character(len=*), intent(in) :: fileName
        double precision, dimension(:, :), intent(in) :: f_eb
        character(len=50) :: fmtString
        integer :: l_sp, l_eb

        write (fmtString, '(A,I2,A)') "(", sp_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file=fileName, status="unknown")
        do l_eb = 1, eb_n
            write (30, fmtString) (f_eb(l_sp, l_eb), l_sp=1, sp_n)
        end do
        close (unit=30)
    end subroutine print_averaged_particle_distro

end module ELECTRON_DISTRIBUTION

module ELECTRON_SCATTERING
    use omp_lib
    use math
    use color_msg
    use phys_const
    use band_structure
    use electron_distribution
    implicit none

    ! 1 -->--Z-->-- 3  k3 = k1 - q
    !        Z
    ! 2 -->--Z-->-- 4  k4 = k2 + q
    !
    ! Transferred wave vector: q = k1 - k3
    ! Transferred energy: hw = e1 - e3
    !
    ! H = V_{1234} \psi^\dag_1 \psi^\dag_2 \psi_4 \psi_3
    ! Spin is conserved for each particle: s_1 = s_3, s_2 = s_4.

    ! Refractive indices of the top and bottom substrates.
    double precision :: refr_ind_top, refr_ind_bottom
    ! Relative dielectric constant.
    double precision :: eps_avg
    ! Which broadening model to use for energy conservation:
    ! 1 = uniform, constant
    ! 2 = uniform, time-dependent
    ! 3 = process-dependent, constant
    integer :: broadening_ee_model
#define BEM_UCO 1
#define BEM_UTD 2
#define BEM_PDC 3
    ! Target value for the broadening [eV], used in the uniform constant model.
    double precision :: broadening_ee_eV
    ! Whether to calculate the Gamma functions using an extra integral over the
    ! frequency mesh.
    logical :: gamma_freq_integral
    ! Energy broadening of the delta functions [eV].
    double precision :: eta_ee
    ! Energy broadening in the denominators of the Lindhard function [eV].
    double precision :: eta_lindhard
    ! Which screening model to use:
    ! 1 = dynamic, 2 = time-dependent Thomas-Fermi, 3 = static, 4 = Thomas-Fermi
    integer :: screening_model
#define SCR_DYN 1
#define SCR_TTF 2
#define SCR_STA 3
#define SCR_STF 4
    ! Thomas-Fermi wave vector [nm^-1].
    double precision :: thomasfermi_q
    ! Process kind index: intra-band (1), inter-band (2).
    integer, parameter :: iv_n = 2
    integer, parameter :: iv_a = 1  ! intra-band
    integer, parameter :: iv_e = 2  ! inter-band
    ! Whether to exclude Auger processes, i.e. two-body scattering proceses
    ! where one electron changes band and the other does not, leading to
    ! a change in the number of carriers (electron/holes) in the two bands.
    logical :: suppress_auger
    ! Maximum energy exchanged in a scattering event.
    double precision :: w_scatt_max
    
    ! Wave vector mesh.
    ! Indices of the scattering wave vectors in the PC (or BZ) mesh.
    integer, dimension(:), allocatable :: lq
    ! Number of scattering wave vectors.
    integer :: lq_n
    ! Maximum transferred wave vector [nm^-1].
    double precision :: lq_qmax
    ! Norm of the scattering wave vectors.
    double precision, dimension(:), allocatable :: norm_lq
    ! Indices of the opposite vectors in the mesh.
    integer, dimension(:), allocatable :: minus_lq

    ! Frequency mesh, where several quantities are calculated [eV].
    double precision, dimension(:), allocatable :: ww
    ! Number of frequencies in the mesh.
    integer :: ww_n
    ! Absolute maximum value of the frequency mesh [eV].
    double precision :: ww_max
    ! Mesh step [eV].
    double precision :: ww_dw
    ! Indices of the opposite frequencies in the mesh.
    integer, dimension(:), allocatable :: minus_ww

    ! Quantities calculated on the wave vector and frequency meshes.
    ! Bare potential.
    complex(kind=8), dimension(:, :), allocatable :: pot_bare
    ! Lindhard function.
    complex(kind=8), dimension(:, :), allocatable :: chi0
    ! Screened potential.
    complex(kind=8), dimension(:, :), allocatable :: pot_scr
    ! Spectral density of particle-hole transitions.
    double precision, dimension(:, :, :), allocatable :: spdeh
    ! Scattering rates.
    double precision, dimension(:, :, :), allocatable :: gamma_ee_out
    double precision, dimension(:, :, :), allocatable :: gamma_ee_in

    ! Angular average of the Lindhard function.
    integer, parameter :: uu_n = 31
    double precision, dimension(:), allocatable :: uu
    complex(kind=8), dimension(:, :), allocatable :: chi0az

    private

    public :: refr_ind_top, refr_ind_bottom
    public :: thomasfermi_q, eta_ee
    public :: init_electron_scattering, calculate_ee_rate
    public :: analyze_spectral_density_support
    public :: print_lindhard, print_spectral_density
    public :: optical_eps

contains

    subroutine init_electron_scattering
        implicit none

        character(len=50) :: fmt_string
        integer :: l_lq, l_ww, l_uu
        ! User-friendly parameters for the input file.
        double precision :: relative_eps_top, relative_eps_bottom
        double precision :: maximum_scatt_wavevector_nm_m_1
        double precision :: maximum_scatt_energy_eV
        double precision :: thomasfermi_wv_nm
        double precision :: eta_lindhard_factor
        double precision :: uumax, uumin
        namelist /coulomb_scattering_params/ relative_eps_top, relative_eps_bottom, &
            maximum_scatt_wavevector_nm_m_1, maximum_scatt_energy_eV, &
            suppress_auger, &
            screening_model, thomasfermi_wv_nm, eta_lindhard_factor, &
            broadening_ee_model, broadening_ee_eV, gamma_freq_integral
        namelist /electron_scattering_values/ eta_ee, ww_max, ww_n, lq_n

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=coulomb_scattering_params)
        close (37)

        ! Set some parameters.
        refr_ind_top = sqrt(relative_eps_top)
        refr_ind_bottom = sqrt(relative_eps_bottom)
        eps_avg = 0.5 * (relative_eps_top + relative_eps_bottom)
        lq_qmax = maximum_scatt_wavevector_nm_m_1
        w_scatt_max = maximum_scatt_energy_eV

        select case (screening_model)
            case (SCR_DYN)
                LOGME trim(msg_aqua("init_electron_scattering: Using dynamic screening."))
            case (SCR_TTF)
                LOGME trim(msg_aqua("init_electron_scattering: Using time-dependent Thomas-Fermi screening."))
            case (SCR_STA)
                LOGME trim(msg_aqua("init_electron_scattering: Using static screening."))
            case (SCR_STF)
                LOGME trim(msg_aqua("init_electron_scattering: Using constant Thomas-Fermi screening."))
            case default
                LOGME trim(msg_red("init_electron_scattering: Undefined screening model."))
                stop
        end select

        if (gamma_freq_integral) then
            LOGME trim(msg_aqua("init_electron_scattering: Calculate Gammas with extra frequency integral."))
        else
            LOGME trim(msg_aqua("init_electron_scattering: Calculate Gammas with interpolation."))
        end if

        ! Thomas-Fermi screening length.  Initialize to the value given by the
        ! user, but it might be modified by the program, or not used,
        ! depending on the screening model.
        thomasfermi_q = thomasfermi_wv_nm

        ! Set-up the calculation of the broadening in the delta functions
        ! of energy conservation.
        select case (broadening_ee_model)
            case (BEM_UCO)
                LOGME trim(msg_aqua("init_electron_scattering: Using uniform constant broadening."))
                eta_ee = broadening_ee_eV
            case (BEM_UTD)
                LOGME trim(msg_aqua("init_electron_scattering: Using uniform time-dependent broadening."))
                eta_ee = -1.0
            case (BEM_PDC)
                LOGME trim(msg_aqua("init_electron_scattering: Using process-dependent broadening."))
                eta_ee = -1.0
            case default
                LOGME trim(msg_red("init_electron_scattering: Undefined broadening model."))
                stop
        end select

        ! Constant and uniform broadening for the Lindhard function.
        eta_lindhard = eta_lindhard_factor * de_avg

        call init_scattering_mesh
        call init_frequency_mesh

        allocate (pot_bare(lq_n, ww_n))
        allocate (chi0(lq_n, ww_n))
        allocate (pot_scr(lq_n, ww_n))
        allocate (spdeh(lq_n, ww_n, iv_n))
        allocate (gamma_ee_out(sp_n, ba_n, pc_n))
        allocate (gamma_ee_in(sp_n, ba_n, pc_n))

        ! The bare potential does not depend on the particle distribution and
        ! can be calculated only once.
        call calculate_bare_potential

        ! Angular average of the Lindhard function.
        if (screening_model .eq. SCR_DYN .or. screening_model .eq. SCR_STA) then
            allocate(chi0az(uu_n,ww_n))
            ! One-dimensional wave vector mesh.
            allocate(uu(uu_n))
            uumax = maxval(norm_lq)
            uumin = minval(norm_lq)
            do l_uu = 1, uu_n
                uu(l_uu) = uumin + (l_uu - 1.0) * (uumax - uumin) / (uu_n - 1.0)
            end do

            open (unit=30, file="out-geo-lindhard-wv-mesh.csv", status="unknown")
            do l_uu = 1, uu_n
                write (30, '(F15.8)') uu(l_uu)
            end do
            close (unit=30)
        end if

        ! Save derived values.
        open (unit=VALUNIT, file=VALFILE, status="old", position="append")
        write (VALUNIT, nml=electron_scattering_values)
        close (VALUNIT)

        ! Save the arrays.
        ! Scattering wave vectors.  Print also their opposites as a check.
        write (fmt_string, "(A,I1,A)") "(", 2*ndim - 1, "(F15.8,', '),F15.8)"
        open (unit=30, file="out-geo-coulomb-scatt-mesh.csv", status="unknown")
        do l_lq = 1, lq_n
            write (30, fmt_string) wv_bz(:, lq(l_lq)), wv_bz(:, lq(minus_lq(l_lq)))
        end do
        close (unit=30)
        ! Frequency interpolation mesh.
        open (unit=30, file="out-geo-frequency-interp-mesh.csv", status="unknown")
        do l_ww = 1, ww_n
            write (30, '(F15.8)') ww(l_ww)
        end do
        close (unit=30)

    end subroutine init_electron_scattering

    subroutine init_frequency_mesh
        ! Initialize the frequency mesh.
        implicit none
        integer :: l_ww
        double precision :: ww_min

#if 0
        integer :: l_1_pc, l_lq, l_3_pc, l_1_ba, l_3_ba
        integer :: ii_3_1, ii_3_2
        integer, dimension(ndim) :: ii_1, ii_q, ii_3_tot
        double precision :: ww_diff, en_1, en_3

        ! Find the maximum energy that can be exchanged in a scattering event
        ! in the same band.
        ww_max = 0.0
        ! Initial wave vector.
        do l_1_pc = 1, pc_n
            ii_1 = ip_pc(:, l_1_pc)
            ! Wave vector difference.
            do l_lq = 1, lq_n
                ii_q = ip_pc(:, lq(l_lq))
                ! Calculate the final wave vector.
                ii_3_tot = ii_1 - ii_q
                call ip_fold(ii_3_1, ii_3_2, ii_3_tot)
                l_3_pc = l_from_ip(ii_3_1, ii_3_2, "init_frequency_mesh")
                ! Initial band.
                do l_1_ba = 1, ba_n
                    ! Final band.
                    l_3_ba = l_1_ba
                    en_1 = en_bz(l_1_ba, l_1_pc)
                    en_3 = en_bz(l_3_ba, l_3_pc)
                    ww_diff = abs(en_3 - en_1)
                    if (ww_diff .gt. ww_max) then
                        ww_max = ww_diff
                    end if
                end do
            end do
        end do
        ww_max = 1.1 * ww_max  ! a little padding
#endif

        ! Set the maximum value of the frequency mesh for the screening,
        ! with a little padding.
        ww_max = en_bz_max - en_bz_min + 0.001

        ! Set the size of the frequency mesh based on the energy uncertainty.
        ! ww_n = 2*int(ceiling(ww_max/(2.0*eta_ee)))

        ! The mesh is symmetric around zero, includes negative frequencies.
        ww_min = -ww_max

        ! Set the size of the frequency mesh based on the average distance
        ! between energy levels.
        ww_n = 2*int(ceiling(ww_max/(2.0*de_avg)))

        ! Initialize the mesh.
        allocate (ww(ww_n))
        do l_ww = 1, ww_n
            ww(l_ww) = ww_min + (ww_max - ww_min)*(l_ww - 1.0)/(ww_n - 1.0)
        end do
        ! Calculate the mesh step.
        ww_dw = ww(2) - ww(1)
        ! Indices of opposite frequencies.  Since the mesh is symmetric, it is
        ! just the index running backwards.
        allocate (minus_ww(ww_n))
        do l_ww = 1, ww_n
            minus_ww(l_ww) = ww_n - l_ww + 1
        end do
    end subroutine init_frequency_mesh

    subroutine init_scattering_mesh
        ! Initialize the mesh of the scattering wave vectors, i.e. the allowed
        ! values of transferred wave vector in a scattering event.
        implicit none
        logical :: found_opposite
        integer :: l_pc, l_lq, l_m_lq
        integer, dimension(ndim) :: ii, ii_m, ii_s, ii_f
        ! Count how many wave vectors in the BZ mesh are smaller than the
        ! maximum.
        l_lq = 0
        do l_pc = 1, pc_n
            if (keep_scatt_vec(l_pc)) then
                l_lq = l_lq + 1
            end if
        end do
        ! Set the size of the scattering wave vectors mesh.
        lq_n = l_lq
        ! Allocate an array which is large enough.
        allocate (lq(lq_n))
        ! Copy the indices into the array.
        l_lq = 0
        do l_pc = 1, pc_n
            if (keep_scatt_vec(l_pc)) then
                l_lq = l_lq + 1
                lq(l_lq) = l_pc
            end if
        end do
        ! Calculate the norm of the scattering wave vectors.
        allocate (norm_lq(lq_n))
        do l_lq = 1, lq_n
            norm_lq(l_lq) = sqrt(wv_bz(1, lq(l_lq))**2 + wv_bz(2, lq(l_lq))**2)
        end do
        ! Initialize the index of the opposite wave vector.
        allocate (minus_lq(lq_n))
        ! Perform a search for all wave vectors in the mesh.
        do l_lq = 1, lq_n
            ii = ip_pc(:, lq(l_lq))
            found_opposite = .false.
            ! Test all wave vectors in the mesh.
            do l_m_lq = 1, lq_n
                ii_m = ip_pc(:, lq(l_m_lq))
                ! Sum the two vectors and fold the result.
                ! If we obtain the Gamma point, the two vectors are opposite.
                ii_s = ii + ii_m
                call ip_fold (ii_f(1), ii_f(2), ii_s)
                if ((ii_f(1) .eq. ip_pc(1,l_pc_gamma)) .and. (ii_f(2) .eq. ip_pc(2,l_pc_gamma))) then
                    minus_lq(l_lq) = l_m_lq
                    found_opposite = .true.
                    ! No need to search any longer for this vector.
                    exit
                end if
            end do
            ! If we have not found the opposite of this vector, raise an error.
            if (.not. found_opposite) then
                LOGME trim(msg_red("init_scattering_mesh: Opposite not found."))
                stop
            end if
        end do
    end subroutine init_scattering_mesh

    function keep_scatt_vec(l_pc)
        ! Returns true if the wave vector corresponding to the index l_pc should
        ! be included in the scattering wave vector mesh.
        implicit none
        logical :: keep_scatt_vec
        integer, intent(in) :: l_pc
        double precision :: wv_length

        wv_length = sqrt(wv_bz(1, l_pc)**2 + wv_bz(2, l_pc)**2)
        if ((wv_length .lt. lq_qmax) .and. (l_pc .ne. l_pc_gamma)) then
            ! Keep the vector if its norm is smaller than a fixed cutoff and if
            ! it is not the Gamma point, i.e. the origin.
            keep_scatt_vec = .true.
        else
            keep_scatt_vec = .false.
        end if
    end function keep_scatt_vec

    subroutine coeff_interp(c_l, c_r, l_l_ww, l_r_ww, en)
        ! Coefficients to calculate an interpolation on the frequency mesh.
        implicit none
        ! Weights of the left and right approximation.
        double precision, intent(out) :: c_l, c_r
        ! Indices of the left and right approximation on the mesh.
        integer, intent(out) :: l_l_ww, l_r_ww
        ! Desired value of the frequency [eV].
        double precision, intent(in) :: en
        ! Index on the frequency mesh.
        integer :: l_ww

        if (en .le. ww(1)) then
#if 0
            ! Extrapolation below.
            l_l_ww = 1
            l_r_ww = 1
            cwL = 0.5
            cwR = 0.5
#else
            LOGME trim(msg_red("coeff_interp: Extrapolation below."))
            stop
#endif
        else if (en .ge. ww(ww_n)) then
#if 0
            ! Extrapolation above.
            l_l_ww = ww_n
            l_r_ww = ww_n
            cwL = 0.5
            cwR = 0.5
#else
            LOGME trim(msg_red("coeff_interp: Extrapolation above."))
            stop
#endif
        else
            ! Find the indices in the mesh which bracket the desired value.
            do l_ww = ww_n - 1, 1, -1
                if (en .ge. ww(l_ww)) then
                    exit
                end if
            end do
            l_l_ww = l_ww
            l_r_ww = l_ww + 1
            ! Coefficients to calculate the linear interpolation of any
            ! function calculated on the frequency mesh.
            c_l = (ww(l_r_ww) - en)/(ww(l_r_ww) - ww(l_l_ww))
            c_r = (en - ww(l_l_ww))/(ww(l_r_ww) - ww(l_l_ww))
        end if
    end subroutine coeff_interp

    subroutine calculate_eta_ee_weighted(f_pc)
        ! Calculate an average broadening, using the occupation as weight.
        implicit none
        ! Electron distribution.
        double precision, dimension(:, :, :), intent(in) :: f_pc
        integer :: l_ba, l_pc
        double precision :: de_tot, w_tot, w

        ! Calculate an effective average level distance, weighted over the
        ! deviation of the level occupation with respect to equilibrium.
        de_tot = 0.0
        w_tot = 0.0
        do l_ba = 1, ba_n
            do l_pc = 1, pc_n
                w = abs(sum(f_pc(:, l_ba, l_pc) - f_eq_pc(:, l_ba, l_pc)))
                de_tot = de_tot + de_pc_ba(l_pc, l_ba) * w
                w_tot = w_tot + w
            end do
        end do
        ! The average level distance is the weighted total divided by the
        ! sum of all weights, renormalized to the target broadening.
        eta_ee = (broadening_ee_eV / de_max) * (de_tot / w_tot)
    end subroutine calculate_eta_ee_weighted

    function delta_ee(l_o_ba, l_o_pc, l_i_ba, l_i_pc, en_d) result(d)
        ! δ(εₒᵤₜ − εᵢₙ − ℏω)
        ! Gaussian approximation to the Dirac delta function of energy
        ! conservation.  The approximation has to be Gaussian to work:
        ! See C. Illg et al., J Theor Appl Phys 10, 1 (2016)
        ! DOI: 10.1007/s40094-015-0193-5
        implicit none
        ! Result [1/eV].
        double precision :: d
        ! Indices of the scattering electron's states.
        integer, intent(in) :: l_o_ba, l_o_pc, l_i_ba, l_i_pc
        ! Positive/negative transferred energy.
        double precision, intent(in) :: en_d
        double precision :: en, en_sq, eta_loc_ee, eta_loc_ee_sq
        ! Minimum allowed broadening [eV].
        double precision, parameter :: eta_min = 0.005

        ! Broadening.
        select case (broadening_ee_model)
            case (BEM_PDC)
                ! Process-dependent value.
                ! Normalize the joint level spacing using the target value.
                eta_loc_ee = broadening_ee_eV * jde_bz(l_o_ba,l_o_pc,l_i_ba,l_i_pc) / jde_avg
                ! Correct for very small broadening, to be on the safe side.
                eta_loc_ee = max(eta_loc_ee, eta_min)
            case default
                ! Use the uniform value, previously calculated.
                eta_loc_ee = eta_ee
        end select
        eta_loc_ee_sq = eta_loc_ee * eta_loc_ee
        ! Energy combination that nullifies at the delta peak.
        ! δ(εₒᵤₜ − εᵢₙ − ℏω)
        en = en_bz(l_o_ba, l_o_pc) - en_bz(l_i_ba, l_i_pc) - en_d
        en_sq = en * en
        ! Calculate the approximation to the delta function.
        ! Avoid very large values of the argument, which produce numerical
        ! exceptions in the exponential.
        if (en_sq .lt. 25.0*eta_loc_ee_sq) then
            d = exp(-en_sq/eta_loc_ee_sq)/pi_sqrt/eta_loc_ee
        else
            d = 0.0
        endif
    end function delta_ee

    function density_matrix_element(l_i_ba, l_o_ba, l_i_pc, l_o_pc) result( m )
        ! The matrix element <out|n_q|in> of the Fourier transform n_q of the
        ! density operator is given by
        ! d_{s_in,s_out} d_{k_out,k_in-q} F_{b_in,b_out}(k_in,k_out),
        ! where d_{,} is the Kronecker delta, s, b, and k label spin, band, and
        ! wave vector, respectively.
        ! This subroutine calculates the coefficient F_{b_in,b_out}(k_in,k_out).
        implicit none
        complex(kind=8) :: m
        ! Band indices.
        integer, intent(in) :: l_i_ba, l_o_ba
        ! Wave vector indices.
        integer, intent(in) :: l_i_pc, l_o_pc
        double precision :: band_sign

        if (l_i_ba .eq. l_o_ba) then
            band_sign = 1.0
        else
            band_sign = -1.0
        end if
        m = 0.5*(band_sign + conjg(fconn_ph_bz(l_i_pc)) * fconn_ph_bz(l_o_pc))
    end function density_matrix_element

    subroutine calculate_bare_potential
        ! Calculate the array with the bare Coulomb potential.
        ! The bare potential does not usually depend on frequency, although
        ! it might if the dielectric response of the substrate depends on
        ! frequency, as in the case of hBN.
        implicit none
        integer :: l_lq, l_ww
        do l_ww = 1, ww_n
            do l_lq = 1, lq_n
                pot_bare(l_lq, l_ww) = 2.0*pi*e_sq/(eps_avg*norm_lq(l_lq))
            end do
        end do
    end subroutine calculate_bare_potential

    subroutine calculate_lindhard(f_pc)
        ! Calculate the array with the Lindhard function.
        implicit none
        ! Electron distribution.
        double precision, dimension(:, :, :), intent(in) :: f_pc
        ! Indices of the final matrix.
        integer :: l_ww, l_lq
        ! Indices for the summations.
        integer :: l_i_pc, l_o_pc, l_i_ba, l_o_ba, l_sp
        ! Integer projections of wave vectors.
        integer :: ii_o_1, ii_o_2
        integer, dimension(ndim) :: ii_i, ii_q, ii_o_tot
        ! Coefficients.
        double precision :: en_i, en_o
        complex(kind=8) :: g, cx
        ! Variable where the sum is accumulated at each wave vector and
        ! frequency.
        complex(kind=8) :: total
        ! Matrix element of the density operator between states <4| and |2>.
        complex(kind=8) :: dme_io
        character(len=50), parameter :: errString = "calculate_lindhard"

!$OMP PARALLEL DO PRIVATE(l_lq, ii_q, l_ww, total, l_i_pc, ii_i, ii_o_tot, ii_o_1 ) &
!$OMP & PRIVATE(ii_o_2, l_o_pc, l_i_ba, en_i, l_o_ba, en_o, g, dme_io, cx, l_sp)
        ! External.  Wave vector difference.
        do l_lq = 1, lq_n
            ! External.  Frequency.
            do l_ww = 1, ww_n
                ! Copy the integer projections to local variable, to avoid
                ! accessing the array many times in the loops.
                ! CHECK: Does this make the algorithm faster?
                ii_q = ip_pc(:, lq(l_lq))
                ! Initialize the result of the summation.
                total = 0.0
                ! Summed over.  Initial wave vector.
                do l_i_pc = 1, pc_n
                    ii_i = ip_pc(:, l_i_pc)
                    ! Calculate the final wave vector.
                    ii_o_tot = ii_i + ii_q
                    call ip_fold(ii_o_1, ii_o_2, ii_o_tot)
                    l_o_pc = l_from_ip(ii_o_1, ii_o_2, errString)
                    ! Summed over.  Initial band.
                    do l_i_ba = 1, ba_n
                        en_i = en_bz(l_i_ba, l_i_pc)
                        ! Summed over.  Final band.
                        do l_o_ba = 1, ba_n
                            en_o = en_bz(l_o_ba, l_o_pc)
                            ! The denominator.
                            g = cmplx(real(en_i - en_o + ww(l_ww), kind=4), real(eta_lindhard, kind=4))
                            ! Matrix element of the density operator.
                            dme_io = density_matrix_element(l_i_ba, l_o_ba, l_i_pc, l_o_pc)
                            cx = dme_io*conjg(dme_io)/g
                            ! Summed over.  Spin, conserved.
                            do l_sp = 1, sp_n
                                total = total + cx*(f_pc(l_sp, l_i_ba, l_i_pc) - f_pc(l_sp, l_o_ba, l_o_pc))
                            end do
                        end do
                    end do
                end do
                chi0(l_lq, l_ww) = total*length_m_2
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine calculate_lindhard

    subroutine calculate_thomasfermi_q(f_pc)
        implicit none
        ! Electron distribution.
        double precision, dimension(:, :, :), intent(in) :: f_pc

        integer :: l_sp, l_eb, l_lk
        ! Effective DOS [eV^-1 nm^-2] that generalizes dn/dmu in the
        ! Thomas-Fermi formula to non-Fermi-Dirac distributions.
        double precision :: edos

        ! Average occupation in each energy bin.
        double precision, dimension(:), allocatable :: f_av

        ! Out-of-equilibrium there is no chemical potential, which is the
        ! cornerston of Thomas-Fermi theory.  However, the bottom line is that
        ! we want to obtain an effective DOS by summing the DOSes at those
        ! energy bins where the Fermi-Dirac distribution has a sizable step.

        edos = 0.0
        allocate( f_av(eb_n) )
        do l_sp = 1, sp_n
            ! Calculate the average occupation in an energy shell.
            f_av = 0.0
            do l_eb = 1, eb_n
                do l_lk = 1, lk_n_eb(l_eb)
                    f_av(l_eb) = f_av(l_eb) + f_pc(l_sp, lk_pc_eb(2, l_lk, l_eb), lk_pc_eb(1, l_lk, l_eb))
                end do
                f_av(l_eb) = f_av(l_eb) / lk_n_eb(l_eb)
            end do
            ! Average the DOS, by summing over the energy shells, using the
            ! gradient of the average occupation as weight.
            do l_eb = 2, eb_n
                edos = edos + 0.5 * (dos_eb(l_eb) + dos_eb(l_eb - 1)) * abs(f_av(l_eb) - f_av(l_eb - 1)) / eb_en_w 
            end do
        end do
        deallocate( f_av )

        ! q_TF = dn/dmu * (2 pi e^2) / bar_eps
        thomasfermi_q = edos * 2.0 * pi * e_sq / eps_avg
    end subroutine calculate_thomasfermi_q

    subroutine calculate_screened_potential
        ! Calculate the array with the RPA screened potential.
        implicit none
        integer :: iwn, iwp, l_lq
        double precision :: thomasfermi_eps

        if (screening_model .eq. SCR_DYN) then
                ! Dynamical screening in the RPA.
                ! Divide by a factor which depends on frequency and wave vector.
                pot_scr = pot_bare/(1 - pot_bare*chi0)
        else if (screening_model .eq. SCR_TTF .or. screening_model .eq. SCR_STF) then
                ! Thomas-Fermi screening.
                ! Divide by a factor that depends on wave vector only.
                do l_lq = 1, lq_n
                    thomasfermi_eps = 1.0 + thomasfermi_q / norm_lq(l_lq)
                    pot_scr(l_lq,:) = pot_bare(l_lq,:) / thomasfermi_eps
                end do
        else if (screening_model .eq. SCR_STA) then
                ! Static screening, taking the limit of the RPA.
                ! Divide by a factor that depends on wave vector only.
                iwn = ww_n / 2  ! negative freq close to zero
                iwp = iwn + 1  ! positive freq close to zero
                do l_lq = 1, lq_n
                    pot_scr(l_lq,:iwn) = pot_bare(l_lq,:iwn)/ &
                        (1-pot_bare(l_lq,:iwn)*chi0(l_lq,iwn))
                    pot_scr(l_lq,iwp:) = pot_bare(l_lq,iwp:)/ &
                        (1-pot_bare(l_lq,iwp:)*chi0(l_lq,iwp))
                end do
        end if
    end subroutine calculate_screened_potential

    subroutine optical_eps(eps_re_w0, w0, w0_n, f_pc)
        ! Calculate the value of the REAL part of thedielectric function at q=0,
        ! by averaging the values at the closest wave vectors.
        ! Perform the calculation on a set of frequencies.
        ! Dielectric function at q = 0 [dimensionless] and several frequencies.
        double precision, dimension(:), intent(out) :: eps_re_w0
        ! Frequencies to use [eV].
        double precision, dimension(:), intent(in) :: w0
        ! Length of the w0 and eps_re_w0 arrays.
        integer, intent(in) :: w0_n
        double precision, dimension(:,:,:), intent(in) :: f_pc
        ! Dielectric function, which is complex in general.
        complex(kind=8), dimension(:), allocatable :: eps_w0

        logical :: extra_used
        integer, parameter :: cl_n = 6
        integer :: l_cl, l_w0, l_ww, l_l_ww, l_r_ww
        double precision :: a_l, a_r
        integer, dimension(ndim,cl_n) :: ii_cl
        integer, dimension(cl_n) :: lq_cl
        character(len=50), parameter :: errString = "optical_eps"

        ! Integer coordinates of the wave vectors around the Gammma point.
        ii_cl(:,1) = (1,0)
        ii_cl(:,2) = (1,1)
        ii_cl(:,3) = (0,1)
        ii_cl(:,4) = (-1,0)
        ii_cl(:,5) = (-1,-1)
        ii_cl(:,6) = (0,-1)

        ! Find the corresponding indices in the scattering wave vector mesh.
        do l_cl = 1, cl_n
            lq_cl(l_cl) = findloc(lq, &
                & l_from_ip(ii_cl(1,l_cl), ii_cl(2,l_cl), errString), &
                & 1)
        end do

        ! Since this function is not called often, we resort to calculating
        ! the screened potential on the entire wave vector mesh.
        ! Hugely inefficient - edit if this function is used more often.
        ! This shortcut also allows us to re-use the logic that chooses the
        ! screening method, otherwise we should factor that out.
        if (screening_model .eq. SCR_DYN .or. screening_model .eq. SCR_STA) then
            call calculate_lindhard(f_pc)
        else if (screening_model .eq. SCR_TTF) then 
            call calculate_thomasfermi_q(f_pc)
        end if

        ! In the main procedure, we calculate the screened potential, not the
        ! dielectric function.  But of course we can retrieve the dielectric
        ! function since W_sc = V_0 / eps  =>  eps = V_0 / W_sc
        call calculate_screened_potential

        allocate( eps_w0(w0_n) )
        extra_used = .false.
        ! Calculate the dielectric function for all input frequencies.
        ! Average over the wave vectors around the Gamma point.
        do l_w0 = 1, w0_n
            ! Initialize the value to zero, so we can sum the contributions.
            eps_w0(l_w0) = 0.0
            if (w0(l_w0) .le. ww(1)) then
                ! Extrapolation by using the smaller frequency in the mesh.
                do l_cl = 1, cl_n
                    eps_w0(l_w0) = eps_w0(l_w0) + pot_bare(lq_cl(l_cl),1) / pot_scr(lq_cl(l_cl),1)
                end do
                extra_used = .true.
            else if (w0(l_w0) .ge. ww(ww_n)) then
                ! Extrapolation by using the largest frequency in the mesh.
                do l_cl = 1, cl_n
                    eps_w0(l_w0) = eps_w0(l_w0) + pot_bare(lq_cl(l_cl),ww_n) / pot_scr(lq_cl(l_cl),ww_n)
                end do
                extra_used = .true.
            else
                ! Bracket the frequency in the mesh.
                do l_ww = 1, ww_n - 1
                    if (ww(l_ww+1) .gt. w0(l_w0)) then
                        l_l_ww = l_ww  ! index of the left frequency
                        l_r_ww = l_l_ww + 1
                        exit
                    end if
                end do
                ! Interpolation coefficients.
                a_l = (ww(l_r_ww) - w0(l_w0)) / (ww(l_r_ww) - ww(l_l_ww))
                a_r = 1.0 - a_r
                ! Interpolate and average.
                do l_cl = 1, cl_n
                    eps_w0(l_w0) = eps_w0(l_w0) &
                        & + a_l * pot_bare(lq_cl(l_cl),l_l_ww) / pot_scr(lq_cl(l_cl),l_l_ww) &
                        & + a_r * pot_bare(lq_cl(l_cl),l_r_ww) / pot_scr(lq_cl(l_cl),l_r_ww)
                end do
            end if
            ! Complete the average procedure by diving by the number of wave
            ! vectors around the Gamma point.
            eps_w0(l_w0) = eps_w0(l_w0) / cl_n
        end do

        ! Take the real part and copy to the array with the results.
        eps_re_w0(:) = real(eps_w0(:))
        deallocate(eps_w0)

        ! Warn if extrapolation has been used.
        ! Probably no need to stop the execution.
        if (extra_used) then
            LOGME trim(msg_aqua("optical_eps: Extrapolation used."))
        end if
    end subroutine optical_eps


    subroutine analyze_spectral_density_support
        ! Calculate how many scattering events fall in each frequency bin.
        implicit none

        ! Indices of the scattering states.
        integer :: iv, l_lq, l_2_pc, l_4_pc, l_2_ba, l_4_ba, l_sp
        ! Index in the frequency mesh.
        integer :: l_ww
        ! Index of the scattering event.
        integer :: i_ev, ev_n
        ! Energy of the scattering event, i.e. difference between the energies
        ! of the outgoing and the incoming state.
        double precision :: en_ev
        ! Integer projections of wave vectors.
        integer :: ii_4_1, ii_4_2
        integer, dimension(ndim) :: ii_2, ii_q, ii_4_tot
        ! Number of scattering events that fall in each frequency bin.
        integer, dimension(:,:), allocatable :: events_ww
        ! List of all the energies of the scattering events for all wave vectors.
        double precision, dimension(:,:), allocatable :: energies_events
        character(len=50) :: fmtString
        character(len=50), parameter :: errString = "analyze_spectral_density_support"

        ! Allocate the array with the count.
        allocate( events_ww(lq_n,ww_n))

        ! Allocate the array with all the energies.
        ev_n = pc_n*ba_n*ba_n*sp_n
        allocate( energies_events(ev_n, lq_n) )
        call print_memory_usage(energies_events, "energies_events")
        energies_events = -100.0

        ! External.  Wave vector difference.
        do l_lq = 1, lq_n
            ! Initialize the counts for the current exchanged wave vector.
            events_ww(:,l_lq) = 0
            i_ev = 0
            ! Copy the integer projections to local variable.
            ii_q = ip_pc(:, lq(l_lq))
            ! Summed over.  Initial wave vector.
            do l_2_pc = 1, pc_n
                ii_2 = ip_pc(:, l_2_pc)
                ! Calculate the final wave vector.
                ii_4_tot = ii_2 + ii_q
                call ip_fold(ii_4_1, ii_4_2, ii_4_tot)
                l_4_pc = l_from_ip(ii_4_1, ii_4_2, errString)
                ! Summed over.  Initial band.
                do l_2_ba = 1, ba_n
                    ! en_2 = en_bz(l_2_ba, l_2_pc)
                    ! Summed over.  Final band.
                    do l_4_ba = 1, ba_n
                        if (l_4_ba .eq. l_2_ba) then
                            iv = iv_a  ! intra-band
                        else
                            iv = iv_e  ! inter-band
                        end if
                        ! Summed over.  Spin, conserved.
                        do l_sp = 1, sp_n
                            ! The energy of the scattering event.
                            en_ev = en_bz(l_2_ba, l_2_pc) - en_bz(l_4_ba, l_4_pc)
                            ! Save the energy to the list.
                            i_ev = i_ev + 1
                            if (i_ev .gt. ev_n) then
                                LOGME trim(msg_red("analyze_spectral_density_support: Event index overflow."))
                                stop
                            end if
                            energies_events(i_ev,l_lq) = en_ev
                            ! The index of the energy in the frequency mesh.
                            l_ww =  int( (en_ev - ww(1)) / ww_dw) + 1
                            ! Increment the counter.
                            events_ww(l_lq,l_ww) = events_ww(l_ww,l_lq) + 1
                        end do
                    end do
                end do
            end do
        ! External.  Wave vector difference.
        end do

        LOGME trim("analyze_spectral_density_support: Sorting energies.")
!$OMP PARALLEL DO PRIVATE(l_lq)
        do l_lq = 1, lq_n
            ! Sort the list of energies.
            call quicksort(energies_events(:,l_lq),1,ev_n)
        end do
!$OMP END PARALLEL DO

        ! Print the energy histogram of events for each wave vector.
        write (fmtString, '(A,I04,A)') "(", ww_n - 1, "(I5, ', '),I5)"
        open (unit=30, file="out-test-spectral-hist.csv", status="unknown")
        do l_lq = 1, lq_n
            write (30, fmtString) events_ww(l_lq,:)
        end do

        ! Print all energies.
        write (fmtString, '(A,I06,A)') "(", ev_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-test-spectral-list.csv", status="unknown")
        do l_lq = 1, lq_n
            write (30, fmtString) energies_events(:,l_lq)
        end do

        deallocate( energies_events )
        deallocate( events_ww )

    end subroutine analyze_spectral_density_support


    subroutine calculate_spectral_density(f_pc, h_pc)
        ! 𝒮(𝐪, ℏ𝜔) = 𝟏/𝑁 ∑_𝐤 f_𝐤 (𝟏 - f_𝐤₊𝐪) δ(𝜀_𝐤₊𝐪 - 𝜀_𝐤 - ℏ𝜔)
        implicit none
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc

        ! Process kind index: intra-band (1) or inter-band (2).
        integer :: iv
        ! Indices of the final matrix.
        integer :: l_ww, l_lq
        ! Indices for the summations.
        integer :: l_2_pc, l_4_pc, l_2_ba, l_4_ba, l_sp
        ! Integer projections of wave vectors.
        integer :: ii_4_1, ii_4_2
        integer, dimension(ndim) :: ii_2, ii_q, ii_4_tot
        ! Coefficients.
        double precision :: a0, cx
        ! Variable where the sum is accumulated at each wave vector and
        ! frequency, separately for intra- and inter-band processes.
        double precision, dimension(iv_n) :: total
        ! Matrix element <4| ... |2> of the density operator.
        complex(kind=8) :: dme_42
        character(len=50), parameter :: errString = "calculate_spectral_density"

!$OMP PARALLEL DO PRIVATE(l_lq, ii_q, l_ww, total, l_2_pc, ii_2, ii_4_tot, ii_4_1) &
!$OMP & PRIVATE(ii_4_2, l_4_pc, l_2_ba, l_4_ba, a0, dme_42, cx, l_sp)

        ! External.  Wave vector difference.
        do l_lq = 1, lq_n
            ! External.  Frequency.
            do l_ww = 1, ww_n
                ! Copy the integer projections to local variable.
                ii_q = ip_pc(:, lq(l_lq))
                ! Initialize the result of the summation.
                total(:) = 0.0
                ! Summed over.  Initial wave vector.
                do l_2_pc = 1, pc_n
                    ii_2 = ip_pc(:, l_2_pc)
                    ! Calculate the final wave vector.
                    ii_4_tot = ii_2 + ii_q
                    call ip_fold(ii_4_1, ii_4_2, ii_4_tot)
                    l_4_pc = l_from_ip(ii_4_1, ii_4_2, errString)
                    ! Summed over.  Initial band.
                    do l_2_ba = 1, ba_n
                        ! Summed over.  Final band.
                        do l_4_ba = 1, ba_n
                            if (l_4_ba .eq. l_2_ba) then
                                iv = iv_a  ! intra-band
                            else
                                iv = iv_e  ! inter-band
                            end if
                            ! The Dirac delta function of energy conservation.
                            ! a0 = delta_ee(en_2 - en_4 + ww(l_ww))
                            a0 = delta_ee(l_2_ba, l_2_pc, l_4_ba, l_4_pc, - ww(l_ww))
                            ! Matrix element of the density operator.
                            dme_42 = density_matrix_element(l_2_ba, l_4_ba, l_2_pc, l_4_pc)
                            cx = a0 * real(dme_42 * conjg(dme_42))
                            ! Summed over.  Spin, conserved.
                            do l_sp = 1, sp_n
                                total(iv) = total(iv) + cx * f_pc(l_sp, l_2_ba, l_2_pc) * h_pc(l_sp, l_4_ba, l_4_pc)
                            end do
                        end do
                    end do
                end do
                spdeh(l_lq, l_ww, :) = total(:) / pc_n
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine calculate_spectral_density

    subroutine calculate_ee_gammas(f_pc, h_pc)
        ! Calculate the scattering rate due to Coulomb interactions.
        ! Γ⁽ᵒᵘᵗ⁾_𝐤 = (𝟐π / ℏ) (𝟏 / 𝑁) ∑_𝐪 ≠ 𝟎 (ℬ² / (𝟐π)⁴) |W(𝐪, 𝜀_𝐤 - 𝜀_𝐤₋𝐪)|² (𝟏 - f_𝐤₋𝐪) 𝒮(𝐪, 𝜀_𝐤 - 𝜀_𝐤₋𝐪)
        !
        ! Γ⁽ⁱⁿ⁾_𝐤 = (𝟐π / ℏ) (𝟏 / 𝑁) ∑_𝐪 ≠ 𝟎 (ℬ² / (𝟐π)⁴) |W(𝐪, 𝜀_𝐤 - 𝜀_𝐤₋𝐪)|² f_𝐤₋𝐪 𝒮(-𝐪, -𝜀_𝐤 + 𝜀_𝐤₋𝐪)
        !
        implicit none
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc

        ! Indices of the distribution and the summations.
        integer :: l_1_pc, l_3_pc, l_lq, l_m_lq, l_1_ba, l_3_ba, l_sp
        ! Indices for the interpolation on the frequency mesh.
        integer :: l_l_ww, l_r_ww
        ! Integer projections of wave vectors.
        integer :: ii_3_1, ii_3_2
        integer, dimension(ndim) :: ii_1, ii_q, ii_3_tot
        ! Coefficients.
        double precision :: en_1, en_3
        double precision :: c_l, c_r, c0
        integer, dimension(iv_n) :: cv
        double precision:: w2s_l, w2s_r, w2s_out, w2s_in, w13_in, w13_out
        ! Variables where the sum is accumulated at each wave vector, band,
        ! and spin, for outgoing and incoming processes.
        double precision :: total_out, total_in
        ! Matrix element <3| ... |1> of the density operator.
        complex(kind=8) :: dme_31
        double precision :: dme_sq
        character(len=50), parameter :: errString = "calculate_ee_gammas"

        c0 = (2*pi/hbar) * area_bz**2 / (2*pi)**4 / pc_n

!$OMP PARALLEL DO PRIVATE(l_1_pc, l_1_ba, l_sp, total_out, total_in, l_lq) &
!$OMP & PRIVATE(ii_1, ii_q, l_m_lq, ii_3_tot, ii_3_1, ii_3_2, l_3_pc) &
!$OMP & PRIVATE(en_1, l_3_ba, cv, en_3, c_l, c_r, l_l_ww, l_r_ww) &
!$OMP & PRIVATE(w2s_l, w2s_r, w2s_out, w2s_in, dme_31, dme_sq, w13_out, w13_in)
        ! External.  Initial wave vector.
        do l_1_pc = 1, pc_n
        ! External.  Initial band.
        do l_1_ba = 1, ba_n
        ! External.  Spin, conserved.
        do l_sp = 1, sp_n
            ! Initialize the results of the summations.
            total_out = 0.0
            total_in = 0.0
            ! Summed over.  Wave vector difference.
            do l_lq = 1, lq_n
                ii_1 = ip_pc(:, l_1_pc)
                ii_q = ip_pc(:, lq(l_lq))
                ! Get the index of the opposite wave vector.
                l_m_lq = minus_lq(l_lq)
                ! Calculate the final wave vector.
                ii_3_tot = ii_1 - ii_q
                call ip_fold(ii_3_1, ii_3_2, ii_3_tot)
                l_3_pc = l_from_ip(ii_3_1, ii_3_2, errString)
                en_1 = en_bz(l_1_ba, l_1_pc)
                ! Summed over.  Final band.
                do l_3_ba = 1, ba_n
                    !-----------------------------------------------------------
                    ! Analyze processes to include or exclude.
                    if (suppress_auger) then
                        if (l_1_ba .eq. l_3_ba) then
                            ! This process is intra-band, must keep
                            ! intra-band only in the spectral term.
                            cv(iv_a) = 1
                            cv(iv_e) = 0
                        else
                            ! This process is inter-band, must keep
                            ! inter-band only in the spectral term.
                            cv(iv_a) = 0
                            cv(iv_e) = 1
                        end if
                    else
                        ! Keep all processes in the spectral term.
                        cv(:) = 1
                    end if
                    !-----------------------------------------------------------
                    en_3 = en_bz(l_3_ba, l_3_pc)
                    ! Neglect the contribution of scattering events that
                    ! exchange an energy larger than the cutoff.
                    if (abs(en_3 - en_1) .lt. w_scatt_max) then
                        ! Interpolate the factors of the integrand on
                        ! the frequency mesh.
                        ! Outgoing.
                        call coeff_interp(c_l, c_r, l_l_ww, l_r_ww, en_1 - en_3)
                        w2s_l = real(pot_scr(l_lq, l_l_ww)*conjg(pot_scr(l_lq, l_l_ww)), kind=8) * sum(spdeh(l_lq, l_l_ww, :) * cv(:))
                        w2s_r = real(pot_scr(l_lq, l_r_ww)*conjg(pot_scr(l_lq, l_r_ww)), kind=8) * sum(spdeh(l_lq, l_r_ww, :) * cv(:))
                        w2s_out = c_l*w2s_l + c_r*w2s_r
                        ! Incoming.  Opposite wave vector and frequency.
                        ! Remember that |W(-q,-e)|^2 = |W(q,e)|^2.
                        call coeff_interp(c_l, c_r, l_l_ww, l_r_ww, en_3 - en_1)
                        w2s_l = real(pot_scr(l_m_lq, l_l_ww)*conjg(pot_scr(l_m_lq, l_l_ww)), kind=8) * sum(spdeh(l_m_lq, l_l_ww, :) * cv(:))
                        w2s_r = real(pot_scr(l_m_lq, l_r_ww)*conjg(pot_scr(l_m_lq, l_r_ww)), kind=8) * sum(spdeh(l_m_lq, l_r_ww, :) * cv(:))
                        w2s_in = c_l*w2s_l + c_r*w2s_r
                        ! Matrix element of the density operator.
                        ! Its squared modulus is the same for both
                        ! in and out processes.
                        dme_31 = density_matrix_element(l_1_ba, l_3_ba, l_1_pc, l_3_pc)
                        dme_sq = real(dme_31 * conjg(dme_31))
                        ! Coefficients of the single-particle rate equation.
                        w13_out = dme_sq * w2s_out
                        w13_in = dme_sq * w2s_in
                        ! Multiply by the distribution function.
                        total_out = total_out + h_pc(l_sp, l_3_ba, l_3_pc) * w13_out
                        total_in = total_in + f_pc(l_sp, l_3_ba, l_3_pc) * w13_in
                    end if
                end do
            end do
            ! Multiply by a coefficient and save the result of the integration.
            gamma_ee_out(l_sp, l_1_ba, l_1_pc) = c0 * total_out
            gamma_ee_in(l_sp, l_1_ba, l_1_pc) = c0 * total_in
        end do
        end do
        end do
!$OMP END PARALLEL DO
    end subroutine calculate_ee_gammas

    subroutine calculate_ee_gammas_2(f_pc, h_pc)
        ! Calculate the scattering rate due to Coulomb interactions.
        ! Introduce an auxiliary broadened delta for the energy conservation.
        implicit none
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc

        ! Indices of the distribution and the summations.
        integer :: l_1_pc, l_3_pc, l_lq, l_m_lq, l_1_ba, l_3_ba, l_sp
        integer :: l_ww, l_m_ww
        ! Integer projections of wave vectors.
        integer :: ii_3_1, ii_3_2
        integer, dimension(ndim) :: ii_1, ii_q, ii_3_tot
        ! Coefficients.
        double precision :: en_1, en_3
        double precision :: c0, a0
        integer, dimension(iv_n) :: cv
        double precision:: w2s_out, w2s_in, w13_in, w13_out
        ! Variables where the sum is accumulated at each wave vector, band,
        ! and spin, for outgoing and incoming processes.
        double precision :: total_out, total_in
        ! Matrix element <3| ... |1> of the density operator.
        complex(kind=8) :: dme_31
        double precision :: dme_sq
        character(len=50), parameter :: errString = "calculate_ee_gammas_2"

        ! Factor including the wave vector summation.
        c0 = (2*pi/hbar) * area_bz**2 / (2*pi)**4 / pc_n
        ! Factor including the transferred energy integral.
        c0 = c0 * ww_dw

!$OMP PARALLEL DO PRIVATE(l_1_pc, l_1_ba, l_sp, total_out, total_in) &
!$OMP & PRIVATE(l_ww, l_m_ww, l_lq) &
!$OMP & PRIVATE(ii_1, ii_q, l_m_lq, ii_3_tot, ii_3_1, ii_3_2, l_3_pc) &
!$OMP & PRIVATE(en_1, l_3_ba, cv, en_3) &
!$OMP & PRIVATE(w2s_out, w2s_in, dme_31, dme_sq, w13_out, w13_in, a0)
        ! External.  Initial wave vector.
        do l_1_pc = 1, pc_n
        ! External.  Initial band.
        do l_1_ba = 1, ba_n
        ! External.  Spin, conserved.
        do l_sp = 1, sp_n
            ! Initialize the results of the summations.
            total_out = 0.0
            total_in = 0.0
            ! Summed over.  Energy difference.
            do l_ww = 1, ww_n
            ! Summed over.  Wave vector difference.
            l_m_ww = minus_ww(l_ww)
            do l_lq = 1, lq_n
                ii_1 = ip_pc(:, l_1_pc)
                ii_q = ip_pc(:, lq(l_lq))
                ! Get the index of the opposite wave vector.
                l_m_lq = minus_lq(l_lq)
                ! Calculate the final wave vector.
                ii_3_tot = ii_1 - ii_q
                call ip_fold(ii_3_1, ii_3_2, ii_3_tot)
                l_3_pc = l_from_ip(ii_3_1, ii_3_2, errString)
                en_1 = en_bz(l_1_ba, l_1_pc)
                ! Summed over.  Final band.
                do l_3_ba = 1, ba_n
                    !-----------------------------------------------------------
                    ! Analyze processes to include or exclude.
                    if (suppress_auger) then
                        if (l_1_ba .eq. l_3_ba) then
                            ! This process is intra-band, must keep
                            ! intra-band only in the spectral term.
                            cv(iv_a) = 1
                            cv(iv_e) = 0
                        else
                            ! This process is inter-band, must keep
                            ! inter-band only in the spectral term.
                            cv(iv_a) = 0
                            cv(iv_e) = 1
                        end if
                    else
                        ! Keep all processes in the spectral term.
                        cv(:) = 1
                    end if
                    !-----------------------------------------------------------
                    en_3 = en_bz(l_3_ba, l_3_pc)
                    ! Neglect the contribution of scattering events that
                    ! exchange an energy larger than the cutoff.
                    if (abs(en_3 - en_1) .lt. w_scatt_max) then
                        ! The factors of the integrand at the exchanged energy.
                        ! Outgoing.
                        w2s_out = real(pot_scr(l_lq, l_ww)*conjg(pot_scr(l_lq, l_ww)), kind=8) * sum(spdeh(l_lq, l_ww, :) * cv(:))
                        ! Incoming.  Opposite wave vector and frequency.
                        ! Remember that |W(-q,-e)|^2 = |W(q,e)|^2.
                        w2s_in = real(pot_scr(l_m_lq, l_m_ww)*conjg(pot_scr(l_m_lq, l_m_ww)), kind=8) * sum(spdeh(l_m_lq, l_m_ww, :) * cv(:))
                        ! Matrix element of the density operator.
                        ! Its squared modulus is the same for both
                        ! in and out processes.
                        dme_31 = density_matrix_element(l_1_ba, l_3_ba, l_1_pc, l_3_pc)
                        dme_sq = real(dme_31 * conjg(dme_31))
                        ! Coefficients of the single-particle rate equation.
                        w13_out = dme_sq * w2s_out
                        w13_in = dme_sq * w2s_in
                        ! Delta of energy conservation.  w = e1 - e3.
                        a0 = delta_ee(l_1_ba, l_1_pc, l_3_ba, l_3_pc, ww(l_ww))
                        ! Multiply by the distribution function.
                        total_out = total_out + h_pc(l_sp, l_3_ba, l_3_pc) * w13_out * a0
                        total_in = total_in + f_pc(l_sp, l_3_ba, l_3_pc) * w13_in * a0
                    end if
                end do
            end do
            end do
            ! Multiply by a coefficient and save the result of the integration.
            gamma_ee_out(l_sp, l_1_ba, l_1_pc) = c0 * total_out
            gamma_ee_in(l_sp, l_1_ba, l_1_pc) = c0 * total_in
        end do
        end do
        end do
!$OMP END PARALLEL DO
    end subroutine calculate_ee_gammas_2

    subroutine calculate_ee_rate(f1_pc, f_pc, h_pc)
        ! Calculate the scattering rate due to Coulomb interactions.
        ! dfₖ/dt = −fₖ Γᵒᵘᵗₖ + (1 − fₖ) Γⁱⁿₖ
        implicit none
        ! Time-derivative of the electron distribution.
        double precision, dimension(:, :, :), intent(out) :: f1_pc
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc

#if 0
        ! Check sanity of electron distribution.
        double precision, dimension(:, :, :), allocatable :: f_test
        allocate( f_test(sp_n, ba_n, pc_n) )

        f_test = 0.0
        if (any(f_pc .lt. f_test)) then
            LOGME trim(msg_red("calculate_ee_rate: Negative values in electron distribution."))
            print *, minval(f_pc)
            stop
        end if
        f_test = 1.0
        if (any(f_pc .gt. f_test)) then
            LOGME trim(msg_red("calculate_ee_rate: Large values in electron distribution."))
            print *, maxval(f_pc)
            stop
        end if
#endif

        ! Calculate distribution-dependent quantities.
        if (screening_model .eq. SCR_DYN .or. screening_model .eq. SCR_STA) then
            call calculate_lindhard(f_pc)
        else if (screening_model .eq. SCR_TTF) then 
            call calculate_thomasfermi_q(f_pc)
        end if
        call calculate_screened_potential
        if (broadening_ee_model .eq. BEM_UTD) then
            call calculate_eta_ee_weighted(f_pc)
        end if
        call calculate_spectral_density(f_pc, h_pc)
        if (gamma_freq_integral) then
            ! Extra frequency integral to smooth the spectral function.
            call calculate_ee_gammas_2(f_pc, h_pc)
        else
            ! Interpolate the spectral function on the frequency mesh.
            call calculate_ee_gammas(f_pc, h_pc)
        end if

        ! Total scattering rate at given wave vector.
        f1_pc = -f_pc*gamma_ee_out + h_pc*gamma_ee_in

        ! call debug_rate(f1_pc, f_pc, h_pc)  ! debug
    end subroutine calculate_ee_rate


    subroutine debug_rate(f1_pc, f_pc, h_pc)
        implicit none
        double precision, dimension(:, :, :), intent(in) :: f1_pc, f_pc, h_pc
        integer :: l_pc, l_ba, l_sp
        character(len=50) :: fmtStr

        write (fmtStr, '(A,I2,A)') "(", sp_n * ba_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-debug-f.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmtStr) (( f_pc(l_sp,l_ba,l_pc), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        write (fmtStr, '(A,I2,A)') "(", sp_n * ba_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-debug-h.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmtStr) (( h_pc(l_sp,l_ba,l_pc), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        write (fmtStr, '(A,I2,A)') "(", sp_n * ba_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-debug-f1.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmtStr) (( f1_pc(l_sp,l_ba,l_pc), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        write (fmtStr, '(A,I2,A)') "(", sp_n * ba_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-debug-gamma-out.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmtStr) (( gamma_ee_out(l_sp,l_ba,l_pc), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        write (fmtStr, '(A,I2,A)') "(", sp_n * ba_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-debug-gamma-in.csv", status="unknown")
        do l_pc = 1, pc_n
            write (30, fmtStr) (( gamma_ee_in(l_sp,l_ba,l_pc), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        stop
    end subroutine debug_rate

    subroutine print_lindhard(f_pc, fileBase)
        implicit none
        character(len=*), intent(in) :: fileBase
        ! Electron distribution.
        double precision, dimension(:, :, :), intent(in) :: f_pc
        integer :: l_ww, l_lq, l_uu
        double precision :: qeta
        complex(kind=8) :: total
        character(len=50) :: fileName, fmtString

        if (screening_model .ne. SCR_DYN .and. screening_model .ne. SCR_STA) return

        ! Calculate the Lindhard function.
        ! This might be redundant, however, it is needed at the first time.
        call calculate_lindhard(f_pc)

        ! Calculate the angular average of the Lindhard function.
        ! For each frequency, calculate the average over the angle.
        ! The Lindhard function is represented on a mesh.  Smear the points of
        ! the mesh with a Gaussian function to perform the integration.
        ! The broadening of the smearing is proportional to the distance between
        ! q points wv_dk1.
        ! The integration tile is (area_bz / pc_n).
        qeta = 0.3 * wv_dk1
        do l_uu = 1, uu_n
            do l_ww = 1, ww_n
                total = 0.0
                do l_lq = 1, lq_n
                    total = total + chi0(l_lq,l_ww) * exp(-(norm_lq(l_lq)- uu(l_uu))**2 / (2.0 * qeta**2))
                end do
                chi0az(l_uu,l_ww) = total * (area_bz / pc_n) / (sqrt(2.0 * pi) * qeta) / (2.0 * pi * uu(l_uu))
            end do
        end do

        write (fmtString, '(A,I04,A)') "(", uu_n - 1, "(F15.8, ', '),F15.8)"

        write (fileName, '(A,A)') trim(fileBase), "-re.csv"
        open (unit=30, file=fileName, status="unknown")
        do l_ww = 1, ww_n
            write (30, fmtString) (real(chi0az(l_uu, l_ww)), l_uu = 1, uu_n)
        end do
        close (unit=30)

        write (fileName, '(A,A)') trim(fileBase), "-im.csv"
        open (unit=30, file=fileName, status="unknown")
        do l_ww = 1, ww_n
            write (30, fmtString) (imag(chi0az(l_uu, l_ww)), l_uu = 1, uu_n)
        end do
        close (unit=30)

#if 0
        write (fmtString, '(A,I04,A)') "(", lq_n - 1, "(F15.8, ', '),F15.8)"

        write (fileName, '(A,A)') trim(fileBase), "-chi0-re.csv"
        open (unit=30, file=fileName, status="unknown")
        do l_ww = 1, ww_n
            write (30, fmtString) (real(chi0(l_lq, l_ww)), l_lq = 1, lq_n)
        end do
        close (unit=30)

        write (fileName, '(A,A)') trim(fileBase), "-chi0-im.csv"
        open (unit=30, file=fileName, status="unknown")
        do l_ww = 1, ww_n
            write (30, fmtString) (imag(chi0(l_lq, l_ww)), l_lq = 1, lq_n)
        end do
        close (unit=30)
#endif

    end subroutine print_lindhard


    subroutine print_spectral_density(f_pc, fileBase)
        implicit none
        character(len=*), intent(in) :: fileBase
        ! Electron distribution.
        double precision, dimension(:, :, :), intent(in) :: f_pc
        integer :: l_p
        double precision, dimension(:,:,:), allocatable :: h_pc
        double precision, dimension(:,:), allocatable :: spdeh_path
        character(len=50) :: fileName, fmtString

        allocate (h_pc(sp_n, ba_n, pc_n))
        h_pc = 1.0 - f_pc
        call calculate_spectral_density(f_pc, h_pc)
        deallocate (h_pc)

        allocate (spdeh_path(l_pc_path_n, ww_n))

        ! Extract the spectral density on a high-symmetry path in the BZ.
        do l_p = 1, l_pc_path_n
             spdeh_path(l_p,:) = sum(spdeh(l_pc_path(l_p),:,:),2)
        end do

        write (fmtString, '(A,I04,A)') "(", ww_n - 1, "(F15.8, ', '),F15.8)"

        write (fileName, '(A,A)') trim(fileBase), ".csv"
        open (unit=30, file=fileName, status="unknown")
        do l_p = 1, l_pc_path_n
            write (30, fmtString) spdeh_path(l_p,:)
        end do
        close (unit=30)

        deallocate (spdeh_path)

    end subroutine print_spectral_density

end module ELECTRON_SCATTERING

