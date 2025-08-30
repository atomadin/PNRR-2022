module TIME_EVOLUTION
    use phys_const
    use time_monitor
    use color_msg
    use band_structure
    use electron_distribution
    use electron_scattering
    use phonon_modes
    use phonon_distribution
    use electron_phonon_scattering
    use differential_transmission

    ! Mesh of save times.
    double precision, dimension(:), allocatable :: tt
    ! Number of save times.
    integer :: tt_n
    ! Initial and final times.
    double precision :: tt_min, tt_max
    ! Number of integration steps between save time.
    integer :: st_n
    ! Time step size.
    double precision :: st_w
    ! Times when a snapshot of the distribution is saved.
    ! The actual number of times is chosen by the user, this is the maximum.
    integer, parameter :: xt_n = 20
    double precision, dimension(xt_n) :: xt
    ! Whether to print the distributions at all times.
    logical :: snapshot_all_times

    ! Electron and hole distributions at the current time in the evolution.
    double precision, dimension(:, :, :), allocatable :: f_t_pc, h_t_pc
    ! Phonon distribution at the current time in the evolution.
    double precision, dimension(:, :), allocatable :: u_t_pc

    ! Observables.
    ! Particle energy density distribution
    double precision, dimension(:, :), allocatable :: rf_t_eb, rh_t_eb
    ! Particle averaged distribution.
    double precision, dimension(:, :), allocatable :: f_t_eb, h_t_eb
    ! Observables in time.
    ! Particle density in time [nm^-2].
    double precision, dimension(:, :, :), allocatable :: rf_tt, rh_tt
    ! Phonon average occupation in time [dimensionless].
    double precision, dimension(:, :), allocatable :: ua_tt
    ! Electron energy density in time [eV nm^-2].
    double precision, dimension(:), allocatable :: re_tt
    ! Phonon energy density in time [eV nm^-2].
    double precision, dimension(:,:), allocatable :: ue_tt
    ! Particle current density in time [fs^-1 nm^-1].
    double precision, dimension(:, :, :, :), allocatable :: rj_tt
    ! Momentum density in time [eV fs nm^-3].
    double precision, dimension(:, :, :, :), allocatable :: rp_tt
    ! Differential transmission.
    double precision, dimension(:,:), allocatable :: dt_wp_tt
    ! Thomas-Fermi wave vector [nm^-1].
    double precision, dimension(:), allocatable :: tf_q_tt
    ! Electron energy broadening [eV].
    double precision, dimension(:), allocatable :: eta_ee_tt

    private

    public :: init_time_evolution, calculate_time_evolution
    public :: analyze_energy_conservation

contains

    subroutine init_time_evolution
        implicit none

        character(len=200) :: msg1, msg2
        integer :: l_tt, l_xt
        ! User-friendly parameters for the input file.
        integer :: num_save_times
        integer, dimension(xt_n) :: snapshot_times_fs
        double precision :: time_max_fs, desired_timestep_fs
        namelist /time_evolution_params/ time_max_fs, desired_timestep_fs, &
            & num_save_times, snapshot_all_times, snapshot_times_fs
        namelist /time_evolution_values/ st_w

        ! Initialize the array of snapshot times to negative values.
        snapshot_times_fs = -1.0

        ! Read parameters from input file.
        open (unit=37, file="input.txt", status="old")
        read (37, nml=time_evolution_params)
        close (37)

        ! Set parameters.
        tt_min = 0.000
        tt_max = time_max_fs
        tt_n = num_save_times
        st_n = ceiling(tt_max / (tt_n - 1.0) / desired_timestep_fs)
        write(msg1, "(I8)") st_n
        write(msg2, "(2A)") "init_time_evolution: Number of steps between save times: ", trim(adjustl(msg1))
        LOGME trim(msg_aqua(msg2))

        ! Set the first snapshot time to be the initial time.
        xt(1) = tt_min
        xt(2:xt_n) = snapshot_times_fs(1:(xt_n-1))

        ! Initialize mesh of save times.
        allocate (tt(tt_n))
        tt = (/(tt_min + (tt_max - tt_min)*(l_tt - 1.0)/(tt_n - 1.0), l_tt=1, tt_n)/)
        ! Time step size.
        st_w = (tt(2) - tt(1))/st_n

        ! Electron and hole distributions.
        allocate (f_t_pc(sp_n, ba_n, pc_n))
        allocate (h_t_pc(sp_n, ba_n, pc_n))

        ! Phonon distribution.
        allocate (u_t_pc(ph_n, pc_n))

        ! Observables.
        allocate (rf_t_eb(sp_n, eb_n))
        allocate (rh_t_eb(sp_n, eb_n))
        allocate (f_t_eb(sp_n, eb_n))
        allocate (h_t_eb(sp_n, eb_n))

        ! Observables in time.
        allocate (rf_tt(sp_n, ba_n, tt_n))
        allocate (rh_tt(sp_n, ba_n, tt_n))
        allocate (ua_tt(ph_n, tt_n))
        allocate (ue_tt(ph_n, tt_n))
        allocate (re_tt(tt_n))
        allocate (rj_tt(2, sp_n, ba_n, tt_n))
        allocate (rp_tt(2, sp_n, ba_n, tt_n))

        ! Checks and monitors.
        allocate (tf_q_tt(tt_n))  ! Thomas-Fermi wave vector [nm^-1]
        allocate (eta_ee_tt(tt_n))  ! electron energy broadening [eV]

        ! To calculate the differential transmission we must set
        ! the parameters of the equilibrium distribution.
        call init_differential_transmission(mu_eq, tempK_eq)
        allocate (dt_wp_tt(wp_n, tt_n))

        ! Init time monitors.
        call time_monitor_set(1, "electrons")
        call time_monitor_set(2, "phonons")
        call time_monitor_set(3, "observables")
        call time_monitor_set(4, "total")

        ! Save derived values.
        open (unit=VALUNIT, file=VALFILE, status="old", position="append")
        write (VALUNIT, nml=time_evolution_values)
        close (VALUNIT)

        ! Print snapshot times.
        open (unit=30, file="out-snap-times.csv", status="unknown")
        do l_xt = 1, xt_n
            if (sign(1.0d0,xt(l_xt)) .gt. 0.0) then
                write (30, "(F15.8)") xt(l_xt)
            else
                ! Stop printing once the first negative value is found.
                exit
            end if
        end do
        close (unit=30)

        ! Print save time mesh.
        open (unit=30, file="out-evol-times.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, "(F15.8)") tt(l_tt)
        end do
        close (unit=30)
    end subroutine init_time_evolution

    subroutine calculate_rate_electron(f1_tot_pc, u_pc, f_pc, h_pc)
        implicit none
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        ! Phonon distribution.
        double precision, dimension(:, :), intent(in) :: u_pc

        ! Total time-derivatives of the electron distribution.
        double precision, dimension(:, :, :), intent(out) :: f1_tot_pc

        ! Time-derivatives of the electron distribution due to several processes.
        double precision, dimension(:, :, :), allocatable :: f1_ee_pc, f1_ep_pc

        allocate (f1_ee_pc(sp_n, ba_n, pc_n))
        allocate (f1_ep_pc(sp_n, ba_n, pc_n))

        ! Electron-electron scattering.
        call time_monitor_start(1)
#if WITHCOULOMB
        call calculate_ee_rate(f1_ee_pc, f_pc, h_pc)
#else
        if (.false.) print*, shape(f_pc), shape(h_pc)
        f1_ee_pc = 0.0
#endif
        call time_monitor_pause(1)

        ! Electron-phonon scattering.
        call time_monitor_start(2)
#if WITHPHONONS
        call calculate_ebp_rate(f1_ep_pc, u_pc, f_pc, h_pc)
#else
        if (.false.) print *, shape(u_pc), shape(f_pc), shape(h_pc)
        f1_ep_pc = 0.0
#endif
        call time_monitor_pause(2)

        ! Sum all the rates.
        f1_tot_pc = f1_ee_pc + f1_ep_pc

        deallocate (f1_ee_pc)
        deallocate (f1_ep_pc)
    end subroutine calculate_rate_electron

    subroutine calculate_rate_phonon(u1_tot_pc, u_pc, f_pc, h_pc)
        implicit none
        ! Electron and hole distributions.
        double precision, dimension(:, :, :), intent(in) :: f_pc, h_pc
        ! Phonon distributions.
        double precision, dimension(:, :), intent(in) :: u_pc

        ! Total time-derivatives of the phonon distribution.
        double precision, dimension(:, :), intent(out) :: u1_tot_pc

        ! Time-derivative of the phonon distribution due to e-h scattering.
        double precision, dimension(:, :), allocatable :: u1_pbe_pc
        ! Time-derivative of the phonon distribution due to relaxation to the
        ! substrate temperature by nonlinear processes.  Phenomenological.
        double precision, dimension(:,:), allocatable :: u1_sbs_pc

        allocate (u1_pbe_pc(ph_n, pc_n))
        allocate (u1_sbs_pc(ph_n, pc_n))

        ! Electron-phonon scattering.
        call time_monitor_start(2)
#if WITHPHONONS
        call calculate_pbe_rate(u1_pbe_pc, u_pc, f_pc, h_pc)
#else
        u1_pbe_pc = 0.0
#endif
        call time_monitor_pause(2)

        ! Phenomenological relaxation to the substrate temperature.
#if WITHSUBSTRATE
        call calculate_sbs_rate(u1_sbs_pc, u_pc)
#else
        u1_sbs_pc = 0.0
#endif

        ! Sum all the rates.
        u1_tot_pc = u1_pbe_pc + u1_sbs_pc
    end subroutine calculate_rate_phonon

    subroutine calculate_one_step(u_0_pc, f_0_pc, h_0_pc, t_0)
        ! The Runge-Kutta integration step.  This is generic and should not be
        ! modified.  Taken from Numerical Recipes.
        implicit none
        double precision, intent(in) :: t_0
        ! Electron and hold distributions at the beginning of each integration step.
        double precision, dimension(:, :, :), intent(in) :: f_0_pc, h_0_pc
        double precision, dimension(:, :), intent(in) :: u_0_pc
        ! Intermediate times for the Runge-Kutta step.
        double precision :: st_w_2, st_w_6, t_m
        ! Time-derivatives of the electron distribution.
        double precision, dimension(:, :, :), allocatable :: f1_0_pc, f1_1_pc, f1_2_pc
        ! Time-derivatives of the phonon distribution.
        double precision, dimension(:, :), allocatable :: u1_0_pc, u1_1_pc, u1_2_pc

        allocate (f1_0_pc(sp_n, ba_n, pc_n))
        allocate (f1_1_pc(sp_n, ba_n, pc_n))
        allocate (f1_2_pc(sp_n, ba_n, pc_n))
        allocate (u1_0_pc(ph_n, pc_n))
        allocate (u1_1_pc(ph_n, pc_n))
        allocate (u1_2_pc(ph_n, pc_n))

        st_w_2 = st_w/2.0
        st_w_6 = st_w/6.0
        t_m = t_0 + st_w_2

        ! Integrate from t_0 to t_m, with the derivative at t_0.
        call calculate_rate_electron(f1_0_pc, u_0_pc, f_0_pc, h_0_pc)
        call calculate_rate_phonon(u1_0_pc, u_0_pc, f_0_pc, h_0_pc)
        f_t_pc = f_0_pc + st_w_2*f1_0_pc
        h_t_pc = 1.0 - f_t_pc
        u_t_pc = u_0_pc + st_w_2*u1_0_pc

        ! Integrate again from t_0 to t_m, with the derivative at t_m.
        call calculate_rate_electron(f1_1_pc, u_t_pc, f_t_pc, h_t_pc)
        call calculate_rate_phonon(u1_1_pc, u_t_pc, f_t_pc, h_t_pc)
        f_t_pc = f_0_pc + st_w_2*f1_1_pc
        h_t_pc = 1.0 - f_t_pc
        u_t_pc = u_0_pc + st_w_2*u1_1_pc

        ! Integrate from t_0 to t, with the new derivative at t_m.
        call calculate_rate_electron(f1_2_pc, u_t_pc, f_t_pc, h_t_pc)
        call calculate_rate_phonon(u1_2_pc, u_t_pc, f_t_pc, h_t_pc)
        f_t_pc = f_0_pc + st_w*f1_2_pc
        h_t_pc = 1.0 - f_t_pc
        u_t_pc = u_0_pc + st_w*u1_2_pc
        ! Sum the two derivatives at t_m.
        f1_1_pc = f1_1_pc + f1_2_pc
        u1_1_pc = u1_1_pc + u1_2_pc

        ! Calculate the derivative at t.
        call calculate_rate_electron(f1_2_pc, u_t_pc, f_t_pc, h_t_pc)
        call calculate_rate_phonon(u1_2_pc, u_t_pc, f_t_pc, h_t_pc)
        ! Integrate from t_0 to t, with a combination of derivatives.
        f_t_pc = f_0_pc + st_w_6*(f1_0_pc + 2.0*f1_1_pc + f1_2_pc)
        h_t_pc = 1.0 - f_t_pc
        u_t_pc = u_0_pc + st_w_6*(u1_0_pc + 2.0*u1_1_pc + u1_2_pc)

        deallocate (f1_0_pc)
        deallocate (f1_1_pc)
        deallocate (f1_2_pc)
        deallocate (u1_0_pc)
        deallocate (u1_1_pc)
        deallocate (u1_2_pc)
    end subroutine calculate_one_step

    subroutine electron_distro_extremes(f_pc, h_pc)
        implicit none
        double precision, dimension(:, :, :), intent(inout) :: f_pc, h_pc
        integer :: l_pc, l_ba, l_sp

        do l_pc = 1, pc_n
            do l_ba = 1, ba_n
                do l_sp = 1, sp_n
                    f_pc(l_sp, l_ba, l_pc) = min(1.0,max(0.0, f_pc(l_sp, l_ba, l_pc)))
                    h_pc(l_sp, l_ba, l_pc) = min(1.0,max(0.0, h_pc(l_sp, l_ba, l_pc)))
                end do
            end do
        end do
    end subroutine electron_distro_extremes

    subroutine calculate_time_evolution
        implicit none
        character(len=200) :: msg1, msg2, msg3
        integer :: l_tt, l_st, c_xt
        double precision :: t_0
        ! Electron, hole, phonon distributions at the beginning of each
        ! integration step.
        double precision, dimension(:, :, :), allocatable :: f_0_pc, h_0_pc
        double precision, dimension(:, :), allocatable :: u_0_pc

        allocate (f_0_pc(sp_n, ba_n, pc_n))
        allocate (h_0_pc(sp_n, ba_n, pc_n))
        allocate (u_0_pc(ph_n, pc_n))

        ! Initialize the distributions.
        f_t_pc = f_i_pc
        h_t_pc = h_i_pc
        u_t_pc = u_i_pc

        ! Counter of performed snapshots printouts.
        c_xt = 0
        ! Evolve *until* tt(l_tt).
        call time_monitor_start(4)
        do l_tt = 1, tt_n
            ! Run evolution, unless it is the first save time.
            if (l_tt .gt. 1) then
                do l_st = 1, st_n
                    ! The start time for this step.
                    t_0 = tt(l_tt - 1) + (l_st - 1)*st_w
                    f_0_pc = f_t_pc
                    h_0_pc = h_t_pc
                    u_0_pc = u_t_pc
                    ! Print a snapshot of the distribution if necessary.
                    call print_snapshot_if_time(c_xt, t_0)
                    ! Evolve from f_0_pc at t0 to the new f_t_pc after st_w.
                    call calculate_one_step(u_0_pc, f_0_pc, h_0_pc, t_0)
                    call electron_distro_extremes(f_t_pc, h_t_pc)
#if STOPAFTER60S
                    ! Stop the evolution after a certain time.
                    if (time_monitor_read(4) .gt. 60) then
                        LOGME trim(msg_blue("calculate_time_evolution: Stop after 60s."))
                        call print_snapshot_at_stop
                        stop
                    end if
#endif
                end do
            end if

            ! Calculate observables.
            call time_monitor_start(3)
            call calculate_obs_at_save_time(l_tt)
            call time_monitor_pause(3)

            ! Periodically log the progress.
            if ((l_tt .gt. 1) .and. ((mod(l_tt,ceiling(tt_n/10.0)) .eq. 0) .or. (l_tt .eq. tt_n))) then
                ! Log the progress.
                write(msg1, '(I8)') l_tt
                write(msg2, '(I8)') tt_n
                write(msg3, '(4A)') "evolution_run: Done evolving until ", trim(adjustl(msg1)), " of ", trim(adjustl(msg2))
                LOGME trim(msg_blue(msg3))
                ! Print time monitors.
                call time_monitor_print(1)
                call time_monitor_print(2)
                call time_monitor_print(3)
            end if

            ! End of the evolution between two save times.
        end do

        deallocate (f_0_pc)
        deallocate (h_0_pc)
        deallocate (u_0_pc)

        ! Print the evolution of the observables on the time mesh.
        call print_evolution_observables
    end subroutine calculate_time_evolution

    subroutine calculate_obs_at_save_time(l_tt)
        ! Calculate and print observables at this stroboscopic time.
        implicit none
        integer, intent(in) :: l_tt
        ! Quantities to be saved and printed at the end of the evolution.
        call calculate_electron_density(rf_tt(:, :, l_tt), rh_tt(:, :, l_tt), f_t_pc, h_t_pc)
        call calculate_phonon_avg_occup(ua_tt(:, l_tt), u_t_pc)
        call calculate_energy_density(re_tt(l_tt), f_t_pc, h_t_pc)
        call calculate_phonon_energy_density(ue_tt(:,l_tt), u_t_pc)
        call calculate_current_density(rj_tt(:, :, :, l_tt), f_t_pc)
        call calculate_momentum_density(rp_tt(:, :, :, l_tt), f_t_pc)
        call calculate_diff_transm(dt_wp_tt(:, l_tt), f_t_pc)
        ! Checks and monitors.
        tf_q_tt(l_tt) = thomasfermi_q
        eta_ee_tt(l_tt) = eta_ee
    end subroutine calculate_obs_at_save_time

    subroutine print_snapshot_if_time(c_xt, t0)
        implicit none
        ! Notice: there might be a difference between the snapshot time
        ! requested by the user and the actual used.  However, the difference
        ! is smaller than the integration timestep and is thus negligible.

        ! Counter of performed snapshots.
        integer, intent(inout) :: c_xt
        ! The current time.
        double precision, intent(in) :: t0
        ! Index in the array of target times.
        integer :: l_xt
        character(len=50) :: fileName

        ! Check that we have not reached the maximum number of snapshots, and
        ! that there is a valid next target time.
        if ((c_xt .lt. xt_n) .and. (sign(1.0d0,xt(c_xt+1)) .gt. 0.0)) then
            l_xt = c_xt+1
            ! Check that we are going to cross the target snapshot time.
            ! within this step.
            if ((t0+st_w) .gt. xt(l_xt)) then
                ! Quantities to be printed at this time.
                call calculate_el_distro_energy(rf_t_eb, rh_t_eb, f_t_pc, h_t_pc)
                call calculate_distro_average(f_t_eb, h_t_eb, f_t_pc, h_t_pc)
                ! Printouts.
                write (fileName, '(A,I4.4,A)') "out-snap-", l_xt, "-el-distro-k.csv"
                call print_particle_distro(f_t_pc, fileName)
                write (fileName, '(A,I4.4,A)') "out-snap-", l_xt, "-ph-distro-k.csv"
                call print_phonon_distro(u_t_pc, fileName)
                write (fileName, '(A,I4.4,A)') "out-snap-", l_xt, "-el-distro-e.csv"
                call print_particle_energy_distro(rf_t_eb, fileName)
                write (fileName, '(A,I4.4,A)') "out-snap-", l_xt, "-el-avg-e.csv"
                call print_averaged_particle_distro(f_t_eb, fileName)
                write (fileName, '(A,I4.4,A)') "out-snap-", l_xt, "-lindhard"  ! base only
                call print_lindhard(f_t_pc, fileName)
                write (fileName, '(A,I4.4,A)') "out-snap-", l_xt, "-spectral"  ! base only
                call print_spectral_density(f_t_pc, fileName)
                ! Increase the counter of performed snapshots.
                c_xt = c_xt + 1
            end if
        end if
    end subroutine print_snapshot_if_time

    subroutine print_snapshot_at_stop
        implicit none
        call print_ebp_gammas("out-stop-ebp-gammas")
    end subroutine print_snapshot_at_stop

    subroutine print_evolution_observables
        ! Print the observables that have been calculated during the
        ! time-evolution.
        implicit none
        integer :: l_sp, l_ba, l_ph, l_tt, l_wp, col_n
        character(len=50) :: fmtSB, fmtStr

        ! Format string for spin- and band-resolved quantities.
        col_n = ba_n*sp_n
        write (fmtSB, '(A,I2,A)') "(", col_n - 1, "(F15.8, ', '),F15.8)"

        ! Electron (hole) density in conduction (valence) band.
        open (unit=30, file="out-evol-el-dens.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtSB) ( (/ rh_tt(l_sp, 1, l_tt), rf_tt(l_sp, 2, l_tt) /), l_sp=1, sp_n)
        end do
        close (unit=30)

        ! Phonon average occupation, mode-resolved.
        write (fmtStr, '(A,I2,A)') "(", ph_n - 1, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-evol-ph-avg-occ.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtStr) (ua_tt(l_ph, l_tt), l_ph=1, ph_n)
        end do
        close (unit=30)

        ! Energy in the electron and phonon subsystems.
        write (fmtStr, '(A,I2,A)') "(", 5, "(F15.8, ', '),F15.8)"
        open (unit=30, file="out-evol-energy.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtStr) re_tt(l_tt), sum(ue_tt(:,l_tt)), (ue_tt(l_ph,l_tt), l_ph = 1, ph_n)
        end do
        close (unit=30)

        ! Electron current and momentum density.
        open (unit=30, file="out-evol-el-current.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtSB) ((sqrt(rj_tt(1, l_sp, l_ba, l_tt)**2 + rj_tt(2, l_sp, l_ba, l_tt)**2), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        open (unit=30, file="out-evol-el-current-dir.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtSB) ((atan2(rj_tt(2, l_sp, l_ba, l_tt), rj_tt(1, l_sp, l_ba, l_tt)), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        open (unit=30, file="out-evol-el-momentum.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtSB) ((sqrt(rp_tt(1, l_sp, l_ba, l_tt)**2 + rp_tt(2, l_sp, l_ba, l_tt)**2), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        open (unit=30, file="out-evol-el-momentum-dir.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtSB) ((atan2(rp_tt(2, l_sp, l_ba, l_tt), rp_tt(1, l_sp, l_ba, l_tt)), l_ba=1, ba_n), l_sp=1, sp_n)
        end do
        close (unit=30)

        ! Differential transmission.
        write (fmtStr, '(A,I2,A)') "(", wp_n - 1, "(E12.4, ', '),E12.4)"
        open (unit=30, file="out-evol-diff-transm.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtStr) (dt_wp_tt(l_wp,l_tt), l_wp = 1, wp_n)
        end do
        close (unit=30)

        ! Thomas-Fermi wave vector.
        write (fmtStr, '(A,I2,A)') "(", wp_n - 1, "(E12.4, ', '),E12.4)"
        open (unit=30, file="out-evol-thomasfermi.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtStr) tf_q_tt(l_tt)
        end do
        close (unit=30)

        ! Electron energy broadening.
        write (fmtStr, '(A,I2,A)') "(", wp_n - 1, "(E12.4, ', '),E12.4)"
        open (unit=30, file="out-evol-eta-ee.csv", status="unknown")
        do l_tt = 1, tt_n
            write (30, fmtStr) eta_ee_tt(l_tt)
        end do
        close (unit=30)

    end subroutine print_evolution_observables

    subroutine analyze_energy_conservation
        ! Calculate a few steps of the evolution and print the values of the
        ! electron energy.
        implicit none

        integer :: i_step
        integer, parameter :: step_n = 10
        double precision, parameter :: step_w = 1.0  ! fs
        double precision :: t_0
        ! Electron, hole, phonon distributions at the beginning of each
        ! integration step.
        double precision, dimension(:, :, :), allocatable :: f_0_pc, h_0_pc
        double precision, dimension(:, :), allocatable :: u_0_pc
        ! Values of the energy density.
        double precision, dimension(:), allocatable :: en_dens

        allocate (f_0_pc(sp_n, ba_n, pc_n))
        allocate (h_0_pc(sp_n, ba_n, pc_n))
        allocate (u_0_pc(ph_n, pc_n))
        allocate (en_dens(step_n))

        ! Initialize the distributions.
        f_t_pc = f_i_pc
        h_t_pc = h_i_pc
        u_t_pc = u_i_pc
    
        LOGME trim("analyze_energy_conservation: Running a few steps.")
        do i_step = 1, step_n
            write (*, '(A, I)') "Step: ", i_step
            ! The time for this step.
            t_0 = (i_step - 1) * step_w
            f_0_pc = f_t_pc
            h_0_pc = h_t_pc
            u_0_pc = u_t_pc
            call calculate_one_step(u_0_pc, f_0_pc, h_0_pc, t_0)
            call electron_distro_extremes(f_t_pc, h_t_pc)
            call calculate_energy_density(en_dens(i_step), f_t_pc, h_t_pc)
        end do

        ! Energy in the electron and phonon subsystems.
        do i_step = 1, step_n
            write (*, '(A,F8.3)') "E_e[%] = ", en_dens(i_step) / en_dens(1) * 100.0
        end do

        deallocate (en_dens)
        deallocate (f_0_pc)
        deallocate (h_0_pc)
        deallocate (u_0_pc)

    end subroutine analyze_energy_conservation

end module TIME_EVOLUTION

