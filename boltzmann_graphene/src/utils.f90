module math
    implicit none
    double precision, parameter :: pi = 3.141592653589793  ! pi
    double precision, parameter :: pi_sqrt = 1.7724538509055159  ! square root of pi
    double precision, parameter :: s3 = 1.7320508075688772  ! square root of 3

    public

contains

    subroutine bubble_sort(arr,n)
        implicit none
        double precision, dimension(:), intent(inout) :: arr
        integer, intent(in) :: n
        integer :: i, j
        double precision :: temp
        logical :: swapped

        do i = 1, n-1
            swapped = .false.
            do j = 1, n - i
                if (arr(j) > arr(j+1)) then
                    temp = arr(j)
                    arr(j) = arr(j+1)
                    arr(j+1) = temp
                    swapped = .true.
                end if
            end do
            if (.not. swapped) exit
        end do
    end subroutine bubble_sort

    recursive subroutine quicksort(arr, left, right)
        implicit none
        integer, intent(in) :: left, right
        double precision, intent(inout) :: arr(:)
        integer :: i, j
        double precision :: pivot, temp

        if (left >= right) return

        pivot = arr((left + right) / 2)
        i = left
        j = right

        do
            do while (arr(i) < pivot)
                i = i + 1
            end do
            do while (arr(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do

        if (left < j) call quicksort(arr, left, j)
        if (i < right) call quicksort(arr, i, right)
    end subroutine quicksort

end module math

module phys_const
    use math
    
    implicit none
    ! Constants with 3 significant digits.
    ! Reduced Planck's constant [eV fs].
    double precision, parameter :: hbar = 0.658
    ! Planck's constant [eV fs].
    double precision, parameter :: h_planck = hbar * 2.0 * pi
    ! Speed of light [nm/fs].
    double precision, parameter :: c_light = 300.0
    ! Planck's constant times speed of light [eV nm].
    double precision, parameter :: hc1240 = h_planck * c_light
    ! Boltzmann constant [eV / K].
    double precision, parameter :: k_boltz = 1.0/11600.0
    ! Electron charge squared [eV nm].
    double precision, parameter :: e_sq = 1.43
    ! One Joule in electronvolts.
    double precision, parameter :: Joule_eV = 6.241509d18

contains

    pure function bose_einstein(en, mu, temp) result( n )
        ! Bose-Einstein distribution.
        implicit none
        double precision :: n
        double precision, intent(in) :: en
        double precision, intent(in) :: mu
        double precision, intent(in) :: temp

        n = 1.0/(exp((en - mu)/temp) - 1.0)
    end function bose_einstein

    pure function fermi_dirac(en, mu, temp) result( f )
        ! Fermi-Dirac distribution.
        implicit none
        double precision :: f
        double precision, intent(in) :: en
        double precision, intent(in) :: mu
        double precision, intent(in) :: temp
        double precision :: x
        x = (en - mu) / temp
        if (x .lt. -20.0) then
            f = 1.0
        else if (x .gt. 20.0) then
            f = 0.0
        else
            f = 1.0/(exp(x) + 1.0)
        end if
    end function fermi_dirac

end module phys_const

module color_msg
    implicit none

    character(len=5), parameter :: c_red = achar(27)//"[31m"
    character(len=5), parameter :: c_green = achar(27)//"[32m"
    character(len=5), parameter :: c_blue = achar(27)//"[34m"
    character(len=5), parameter :: c_purple = achar(27)//"[35m"
    character(len=5), parameter :: c_aqua = achar(27)//"[36m"
    character(len=5), parameter :: c_peach = achar(27)//"[91m"
    character(len=5), parameter :: c_pink = achar(27)//"[95m"
    character(len=4), parameter :: c_reset = achar(27)//"[0m"

    private

    public :: msg_aqua, msg_blue, msg_red
    public :: print_memory_usage

contains

    function msg_aqua(m) result (mc)
        implicit none
        character(len=*), intent(in) :: m
        character(len=200) :: mc
        mc = c_aqua // trim(m) // c_reset
    end function

    function msg_blue(m) result (mc)
        implicit none
        character(len=*), intent(in) :: m
        character(len=200) :: mc
        mc = c_blue // trim(m) // c_reset
    end function

    function msg_red(m) result (mc)
        implicit none
        character(len=*), intent(in) :: m
        character(len=200) :: mc
        mc = c_red // trim(m) // c_reset
    end function

    subroutine print_memory_usage(arr, arrname)
        implicit none
        double precision, dimension(:,:), intent(in) :: arr
        character(len=*), intent(in) :: arrname
        integer :: mbytes
        character(len=50) :: msg

        ! Compute and print memory usage
        mbytes = size(arr) * storage_size(arr(1,1)) / 8388608
        write(msg, '(A,A,A,I5)') "Mb used by array ", arrname, ": ", mbytes
        LOGME trim(msg_aqua(msg))
    end subroutine print_memory_usage

end module color_msg


module time_monitor
    use color_msg
    ! A module that contains a few time monitors.
    ! Each time monitor accumulates elapsed seconds in its counter.
    implicit none
    integer, parameter :: mon_n = 7
    integer, parameter :: t_day = 24*3600  ! seconds in a days
    character(len=50), dimension(mon_n) :: mon_msg
    logical, dimension(mon_n) :: mon_running
    integer, dimension(mon_n) :: mon_t0, mon_dt


    private

    public :: time_monitor_set, time_monitor_print, &
        & time_monitor_start, time_monitor_read, time_monitor_pause

contains

    function get_mseconds() result (t)
        implicit none
        integer :: t
        character(len=50) :: cdate, ctime, czone
        integer, dimension(8) :: ival

        call date_and_time(cdate, ctime, czone, ival)
        !   hours               minutes           seconds        mseconds
        t = ival(5)*3600*1000 + ival(6)*60*1000 + ival(7)*1000 + ival(8)
    end function get_mseconds

    subroutine time_monitor_start(i)
        implicit none
        integer, intent(in) :: i
        mon_t0(i) = get_mseconds()
        mon_running(i) = .true.
    end subroutine time_monitor_start

    function time_monitor_read(i) result (dt)
        implicit none
        integer, intent(in) :: i
        integer :: t, dt
        t = get_mseconds()
        dt = int(ceiling( (t - mon_t0(i) + mon_dt(i)) / 1000.0 ))
    end function time_monitor_read

    subroutine time_monitor_pause(i)
        implicit none
        integer, intent(in) :: i
        integer :: t
        t = get_mseconds()
        if (t .lt. mon_t0(i)) then
            ! Assume that the previous measurement was before midnight and this
            ! measurement is after midnight.  Ignore the case that more than
            ! one day elapses between the two measurements.
            mon_dt(i) = mon_dt(i) + t + t_day - mon_t0(i)
        else
            mon_dt(i) = mon_dt(i) + t - mon_t0(i)
        end if
        mon_t0(i) = -1
        mon_running(i) = .false.
    end subroutine time_monitor_pause

    subroutine time_monitor_set(i, msg)
        implicit none
        integer, intent(in) :: i
        character(len=*) :: msg
        mon_msg(i) = trim(msg)
        mon_dt(i) = 0
        mon_t0(i) = -1
        mon_running(i) = .false.
    end subroutine time_monitor_set

    subroutine time_monitor_print(i)
        implicit none
        integer, intent(in) :: i
        character(len=200) :: msg1, msg2

        write(msg1, '(I8)') int(ceiling(mon_dt(i) / 1000.0))
        write(msg2, '(4A)') "time_monitor_print: ", trim(mon_msg(i)), ", elapsed seconds: ", trim(adjustl(msg1))
        LOGME trim(msg_aqua(msg2))
    end subroutine time_monitor_print
    
end module time_monitor
