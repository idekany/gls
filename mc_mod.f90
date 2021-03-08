module mc_mod

    implicit none
    integer, parameter :: dp=kind(0.d0)
    save

contains

    subroutine gauss(garr,sig,n)
    ! This funtion generates an array of Gaussian random variables
    ! with zero norm and sig**2 variance.

    implicit none

    integer, intent(in) :: n
    real(dp),intent(in) :: sig

    real(dp),dimension(n),intent(out) :: garr

    integer :: i
    real(dp) :: rnd1, rnd2, p1, p2, rr, fac

    do i=1,n
        rr=2.0d0
        do while (rr>1.0d0)
            call random_number(rnd1)
            call random_number(rnd2)
            p1=2.0d0*rnd1-1.0d0
            P2=2.0d0*rnd2-1.0d0
            rr=p1*p1+p2*p2
            enddo

        fac=DSQRT(-2.0d0*DLOG(rr)/rr)
        garr(i)=p1*fac
          enddo

    garr=sig*garr

    return
    end subroutine gauss

    ! ========================================================================

    subroutine ssort (x, iy)
    implicit none
    !
    !    This code performs a "bubble sort".
    !
    !    Sort an array and make the same interchanges in
    !    an auxiliary array.  The array is sorted in
    !    decreasing order.
    !
    !   Description of parameters
    !      X - array of values to be sorted   (usually abscissas)
    !      IY - array to be carried with X (all swaps of X elements are
    !          matched in IY .  After the sort IY(J) contains the original
    !          position of the value X(J) in the unsorted X array.
    !      N - number of values in array X to be sorted
    !
    real(dp), dimension(:), intent(inout) :: x
    integer, dimension(:), intent(inout):: iy

    integer :: n,i,j,jmax,itemp
    real(dp) :: temp

    n=size(x)

    if(n/=size(iy)) then
        write(0,*) "ERROR (ssort): input arrays of different sizes."
    endif

    jmax=n-1
    do i=1,n-1
        temp=1.e38
        do j=1,jmax
            if ( x(j) > x(j+1) ) cycle
            temp=x(j)
            x(j)=x(j+1)
            x(j+1)=temp
            itemp=iy(j)
            iy(j)=iy(j+1)
            iy(j+1)=itemp
            enddo
        if(temp==1.e38) exit
        jmax=jmax-1
    enddo
    return
    end subroutine ssort

    ! ========================================================================

    subroutine interx1(x,y,xd,yd)

    !    ----
    !    This routine interpolates the sequence:
    !    X(I),Y(I), I=1,2,...,N1
    !
    !    in X(I), into the value xd
    !
    !    output is given as:
    !    yd
    !    ----

    implicit none

    real(dp),dimension(:),intent(in) :: x,y
    real(dp),intent(in) :: xd

    real(dp),intent(out) :: yd

    integer :: iq,j,jj,j1,j2,j3,n1
    real(dp) :: aa,bb,cc,d,dev

    n1=size(x)
    if(n1/=size(y)) then
        write(0,*) "ERROR (interx1): input arrays of different sizes."
    endif

    iq=2
    !    IQ = 2 --> QUADRATIC INTERPOLATION
    !       = 1 --> LINEAR INTERPOLATION

    !    Find the proper array indices for the interpolation regime:

    dev=1.0d30

    !    in this cycle we find the x value
    !    which is the closest one to the actual xd
    do j=1,n1
        d=dabs(xd-x(j))
        if(d.gt.dev) cycle
        dev=d
        jj = j
        enddo

    if (jj.eq.1) then
          j1=1
          j2=2
          j3=3
    else if (jj.eq.n1) then
          j1=n1-2
          j2=n1-1
          j3=n1
    else
          j1=jj-1
          j2=jj
          j3=jj+1
    endif

    if(IQ==1) then
        ! LINEAR INTERPOLATION
        YD=Y(J1)+(xd-X(J1))*(Y(J2)-Y(J1))/(X(J2)-X(J1))
    else if (IQ==2) then
        ! QUADRATIC INTERPOLATION
        AA=(xd-X(J2))*(xd-X(J3))/(X(J1)-X(J2))/(X(J1)-X(J3))
        BB=(xd-X(J1))*(xd-X(J3))/(X(J2)-X(J1))/(X(J2)-X(J3))
        CC=(xd-X(J1))*(xd-X(J2))/(X(J3)-X(J1))/(X(J3)-X(J2))
        YD=Y(J1)*AA+Y(J2)*BB+Y(J3)*CC
        endif

    return
    end subroutine interx1

    ! ========================================================================

    subroutine shuffle(a)

    implicit none

    real(dp), dimension(:), intent(inout) :: a

    integer :: i, randpos, n
    real(dp) :: temp,r

    n=size(a)

    do i = n, 2, -1
        call random_number(r)
        randpos = int(r * i) + 1
        temp = a(randpos)
        a(randpos) = a(i)
        a(i) = temp
        end do

    end subroutine shuffle

    ! ========================================================================

    subroutine init_random_seed()
      use iso_fortran_env, only: int32
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid
      integer(int32) :: t

      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
         ! Fall back to the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(t)
         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int32 * 24 * 60 * 60 * 1000 &
                 + dt(2) * 31_int32 * 24 * 60 * 60 * 1000 &
                 + dt(3) * 24_int32 * 60 * 60 * 1000 &
                 + dt(5) * 60 * 60 * 1000 &
                 + dt(6) * 60 * 1000 + dt(7) * 1000 &
                 + dt(8)
         end if
         pid = getpid()
         t = ieor(t, int(pid, kind(t)))
         do i = 1, n
            seed(i) = lcg(t)
         end do
      end if
      call random_seed(put=seed)
    contains
      ! This simple PRNG might not be good enough for real work, but is
      ! sufficient for seeding a better PRNG.
      function lcg(s)
        integer :: lcg
        integer(int32) :: s
        if (s == 0) then
           s = 104729
        else
           s = mod(s, 4294967296_int32)
        end if
        s = mod(s * 279470273_int32, 4294967291_int32)
        lcg = int(mod(s, int(huge(0), int32)), kind(0))
      end function lcg
    end subroutine init_random_seed

end module mc_mod
