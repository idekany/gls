module gls_mod

    implicit none
    integer, parameter :: dp=kind(0.d0)
    real(dp), parameter :: twopi=2.0_dp*acos(-1.0_dp)
    save

contains

    subroutine readfile(name,nmax,wexp,deferr,wsum,t,flx,flxerr,ww)
    use omp_lib
    implicit none

    character, intent(in) :: name*80
    integer, intent(in) :: nmax
    integer, intent(inout) :: wexp
    real(dp), intent(in) :: deferr
    real(dp), intent(out) :: wsum
    real(dp), dimension(:), allocatable, intent(out) :: flx, flxerr, ww, t

    character :: line*1024
    logical :: isinputfile, iforce
    integer :: io=0, ae=0, ioc=0, i, ntok, n
    real(dp) :: tmin

    iforce=.false.
    tmin=1.0d30

    inquire(file=name, exist=isinputfile)
    if (isinputfile.eqv.(.false.)) then
        write(0,*) "ERROR: input file does not exist."
        stop
        endif

    ! >>> Read time-series and calculate some sums.

    open(unit=1, file=name, status='old' )
    i=1
    do
        read(1,'(a)',IOSTAT=io) line
        if (io>0) then
            print *,name
            write(0,*) "I/O ERROR ",io,": unreadable input."
            close(1)
            stop
        endif
        if (io<0) exit
        if (index(line,'#')/=0) cycle                  ! skip comment line
        if (lnblnk(line)==0) cycle                     ! skip empty line
        N=i
        i=i+1 
    enddo

    allocate(flx(n), STAT = ae)
    allocate(flxerr(n), STAT = ae)
    allocate(ww(n), STAT = ae)
    allocate(t(n), STAT = ae)
    
    if(ae/=0) then
      write(0,*) "ALLOCATION ERROR ",ae
    endif
        
    !print *,N," lines read."
    
    rewind(1,IOSTAT=io)
    if(io/=0) then
      write(0,*) "I/O ERROR ",io,": rewind failed."
    endif

    if (index(name,'#') .gt. 0)  read(1,*)     ! skip header
    wsum = 0.0_dp
    i=1
    do
        read(1,'(a)',IOSTAT=io) line                   ! read line from datafile
        if (io>0) then
            write(0,*) "I/O ERROR: unreadable input."
            close(1)
            stop
        endif
        if (io<0) exit                                 ! End of file, exit loop.
        ! print *,line
        if (index(line,'#')/=0) cycle                  ! skip comment line
        if (lnblnk(line)==0) cycle                     ! skip empty line
        ntok=ntokens(line)
    
        if ( ntok == 3 ) then                          ! error values present (3rd column)
            read(line,*,IOSTAT=ioc) t(i),flx(i),flxerr(i)
            !write(6,*) i,t(i),flx(i),flxerr(i)
            if (flxerr(i)==0.0_dp) flxerr(i)=deferr
            if (ioc>0) then
                write(0,*) "I/O ERROR: unreadable input."
                close(1)
                stop
            endif
        else if ( ntok == 2 ) then                      ! no error values, will use variance
            read(line,*,IOSTAT=ioc) t(i),flx(i)
            if (ioc>0) then
                write(0,*) "I/O ERROR: unreadable input."
                close(1)
                stop
            endif
            flxerr(i) = 1.0_dp
            if ( wexp == 2 ) then
                write(0,*) "WARNING: no error values in input file!"
                write(0,*) "  Forcing wexp = 0"
                wexp=0
                iforce=.true.
            endif
        else
            write(0,*) "I/O ERROR: number of input columns must be 2 or 3."
            close(1)
            stop
        endif
    
        N = i
        tmin = MIN(tmin,t(i))
        ww(i) = (1.0_dp/FLXerr(i))**wexp
        wsum  = wsum + ww(i)
        i=i+1
        if (i>nmax) then
            write(0,*) "WARNING: input has reached maximum number of points (",nmax,")!"
            exit
        endif
    
    end do
    close(1)

    t=t-tmin   ! Shift times by tmin.
    
    ! >>> Shift the time values, compute time base and normalize weights:
    ww = ww / wsum                 ! normalize weights

    if (iforce) then    ! wexp = 0 was forced (2-column input line was found)
        wsum=n
        ww=1.0_dp/DFLOAT(n)
        flxerr=1.0_dp
    endif

    contains

      integer function ntokens(line)
      character,intent(in):: line*(*)
      integer i, n, toks

      i = 1;
      n = len_trim(line)
      toks = 0
      ntokens = 0
      do while(i <= n)
         do while(line(i:i) == ' ') 
           i = i + 1
           if (n < i) return
         enddo
         toks = toks + 1
         ntokens = toks
         do
           i = i + 1
           if (n < i) return
           if (line(i:i) == ' ') exit
         enddo
      enddo
      end function ntokens 

    end subroutine readfile

! ========================================================================

    subroutine glsspec(nthread,name,isout,Fbeg,Fend,ofac,t,flx,ww,yy,psin,powsin,amp,ph,cbest,rms)

    implicit none

    character, intent(in) :: name*80
    integer, intent(in) :: isout, nthread
    real(dp), intent(in) :: ofac,fbeg,fend
    real(dp), dimension(:), intent(in) :: t,flx,ww
    real(dp), intent(out) :: psin,powsin,rms,yy

    integer :: i,n,nstep,jj,ix
    real(dp),dimension(size(t)) :: wy,v
    real(dp), allocatable, dimension(:) :: freq, pow, sw, powls
    real(dp) :: tbase,step,flxmean, &
                a,b,c,ph,amp,cbest,pls
    
    !$call omp_set_num_threads(nthread)
    
    ! Number of epochs:
    n=size(t)

    ! Time base:
    tbase=MAXVAL(t)
    
    ! Frequency step:
    step = 1.0_dp/tbase/ofac
    nstep= int(dabs(Fend-Fbeg)/step)+1
    
    !print *,nstep
    allocate(freq(nstep))
    allocate(pow(nstep))
    allocate(powls(nstep))
    if(isout==1) allocate(sw(nstep))

    !----

    flxmean=0.0_dp
    do i=1,n
        flxmean = flxmean + flx(i)*ww(i)     ! weighted mean
    end do

    yy=0.0_dp
    do i=1,n
        wy(i) = flx(i)-flxmean            ! centering flux
        yy    = yy + wy(i)**2*ww(i)       ! sum for chi2 above the mean
        wy(i) = wy(i)*ww(i)               ! attach weights to centered data (saves many multiplications)
    end do

    rms=sqrt(yy*n/(n-1))

    ! >>> Compute GLS spectrum.

    !$omp parallel do private(v,a,b,c)
    
    do jj=1,nstep
        freq(jj)=fbeg+(jj-1)*step
        v = t * twopi * freq(jj)
        call sinefit(v,wy,ww,yy,a,b,c,pow(jj),powls(jj))
        if(isout==1) call specwin(v,sw(jj))
    enddo
    !$omp end parallel do

    powsin=MAXVAL(pow)
    ix=MAXLOC(pow,1)
    v=t*twopi*freq(ix)
    call sinefit(v,wy,ww,yy,a,b,c,powsin,pls)
    ph = MOD(DATAN2(a,b)+twopi, twopi)
    amp = DSQRT(a**2+b**2)
    cbest = flxmean+c
    psin = 1.0_dp/freq(ix)
    
    if(isout==1) then    ! open output file for spectrum
        open(unit=10, file=TRIM(name)//'.gls')
        write(10,'(4A18)')  "#frequency","pow","window","powLS"    
        do jj=1,nstep
            write(10,*) freq(jj), pow(jj), sw(jj), powls(jj)
        enddo
        close(10)
    endif
    
    deallocate(freq)
    deallocate(pow)
    if(allocated(sw)) deallocate(sw)
    deallocate(powls)

    end subroutine glsspec


! ========================================================================

    pure subroutine sinefit(v,z,ww,yy,A,B,off,pow,powls)
    ! fit sine wave y=A*cosx+B*sinx+off
    ! A,B,off - fit parameter
    implicit none
    
    real(dp), dimension(:), intent(in) :: v, z, ww
    real(dp), intent(in) :: yy
    
    real(dp), intent(out) :: a,b,off,pow,powls
    
    integer :: i,n
    real(dp) :: cc,ss,cs,c,s,yc,ys,d,cosx,sinx

    n=size(v)

    CC = 0.
    CS = 0.
    C  = 0.
    S  = 0.
    YC = 0.
    YS = 0.
    do i=1,N                        ! calculate the weighted sums
       cosx = cos(v(i))
       sinx = sin(v(i))
       CC = CC + ww(i) * cosx**2
       CS = CS + ww(i) * cosx*sinx
       C  = C  + ww(i) * cosx
       S  = S  + ww(i) * sinx
       YC = YC + z(i) * cosx
       YS = YS + z(i) * sinx
    enddo
    SS = 1. - CC
    D  = CC*SS - CS*CS
    powLS = (SS*YC**2/D+CC*YS**2/D-2*CS*YC*YS/D) / YY         ! Lomb-Scargle power
    CC = CC - C*C
    SS = SS - S*S
    CS = CS - C*S
    D  = CC*SS - CS*CS
    A = (YC*SS-YS*CS) / D
    B = (YS*CC-YC*CS) / D
    off = -A*C-B*S
    !      pow= (A*YC+B*YS)/ YY
    pow = (SS*YC**2/D+CC*YS**2/D-2*CS*YC*YS/D) / YY           ! GLS power
    end subroutine sinefit
      
! ========================================================================

    subroutine specwin(v,WS)
    implicit none
    
    real(dp), dimension(:), intent(in) :: v
    
    real(dp), intent(out) :: ws
    
    integer :: i,n
    real(dp) :: wc
    
    n=size(v)
    
    wc = 0.0_dp
    ws = 0.0_dp
    do i=1,n
       wc = wc + cos(v(i))
       ws = ws + sin(v(i))
    enddo
    ws = (wc**2.0_dp+ws**2.0_dp) / n**2.0_dp
    return
    end subroutine specwin

end module gls_mod
