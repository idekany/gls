! -----------------------------------------------------------------------
! This code computes the generalized Lomb-Scargle (GLS) periodogram
! for one time-series read from an ASCII input file.
! The following paper describes the method:
!        Zechmeister & Kuerster, 2009, A&A, 496, 577)
! -----------------------------------------------------------------------

PROGRAM gls

use gls_mod, only: dp, twopi, readfile, glsspec, sinefit, specwin
use mc_mod, only: gauss, ssort, interx1, shuffle, init_random_seed

implicit none

integer, parameter :: nmax=1000000, nsimmax=10000
integer :: i,N,wexp=2,io=0,index,iargc,en
! wexp = weightening exponent (default Chi2)
! weightening of errors: wexp=0 (variance), wexp=2 (Chi2)
integer :: numarg,nstep,istep,iter,nthread=1
integer :: iverb=1,isout=1,ilog=0,isim=0,nsim=100,iboot=0
logical :: isparfile, longhelp=.false.
character :: TV*40, name*80, fittype*8, logfile*80
integer, external :: ntokens

real(dp) :: tbase=0,FLXmean=0,wsum,YY=0, &
rmssin=0,chisin=0,wrmssin=0,powSin=0,dFLXsin, &
fbeg,fend,fendmax,step,ofac,PSin,Amp,CBest,ph, &
w,M,FAP,prob,probMC=1,deferr,rms0=1,rms,PSin0,PowSin0=0
        
real(dp), dimension(:), allocatable :: flx,flxerr,t,ww
real(dp), allocatable, dimension(:) :: powMC,pp
integer, allocatable, dimension(:) :: indx

fittype="unknown"
logfile='gls.log'
fendMax=50.0d0
ofac=20.0d0
deferr=0.001d0

!========================= USER INTERFACE ========================
   
! >>> Read parameter file.

inquire(file="gls.par", exist=isparfile)
if (isparfile) then    ! (else defaults will be used)
   open(unit=2, file= "gls.par", status= "OLD")
   read(2,*) ofac                ! default oversampling factor = 20
   read(2,*) wexp                ! default wexp = 0 (variance)
   read(2,*) Fbeg                ! default Fbeg = 1./tbase/10.
   read(2,*) Fend                ! default Fend = N*Fbeg*10. (average Nyquist frequency)
   read(2,*) deferr              ! default error (if error value is 0 in input file)
   read(2,*) isim                ! estimate significance by MC?
   read(2,*) nsim                ! number of realizations in MC
   read(2,*) iboot               ! use bootstrapping instead of white noise?
   read(2,*) iverb               ! verbose mode?
   read(2,*) isout               ! write out spectrum?
   read(2,*) ilog                ! write logfile?
   read(2,*) nthread             ! number of threads
   close(2)
endif

ofac=DABS(ofac)
nsim=ABS(nsim)
nthread=ABS(nthread)
if (nsim>nsimmax) then
  write(0,*) "WARNING: nsim > ",nsimmax," !"
  write(0,*) "  nsim = ",nsimmax," has been set."
  nsim=nsimmax
end if

if (fbeg>fend) then
    write(0,*) "WARNING: minimum frequency is higher than maximum frequency!"
    write(0,*) "  Minimum and maximum frequencies have been swapped."
    fbeg = fbeg + fend
    fend = fbeg - fend
    fbeg = fbeg - fend
endif

if (fend>fendmax) then
    fend=fendmax
    write(0,*) "WARNING: maximum frequency is larger than ",fendmax," !"
    write(0,*) "  Maximum frequency have been set to ",fendmax," ."
endif

if (.NOT.(wexp==0 .or. wexp==2)) then
    write(0,*) "WARNING: weightening exponent has an illegal value!"
    write(0,*) "  Weightening exponent has been set to 0."
endif

! >>> Read command line arguments.

numarg=iargc()

if (numarg==0) then
    call print_help(longhelp)
endif

call getarg(1,name)

if(name=='-help') then
    longhelp=.true.
    call print_help(longhelp)
endif

if (numarg>1) then    ! if the are any command line arguments
    do i=1,numarg-1
        call getarg(i,name)             ! testing i-th command line argument
        if (name=='-e2') wexp=2          ! i.e. chi2 (use error bars)
        if (name=='-e1') wexp=1          ! i.e. chi2 (use error bars)
        if (name=='-ne') wexp=0          ! i.e. variance (no error bars)
        if (name=='-s')  isout=1         ! write out spectrum
        if (name=='-ns') isout=0         ! don't write out spectrum
        if (name=='-v')  iverb=1
        enddo
endif

if (wexp==2) fittype = "chi2"
if (wexp==0) fittype = "variance"

! >>> Read input file name (last command line argument):
call getarg(numarg,name)


!========================= DATA INPUT ========================

call readfile(name,nmax,wexp,deferr,wsum,t,flx,flxerr,ww)

n=size(t)

!========================= GLS SPECTRUM ========================

! >>> Check if parameter file exists, if no, set Fbeg, Fend to defaults.

tbase=MAXVAL(t)

if(isparfile.eqv.(.false.)) then    ! no parameter file -> no frequency range specified
    Fbeg = 1./tbase/10.              ! default value
    Fend = N*Fbeg*10.                ! default value (average Nyquist frequency)
    !if (Fbeg.lt.1./tbase/10.) Fbeg = 1./tbase/10.              ! default value      
endif
! Set Fbeg to small positive value:
step = 1.0d0/tbase/ofac
if(Fbeg.lt.step) Fbeg=1./tbase/10.

! >>> Compute GLS Spectrum:    
call glsspec(nthread,name,isout,Fbeg,Fend,ofac,t,flx,ww,yy,psin,powsin,amp,ph,cbest,rms)

powSin0=powSin
PSin0=PSin
rms0=rms
do i=1,N
    dFLXsin = FLX(i)-(Amp*sin(t(i)*twopi/PSin0+ph)+CBest)
    rmssin = rmssin + dFLXsin**2
    chisin = chisin + (dFLXsin/FLXerr(i))**2
end do
rmssin = sqrt(rmssin/(N-4))           ! unbiased rms
wrmssin = sqrt(chisin/wsum*N/(N-4))   ! weighted rms

! >>> Calculate analytical estimator of false alarm probability.

prob= (1.0d0-powSin0)**((dfloat(N)-3.0d0)/2.0d0)        ! spectral false alarm probability
M   = tbase * dabs(Fend-Fbeg)         ! number of independent frequencies
FAP = M * prob
if (FAP>0.01d0) FAP = 1.0d0 - (1.0d0-prob)**M


!========================= MONTE CARLO FAP ESTIMATION ========================

if(isim==1) then    ! Perform MC estimation of FAP.

    allocate(powmc(nsim))
    allocate(pp(nsim))
    allocate(indx(nsim))
    
    ww=1.0_dp/DFLOAT(n)  ! No error bars in this case..
    do iter=1,nsim   ! Loop for MC (will become iterated if isim=1)
        
        if(iboot.eq.0) then
            call gauss(flx,rms0,n)
            !open(7,file="tmp")
            !do i=1,n
            !  write(7,*) t(i),flx(i)
            !enddo
            !close(7)
        else
            call shuffle(flx)
        end if
        
        call glsspec(nthread,name,isout,Fbeg,Fend,ofac,t,flx,ww,yy,psin,powsin,amp,ph,cbest,rms)
        
        powMC(iter)=powsin
        pp(iter)=1.0d0-DFLOAT(nsim-iter)/DFLOAT(nsim)
        
    end do   !  End loop for MC.
    call ssort(powMC,indx)
    if(powSin0.le.powMC(1)) then
        call interx1(powMC,pp,powsin0,probmc)
    else
        probmc=pp(1)
    end if
    
    open(7,file="tmp1")
    do i=1,nsim
        write(7,*) pp(i),powmc(i)
        end do
    close(7)
    
    deallocate(powmc)
    deallocate(pp)
    deallocate(indx)

end if

if (iverb.eq.1) then
    print *, "Number of data points = ", N
    print *, "Total time base =", tbase," days"
    print *, "Number of independent frequencies = ", 1./tbase, " 1/d"
    print *, fittype(:lnblnk(fittype))//"_0 = ", YY*wsum
    print *, "rms_0 = ",rms
    write(6,*)"Number of frequency steps = ",int((Fend-Fbeg)/step+1)
    endif


if (iverb.eq.1) then
    print *
    write (*,'(1x,a,i3,a,i2)') 'Used parameters: ofac =',int(ofac),' , wexp =',int(wexp)
    print *, " ",fittype(:lnblnk(fittype))//" for constant:",YY*wsum
    print *, " RMS for constant: ",rms0
    if (wexp .eq. 2) print *, " mean weighted internal error:", sqrt(N/wsum)
    print *
    print *, "Best sine fit values for ", name(:lnblnk(name))
    print *, " best sine period [d] =       ", PSin0," +/-",sqrt(6./N)*2/twopi/tbase*wrmssin/Amp*PSin0**2
    print *, " amplitude =                  ", Amp," +/-",sqrt(2./N)*wrmssin
    write(6,*)
    print *, fittype(:lnblnk(fittype))," of residuals:",(1.-powSin0)*YY*wsum,chisin ! variance after fit
    print *, " weighted rms of residuals:  ", wrmssin
    print *, " rms of residuals:           ", rmssin
    print *, " highest power p[0..1]:      ", powSin0
    print *, "               P[0..(N-1)/2]:", powSin0*(N-1)/2
    print *, " Analytical FAP(spectral)=(1-p)^((N-3)/2):", prob
    print *, " Analytical FAP:                        ", FAP
    if(isim.eq.1) then
        print *, " FAP (MC):", probMC
        print *, " Number of realizations:", nsim
        end if
    print *, "----------------------------------------------------------------"
    print *, "RESULTS:"
    print *
    endif

! >>> Write results to STDOUT.
if(isim.eq.0) then
  write(*,'(a,1x,i5,1x,f9.6,1x,f7.3,1x,f6.2,1x,i7,2x,e15.3e4,1x,f8.6)') &
  TRIM(name),N,1.0d0/PSin0,powSin0,Fend,int(M),FAP,1.0d0/tbase
else
  write(*,'(a,1x,i5,1x,f9.6,1x,f7.3,1x,f6.2,1x,i7,2x,2(e15.3e4,1x),f8.6)') &
  TRIM(name),N,1.0d0/PSin0,powSin0,Fend,int(M),FAP,probmc,1.0d0/tbase
endif
  
if (allocated(flx)) deallocate(flx)
if (allocated(flxerr)) deallocate(flxerr)
if (allocated(ww)) deallocate(ww)
if (allocated(t)) deallocate(t)


!  TRIM(name),N,1.0d0/PSin0,powSin0*(N-1)/2,Fend,int(M),FAP

! >>>  Add results to logfile.
if (ilog==1) then
    open(unit=3, file=TRIM(logfile) ,access='APPEND')
    write(3,*)"file='",name(:lnblnk(name)),"'"
    write(3,*)"PSin=",PSin0,";  Amp=",Amp,";"
    close(3)
endif

contains

subroutine print_help(longhelp)
    implicit none
    logical, intent(in) :: longhelp
    
    if (longhelp) then
        print *,"Usage:      gls [<flag>] <inputfile>"
        print *,"Where:"
        print *," <inputfile> contains the input time series in columns:"
        print *,"             time, flux [, flux error]"
        print *," <flag> may be one or more of the following options:"
        print *,"        -e2 : use chi^2 weightening (= use error bars)"
        print *,"        -ne : use variance weightening (= don't use error bars)"
        print *,"        -s  : write out spectrum (output file is <inputfile>.gls)"
        print *,"        -ns : don't write out spectrum"
        print *,"        -v  : verbose mode (writes out information on STDOUT)"
        print *,"Optionally, a parameter file called 'gls.par' can be provided."
        print *,"The parameter file overrides default settings, flags override both."
        print *,"Normal output is one line to STDOUT with the following values:"
        print *,"  <filename>, number of points in time series, main frequency [d^-1],"
        print *,"  GLS power of main frequency, high frequency limit, number of independent frequencies,"
        print *,"  analytical false alarm probability (FAP), [Monte Carlo FAP,] nominal frequency resolution"
        print *,""
        stop
    else
        print *,"Usage:      gls [<flag>] <inputfile>"
	       print *,"Type 'gls -help' for more info."
	       stop
    endif
    
end subroutine print_help
    
end program gls 

