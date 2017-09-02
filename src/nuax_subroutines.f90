module nuax_subroutines

contains

!************************************************************************

  subroutine gnutrino(u2,tt2)

!  Neutrino rates from Itoh, Adachi, Nakagawa, Kohyama, & Munkata, 1989, 
!  ApJ, 339, 354. (+Errata 1990, ApJ, 360, 741)
!    (Photo, Pair)
! Plasmon neutrino rates from Itoh et al. 1996, ApJS, 102, 411
! Recombination neutrinos from Kohyama, Itoh, Obama, & Mutoh 1993, ApJ,
!  415, 267
!  Bremsstrahlung from:
! 1. Liquid Metal: Itoh & Kohyama, 1983, ApJ, 275, 858.
! 2. Crystal Lattice: Itoh etal., 1984, ApJ, 279, 413.
! 3. Low T Quantum Corr.: Itoh etal., 1984, ApJ, 280, 787.
! 4. Phonon Contr.: Itoh etal., 1984, ApJ, 285, 304.
! 5. Partially Degen. e: Munkata, Kohyama & Itoh, 1987, ApJ, 316, 708.
! Prepared for WDXD by PAB: Oct. 1990
! Last messed with on 11 09 05 (new plasmon neutrino rates ABK)

    use misc_const
    use neut
    use temp
    use jee
    use comp
    use cc
    use rhomooee
    use phys
    use flags

    implicit double precision (a-h,o-z)

    real*8 :: mue, qplas
    real*8, pointer :: tx2, e2
    real*8, parameter :: gamc = 3577200., gamo = 5778200.
    real*8, parameter :: gamhe = 573264., gamh = 227500.
    real*8 :: cv,cvp,ca,cap,qphot
    real*8 :: aavg,zbar,t,qrec

!** initialize ***
    tx2 => t2
    e2 => ea2
    t=10.d0**tt2
    rho=10.d0**u2
    do i=1,5
       en(i) = 0.0
    end do
    fn = 0.0
    qbrem = 0.0
    qplas = 0.0
    qphot = 0.0
    qpair = 0.0
    qrec  = 0.0

! Escape clause (enter when T < 10^7K or Rho < 1 g/cc)
! With this, I don't need an escape clause for flam

    if(tt2.lt.7.0 .or. u2.lt.0.0)then
       goto 86
    endif
!** compute gamma and mue for mixture ***
!    if(je.eq.1)then
    if (x(1) .eq. 1.d0) then
!** assume pure H for now. (Note: for H, I shouldn't be here anyway!) ***
       ce = 1.
       mue = 1.
       gamma = gamh * 10.**(u2/3. - tt2)
       zbar = 1.0
       aavg = 1.0
    else
!       if(je.eq.2)then
       if (x(2) .eq. 1.d0) then
!** pure helium ***
          omue = 1./(x(2)/2.)
          gamma = gamhe * 10.**(u2/3. - tt2)
          zbar = 2.0
          aavg = 4.0
       else
          xc2 = x(3)
          xo2 = x(4)
!** helium, carbon, oxygen, and iron (the latter for tests) ***
          gammahe = gamhe * 10.**(u2/3. - tt2)
          gammac  = gamc  * 10.**(u2/3. - tt2)
          gammao  = gamo  * 10.**(u2/3. - tt2)
          gamma = xc2*gammac + xo2*gammao
          omue = xc2*6./12. + xo2*8./16.
          zbar = xc2*6.0 + xo2*8.0
          aavg = xc2*12.0 + xo2*16.0
       endif
    endif

! zbar is the mean charge on a nucleus, for bremm and pair.  
! flam was x in neu 

    flam=t/5.9302d+9
    mue = 1.d0 / omue
    rhomooe = rho/mue
    rhomue23 = (rhomooe)**(2./3.)
    sqfac = 1. + 1.018*rhomue23
    tfermi = 5.9302d+9*(dsqrt(sqfac) - 1.)
    degtemp = t/tfermi
    y = (rhomooe)**(1./3.)/(1000.*flam)
    s2theta = 0.217
    cv=0.5d0+2.0d0*s2theta
    cvp=1.d0-cv
    ca=0.5d0
    cap=1.d0-ca
    if(y .le. 100.)then
       call plasmon(tt2,rhomooe,qplas)
    endif
    if(y .le. 100.)then
       call photon(flam,tt2,rho,mue,cv,cvp,ca,cap,qphot)
    endif
    if(flam.ge.0.10)then
       call pairn(flam,rho,mue,cv,cvp,ca,cap,qpair)
    endif
    call recombine(aavg,zbar,t,rho,cv,cvp,ca,cap,qrec)   
    if(u2.gt.4.0 .and. degtemp.gt.0.3)then
       call brempde(u2,tt2,gamma,qbrempde)
       qbrem = qbrempde
       if(qbrem.lt.0.0)then
          qbrem = 0.0
       endif
!** skip over other neutrino rates if we came here. ***
       go to 50
    endif
    if(u2.ge.1.0)then
       if(gamma.gt.1.0 .and. gamma.lt.171.0)then
          if(x(2).gt.0.00001)then
             call bremliqh(u2,tt2,gamma,qbremlh)
          else 
             qbremlh = 0.0
          endif
          if(x(3).gt.0.00001)then
             call bremliqc(u2,tt2,gamma,qbremlc)
          else 
             qbremlc = 0.0
          endif
          if(x(4).gt.0.00001)then
             call bremliqo(u2,tt2,gamma,qbremlo)
          else 
             qbremlo = 0.0
          endif
          qbrem = qbremlh*x(2) + qbremlc*x(3) + qbremlo*x(4)
       elseif(gamma.ge.171.0)then
          if(x(3).gt.0.00001)then
             call bremsolc(u2,tt2,gamma,qbremsc)
          else 
             qbremsc = 0.0
          endif
          if(x(4).gt.0.00001)then
             call bremsolo(u2,tt2,gamma,qbremso)
          else 
             qbremso = 0.0
          endif
          qbrem = qbremsc*x(3) + qbremso*x(4)
          if(qbrem.lt.0.0)then
             if (verbose) write(*,*) 'Ooops!',qbrem
             qbrem = 0.0
          endif
       endif
    endif

!** convert to lumonisities ***
50  eplas=qplas/rho
    en(1) = eplas
    ephot=qphot/rho
    en(2) = ephot
    epair=qpair/rho
    en(3) = epair
    erec=qrec/rho
    en(4) = erec

!** output of Bremsstrahlung routines already is divided by rho ***
    ebrem=qbrem
    en(5) = ebrem
!    etot=eplas+ephot+epair+erec+ebrem
    etot=ebrem
    fn=etot
    nucounter = nucounter + 1   
86  return
  end subroutine gnutrino

!************************************************************************
! Plasmon Neutrinos according to Itho et al. 1996, APJS, 102, 411
! According to Itoh et al. 1996, ApJ, 470, 1015, those rates are 
! numerically more accurate than the ones from 1992 (which
! are the same as the ones from 1989 except for ultra-relativistic 
! electrons).
! Replaces Paul's plasmon routine.
! ABK 2005
!**********************************************************************

  subroutine plasmon(t2, rhomooe, qplas)

    implicit none
    real*8 :: t2, logrhomooe, qplas, rhomooe

    logrhomooe = dlog10(dble(rhomooe))
    call formulae(t2, logrhomooe, qplas)
    return
  end subroutine plasmon

!**********************************************************************
! Plasma process fitting formulae as written by Itoh et al. 1996
! Given in AAS CD-ROM series, volume 7, 1996 December.
!**********************************************************************
 
  subroutine formulae ( t, r, qn)      
       
    implicit real*8 (a-z)
    
    real*8 :: r
       
    tn      = 1.d1**t
    rn      = 1.d1**r
    rdouble = dlog10(2.0d0*rn)
     
    x = (  1.75d1 + rdouble - 3.d0*t )/6.d0
    y = ( -2.45d1 + rdouble + 3.d0*t )/6.d0
    
    if ( dabs(dble(x)).gt.7.d-1.or.y.lt.0.d0 ) then
       fxy = 1.d0
    elseif(y-1.6d0+1.25d0*x.lt.0.d0) then
       fxy = 1.05d0 +(0.39d0 -1.25d0*x - 0.35d0*dsin(4.5d0*x) &
            -0.3d0*dexp(-1.d0*(4.5d0*x + 0.9d0)**2.d0)) &
            *dexp(-1.d0*((y-1.6d0+1.25d0*x)/(0.57d0-0.25d0*x))**2.d0)       
    else  
       fxy = 1.05d0 +(0.39d0 - 1.25d0*x - 0.35d0*dsin(4.5d0*x) &
            -0.3d0*dexp(-1.d0*(4.5d0*x + 0.9d0)**2.d0)) 
    endif
    
    lamb = (1.686d-10)*tn
    
    gamd = 1.1095d11*rn/((tn**2.d0)*((1.d0+(1.019d-6*rn)**(2.d0/3.d0))**0.5d0)) 
    ft=2.4d0+0.6d0*(gamd**0.25d0)+0.51d0*(gamd**0.5d0) + 1.25d0*(gamd**0.75d0) 
    fl=(8.6d0*gamd+1.35d0*(gamd**(7.d0/4.d0)))/(2.25d2-1.7d1*(gamd**0.5d0)+gamd)
 
    if(dabs(gamd**0.5d0).gt.700) then 
       qap=0.d0
    else 
       qap = fxy*(fl+ft)*dexp(-1.d0*gamd**0.5d0)*(gamd**3.d0)*(lamb**9.d0)*3.d21
    endif
      
    cvd  = ( 0.5d0 + 2.d0*0.2319d0)**2.d0 + n*(0.5d0 - 2.d0*0.2319d0)**2.d0
 
    if ( qap*cvd.le.0.d0 ) then 
       q = 9.9999999d-99
    else
       q    = dlog10(dble(qap)*cvd)
       qn   = qap*cvd 
    endif
    return 

  end subroutine formulae    

!****************************************************************************

  subroutine photon(flam,t2,rho,mue,cv,cvp,ca,cap,qphot)

    use misc_const, only : pi

! Photoneutrinos
! a(i) terms are different in the newest paper. Pain to crank out too.
! Coefficients valid for 10^7 to 10^8K, then 10^8 to 10^9K, and
! finally above 10^9K. 

    implicit double precision (a-h,o-z)

    real*8, dimension(3) :: a,b !a(0:2)
    real*8, dimension(3,5) :: d
    real*8, dimension(3,7) :: c
    real*8 :: cee,xi,flam,f
    real*8 :: mue, rho, cv, cvp, ca, cap, qphot, t2

    xi=(((rho/mue)/1.d9)**(1./3.))/flam

! b(i) terms same as before. Load 'em up.

    b(1)=6.290d-3
    b(2)=7.483d-3
    b(3)=3.061d-4

! Calculate a(i) numerical factors.
! and tabulated coefficients for sines and cosines.
! tabulated stuff is in form c(i,j), starting with 
! i=1, j=1 instead of 0,0

    half = 1.d0 / 2.d0
    facpi = 5.d0 * pi / 3.d0 

!*** 10^7 to 10^8 K ***
    if(t2.ge.7.0 .and. t2.lt.8.0)then
       tau = t2 - 7.0
       cee = 0.5654 + tau
       c(1,1) = 1.008d+11 * half
       c(1,2) = 0.0
       c(1,3) = 0.0
       c(1,4) = 0.0
       c(1,5) = 0.0
       c(1,6) = 0.0
       c(1,7) = 0.0
       c(2,1) =  8.156d+10 * half
       c(2,2) =  9.728d+8
       c(2,3) = -3.806d+9
       c(2,4) = -4.384d+9
       c(2,5) = -5.774d+9
       c(2,6) = -5.249d+9
       c(2,7) = -5.153d+9 * half
       c(3,1) =  1.067d+11 * half
       c(3,2) = -9.782d+9
       c(3,3) = -7.193d+9
       c(3,4) = -6.936d+9
       c(3,5) = -6.893d+9
       c(3,6) = -7.041d+9
       c(3,7) = -7.193d+9 * half
       d(1,1) = 0.0
       d(1,2) = 0.0
       d(1,3) = 0.0
       d(1,4) = 0.0
       d(1,5) = 0.0
       d(2,1) = -1.879d+10
       d(2,2) = -9.667d+9
       d(2,3) = -5.602d+9
       d(2,4) = -3.370d+9
       d(2,5) = -1.825d+9
       d(3,1) = -2.919d+10
       d(3,2) = -1.185d+10
       d(3,3) = -7.270d+9
       d(3,4) = -4.222d+9
       d(3,5) = -1.560d+9

!*** 10^8 to 10^9 K ***
    elseif(t2.ge.8.0 .and. t2.lt.9.0)then
       tau = t2 - 8.0
       cee = 1.5654
       c(1,1) =  9.889d+10 * half
       c(1,2) = -4.524d+8
       c(1,3) = -6.088d+6
       c(1,4) =  4.269d+7
       c(1,5) =  5.172d+7
       c(1,6) =  4.910d+7
       c(1,7) =  4.388d+7 * half
       c(2,1) =  1.813d+11 * half
       c(2,2) = -7.556d+9
       c(2,3) = -3.304d+9
       c(2,4) = -1.031d+9
       c(2,5) = -1.764d+9
       c(2,6) = -1.851d+9
       c(2,7) = -1.928d+9 * half
       c(3,1) =  9.750d+10 * half
       c(3,2) =  3.484d+9
       c(3,3) =  5.199d+9
       c(3,4) = -1.695d+9
       c(3,5) = -2.865d+9
       c(3,6) = -3.395d+9
       c(3,7) = -3.418d+9 * half
       d(1,1) = -1.135d+8
       d(1,2) =  1.256d+8
       d(1,3) =  5.149d+7
       d(1,4) =  3.436d+7
       d(1,5) =  1.005d+7
       d(2,1) =  1.652d+9
       d(2,2) = -3.119d+9
       d(2,3) = -1.839d+9
       d(2,4) = -1.458d+9
       d(2,5) = -8.956d+8
       d(3,1) = -1.549d+10
       d(3,2) = -9.338d+9
       d(3,3) = -5.899d+9
       d(3,4) = -3.035d+9
       d(3,5) = -1.598d+9

!*** Temperature is greater than 10^9 K ***
    elseif(t2.ge.9.0)then
       tau = t2 - 9.0
       cee = 1.5654
       c(1,1) =  9.581d+10 * half
       c(1,2) =  4.107d+8
       c(1,3) =  2.305d+8
       c(1,4) =  2.236d+8
       c(1,5) =  1.580d+8
       c(1,6) =  2.165d+8
       c(1,7) =  1.721d+8 * half
       c(2,1) =  1.459d+12 * half
       c(2,2) =  1.314d+11
       c(2,3) = -1.169d+11
       c(2,4) = -1.765d+11
       c(2,5) = -1.867d+11
       c(2,6) = -1.983d+11
       c(2,7) = -1.896d+11 * half
       c(3,1) =  2.424d+11 * half
       c(3,2) = -3.669d+9
       c(3,3) = -8.691d+9
       c(3,4) = -7.967d+9
       c(3,5) = -7.932d+9
       c(3,6) = -7.987d+9
       c(3,7) = -8.333d+9 * half
       d(1,1) =  4.724d+8
       d(1,2) =  2.976d+8
       d(1,3) =  2.242d+8
       d(1,4) =  7.937d+7
       d(1,5) =  4.859d+7
       d(2,1) = -7.094d+11
       d(2,2) = -3.697d+11
       d(2,3) = -2.189d+11
       d(2,4) = -1.273d+11
       d(2,5) = -5.705d+10
       d(3,1) = -2.254d+10
       d(3,2) = -1.551d+10
       d(3,3) = -7.793d+9
       d(3,4) = -4.489d+9
       d(3,5) = -2.185d+9
    endif

! set up the variables (sines and cosines) for the loops

    su1 = dsin(facpi*1.d0*tau)
    su2 = dsin(facpi*2.d0*tau)
    su3 = dsin(facpi*3.d0*tau)
    su4 = dsin(facpi*4.d0*tau)
    su5 = dsin(facpi*5.d0*tau)
    cu1 = dcos(facpi*1.d0*tau)
    cu2 = dcos(facpi*2.d0*tau)
    cu3 = dcos(facpi*3.d0*tau)
    cu4 = dcos(facpi*4.d0*tau)
    cu5 = dcos(facpi*5.d0*tau)
    cu  = dcos(facpi*10.d0*tau)

!  Sum up the c's and d's to make a's

    csum1=c(1,1)+c(1,2)*cu1 + c(1,3)*cu2 + c(1,4)*cu3 + c(1,5)*cu4 + c(1,6)*cu5
    dsum1 = d(1,1)*su1 + d(1,2)*su2 + d(1,3)*su3 + d(1,4) * su4 + d(1,5) * su5 
    a(1)  = csum1 + dsum1 + c(1,7) * cu
    csum2=c(2,1)+c(2,2)*cu1 + c(2,3)*cu2 + c(2,4)*cu3 + c(2,5)*cu4 + c(2,6)*cu5
    dsum2 = d(2,1)*su1+d(2,2) * su2 + d(2,3) *su3 + d(2,4) * su4 + d(2,5) * su5 
    a(2)  = csum2 + dsum2 + c(2,7) * cu
    csum3=c(3,1)+c(3,2)*cu1 + c(3,3)*cu2 + c(3,4)*cu3 + c(3,5)*cu4 + c(3,6)*cu5
    dsum3 = d(3,1)*su1+d(3,2) * su2 + d(3,3) *su3 + d(3,4) * su4 + d(3,5) * su5 
    a(3)  = csum3 + dsum3 + c(3,7) * cu

! There, now we have the a(i) coefficients! Lot's o'work.

    call qcalc(a,b,cee,xi,flam,f)
    n=2
    temp=1.875d8*flam+1.653d8*flam**2+8.499d8*flam**3-1.604d8*flam**4
    smq=0.666d0*((1.+2.045d0*flam)**(-2.066d0))/(1.d0+rho/mue/temp)
    temp=(cv**2-ca**2)+n*(cvp**2-cap**2)
    temp=temp/((cv**2+ca**2)+n*(cvp**2+cap**2))
    bigq=0.5*((cv**2+ca**2)+n*(cvp**2+cap**2))*(1.d0-temp*smq)
    qphot=bigq*f*(rho/mue)*(flam**5)

    return
  end subroutine photon

!************************************************************************

  subroutine pairn(flam,rho,mue,cv,cvp,ca,cap,qpair)

! Pair production neutrinos
! Munkata etal version is unchanged except for log T > 10, which
! we better not hit, although I've included the appropriate coefficients.

    implicit double precision (a-h,o-z)
    
    real*8, dimension(3) :: a, b 
    real*8 ::  mue,c,xi,flam,f
  
    xi=(((rho/mue)/1e9)**(1./3.))/flam
    a(1)=6.002e19
    a(2)=2.084e20
    a(3)=1.872e21
    n = 2

! pick out appropriate b's and c's according to temp
! flam = 1.6863 is T = 10^10 K

    if(flam.lt.1.6863)then
       b(1)=9.383e-1
       b(2)=-4.141e-1
       b(3)=5.829e-2
       c=5.5924
    else
       b(1)= 1.2383
       b(2)=-8.141e-1
       b(3)=0.0
       c=4.9924
    endif
    call qcalc(a,b,c,xi,flam,f)
    g = 1.d0-13.04d0*(flam**2) + 133.5d0*(flam**4) + 1534d0*(flam**6) + &
         918.6d0*(flam**8)
    gef=g*dexp(-2.d0/flam)*f
    temp=1.d0/(10.7480d0*flam**2+0.3967d0*flam**0.5d0+1.0050d0)
    smq=temp*(1.d0+(rho/mue)/(7.692d7*flam**3+9.715d6*flam**0.5d0))**(-0.3)
    cvm=(cv**2-ca**2)+n*(cvp**2-cap**2)
    cvp=(cv**2+ca**2)+n*(cvp**2+cap**2)
    qpair=0.5d0*cvp*(1.d0+(cvm/cvp)*smq)*gef
    return

  end subroutine pairn

!**********************************************************************

  subroutine qcalc(a,b,c,xi,flam,f)

! Also handles rates for pair, photo-, and plasma rates in Munkata etal.
! (1985) paper.

    implicit double precision (a-h,o-z)
    
    real*8, dimension(3) :: a,b
    real*8 :: c, xi, flam, f

    top=(a(1)+a(2)*xi+a(3)*(xi**2.))*dexp(-dble(c)*xi)	
    bottom=(xi**3.)+(b(1)/flam)+(b(2)/(flam**2))+(b(3)/(flam**3))
    f=top/bottom
    
    return

  end subroutine qcalc
  
!************************************************************************

  subroutine recombine(aavg,zbar,t,rho,cv,cvp,ca,cap,qrec)       

! Recombination neutrino rates of Kohyama, Itoh, Obama, & Mutoh,
! 1993, ApJ, 415, 267.
! These are much lower than Beaudet, Petrosian, & Salpeter rates,
! and are not believed to be important for white dwarf stars.

    use comp

!    implicit double precision(a-h,o-z)
    implicit none

    real*8 :: n,qfac1,qfac2,qfac3,aavg,zbar
    real*8 :: gnoo,gnoofac, zetfac, rjayarg
    real*8 :: t,rho,cv,cvp,ca,cap,qrec

    n = 2.0
    zetfac = 1.579d+5 * zbar * zbar / t

! gnoofac for case where A/Z = 2 
! no H version

    gnoofac = 5.526d+7*rho/(t**1.5)

!*** H version ***
    if(x(1) .ne. 0.0)then
       gnoofac = gnoofac * (1.d0 + (0.992*x(1)/1.008))
    endif

!*** now determine gnoo from gnoofac via inverse of F-D integral ***
    gnoo = finv(gnoofac)

!*** if nu falls outside the range, exit ***
    if(gnoo.lt.-20.0 .or. gnoo.gt.10.0)then
       qrec = 0.d0
       goto 96
    endif

!*** now compute Q from pieces ***
    qfac1 = (cv*cv + 1.5*ca*ca) + (n*(cvp*cvp +1.5*cap*cap))
    qfac2 = 2.649d-18 * zbar**14 * rho / aavg
    qfac3 = 1.d0 / (1.d0 + dexp(zetfac+gnoo))

!*** determine Q for recombination neutrinos. Qo for comparison to paper. ***
    call rjay(zetfac,gnoo,rjayarg)
    qrec = qfac1*qfac2*qfac3*rjayarg
    
96  return      
    
  end subroutine recombine

!******************************************************************************

  function finv(f)       

! Compute Fermi-Dirac integral F1/2 to accuracy of 10^-4 using 
! rational function approximation of Anita (1993, ApJS, 84, 101).

    implicit double precision(a-h,o-z)
    real*8, dimension(3) :: a1, b1, a2, b2
    real*8, parameter :: an = 0.5d0
    integer, parameter :: m1 = 2, k1 = 2, m2 = 2, k2 = 2

! low precision formulae

    data a1/4.4593646d+01,+1.1288764d+01,1.0000000d+00/
    data b1/3.9519346d+01,-5.7517464d+00,2.6594291d-01/
    data a2/3.4873722d+01,-2.6922515d+01,1.0000000d+00/
    data b2/2.6612832d+01,-2.0452930d+01,1.1808945d+01/

! small eta rational function       

    if(f .lt. 4.d0)then
       rn = f + a1(m1)
       do i=m1-1,1,-1
          rn = rn*f + a1(i)
       end do
       den = b1(k1+1)
       do i=k1,1,-1
          den = den*f + b1(i)
       end do
       finv = dlog(f * rn / den)
    else

! large eta rational function       

       ff = 1.d0 / f**(1.d0/(1.d0+an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
          rn = rn*ff + a2(i)
       end do
       den = b2(k2+1)
       do i=k2,1,-1
          den = den*ff + b2(i)
       end do
       finv =  rn / (den * ff)
    endif

    return      
     
  end function finv

!****************************************************************************
! Formerly function rjay. For some reason intel fortan compiler just could
! not recognize it as a function and not a variable so changed to a 
! subroutine.

  subroutine rjay(zetfac,gnoo,rjayarg)       

! Function to determine value of J for each element

    implicit double precision(a-h,o-z)

    if(gnoo.ge.-20.0 .and. gnoo.lt.0.d0)then
       zfac = zetfac
       rjnum = (0.0151/zfac)+(0.242/zfac**2.25)+(1.21/zfac**4.55)
       rjden = 1.d0 + 0.0371*dexp(0.906*gnoo)*(1.d0+0.928*zfac)
       rjayarg = rjnum*dexp(gnoo)/rjden
    elseif(gnoo .le. 10.0)then
       den = 1.d0+(-0.0120*gnoo)+(0.029*gnoo*gnoo)+(-0.00104*gnoo**3)
       zfac = zetfac/den
       rjnum = (0.0123/zfac)+(0.266/zfac**2.25)+(1.30/zfac**4.55)
       rjden = 1.d0 + 0.117*dexp(0.897*gnoo)*(1.d0+0.177*zfac)
       rjayarg = rjnum*dexp(gnoo)/rjden
    endif

    return      

  end subroutine rjay

!************************************************************************

  subroutine brempde(dens,t,gamma,qbremde)

! this subroutine returns the neutrino pair bremsstrahlung rxn rates
! for partially degenerate electrons for the liquid metal regime.
! From Munkata, Kohyama, & Itoh 1987, ApJ, 316, 708.
! Valid in the range:  Gamma 1 to 178. (DO NOT use below Gamma = 1)
!                      Temp  10**8 to 10**11
!                      Dens  1     to 10**12
! Z and AT are only places where composition dependence shows up.

    use comp

    implicit double precision (a-z)

    real*8, dimension(6) :: a
    real*8, dimension(5) :: b
    real*8 :: dens, t, gamma, qbremde

    t8 = 10.**(t - 8.0)
    d = 10.**(dens)
    at = 1. /(x(2)*4. + x(3)*12. + x(4)*16.)
    z = (x(2)*2. + x(3)*6. + x(4)*8.)
    en = 2.
    s2theta = 0.217
    cv = 0.5 + 2.*s2theta
    ca = 0.5
    cvp = 1. - cv
    cap = 1. - ca

! load coefficient arrays. NOTE: a(1) = a_0 etc.

    a(1) = 23.5 
    a(2) = 6.83e+04
    a(3) = 7.81e+04
    a(4) = 230.
    a(5) = 6.70e+05
    a(6) = 7.66e+09
    b(1) = 1.47
    b(2) = 0.0329
    b(3) = (7.75e+05*t8**1.5) + (247.*t8**3.85)
    b(4) = 4.07 + (0.0240*t8**1.40)
    b(5) = 4.59e-05*t8**(-0.110)

! crank it out.

    xi = d / ((7.05e+06*t8**1.5) + (5.12e+04*t8**3))
    x2 = xi*xi
    
    f1 = 1./(a(1)+(a(2)/t8**2)+(a(3)/t8**5))
    f2 = 1.26*(1.+(1./xi))/(1.+(b(1)/xi)+(b(2)/x2))
    fsum = f1 + f2
    g1 = (1. + 1.e-9*d)*(a(4)+(a(5)/t8**2))+(a(6)/t8**5)
    g2 = (b(3)/d) + b(4) + b(5)*d**0.656
    gsum = (1./g1) + (1./g2)
 
    first = 0.5*(cv*cv + ca*ca + en*cvp*cvp + en*cap*cap)*fsum
    sec = 0.5*(cv*cv - ca*ca + en*cvp*cvp - en*cap*cap)*gsum
    qbremde = 0.5738*(z*z/at)*(t8**6)*(first - sec)

    return

  end subroutine brempde

!*************************************************************************

  subroutine bremliqh(dens,t,gamma,qbremlh)

! This subroutine returns the bremsstrahlung neutrino rates computed for 
! Liquid Helium using the results of Itoh, Kohyama, Matsumoto, & Seki 
! 1984, ApJ, 280, 787 and Itoh & Kohyama 1983, ApJ, 275, 858
! These formulae include liquid metal and low temp. quantum correction conts.
! 10 12 92 Put in typo correction from erratum

    use bremliqh_data
    use misc_const, only: pi,g

    implicit double precision (a-z)

    real*8, dimension(6) :: a,e,i,p,aq,eq
    real*8, dimension(5) :: alphaq
    real*8, dimension(4) :: b,f,j,q,beta,bq,fq,alpha
    real*8  :: c, d, gg, h, k, l, r, s
    real*8 :: cq,dq,gq,hq
    real*8 :: dens, t, gamma, qbremlh

    do ii=1,6
       a(ii) = lnah(ii)
       e(ii) = lneh(ii)
       i(ii) = lnih(ii)
       p(ii) = lnph(ii)
       aq(ii) = qnah(ii)
       eq(ii) = qneh(ii)
    end do

    do ii=1,5
       alphaq(ii) = alphaqnh(ii)
    end do

    do ii=1,4
       b(ii) = lnbh(ii)
       f(ii) = lnfh(ii)
       j(ii) = lnjh(ii)
       q(ii) = lnqh(ii)
       alpha(ii) = alphalnh(ii)
       beta(ii) = betalnh(ii)
       bq(ii) = qnbh(ii)
       fq(ii) = qnfh(ii)
    end do

    e = lneh
    c = lnch
    d = lndh
    gg = lngh
    h = lnhh
    k = lnkh
    l = lnlh
    r = lnrh
    s = lnsh
    cq = qnch
    dq = qndh
    gq = qngh
    hq = qnhh

! en is the number of neutrino flavors besides electron.
! s2theta is sin(theta)**2, where theta is the Weinberg angle

    at = 4.
    z = 2.
    twopi = 2.*pi
    pisqr = pi**2
    en = 2.
    s2theta = 0.217
    cv = 0.5 + 2.*s2theta
    ca = 0.5
    cvp = 1. - cv
    cap = 1. - ca
    x13 = -1./3.
    xone3 = 1./3.
    x14 = -0.25000000
    x23 = -2./3.
    x24 = -0.50000000
    x34 = -0.75000000
    x43 = 4./3.
    t8 = 10.**t / 1.e+08
    d6 = 10.**dens / 1.e+06 
    alfa = 0.07279 * gamma * (d6**xone3) / ((z*at)**x43)
    
! compute the interpolation coefficients v and w.
! NOTE: u is defined differently here than in BREMSOL*.f
    u = twopi * ( dens - 3. ) / 10.
!*** liquid variables ***
    v = alpha(1) + alpha(2)*gamma**x13 + alpha(3)*gamma**x23 + alpha(4)/gamma
    w = beta(1) + beta(2)*gamma**x13 + beta(3)*gamma**x23 + beta(4)/gamma
!*** quantum correction variables ***
    vq = alphaq(1) + alphaq(2)*gamma**x14 + alphaq(3)*gamma**x24 + &
         alphaq(4)*gamma**x34 + alphaq(5) / gamma

! set up the variables (sines and cosines) for the loops
    su  = dsin(u)
    su2 = su * su 
    cu2 = 1. - su2
    cu  = dsqrt( cu2 )
    s2u = 2. * su * cu
    c2u = cu2 - su2
    s3u = su * ( 3. - 4. * su2 )
    c3u = cu * ( 4. * cu2 - 3. )
    c4u = -8. * cu2 * su2 + 1.
    s4u = su * ( (2.*c3u) + (2.*cu) )
    c5u = cu * ( (16.*cu2*cu2) - (20.*cu2) + 5. )

!  Liquid contributions

    asum = a(1) + a(2) * cu + a(3) * c2u  + a(4) * c3u + a(5) * c4u + a(6) * c5u
    bsum = b(1) * su + b(2) * s2u + b(3) * s3u + b(4) * s4u + c*u + d
    esum = e(1) + e(2) * cu + e(3) * c2u + e(4) * c3u + e(5) * c4u + e(6) * c5u
    fsum = f(1) * su + f(2) * s2u + f(3) * s3u + f(4) * s4u + gg*u + h
    isum = i(1) + i(2) * cu + i(3) * c2u + i(4) * c3u + i(5) * c4u + i(6) * c5u
    jsum = j(1) * su + j(2) * s2u + j(3) * s3u + j(4) * s4u + k*u + l
    psum = p(1) + p(2) * cu + p(3) * c2u + p(4) * c3u + p(5) * c4u + p(6) * c5u
    qsum = q(1) * su + q(2) * s2u + q(3) * s3u + q(4) * s4u + r*u + s

!  Quantum correction contributions
    aqsum = aq(1) + aq(2) * cu + aq(3) * c2u + aq(4) * c3u + aq(5) * c4u + &
         aq(6) * c5u
    bqsum = bq(1) * su + bq(2) * s2u + bq(3) * s3u + bq(4) * s4u + cq*u + dq
    eqsum = eq(1) + eq(2) * cu    + eq(3) * c2u + eq(4) * c3u + eq(5) * c4u + & 
         eq(6) * c5u
    fqsum = fq(1) * su + fq(2) * s2u + fq(3) * s3u + fq(4) * s4u + gq*u + hq
    
! sum up Liquid contributions

    fliq1 = asum + bsum
    fliq160  = esum + fsum
    gliq1 = isum + jsum
    gliq160  = psum + qsum

    fliq    = v * fliq1 + (1. - v) * fliq160
    gliq    = w * gliq1 + (1. - w) * gliq160

! sum up Quantum correction contributions
    rquan1 = aqsum + bqsum
    rquan160 = eqsum + fqsum
    delfof = (1. - vq)*rquan1 + vq * rquan160
    delgog = (1. - vq)*rquan1 + vq * rquan160

! sum up for total neutrino loss contributions in solid phase
    ftot = (1. - alfa * delfof) * fliq
    gtot = (1. - alfa * delgog) * gliq

! compute neutrino loss rate
    first = 0.5*(cv*cv + ca*ca + en*cvp*cvp + en*cap*cap)*ftot
    sec = 0.5*(cv*cv - ca*ca + en*cvp*cvp - en*cap*cap)*gtot
    qbremlh = 0.5738*(z*z/at)*(t8**6)*(first - sec)

    return

  end subroutine bremliqh
  
!************************************************************************

  subroutine bremliqc(dens,t,gamma,qbremlc)

! This subroutine returns the bremsstrahlung neutrino rates computed for 
! Liquid Carbon using the results of Itoh, Kohyama, Matsumoto, & Seki 
! 1984, ApJ, 280, 787 and Itoh & Kohyama 1983, ApJ, 275, 858
! These formulae include liquid metal and low temp. quantum correction conts.

    use bremliqc_data
    use misc_const, only: pi,g

    implicit double precision (a-z)

    real*8, dimension(6) :: a,e,i,p,aq,eq
    real*8, dimension(5) :: alphaq
    real*8, dimension(4) :: b,f,j,q,alpha,beta,bq,fq
    real*8 :: c,d,gg,h,k,l,r,s,cq,dq,gq,hq
    real*8 :: dens, t, gamma, qbremlc
    
    
    do ii=1,6
       a(ii) = lnac(ii)
       e(ii) = lnec(ii)
       i(ii) = lnic(ii)
       p(ii) = lnpc(ii)
       aq(ii) = qnac(ii)
       eq(ii) = qnec(ii)
    end do

    do ii=1,5
       alphaq(ii) = alphaqnc(ii)
    end do

    do ii=1,4
       b(ii) = lnbc(ii)
       f(ii) = lnfc(ii)
       j(ii) = lnjc(ii)
       q(ii) = lnqc(ii)
       alpha(ii) = alphalnc(ii)
       beta(ii) = betalnc(ii)
       bq(ii) = qnbc(ii)
       fq(ii) = qnfc(ii)
    end do

    c = lncc
    d = lndc
    gg = lngc
    h = lnhc
    k = lnkc
    l = lnlc
    r = lnrc
    s = lnsc 
    cq = qncc
    dq = qndc
    gq = qngc
    hq = qnhc

    at = 12.
    z = 6.
    twopi = 2.*pi
    pisqr = pi**2
    en = 2.
    s2theta = 0.217
    cv = 0.5 + 2.*s2theta
    ca = 0.5
    cvp = 1. - cv
    cap = 1. - ca
    x13 = -1./3.
    xone3 = 1./3.
    x14 = -0.25000000
    x23 = -2./3.
    x24 = -0.50000000
    x34 = -0.75000000
    x43 = 4./3.
    t8 = 10.**t / 1.e+08
    d6 = 10.**dens / 1.e+06 
    alfa = 0.07279 * gamma * (d6**xone3) / ((z*at)**x43)

! compute the interpolation coefficients v and w.
! NOTE: u is defined differently here than in BREMSOL*.f

    u = twopi * ( dens - 3. ) / 10.

!*** liquid variables ***
    v = alpha(1) + alpha(2)*gamma**x13 + alpha(3)*gamma**x23 + alpha(4) / gamma
    w = beta(1) + beta(2) * gamma**x13 + beta(3) * gamma**x23 + beta(4) / gamma

!*** quantum correction variables ***
    vq = alphaq(1) + alphaq(2)*gamma**x14 + alphaq(3)*gamma**x24 &
         + alphaq(4)*gamma**x34 + alphaq(5) / gamma

! set up the variables (sines and cosines) for the loops

    su  = dsin(u)
    su2 = su * su 
    cu2 = 1. - su2
    cu  = dsqrt( cu2 )
    s2u = 2. * su * cu
    c2u = cu2 - su2
    s3u = su * ( 3. - 4. * su2 )
    c3u = cu * ( 4. * cu2 - 3. )
    c4u = -8. * cu2 * su2 + 1.
    s4u = su * ( (2.*c3u) + (2.*cu) )
    c5u = cu * ( (16.*cu2*cu2) - (20.*cu2) + 5. )

!  Liquid contribution sums

    asum = a(1) + a(2) * cu + a(3) * c2u + a(4) * c3u + a(5) * c4u + a(6) * c5u
    bsum = b(1) * su + b(2) * s2u + b(3) * s3u + b(4) * s4u + c*u + d
    esum = e(1) + e(2) * cu + e(3) * c2u + e(4) * c3u + e(5) * c4u + e(6) * c5u
    fsum = f(1) * su + f(2) * s2u + f(3) * s3u + f(4) * s4u + gg*u + h
    isum = i(1) + i(2) * cu + i(3) * c2u + i(4) * c3u + i(5) * c4u + i(6) * c5u
    jsum = j(1) * su + j(2) * s2u + j(3) * s3u + j(4) * s4u + k*u + l
    psum = p(1) + p(2) * cu + p(3) * c2u + p(4) * c3u + p(5) * c4u + p(6) * c5u
    qsum = q(1) * su + q(2) * s2u + q(3) * s3u + q(4) * s4u + r*u + s

!  Quantum correction contributions

    aqsum = aq(1) + aq(2) * cu + aq(3) * c2u + aq(4) * c3u + aq(5) * c4u + &
         aq(6) * c5u
    bqsum = bq(1) * su + bq(2) * s2u + bq(3) * s3u + bq(4) * s4u + cq*u + dq
    eqsum = eq(1) + eq(2) * cu    + eq(3) * c2u + eq(4) * c3u + eq(5) * c4u + &
         eq(6) * c5u
    fqsum = fq(1) * su + fq(2) * s2u + fq(3) * s3u + fq(4) * s4u + gq*u + hq
    
! sum up Liquid contributions

    fliq1 = asum + bsum
    fliq160  = esum + fsum
    gliq1 = isum + jsum
    gliq160  = psum + qsum

    fliq    = v * fliq1 + (1. - v) * fliq160
    gliq    = w * gliq1 + (1. - w) * gliq160

! sum up Quantum correction contributions

    rquan1 = aqsum + bqsum
    rquan160 = eqsum + fqsum

    delfof = (1. - vq)*rquan1 + vq * rquan160
    delgog = (1. - vq)*rquan1 + vq * rquan160

! sum up for total neutrino loss contributions in liquid phase
    
    ftot = (1. - alfa * delfof) * fliq
    gtot = (1. - alfa * delgog) * gliq

! compute neutrino loss rate

    first = 0.5*(cv*cv + ca*ca + en*cvp*cvp + en*cap*cap)*ftot
    sec = 0.5*(cv*cv - ca*ca + en*cvp*cvp - en*cap*cap)*gtot
    qbremlc = 0.5738*(z*z/at)*(t8**6)*(first - sec)
    
    return

  end subroutine bremliqc

!************************************************************************

  subroutine bremliqo(dens,t,gamma,qbremlo)

! This subroutine returns the bremsstrahlung neutrino rates computed for 
! Liquid Oxygen using the results of Itoh, Kohyama, Matsumoto, & Seki 
! 1984, ApJ, 280, 787 and Itoh & Kohyama 1983, ApJ, 275, 858
! These formulae include liquid metal and low temp. quantum correction conts.
  
    use bremliqo_data
    use misc_const, only: pi,g

    implicit double precision (a-z)
  
! Because variables used in this routine do not have the same name as those in
! the common block modules, use pointers (sorry).

    real*8, dimension(6) :: a,e,i,p,aq,eq
    real*8, dimension(5) :: alphaq
    real*8, dimension(4) :: b,f,j,q,alpha,beta,bq,fq
    real*8 :: c,d,gg,h,k,l,r,s,cq,dq,gq,hq
    real*8 :: dens, t, gamma, qbremlo

    do ii=1,6
       a(ii) = lnao(ii)
       e(ii) = lneo(ii)
       i(ii) = lnio(ii)
       p(ii) = lnpo(ii)
       aq(ii) = qnao(ii)
       eq(ii) = qneo(ii)
    end do

    do ii=1,5
       alphaq(ii) = alphaqno(ii)
    end do

    do ii=1,4
       b(ii) = lnbo(ii)
       f(ii) = lnfo(ii)
       j(ii) = lnjo(ii)
       q(ii) = lnqo(ii)
       alpha(ii) = alphalno(ii)
       beta(ii) = betalno(ii)
       bq(ii) = qnbo(ii)
       fq(ii) = qnfo(ii)
    end do

    c = lnco
    d = lndo
    gg = lngo
    h = lnho
    k = lnko
    l = lnlo
    r = lnro
    s = lnso
    cq = qnco
    dq = qndo
    gq = qngo
    hq = qnho

    at = 16.
    z = 8.
    twopi = 2.*pi
    pisqr = pi**2
    en = 2.
    s2theta = 0.217
    cv = 0.5 + 2.*s2theta
    ca = 0.5
    cvp = 1. - cv
    cap = 1. - ca
    x13 = -1./3.
    xone3 = 1./3.
    x14 = -0.25000000
    x23 = -2./3.
    x24 = -0.50000000
    x34 = -0.75000000
    x43 = 4./3.
    t8 = 10.**t / 1.e+08
    d6 = 10.**dens / 1.e+06 
    alfa = 0.07279 * gamma * (d6**xone3) / ((z*at)**x43)

! compute the interpolation coefficients v and w.
! NOTE: u is defined differently here than in BREMSOL*.f

    u = twopi * ( dens - 3. ) / 10.

!*** liquid variables ***
    v = alpha(1) + alpha(2)*gamma**x13 + alpha(3)*gamma**x23 + alpha(4)/gamma
    w = beta(1) + beta(2) * gamma**x13 + beta(3) * gamma**x23 + beta(4)/gamma

!*** quantum correction variables ***
    vq = alphaq(1) + alphaq(2)*gamma**x14 + alphaq(3)*gamma**x24 &
         + alphaq(4)*gamma**x34 + alphaq(5) / gamma

! set up the variables (sines and cosines) for the loops

    su  = dsin(u)
    su2 = su * su 
    cu2 = 1. - su2
    cu  = dsqrt( cu2 )
    s2u = 2. * su * cu
    c2u = cu2 - su2
    s3u = su * ( 3. - 4. * su2 )
    c3u = cu * ( 4. * cu2 - 3. )
    c4u = -8. * cu2 * su2 + 1.
    s4u = su * ( (2.*c3u) + (2.*cu) )
    c5u = cu * ( (16.*cu2*cu2) - (20.*cu2) + 5. )

!  Liquid contributions

    asum = a(1) + a(2) * cu + a(3) * c2u + a(4) * c3u + a(5) * c4u + a(6) * c5u
    bsum = b(1) * su + b(2) * s2u + b(3) * s3u + b(4) * s4u + c*u + d
    esum = e(1) + e(2) * cu + e(3) * c2u + e(4) * c3u + e(5) * c4u + e(6) * c5u
    fsum = f(1) * su + f(2) * s2u + f(3) * s3u + f(4) * s4u + gg*u + h
    isum = i(1) + i(2) * cu + i(3) * c2u + i(4) * c3u + i(5) * c4u + i(6) * c5u
    jsum = j(1) * su + j(2) * s2u + j(3) * s3u + j(4) * s4u + k*u + l
    psum = p(1) + p(2) * cu + p(3) * c2u + p(4) * c3u + p(5) * c4u + p(6) * c5u
    qsum = q(1) * su + q(2) * s2u + q(3) * s3u + q(4) * s4u + r*u + s

! Quantum correction contributions

    aqsum = aq(1) + aq(2) * cu + aq(3) * c2u + aq(4) * c3u + aq(5) * c4u + &
         aq(6) * c5u
    bqsum = bq(1) * su + bq(2) * s2u + bq(3) * s3u + bq(4) * s4u + cq*u + dq
    eqsum = eq(1) + eq(2) * cu + eq(3) * c2u + eq(4) * c3u + eq(5) * c4u + &
         eq(6) * c5u
    fqsum = fq(1) * su + fq(2) * s2u + fq(3) * s3u + fq(4) * s4u + gq*u + hq

! sum up Liquid contributions

    fliq1 = asum + bsum
    fliq160  = esum + fsum
    gliq1 = isum + jsum
    gliq160  = psum + qsum

    fliq    = v * fliq1 + (1. - v) * fliq160
    gliq    = w * gliq1 + (1. - w) * gliq160

! sum up Quantum correction contributions

    rquan1 = aqsum + bqsum
    rquan160 = eqsum + fqsum

    delfof = (1. - vq)*rquan1 + vq * rquan160
    delgog = (1. - vq)*rquan1 + vq * rquan160

! sum up for total neutrino loss contributions in liquid phase

    ftot = (1. - alfa * delfof) * fliq
    gtot = (1. - alfa * delgog) * gliq

! compute neutrino loss rate

    first = 0.5*(cv*cv + ca*ca + en*cvp*cvp + en*cap*cap)*ftot
    sec = 0.5*(cv*cv - ca*ca + en*cvp*cvp - en*cap*cap)*gtot
    qbremlo = 0.5738*(z*z/at)*(t8**6)*(first - sec)

    return

  end subroutine bremliqo

!************************************************************************

  subroutine bremsolc(dens,t,gamma,qbremsc)

! This routine returns the bremsstrahlung neutrino rates computed for 
! Solid Carbon using the results of 
! Itoh, Kohyama, Matsumoto, & Seki 1984, ApJ, 285, 304 
! (+Errata 1987, ApJ, 322, 584)
! These formulae include lattice and phonon (lattice vibration) contributions

    use bremsolc_data
    use misc_const, only: pi,g
    use soldatc
    use phondatc

    implicit double precision (a-z)

! Because variables used in this routine do not have the same name as those in
! the common block modules, use pointers (sorry).

    real*8, pointer, dimension(:) :: a,e,i,p,alpha,beta,b, f, j, q
    real*8, pointer :: c,d,gg,h,k,l,r,s
    real*8, pointer, dimension(:) :: ap,ip,alphap,betap,bp,jp
    real*8, pointer :: cp,dp,kp,lp
    real*8 :: dens, t, gamma, qbremsc

! soldatc variables

    a => cb_snac
    b => cb_snbc
    e => cb_snec
    f => cb_snfc
    i => cb_snic
    j => cb_snjc
    p => cb_snpc
    q => cb_snqc
    alpha => cb_alphasnc
    beta => cb_betasnc
    c => cb_sncc
    d => cb_sndc
    gg => cb_sngc
    h => cb_snhc
    k => cb_snkc
    l => cb_snlc
    r => cb_snrc
    s => cb_snsc 

! phondatc variables

    ap => cb_pnac
    bp => cb_pnbc
    ip => cb_pnic
    jp => cb_pnjc
    alphap => cb_alphapnc
    betap => cb_betapnc
    cp => cb_pncc
    dp => cb_pndc
    kp => cb_pnkc
    lp => cb_pnlc

    at = 12.
    z = 6.
    twopi = 2.*pi
    pisqr = pi**2
    en = 2.
    s2theta = 0.217
    cv = 0.5 + 2.*s2theta
    ca = 0.5
    cvp = 1. - cv
    cap = 1. - ca
    x13 = -1./3.
    x23 = -2./3.
    t8 = 10.**t / 1.e+08

! compute the interpolation coefficients v and w

    u = twopi * ( dens - 3. ) / 9.

! lattice variables If gamma>5000., set to 5000. for the lattice results.

    if(gamma.gt.5000.)then
       gamma2 = 5000.
    else
       gamma2 = gamma
    endif
    v = alpha(1) + alpha(2) * gamma2**x13 + alpha(3) * gamma2**x23 + &
         alpha(4) / gamma2
    w = beta(1) + beta(2) * gamma2**x13 + beta(3) * gamma2**x23 + &
         beta(4) / gamma2
    
! phonon variables 

    vp = alphap(1)+alphap(2)*gamma**x13 + alphap(3)*gamma**x23+alphap(4)/gamma
    wp = betap(1) + betap(2)*gamma**x13 + betap(3)*gamma**x23 + betap(4)/gamma

! set up the variables (sines and cosines) for the loops

    su  = dsin(u)
    su2 = su * su 
    cu2 = 1. - su2
    cu  = dsqrt( cu2 )
    s2u = 2. * su * cu
    c2u = cu2 - su2
    s3u = su * ( 3. - 4. * su2 )
    c3u = cu * ( 4. * cu2 - 3. )
    c4u = -8. * cu2 * su2 + 1.

! Lattice contributions

    asum = a(1) + a(2) * cu    + a(3) * c2u + a(4) * c3u + a(5) * c4u
    bsum = b(1) * su + b(2) * s2u + b(3) * s3u + c*u + d
    esum = e(1) + e(2) * cu    + e(3) * c2u + e(4) * c3u + e(5) * c4u
    fsum = f(1) * su + f(2) * s2u + f(3) * s3u + gg*u + h
    isum = i(1) + i(2) * cu    + i(3) * c2u + i(4) * c3u + i(5) * c4u
    jsum = j(1) * su + j(2) * s2u + j(3) * s3u + k*u + l
    psum = p(1) + p(2) * cu    + p(3) * c2u + p(4) * c3u + p(5) * c4u
    qsum = q(1) * su + q(2) * s2u + q(3) * s3u + r*u + s

!  Phonon contributions

    apsum = ap(1) + ap(2)*cu + ap(3)*c2u + ap(4)*c3u + ap(5)*c4u
    bpsum = bp(1)*su + bp(2)*s2u + bp(3)*s3u + cp*u + dp
    ipsum = ip(1) + ip(2)*cu + ip(3)*c2u + ip(4)*c3u + ip(5)*c4u
    jpsum = jp(1)*su + jp(2)*s2u + jp(3)*s3u + kp*u + lp
  
! sum up Lattice contributions

    flat171 = asum + bsum
    flat5k  = esum + fsum
    glat171 = isum + jsum
    glat5k  = psum + qsum
    flat    = (1.-v) * flat171 + v * flat5k
    glat    = (1.-w) * glat171 + w * glat5k

! sum up Phonon contributions
! We just extrapolate this is gamma > 5000.

    fplat171 = vp*(apsum + bpsum)
    gplat171 = wp*(ipsum + jpsum)
    
! sum up for total neutrino loss contributions in solid phase

    ftot = flat + fplat171
    gtot = glat + gplat171

! compute neutrino loss rate

    first = 0.5*(cv*cv + ca*ca + en*cvp*cvp + en*cap*cap)*ftot
    sec = 0.5*(cv*cv - ca*ca + en*cvp*cvp - en*cap*cap)*gtot
    qbremsc = 0.5738*(z*z/at)*(t8**6)*(first - sec)
    
    return

  end subroutine bremsolc

!************************************************************************

  subroutine bremsolo(dens,t,gamma,qbremso)

! This routine returns the bremsstrahlung neutrino rates computed for 
! Solid Oxygen using the results of 
! Itoh, Kohyama, Matsumoto, & Seki 1984, ApJ, 285, 304 
! (+Errata 1987, ApJ, 322, 584)
! Includes lattice and phonon (lattice vibrations) contributions.

  use bremsolo_data
  use misc_const, only: pi,g
  use soldato
  use phondato

  implicit double precision (a-z)

! Because variables used in this routine do not have the same name as those in
! the common block modules, use pointers (sorry).

    real*8, pointer, dimension(:) :: a,e,i,p,alpha,beta,b,f,j,q
    real*8, pointer :: c,d,gg,h,k,l,r,s
    real*8, pointer, dimension(:) :: ap,ip,alphap,betap,bp,jp
    real*8, pointer :: cp,dp,kp,lp
    real*8 :: dens, t, gamma, qbremso

! soldato variables

    a => cb_snao
    b => cb_snbo
    e => cb_sneo
    f => cb_snfo
    i => cb_snio
    j => cb_snjo
    p => cb_snpo
    q => cb_snqo
    alpha => cb_alphasno
    beta => cb_betasno
    c => cb_snco
    d => cb_sndo
    gg => cb_sngo
    h => cb_snho
    k => cb_snko
    l => cb_snlo
    r => cb_snro
    s => cb_snso

! phondato variables

    ap => cb_pnao
    bp => cb_pnbo
    ip => cb_pnio
    jp => cb_pnjo
    alphap => cb_alphapno
    betap => cb_betapno
    cp => cb_pnco
    dp => cb_pndo
    kp => cb_pnko
    lp => cb_pnlo

    at = 16.
    z = 8.
    twopi = 2.*pi
    pisqr = pi**2
    en = 2.
    s2theta = 0.217
    cv = 0.5 + 2.*s2theta
    ca = 0.5
    cvp = 1. - cv
    cap = 1. - ca
    x13 = -1./3.
    x23 = -2./3.
    t8 = 10.**t / 1.e+08

! compute the interpolation coefficients v and w

    u = twopi * ( dens - 3. ) / 9.

! lattice variables If gamma>5000., set to 5000. for lattice results.

    if(gamma.gt.5000.)then
       gamma2 = 5000.
    else
       gamma2 = gamma
    endif
    v = alpha(1) + alpha(2) * gamma2**x13 + alpha(3) * gamma2**x23 + &
         alpha(4) / gamma2
    w = beta(1) + beta(2)*gamma2**x13 + beta(3)*gamma2**x23 + beta(4)/gamma2

! phonon variables

    vp = alphap(1)+ alphap(2)*gamma**x13+ alphap(3)*gamma**x23+ alphap(4)/gamma
    wp = betap(1) + betap(2)*gamma**x13 + betap(3)*gamma**x23 + betap(4) /gamma

! set up the variables (sines and cosines) for the loops

    su  = dsin(u)
    su2 = su * su 
    cu2 = 1. - su2
    cu  = dsqrt( cu2 )
    s2u = 2. * su * cu
    c2u = cu2 - su2
    s3u = su * ( 3. - 4. * su2 )
    c3u = cu * ( 4. * cu2 - 3. )
    c4u = -8. * cu2 * su2 + 1.

!  Lattice contributions

    asum = a(1) + a(2) * cu    + a(3) * c2u + a(4) * c3u + a(5) * c4u
    bsum = b(1) * su + b(2) * s2u + b(3) * s3u + c*u + d
    esum = e(1) + e(2) * cu    + e(3) * c2u + e(4) * c3u + e(5) * c4u
    fsum = f(1) * su + f(2) * s2u + f(3) * s3u + gg*u + h
    isum = i(1) + i(2) * cu    + i(3) * c2u + i(4) * c3u + i(5) * c4u
    jsum = j(1) * su + j(2) * s2u + j(3) * s3u + k*u + l
    psum = p(1) + p(2) * cu    + p(3) * c2u + p(4) * c3u + p(5) * c4u
    qsum = q(1) * su + q(2) * s2u + q(3) * s3u + r*u + s

!  Phonon contributions

    apsum = ap(1) + ap(2)*cu + ap(3)*c2u + ap(4)*c3u + ap(5)*c4u
    bpsum = bp(1)*su + bp(2)*s2u + bp(3)*s3u + cp*u + dp
    ipsum = ip(1) + ip(2)*cu + ip(3)*c2u + ip(4)*c3u + ip(5)*c4u
    jpsum = jp(1)*su + jp(2)*s2u + jp(3)*s3u + kp*u + lp

! sum up Lattice contributions

    flat171 = asum + bsum
    flat5k  = esum + fsum
    glat171 = isum + jsum
    glat5k  = psum + qsum
    flat    = (1.-v) * flat171 + v * flat5k
    glat    = (1.-w) * glat171 + w * glat5k

! sum up Phonon contributions
! We just extrapolate this is gamma > 5000.

    fplat171 = vp*(apsum + bpsum)
    gplat171 = wp*(ipsum + jpsum)

! sum up for total neutrino loss contributions in solid phase
    
    ftot = flat + fplat171
    gtot = glat + gplat171
    
! compute neutrino loss rate

    first = 0.5*(cv*cv + ca*ca + en*cvp*cvp + en*cap*cap)*ftot
    sec = 0.5*(cv*cv - ca*ca + en*cvp*cvp - en*cap*cap)*gtot
    qbremso = 0.5738*(z*z/at)*(t8**6)*(first - sec)
    
    return

  end subroutine bremsolo

!******************************************************************************
! Axion rates according to Nakagawa, Kohyama and Itoh (1987), Ap.J. 322, 291 
! and Nakagawa, Adachi, Kohyama and Itoh (1988), Ap.J. 326, 241 (phonon
! contributions).
! Coupling constant g given in Isern, Hernanz and Garcia-Berro (1992), Ap.J.
! 392, L23.
!******************************************************************************
! June 2004, Agnes.
!******************************************************************************

  subroutine axions(Xc, d, t, epsaxion)

! maxion is actually ma * cos^2(beta), the "mass" of the axion in eV. 
! It is an input parameter, defined in wd40.axions.dat

    use ax
    use misc_const, only: pi,g

    implicit none
      
    real*8 :: T8,T7,rho,rho6,gamma,gax,alphaa,epsaxion, &
         d,t,Fliquid,Flattice,Fphonon,Xc,F12,F16,Zc,Zo,Ac,Ao, &
         gammc,gammo,eps_c,eps_o,Xo

    Xo=1.-Xc

    Zc=6.
    Zo=8.
    Ac=12.
    Ao=16.
    
! First calculate gamma, which is a measure of the degree of crystallization
    T8 = 10**t/1.d8
    T7 = 10**t/1.d7
    rho  = 10**d
    rho6 = 10**d/1.d6
    gammc = 2.275d-1 * Zc**2/T8 *(rho6/Ac)**(1./3.)
    gammo = 2.275d-1 * Zo**2/T8 *(rho6/Ao)**(1./3.)
    gamma = Xc*gammc + Xo*gammo

! For gamma < 178 (non-crystallized), contribution from Fliquid
    if (gamma .le. 178.) then
       call calcFliquid12(rho, gammc, Fliquid)
       F12 = Fliquid
       call calcFliquid16(rho, gammo, Fliquid)
       F16 = Fliquid
    end if

! For gamma > 178 (crystallized), contribution from Flattice and Fhonon
    if (gamma .ge. 178.) then
       call calcFlattice12(rho, gammc, Flattice)
       call calcFphonon12( rho, gammc, Fphonon)
       F12 = Flattice + Fphonon
       call calcFlattice16(rho, gammo, Flattice)
       call calcFphonon16( rho, gammo, Flattice)
       F16 = Flattice + Fphonon
    end if

    gax = 8.5e-11/3.*maxion
    alphaa = gax**2/(4*pi)
    eps_c = 1.08d23 * alphaa * Zc**2/Ac * T7**4 * F12
    eps_o = 1.08d23 * alphaa * Zo**2/Ao * T7**4 * F16
    epsaxion = Xc*eps_c + Xo*eps_o

  end subroutine axions

!************************************************************************* 
! Axion routine that calculates a number called F in the case of crystal-
! lized pure carbon. Called by axions.
!*************************************************************************

  subroutine calcFlattice12(rho, gamma, Flattice)
  
    use misc_const, only: pi,g
    
    implicit none
  
    integer :: n
    real*8 ::  k,l,r,s,rho,gamma,Flattice,rho6,u,w,F180,F5000,logFlattice

! Coefficients for C12
    real*8, dimension(5) :: i,p
    real*8, dimension(4) :: j,q,beta
    data i / 3.565, -1.1495, -0.2968, -0.0738, -0.0694 /
    data p / 3.621, -1.1746, -0.297,  -0.0745, -0.0696 /
    data j     / -0.4141, 0.0072, 0.0007, -0.003  /
    data q     / -0.4223, 0.012,  0.0007, -0.0033 /
    data beta / 0.1906, -5.6458, 32.461, 142.6  /  
      
    k=-0.3664
    l=-1.7456
    r=-0.4079
    s=-1.4663
  
    rho6 = rho/1.e6
    u = 2.*pi*(log10(rho)+5.)/18.
    w = 0
    do n=1,4
       w = w + beta(n) * gamma**(-real(n-1)/3.)
    enddo
    F180 = i(1)/2. + k*u + l
    F5000 = p(1)/2. + r*u + s
    do n=1,4
       F180  =  F180 +  i(n+1)*cos(real(n)*u)  +  j(n)*sin(real(n)*u)
       F5000 = F5000 +  p(n+1)*cos(real(n)*u)  +  q(n)*sin(real(n)*u)
    enddo
  
    logFlattice = w*F180 + (1.- w)*F5000
    Flattice = 10**logFlattice

  end subroutine calcFlattice12

!************************************************************************
! Axion routine that calculates a number called F in the case of crystal-
! lized pure oxygen. Called by axions.
!*************************************************************************
  subroutine calcFlattice16(rho, gamma, Flattice)
      
    use misc_const, only: pi,g
    
    implicit none
    
    integer :: n
    real*8 :: k,l,r,s,rho,gamma,Flattice,rho6,u,w,F180,F5000,logFlattice

! Coefficients for O16
    real*8, dimension(5) :: i,p
    real*8, dimension(4) :: j,q,beta
    data i / 3.3487,-1.0898, -0.2755, -0.067 , -0.063  /
    data p / 3.3263,-1.0877, -0.2682, -0.0686, -0.0621 /
    data j     / -0.3875, 0.0113, 0.0049, -0.0042 /
    data q     / -0.4107, 0.0197, 0.003 , -0.0036 /
    data beta  / 0.2263, -7.7843, 66.147, 13.781 /  
  
    k=-0.3303
    l=-1.6563
    r=-0.4079
    s=-1.1416
  
    rho6 = rho/1.e6
    u = 2.*pi*(log10(rho)+5.)/18.
    w = 0
    do n=1,4
       w = w + beta(n) * gamma**(-real(n-1)/3.)
    enddo
    F180 = i(1)/2. + k*u + l
    F5000 = p(1)/2. + r*u + s
    do n=1,4
       F180  =  F180 +  i(n+1)*cos(real(n)*u)  +  j(n)*sin(real(n)*u)
       F5000 = F5000 +  p(n+1)*cos(real(n)*u)  +  q(n)*sin(real(n)*u)
    enddo
    
    logFlattice = w*F180 + (1.- w)*F5000
    Flattice = 10**logFlattice
  
  end subroutine calcFlattice16
      
!*************************************************************************   
! Axion routine that calculates a number called F in the case of liquid
! pure carbon. Called by axions.
!*************************************************************************

  subroutine calcFliquid12(rho, gamma, Fliquid)

    use misc_const, only: pi,g      
    
    implicit none
    integer :: i
    real*8 :: c,d,gg,h,rho,gamma,Fliquid,rho6,u,v,F1,F160,logFliquid
    
!Coefficients for C12
    real*8, dimension(5) :: a,e
    real*8, dimension(4) :: b,f,alpha
    data a / 2.7337,-0.8648,-0.2367,-0.0715,-0.0477 /
    data e / 3.1029,-1.0355,-0.247, -0.0551,-0.0558 /
    data b     /-0.345, -0.0135, 0.0132, 0.0022 /
    data f     /-0.3332, 0.0271, 0.005, -0.0026 /
    data alpha /-0.3423, 2.2053,-2.1415, 1.3474 /  

    c=-0.1395
    d=-1.3894
    gg=-0.3146
    h=-1.3409
  
    rho6 = rho/1.e6
    u = 2.*pi*(log10(rho)+5.)/18.
    v = 0
    do i=1,4
       v = v + alpha(i) * gamma**(-real(i-1)/3.)
    enddo
    F1 = a(1)/2. + c*u + d
    F160 = e(1)/2. + gg*u + h
    
    do i=1,4
       F1   =   F1 +  a(i+1)*cos(real(i)*u)  +  b(i)*sin(real(i)*u)
       F160 = F160 +  e(i+1)*cos(real(i)*u)  +  f(i)*sin(real(i)*u)
    enddo
    
    logFliquid = v*F1 + (1.-v)*F160
    Fliquid = 10**logFliquid

  end subroutine calcFliquid12

!*************************************************************************
! Axion routine that calculates a number called F in the case of liquid
! pure oxygen. Called by axions.
!*************************************************************************

  subroutine calcFliquid16(rho, gamma, Fliquid)
  
    use misc_const, only: pi,g
    
    implicit none
    integer :: i
    real*8 :: c,d,gg,h,rho,gamma,Fliquid,rho6,u,v,F1,F160,logFliquid
    
! Coefficients for O16
    real*8, dimension(5) :: a,e
    real*8, dimension(4) :: b,f,alpha
    data a / 2.7258,-0.862, -0.2375,-0.0717,-0.0468 /
    data e / 3.1015,-1.0292,-0.2499,-0.0572,-0.0553 /
    data b     /-0.3408,-0.0139, 0.0142, 0.0029 /
    data f     /-0.34  , 0.0218 ,0.0077,-0.003  /
    data alpha /-0.1407, 0.4307, 1.9644,-1.2932 /  
           
    c = -0.1258
    d = -1.4084 
    gg = -0.2959
    h = -1.3475

    rho6 = rho/1.e6
    u = 2.*pi*(log10(rho)+5.)/18.
    v = 0
    do i=1,4
       v = v + alpha(i) * gamma**(-real(i-1)/3.)
    enddo
    F1 = a(1)/2. + c*u + d
    F160 = e(1)/2. + gg*u + h
    do i=1,4
       F1   =   F1 +  a(i+1)*cos(real(i)*u)  +  b(i)*sin(real(i)*u)
       F160 = F160 +  e(i+1)*cos(real(i)*u)  +  f(i)*sin(real(i)*u)
    enddo
    
    logFliquid = v*F1 + (1.-v)*F160
    Fliquid = 10**logFliquid
    
  end subroutine calcFliquid16

!*********************************************************************
! Axion routine that calculates the phonon contribution to the number 
! called F in the case pure carbon. Called by axions.
!*************************************************************************  

  subroutine calcFphonon12(rho, gamma, Fphonon)

    use misc_const, only: pi,g      

    implicit none
    integer :: n
    real*8 :: k,l,rho,gamma,Fphonon,rho6,x,F180,u,logFphonon
    
! Coefficients for C12
    real*8, dimension(5) :: i
    real*8, dimension(4) :: j,gama
    data i    / 3.8289, -1.1987, -0.3269, -0.0939, -0.0787 /
    data j    / 0.5103, -0.0101, -0.0034, -0.0018  /
    data gama / 12.118, -197.1, 1253.6, -2804.8 /  
      
    k=-0.4501
    l=-1.9453
    
    rho6 = rho/1.e6
    u = 2.*pi*(log10(rho)+5.)/18.
    x = 0
    do n=1,4
       x = x + gama(n) * gamma**(-real(n-1)/3.)
    enddo
    F180 = i(1)/2. + k*u + l
    do n=1,4
       F180  =  F180 +  i(n+1)*cos(real(n)*u)  +  j(n)*sin(real(n)*u)
    enddo
    
    logFphonon = x*F180
    Fphonon = 10**logFphonon
    
  end subroutine calcFphonon12

!************************************************************************
! Axion routine that calculates the phonon contribution to the number 
! called F in the case pure carbon. Called by axions.
!*************************************************************************

  subroutine calcFphonon16(rho, gamma, Fphonon)
    
    use misc_const, only: pi,g
    
    implicit none
    
    integer :: n
    real*8 ::  k,l,rho,gamma,Fphonon,rho6,x,F180,u,logFphonon
    
!Coefficients for O16
    real*8, dimension(5) :: i
    real*8, dimension(4) :: j,gama
    data i    /  3.3882, -1.0837, -0.2774, -0.0776, -0.0683 /
    data j    / -0.4482,  0.0079, -0.0036, -0.0035  /
    data gama / 13.343, -219, 1381.9, -3053.7 /  
    
    k=-0.4237
    l=-1.5933
    
    rho6 = rho/1.e6
    u = 2.*pi*(log10(rho)+5.)/18.
    x = 0
    do n=1,4
       x = x + gama(n) * gamma**(-real(n-1)/3.)
    enddo
    F180 = i(1)/2. + k*u + l
    do n=1,4
       F180  =  F180 +  i(n+1)*cos(real(n)*u)  +  j(n)*sin(real(n)*u)
    enddo
    
    logFphonon = x*F180
    Fphonon = 10**logFphonon
  
  end subroutine calcFphonon16

!****************************************************************************

end module nuax_subroutines
