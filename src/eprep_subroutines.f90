module eprep_subroutines

contains

  !************************************************************************

  subroutine eprep

    ! EPREP -- white dwarf carbon stratified prep code subroutine
    !  writes prepped models to tape50

    ! Subroutines
    use utils_subroutines, only : check_time
    use chemprofiles_subroutines, only : comp_env

    ! Common blocks
    use alpha
    use bldx
    use cc
    use comp
    use crash
    use d 
    use dfin
    use dprep
    use flags
    use grad
    use jee
    use kapder
    use kappa, only : opalz
    use misc_const
    use opacfin
    use outfile
    use shells, only : nshint
    use shells2
    use tape50
    use temp
    use terp
    use tfitt
    use wprep
    use wprop
    use xcompp
 

    implicit double precision(a-h,o-z)

    real*8, dimension(7) :: xavec
    integer :: newtime
    integer, pointer :: j

    j => jj
987 format(a80)
1002 format(1x,4e16.9)
1003 format(4e16.9)
1004 format(1x,16i4)
1005 format(16i4)

    xheneg = 0
    icmp = 0
    small = 1.e-4
    nl=nshint
    xc2 = 1.
    xo2 = 0.

    !  call env for the converged model
35  call env

    call check_time(newtime)
    if ((newtime - t0) .gt. xmaxtime) goto 6
    xheneg = 0.0
    ifin = 0
    nel=nlim

    if( iprep .eq. 0 )then
       !*** no reason to continue here ***
       goto 6
    endif
    nel = nel - 1
    n = 0
    do   
34     if (n .ge. nel) exit
       n = n + 1
       index = nel-n+1
       !*** these are placed in reverse order ***
       rrr = r(index)
       amm = mm(index)
       ttt = t(index)
       rho = rl1(index)
       ppp = p(index)

       aa(1,nl+n)=r(index)

       x(1) = xcomp(index,1)   !Hydrogen
       x(2) = xcomp(index,2)   !Helium
       x(3) = xcomp(index,3)   !Carbon
       x(4) = xcomp(index,4)   !Oxygen
       xcomp_orig(nl+n,:) = xcomp(index-1,:)

       if ( n .eq. nel ) then
          x(1) = 0.
          x(2) = 0.
          x(3) = 1.
          x(4) = 0.
       endif

       !*** mass (grams) ***
       aa(2,nl+n)=10.**mm(index)*amsun

       !*** temperature (kelvins) ***
       aa(4,nl+n)=10.**t(index)

       !*** density (g/cm^3) ***
       aa(5,nl+n)=10.**rl1(index)

       !*** pressure (dynes/cm^2) ***
       aa(6,nl+n)=10.**p(index)

       !*** energy generation rate + derivatives (eps, epsr, and epst) ***
       aa(7,nl+n)=0.0
       aa(11,nl+n)=0.0
       aa(12,nl+n)=0.0

       ! use delrad as temp. gradient except in convection zone, then use
       ! deltrue from gradt (or gradt2).

       if(rtg(index).lt.atg(index))then
          tdel=rtg(index)
       else
          tdel=ttg(index)
       endif

       !*** temperature gradient (del) ***
       aa(15,nl+n)=tdel

       !*** adiabatic temperature gradient (delad) ***
       aa(16,nl+n)=atg(index)

       !*** helium profile ***
!       print *, n, xcomp(n,2)
!       aa(17,n+nl) = xcomp(n,2)
!       aa(17,index+n) = xcomp(index+n,2)
       aa(17,nl+n)=xcomp(index,2)
!       do i=1,450
!          print *, i, xcomp(i,2)
!       end do

       !*** radiative luminosity profile ***
       aa(3,nl+n)=alsun*(10.0**(blfin))*(ttg(index)/rtg(index))

       !*** heat capacity Cv ***
       aa(8,nl+n)=ecv(index)

       !*** pressure gradients chi rho and chi T ***
       aa(9,nl+n)=exr(index)
       aa(10,nl+n)=ext(index)

       !*** opacity ***
       tl=t(index)
       pl=p(index)
       rhol = rl1(index)

       ifinal = 1

       if (.not. kap_new ) then

          oplo = opac3(rhol,tl)
          aa(18,nl+n)=10.**oplo

          !*** start 5-point opacity derivative stuff here. ***
          if(x(1).eq.1.0 .or. x(2).eq.1.0 .or. x(3).eq.1.0)then
             if(x(1).eq.1.0)then
                je=1
             endif
             if(x(2).eq.1.0)then
                je=2
             endif
             if(x(3).eq.1.0)then
                je=3
             endif
             call opacder(rhol,tl,opcs,dlklr,dlklt)
             aa(13,nl+n)=dlklr
             aa(14,nl+n)=dlklt
             !*** note: replace previous opacity value if we can get it here. ***
             aa(18,nl+n)=10.**opcs
          else
             !*** 5 point differencing ***
             oplo5=opac3(rhol+0.02,tl)
             oplo1=opac3(rhol+0.01,tl)
             oplo2=opac3(rhol-0.01,tl) 
             oplo6=opac3(rhol-0.02,tl)
             oplo7=opac3(rho,tl+0.02)
             oplo3=opac3(rhol,tl+0.01)
             oplo4=opac3(rhol,tl-0.01)
             oplo8=opac3(rho,tl-0.02)
             !*** opacitry derivatives (kapr(13) and kapt(14)) ***
             aa(13,nl+n)=(oplo6-(8.*oplo2)+(8.*oplo1)-oplo5)/(12.*0.01)
             aa(14,nl+n)=(oplo8-(8.*oplo4)+(8.*oplo3)-oplo7)/(12.*0.01)
          endif
       else
          xavec(:) = 0.0
          xavec(1) = x(1)
          xavec(2) = x(2)
          xavec(3) = x(3)
          call opalz(10**tl,10**rhol, xavec, fkappa, opr, opt)
          !write(*,'(2f8.4,6es13.5,4f8.3)') rhol, tl, aa(18,nl+n), fkappa, aa(13,nl+n), opr, aa(14,nl+n), opt
          !write(*,'(2f8.4,6es13.5,4f8.3)') rhol, tl, aa(18,nl+n)/fkappa, aa(13,nl+n)/opr, aa(14,nl+n)/opt
          aa(18,nl+n)=fkappa
          aa(13,nl+n)=opr
          aa(14,nl+n)=opt
       endif

       ! Check that Xhe is positive. If not, go back and try again. This time 
       ! around, will use a fail safe built into ccomp.
       if ( (aa(17,nl+n).lt.0) .or. (aa(17,nl+n-1).lt.0) ) then
          xheneg=1
          goto 35
       endif
       !*** if He fraction is constant, set Ledoux to zero. ***
       if(aa(17,nl+n).eq.aa(17,nl+n-1))then
          aa(19,nl+n) = 0.0
          goto 34
       endif
       if(aa(17,nl+n-1).eq.0.0 .or. aa(17,nl+n).eq.0.0)then
          aa(19,nl+n) = 0.0
          goto 34
       endif
       if(index.eq.nel)then
          aa(19,nl+n) = 0.0
       else
          !*** call up modified ledoux term subroutine ***
          call ledoux
          dly = ( dlog(aa(17,nl+n)) - dlog(aa(17,nl+n-1)) )
          dlp = ( dlog(aa(6,nl+n)) - dlog(aa(6,nl+n-1)) )
          bledoux=(1./aa(10,nl+n))*bled*(dly/dlp)
          if (ihotspot .eq. 1) goto 6

! Note: there is a minus sign difference between this ledoux term 
! in the core and the one in Brassard et al. (1991, ApJ, 367, 601).
! This was "corrected" by the absolute value below, but I want to be
! able to see if certain configurations are convectively unstable in the
! core, so I am going to remove the dabs() (MHM April 1998).
          aa(19,nl+n)=bledoux
          if(x(1).ne.0.0 .and. x(2).gt.0.99 .and. &
               x(3).ne.0.0 .and. bledoux.lt.0.) then
             aa(19,nl+n)= 0.d0
          endif
       endif
    enddo

    ifinal = 0

    ! now write it out. 

    nshell=nl+nel

    if (evoloutput) then
       if (tfit(ifit).lt.1.d-3) then
          write(50,1030) nshell, nl, nel, amhyhe,amheca, &
               alph(1),alph(2),alph(3),alph(4),stpms
       end if
    endif

1010 format(i5,'  =  ',i4,'  +  ',i4)
1030 format(i5,'  =  ',i4,'  +  ',i4,2(2x,1pe9.3),0p,4(f7.4,x) 1x,1pe10.3)

    ! smooth out del and delad at core-envelope boundary
    ! set mver for appropriate variable
    ! 1. delad; 2. del; 3. chiT; 4. chiRho

    do i=1,nshell
       if(i.eq.nl)then
          !*** delad ***
          diff = aa(16,i+1) - aa(16,i)
          if(diff.gt.0.002)then
             mver = 1
             call smooth(nl,mver)
             if (verbose) write(*,*) 'smoothing delad'
          endif
          !*** del ***
          diff2 = aa(15,i+1) - aa(15,i)
          if(diff2.gt.0.002)then
             mver = 2
             call smooth(nl,mver)
             if (verbose) write(*,*) 'smoothing del'
          endif
          !*** chi T ***
          diff3 = aa(10,i+1) - aa(10,i)
          if(diff3.gt.0.002)then
             mver = 3
             call smooth(nl,mver)
             if (verbose) write(*,*) 'smoothing chi T'
          endif
          !*** chi rho ***
          diff2 = aa(9,i+1) - aa(9,i)
          if(diff4.gt.0.002)then
             mver = 4
             call smooth(nl,mver)
             if (verbose) write(*,*) 'smoothing chi rho'
          endif
       endif
    enddo

    ! Write out model to tape50

    if (evoloutput) then
       if (tfit(ifit).lt.1.d-3) then
          do i=1,20
             write(50,1020) (aa(i,n), n=1,nshell)
          enddo
       endif
    end if

    ! Write out surface logg to screen
    if (tfit(ifit).lt.1.d-3) then
       bigg = 6.67259d-8
       alogg = log10(bigg) + log10(aa(2,nshell)) - 2.*log10(aa(1,nshell))
    endif

1020 format(1p,4e22.15)

    if (evoloutput) then
       if (tfit(ifit).lt.1.d-3) then
          close(50)
       endif
    end if
       
6   return
  end subroutine eprep

  !****************************************************************************

  subroutine env

!     calculates white dwarf envelopes

!     the variables are r the radius,   m the log of the mass within a
!     sphere of radius r,   p the log of the pressure,   t the log of the
!     temperature,   el1 the degeneracy parameter,   ol1 the log of the
!     total opacity,   rl1 the log of the density,   atg the adiabatic
!     gradient,   rtg the radiative gradient,   phs the mixing length,
!     drdtp the derivative of the density(log) with respect to the
!     temperature(log) at constant pressure,   cp the specific heat per
!     gram at constant pressure,   ttg the actual temperature gradient.

!     the control parameters are:
!     nmax,  maximum number of shells
!     nsuff,  indicates if max number of shells is exceeded 
!     nlim1,  shell number at the bottom of the convection zone
!     kind,  shell number at the top of the convection zone 
!     nlim,  total number of shells for a model
!     1/nhp,  fraction of appropriate scale height used as integration step
!     iter,  total number of iterations 
!     nlum,  number of luminosities
!     lmax,  number of isotherms in eos table
!     ml,  number of density points per isotherm (eos)
!     lmam,  number of isotherms for which optical depths are computed

!     the units are cgs throughout
! Subroutines
    use utils_subroutines, only : check_time
    use chemprofiles_subroutines, only : comp_env 

! Common blocks
    use misc_const
    use shells2
    use tenv
    use threshh
    use jee
    use comp
    use xcompp
    use slr
    use cc
    use d 
    use grad
    use surf
    use terp
    use wprop
    use writeitt
    use oldval
    use fin
    use dfin
    use wprep
    use dprep
    use mixl
    use crash
    use contrl
    use flags, only: verbose, evoloutput, kap_new
    use kappa, only : opalz

    implicit double precision(a-h,o-z)

    real*8, dimension(2,3) :: w
    real*8, dimension(7) :: xavec
    integer :: newtime
    real*8, parameter :: dlummax = 0.1d0, dtemax = 0.5d0

    stop=stpms+sm-dlog10(amsun)
    smm=10.**(sm-dlog10(amsun))
    gslast=1e8
    xmass = smm

    if ( ifin .eq. 0 ) then
       rs=10.**(ww(2,is)-10.8425968)
       sl(1)=10.**ww(1,is)
       writeit = .false.
    elseif ( ifin .eq. 1 ) then
       rs = 10.**(rfin-10.8425968)
       sl(1) = 10.**blfin
    endif

    !     compute the surface gravity gs and scale the optical depth tables
    !     accordingly

    gs=27398.*smm/(rs*rs)
    glog=dlog10(gslast/gs)
    do jm=1,3 
       do l=1,lmam
          kmax=ml(jm,l)
          do kk=1,kmax
             od(jm,l,kk)=od(jm,l,kk)+glog
          enddo
       enddo
    enddo

!  set chem comp for gray atm integration
!  specify a luminosity and compute the effective temperature te
    call comp_env
    jj = 1
    nsuff=1
    nlim1=0
! am=(x(1)*az(1) + x(2)*az(2) * x(3)*az(3)) * 1.660531e-24
    am=(x(1)*az(1) + x(2)*az(2) + x(3)*az(3)) * 1.660531e-24
    te=5800.*(sl(jj)/(rs*rs))**0.25

!  if stepped far enough in h-r diagram, then write out envelope to tape6

    if ( ifin .eq. 1 ) then
       dlum = dabs( bold - bm ) 
       dteff = dabs( teold - dlog10(te) ) 
       if ( dlum .gt. dlummax .or. dteff .gt. dtemax ) then
          bold = bm
          teold = dlog10( te ) 
          writeit = .true.
          if (verbose) then
             write(*,*) 'writing envelope to tape6.  teff =',te
             write(*,*) '--------------------------------------------'
          end if
       endif
    endif

    if ( evoloutput .and. writeit ) write(26,103) amhyhe,amheca
103 format(1x,' amhyhe=',1pe14.4,'  amheca=',1pe14.4)

!  integrate gray atmosphere
    call gray       
    call check_time(newtime)
    if ((newtime - t0) .gt. xmaxtime) goto 16
    it=2
    nhp=8
    iter=1

!     start the inward integration with values provided by the gray
!     subroutine,  if there is convection take the geometrical depth
!     as zero

    p(1)=ps
    t(1)=ts
    r(1)=rs*6.96e+10
    mm(1)=dlog10(smm)
    dept=1e-20
    deptx=dept
    kind1=0

    do n=1,nmax 
       ko=n
       tl=t(n)
       pl=p(n)
       xmass=10.**mm(n)
       call comp_env
       ! am=(x(1)*az(1) + x(2)*az(2) * x(3)*az(3)) * 1.660531e-24
       am=(x(1)*az(1) + x(2)*az(2) + x(3)*az(3)) * 1.660531e-24
       do iii = 1,4
          xcomp(n,iii) = x(iii)
       enddo

!  check that our (p,t) point is on the eeos table
!    (inserted 6/13/88)
!  skip over check if teff > thresh (30,000)

       if (te .gt. thresh) then
          goto 987
       end if
       itchk = (tl - tk(je,1)) / (tk(je,2) - tk(je,1)) + 1
       if ( tl .lt. tk(je,1) .or. tl .gt. tk(je,lmax) ) then
          if ( evoloutput ) then
             write(26,*) 'subroutine env' 
             write(26,321) tl, tk(je,1), tk(je,lmax), je
321          format(' off eeos.  t:',f5.2,' tk(1), tk(lmax): ',2f6.2, &
                  '    je: ',i1)
          end if
          if ( tl .lt. tk(je,1) ) then
             itchk = 1
          endif
          if ( tl .gt. tk(je,lmax) ) then
             itchk = lmax
          endif
       endif

       mltmp = ml(je,itchk) 

       !  end eeos check

       if ( evoloutput ) then
          if ( pl .lt. pp(je,itchk,1) .or. pl .gt. pp(je,itchk,mltmp ) ) then
             write(26,*) 'subroutine env' 
             write(26,432) pl, tl, je, itchk, pp(je,itchk,1), &
                  je, itchk, mltmp, pp(je,itchk,mltmp) 
432          format(' off eeos. '/'  p:',f7.3,'  t:',f7.3,/ &
                  '  pp(',i1,',',i2.2,',1): ',f7.3,/ &
                  '  pp(',i1,',',i2.2,',',i2.2,'):',f7.3) 
          endif
       end if

987    continue

       !     given p and t at r call the subroutine state and store the values
       !     of el1,rl1,ol1,atg,rtg,drdtp,cp
       call estatem
       call check_time(newtime)
       if ((newtime - t0) .gt. xmaxtime) goto 16

       if (.not. kap_new ) then
          ol1x = opac3(rl1x,tl)
       else
          xavec(:) = 0.0
          xavec(1) = x(1)
          xavec(2) = x(2)
          xavec(3) = x(3)
          call opalz(10**tl,10**rl1x, xavec, fkappa, opr, opt)
          ol2x = dlog10(fkappa)
          ! write(*,'(2f8.4,6es13.5,4f8.3)') rl1x, tl, ol1x, ol2x, xavec(1:3)
          ol1x = ol2x 
       endif

       rtgx = 7.732749e+9*10.**(ol1x+pl-4.*tl)*sl(jj)/xmass
       ol1(n)=ol1x
       rtg(n)=rtgx
       el1(n)=el1x
       rl1(n)=rl1x
       atg(n)=atgx
       drdtp(n)=drdtpx
       cp(n)=cpx
       gam1=((xt**2)*10.e0**(p(n)-rl1x-t(n))/cv)+xr
       bn2(n)=-drdtpx*(rtgx-atgx)*((g*amsun*10.e0**mm(n)/r(n)**2)**2)*&
            10.e0**(rl1x-p(n)) 
       ac2(n)=gam1*10.e0**(p(n)-rl1x)/r(n)**2
       dpdr=-27398.*(10.**mm(n))*rho*(6.96e+10/r(n))**2
       dmdr=4.*pi*rho*r(n)**2
       if ( ifin .eq. 1 ) then
          ecv(n) = cv
          ext(n) = xt
          exr(n) = xr
       endif
       !   define the integration step as a fraction of the pressure scale 
       !    height

       dr=10.**p(n)/nhp/dpdr

       !     Schwarzschild's test, if radiative equilibrium take rtg
       !     as the actual gradient, if convective define the mixing length
       !     as the smaller of either the local pressure scale height or the 
       !     geometrical depth from the top of the convection zone and call
       !     gradt to evaluate the true gradient

       if ( atg(n) .le. rtg(n))  then 
          nlim1=n 
          kind=n-kind1
          kind1=kind1+1
          phs(n)=dabs(dr*nhp)
          dept=dept+dabs(dr) 
          deptx=dept
          pecx=9./8./(1.601743e+07*cpx*10.**(ol1x+2.5*rl1x+mm(n) &
               -3.*t(n)-0.5*p(n))*(dabs(drdtpx))**0.5 &
               *(phs(n)**2)*(6.96e+10/r(n))**2)
          pecx = pecx/aml/aml
          if(mls.eq.0)then
             !*** call Bohm-Vitense version ***
             call gradt
          else
! Some scaling to get the right convection zone depth
             pecx=(8./9.)*3./4.242440687*pecx  ! original, mysterious scaling
             pecx = pecx_factor * pecx ! new scaling
             pecx = pecx * (8./9.)/dsqrt(2.d0)
             call gradt2
          endif
          ttg(n)=ttgx       
       else
          ttg(n)=rtg(n)
          kind1=0 
          dept=0. 
       endif

       dtdr=ttg(n)*dpdr*10.**(t(n)-p(n))

!     evaluate r,p,t,m at the next shell and quit if the envelope
!     mass fraction exceeds -stpms or the temperature exceeds tk(lmax)
       call intg
       call check_time(newtime)
       if ((newtime - t0) .gt. xmaxtime) goto 16
       je=3
       if ( t(n+1) .gt. tk(je,lmax) ) then
          if ( evoloutput ) then
             write (26,*)' off eos. t(n+1) gt tk(lmax)', &
                  t(n+1),tk(je,lmax),n
          end if
          goto 13
       elseif ( stop .gt. mm(n) ) then
          if ( evoloutput ) then
             write (26,*)' stop gt m(n).  n=',n, ' n -> n-1 '
          end if
          goto 13
       elseif ( n+2 .gt. nmax ) then
          if ( evoloutput ) then
             write (26,*)' too many shells.', n+2
          end if
          goto 11
       endif
       nlim=n+1
    enddo

    ! end of large loop (envelope integration) at statement 3
    ! if we go to 11, then we ran out of memory, and didn't integrate
    ! all the way to the fitting point.
11  nsuff=-1    
13  plim2=10.**p(nlim)
    plim1=plim2
    acu=0.
    n = n-1
    call write4

    if (ihotspot .eq. 1) then
       goto 16
    end if

    if ( writeit ) writeit = .false.

    gslast=gs
16  return
  end subroutine env

!****************************************************************************
!  func1 (for diffusion exponents) Top of H/He transition zone.
!****************************************************************************   
  function func1(x)

    use alpha

    implicit double precision(a-h,o-z)

    x5 = x**alph(1)
    result = ( 2. - x5 ) / ( 2. + 3. * x5 )
    func1 = result

    return

  end function func1

!****************************************************************************
!  func2 (for diffusion exponents) Bottom of H/He transition zone.
!****************************************************************************

  function func2(x)

    use alpha

    implicit double precision(a-h,o-z)

    xx = x**alph(2)
    result =  xx / ( 8. - 3. * xx ) 
    func2 = result

    return

  end function func2

!****************************************************************************
!  func3 (for diffusion exponents) Top of He/C transition zone.
!****************************************************************************

  function func3(x)

    use alpha

    implicit double precision(a-h,o-z)

    x2 = x**alph(3)
    result = ( 2. - x2 ) / ( 2. + 2. * x2 )
    func3 = result

    return

  end function func3

!*****************************************************************************
!  func4 (for diffusion exponents) Bottom of He/C transition zone.
!*****************************************************************************

  function func4(x)

    use alpha

    implicit double precision(a-h,o-z)

    xx = x**alph(4)
    result = xx / ( 6. - 2. * xx )
    func4 = result

    return

  end function func4

!****************************************************************************

  subroutine gray

    !     gray atmosphere stratification,   the pressure(log)  is patm, the
    !     temperature(log) is tatm and the optical depth(log) is tau

    ! Subroutines
    use utils_subroutines, only : check_time

    ! Common blocks
    use misc_const
    use cc
    use vca
    use dvca
    use tenv
    use jee
    use comp
    use slr
    use d 
    use writeitt
    use oldval
    use threshh
    use crash
    use flags, only: verbose, evoloutput, kap_new
    use kappa, only : opalz

    implicit double precision(a-h,o-z)

    real*8, dimension(7) :: xavec
    integer :: newtime
    integer, pointer :: j

    j => jj

139 format(/'-------------------------------------------------------')
140 format(/'Model ',i3,'  gray atmosphere structure (wdxd)') 
142 format(/'       composition (h he c): ',3(1x,f6.3))
150 format(/8h m/msun=,f6.3,2x,8h r/rsun=,1pe10.3,2x,6h grav=,e11.4, &
         2 x,6h teff=,e11.4,2x,8h l/lsun=,e10.3/)
141 format(7x,4hlogt,10x,6hlogrho,10x,6hlogtau,9x,4hlogp,12x,3hatg, &
         12x,3hrtg,9x,8hlogkappa,9x,3heta)
200 format(8(3x,e12.5))
101 format(/26h still radiative when tau=,f8.3) 
102 format(/27h convection sets in at tau=,f6.3)
201 format(' new envelope.  gray atm still radiative at tau =',f8.3, &
         '  teff=',f7.0)
202 format(' new envelope.  gray atm convective at tau =',f6.3, &
         '  teff=',f7.0)
152 format(F10.0,2F8.3)

    if ( writeit .and. evoloutput ) then
       write(26,139)
       write(17,140) modnr
       write(26,140) modnr
       write(26,142) x(1), x(2), x(3)
       write(17,150) smm,rs,gs,te,sl(j)
155    format(F10.0,2F10.3)
!       write(*,155) te, smm, dlog10(gs)
       write(26,150) smm,rs,gs,te,sl(j)
       write(315,152) te,smm,dlog10(gs)
       write(26,141)
    end if

    !    interpolate linearly in tk-od table to find the pressure
    !     corresponding to initial tau and tatm

    tau=-4.
    tatm=dlog10(te*(0.75*0.5773)**0.25)
    patm = dlog10( gs ) + tau + 1.
98  continue
    patm1 = patm

    !  integrate equation of hydrostatic equilibrium using log pressure
    !  as independant variable,  quit when convection sets in or when tau=1
    !  if radiative equilibrium

    !  integration loops back to statement 34 after incrementing pressure 

    xmass=smm
34  tl=tatm
    pl=patm
    !*** trying out stop at tau = 0.1 instead of tau = 10.0 ***
    if ( tau .ge. -1. ) then 
       tauz=10.**tau
       if ( evoloutput ) then
          if ( writeit ) write(26,101) tauz
          write(17,201) tauz, te
          write(20,201) tauz, te
          write(26,201) tauz, te
       end if
       if (verbose) print 201, tauz, te
       ps=patm
       ts=tatm
       goto 306
    endif

    !  check that our (p,t) point is on the eeos table
    !    (inserted 6/13/88)
    !  skip over check if teff > thresh (30,000)

    if (te .gt. thresh) goto 987
    itchk = (tl - tk(je,1)) / (tk(je,2) - tk(je,1)) + 1
    if ( tl .lt. tk(je,1) .or. tl .gt. tk(je,lmax) ) then
       if ( evoloutput ) then
          write(26,*) 'subroutine gray:'
          write(26,321) tl, tk(je,1), tk(je,lmax), je
          write(20,*) 'subroutine gray:' 
          write(20,321) tl, tk(je,1), tk(je,lmax), je
321       format(' off eeos.  t:',f5.2,' tk(je,1), tk(je,lmax): ',2f6.2,' &
               je: ',i1)
       end if
       if ( tl .lt. tk(je,1) ) itchk = 1 
       if ( tl .gt. tk(je,lmax) ) itchk = lmax
    endif

    mltmp = ml(je,itchk)

    if ( evoloutput ) then
       if ( pl .lt. pp(je,itchk,1) .or. &
            pl .gt. pp(je,itchk,mltmp ) ) then
          write(26,*) 'subroutine gray:'
          write(26,432) pl, tl, je, itchk, pp(je,itchk,1), &
               je, itchk, mltmp, pp(je,itchk,mltmp) 
          write(20,*) 'subroutine gray:' 
          write(20,432) pl, tl, je, itchk, pp(je,itchk,1), &
               je, itchk, mltmp, pp(je,itchk,mltmp) 
432       format(' off eeos. '/'  p:',f7.3,'  t:',f7.3,/ &
               '  pp(',i1,',',i2.2,',1): ',f7.3,/ &
               '  pp(',i1,',',i2.2,',',i2.2,'):',f7.3) 
       endif
    end if

    !  end eeos check

987 continue

    !  ccomp called before gray, no need to call here 

    call estatem
    call check_time(newtime)
    if ((newtime - t0) .gt. xmaxtime) goto 306
    if (.not. kap_new ) then
       ol1x = opac3(rl1x,tl)
    else
       xavec(:) = 0.0
       xavec(1) = x(1)
       xavec(2) = x(2)
       xavec(3) = x(3)
       call opalz(10**tl,10**rl1x, xavec, fkappa, opr, opt)
       ol2x = dlog10(fkappa)
       !write(*,'(2f8.4,6es13.5,4f8.3)') rl1x, tl, ol1x, ol2x, xavec(1:3)
       ol1x = ol2x 
    endif

    rtgx = 7.732749e+9*10.**(ol1x+pl-4.*tl)*sl(j)/xmass
    if ( evoloutput .and. writeit ) then
       write(26,200) tl,rl1x,tau,pl,atgx,rtgx,ol1x,el1x
    end if
    diff=rtgx-atgx
    if ( atgx .le. 0.0 ) diff = rtgx - 0.40

    ! return if convective
543 format(//'first shell of gray atm convective'/, &
         'backing off patm trying again: ',f8.3)
    if(diff .gt. 0.) then
       if ( patm .eq. patm1 ) then
          if ( evoloutput ) then
             write(17,543)
             write(20,543)
          end if
          if (verbose) print 543, patm-.3
          patm = patm - .3
          goto 98
       endif
       taux=10.**t1+diff1*(10.**tau-10.**t1)/(diff1-diff) 
       ps=patm
       ts=tatm
       if ( verbose ) write(*,202) taux, te
       if ( evoloutput ) then
          if ( writeit ) write(26,102) taux
          write(17,202) taux, te
          write(20,202) taux, te
       end if
       goto 306
    endif

    ! otherwise continue with the atm integration

    diff1=diff
    t1=tau
    dtaudp=(2.30258*10.**(pl+ol1x))/gs
    tau1=10.**tau+0.05*dtaudp
    tau=dlog10(tau1)
    qtau=0.7104-0.1331*dexp(-3.4488d0*tau1)
    tatm=dlog10(te*((tau1+qtau)*0.75)**0.25)
    patm=patm+0.05
    goto 34

    !  goto 34 ==> loop back to top and continue integration
306 return

  end subroutine gray

  !************************************************************************

  subroutine estatem

    !  compute thermodynamic quantities for the compsition just specified 
    !  by ccomp.  use the additive volume technique (see istatco for ref).
    !  see Fontaine, Graboske & Van Horn 1977, ApJS, 35, 293

    ! Subroutines
    use utils_subroutines, only : check_time
    use wd_eos_mod, only : wd_eos, crad

    ! Common blocks
    use misc_const
    use jee
    use comp
    use cc
    use grad
    use crash
    use flags, only : ios_new

    implicit double precision(a-h,o-z)

    double precision prad, pgas
    real*8, dimension(11) :: res
    real*8, dimension(7) :: xavec

    integer :: newtime


    call check_time(newtime)
    if ((newtime - t0) .gt. xmaxtime) goto 406
    
!  initialize summation variables
    
    rhotmp  = 0.0
    etatmp  = 0.0
    cptmp   = 0.0
    cpdatmp = 0.0
    deltmp  = 0.0
    dpttmp  = 0.0
    dtptmp  = 0.0
    
!  call estate for each composition i for which x(i) ne 0.

    do je = 1,3
       if ( x(je) .gt. 0. ) then
          xx = x(je)
          xje = xx
          call estate
          rhotmp = rhotmp + xje/ rho
          etatmp = etatmp + xje * el1x
          cptmp  = cptmp + xje * cpx 
          cpdatmp = cpdatmp + xje * cpx * atgx 
          dpttmp = dpttmp + xje / ( rho * drdptx )
          dtptmp = dtptmp + xje / ( rho * drdtpx )
       endif
    enddo
    
    if (ios_new) then
       xavec(:) = 0.d0
       xavec(1) = x(1)
       xavec(2) = x(2)
       xavec(3) = x(3)
       xavec(5) = x(4)
! radiation pressure is neglibible by the ZZ Ceti I.S.
       prad = (crad * 10**(4.*tl))/3
       pgas = 10**pl - prad   ! can be negative for hot models ~ 100,000 K
       call wd_eos(10**pl,10**tl,xavec,'p',res)
!       call wd_eos(pgas,10**tl,xavec,'p',res)
       
       rho    = res(1)
       rl1x   = dlog10(rho)
       el1x   = res(9)
       cpx    = res(2)
       atgx   = res(4)
       xr     = res(5)
       xt     = res(6)
       cv     = res(3)
       drdtpx = xt/xr
       drdptx = -drdtpx/xt

       !write(*,'(30e13.5)') pl,dlog10(pgas),dlog10(prad),tl,rho,rho2,cpx,cpx2,atgx,atgx2, xr,xr2,xt,xt2,cv,cv2,drdtpx,drdtpx2,drdptx,drdptx2
       !write(*,'(20e13.5)') pl,tl,rho,rho/rho2,cpx/cpx2,atgx/atgx2,xr/xr2,xt/xt2,cv/cv2,drdtpx/drdtpx2,drdptx/drdptx2
          
    else
       rho    = 1/rhotmp
       rl1x   = dlog10(rho)
       el1x   = etatmp
       cpx    = cptmp
       atgx   = cpdatmp / cpx
       drdptx = 1. / (rho * dpttmp)
       drdtpx = 1. / (rho * dtptmp)
       xr     = -1./ drdptx
       xt     = -drdtpx / drdptx
       cv     = cpx * ( 1. - xt * atgx )
    endif
       
406 return

  end subroutine estatem

 !****************************************************************************

  subroutine estate

    !     3-point lagragian interpolation in the pp-tk grid to find thermo
    !     quantities el1,rl1,ol1,atg,rtg,drdtp,cp [PROFILE: 23% runtime]

    ! Subroutines
    use utils_subroutines, only : check_time

    ! Common blocks
    use misc_const
    use tenv
    use jee
    use comp
    use cc
    use d 
    use grad
    use crash

    implicit double precision(a-h,o-z)

    real*8, dimension(3) :: xtl,xrl,el,rl,at
    real*8, pointer :: sm
    integer, dimension(3) :: k1,k2,k3
    integer :: newtime
    integer, pointer :: j

    sm => smm
    j => jj

    call check_time(newtime)
    if ((newtime - t0) .gt. xmaxtime) goto 506

    ! trap point in t, after checking if = low boundry or less, high
    ! boundry or more

    if ( tl .lt. tk(je,2) ) then
       l1=1
       l2=2
       l3=3
    else if ( tl .gt. tk(je,lmax-1) ) then
       l3=lmax
       l2=lmax-1
       l1=lmax-2
    else

! The last 6 isotherms are at a different
! spacing in log T of 0.25  in the current incarnation of the
! 'envelope' EOS extended to high temperatures.  

       if (tl.lt.7.54d0) then
          l3=int((tl-tk(je,1))/(tk(je,2)-tk(je,1)))+2
          l2=l3-1
          l1=l3-2
       else
          l3=int((tl-tk(je,54))/(tk(je,55)-tk(je,54)))+56
          l2=l3-1
          l1=l3-2
       endif
    endif

! get fractional distance between mesh points in temp lagrange coefficients.

    dtdt=(tl-tk(je,l2))/(tk(je,l3)-tk(je,l2))
    dt1=tk(je,l1)-tl 
    dt2=tk(je,l2)-tl 
    dt3=tk(je,l3)-tl 
    dtdt1=dt2*dt3/((dt1-dt2)*(dt1-dt3))
    dtdt2=dt1*dt3/((dt2-dt1)*(dt2-dt3))
    dtdt3=dt1*dt2/((dt3-dt1)*(dt3-dt2))

    ! Now trap point in p
    ! Note, we come loop back to 8 for l = l1, l2, l3 
    ! This is a *slow* (linear) search thru the grid in pressure
    ! and I wouldn't be surprised if the program spends a large 
    ! fraction of its time right here

    l=l1
8   kmax=ml(je,l) 
    lp=l+1-l1
    if ( pl .lt. pp(je,l,2) ) then
       k1(lp)=1
       k2(lp)=2
       k3(lp)=3
    else
       ! This bit of old style code fails in F90 because the goto 16 fails to
       ! make the do loop iterate the k variable.
!!$16     do k=3,kmax
!!$          if ( pl .gt. pp(je,l,k) ) then
!!$             goto 16
!!$          else
!!$             k3(lp)=k
!!$             k2(lp)=k-1
!!$             k1(lp)=k-2
!!$             go to 14
!!$          endif
!!$       enddo

       ! Replace with the following clumsy code
       k = 3

16     if (k .gt. kmax) goto 15
       if ( pl .gt. pp(je,l,k) ) then
          k = k + 1
          goto 16
       else
          k3(lp)=k
          k2(lp)=k-1
          k1(lp)=k-2
          go to 14
       end if
15     continue

       k3(lp)=kmax
       k2(lp)=kmax-1
       k1(lp)=kmax-2
    endif

14  kl1=k1(lp)
    kl2=k2(lp)
    kl3=k3(lp)
    dp1=pp(je,l,kl1)-pl
    dp2=pp(je,l,kl2)-pl
    dp3=pp(je,l,kl3)-pl
    dpdp1=(pl-pp(je,l,kl2))*(pl-pp(je,l,kl3))/(pp(je,l,kl1)- &
         pp(je,l,kl2))/(pp(je,l,kl1)-pp(je,l,kl3))
    dpdp2=(pl-pp(je,l,kl1))*(pl-pp(je,l,kl3))/(pp(je,l,kl2)- &
         pp(je,l,kl1))/(pp(je,l,kl2)-pp(je,l,kl3))
    dpdp3=(pl-pp(je,l,kl1))*(pl-pp(je,l,kl2))/(pp(je,l,kl3)- &
         pp(je,l,kl2))/(pp(je,l,kl3)-pp(je,l,kl1))

! Interpolate in P for each T point

    el(lp)=eta(je,l,kl1)*dpdp1+eta(je,l,kl2)*dpdp2+eta(je,l,kl3)*dpdp3
    at(lp)=datg(je,l,kl1)*dpdp1+datg(je,l,kl2)*dpdp2+datg(je,l,kl3)*dpdp3
    rl(lp)=rhok(je,l,kl1)*dpdp1+rhok(je,l,kl2)*dpdp2+rhok(je,l,kl3)*dpdp3
    xtl(lp)=xtt(je,l,kl1)*dpdp1+xtt(je,l,kl2)*dpdp2+xtt(je,l,kl3)*dpdp3
    xrl(lp)=xrt(je,l,kl1)*dpdp1+xrt(je,l,kl2)*dpdp2+xrt(je,l,kl3)*dpdp3

    if ( lp .ge. 3 ) goto 19
    l=l+1
    go to 8

! this is the bottom of the loop
! Now interpolate in T for actual pressure.

19  el1x=el(2)+0.5d0*dtdt*(el(3)-el(1))+dtdt*dtdt*(0.5d0*(el(1)+el(3))-el(2))

    atgx=at(2)+0.5d0*dtdt*(at(3)-at(1))+dtdt*dtdt*(0.5d0*(at(1)+at(3))-at(2))
    rl1x=rl(2)+0.5d0*dtdt*(rl(3)-rl(1))+dtdt*dtdt*(0.5d0*(rl(1)+rl(3))-rl(2))
    xt1x=xtl(2)+0.5d0*dtdt*(xtl(3)-xtl(1))+dtdt*dtdt*(0.5d0*(xtl(1)+xtl(3)) &
         -xtl(2))
    xr1x=xrl(2)+0.5d0*dtdt*(xrl(3)-xrl(1))+dtdt*dtdt*(0.5d0*(xrl(1)+xrl(3)) &
         -xrl(2))

! compute final thermo properties

    xt = xt1x
    xr = xr1x
    drdtpx = xt/xr
    drdptx = -1.0d0/xr
    rho  = 10.d0**rl1x
    cpx  = drdtpx/rho/atgx*10.**(pl-tl)
    cv   = cpx*(1.0d0-(xt*atgx))

506 return

  end subroutine estate

  !************************************************************************

  subroutine gradt

    !     solution of the cubic equation characteristic of the mixing length
    !     theory  (u is essentially the ratio of radiative to convective
    !     conductivities)  Bohm-Vitense version  (ML1)

    use misc_const
    use cc
    use grad

    implicit double precision(a-h,o-z)

    u=pecx
    if ( u .lt. 1.e-10 ) then
       ttgx=atgx
       goto 606
    elseif ( u .gt. 1.e+10 ) then
       ttgx=rtgx
       goto 606
    endif

    q=368./729.*u**2
    r=4672./19683.*u**3+4./9.*u*(rtgx-atgx)
    s1=(r+(q**3+r**2)**0.5)**0.3333333
    s2=-q/s1

    if ( u .le. 1. ) then
       super=(s1+s2+19./27.*u)**2-u**2
       if ( super .le. 0. ) super=0.
       ttgx=atgx+super
       goto 606
    else
       super=9.*(s1+s2-8./27.*u)**3/u/8.
       if ( super .le. 0. ) super=0.
       ttgx=rtgx-super
       goto 606
    endif
606 return

  end subroutine gradt

  !************************************************************************

  subroutine gradt2

    !     solution of the cubic equation characteristic of the mixing length
    !     theory  (u is essentially the ratio of radiative to convective
    !     conductivities)  Bohm & Cassinelli version (ML2 or ML3)

    use misc_const
    use cc
    use grad

    implicit double precision(a-h,o-z)

    u=pecx
    if ( u .lt. 1.e-10 ) then
       ttgx=atgx
       goto 706
    elseif ( u .gt. 1.e+10 ) then
       ttgx=rtgx
       goto 706
    endif

    ! ml2 (and ml3) versions use different constants in q and r.

    q=17./81.*u**2
    r=26./729.*u**3+1./6.*u*(rtgx-atgx)
    s1=(r+(q**3+r**2)**0.5)**0.3333333
    s2=-q/s1

    if ( u .le. 1. ) then
       super=(s1+s2+8./9.*u)**2-u**2
       if ( super .le. 0. ) super=0.
       ttgx=atgx+super
       goto 706
    else
       super=3.*(s1+s2-1./9.*u)**3/u
       if ( super .le. 0. ) super=0.
       ttgx=rtgx-super
       goto 706
    endif

706 return

  end subroutine gradt2

  !************************************************************************

  subroutine intg

    !     Runge-Kutta-Gill method of integration.  See Ralston & Wiff
    !     Math Methods for Digital Computers  p.110  (Wiley,1960)

    ! Subroutines
    use utils_subroutines, only : check_time
    use chemprofiles_subroutines, only : comp_env

    ! Common blocks  
    use misc_const
    use cc
    use tenv
    use shells2
    use d 
    use slr
    use jee
    use comp
    use grad
    use threshh
    use mixl
    use crash
    use flags, only : kap_new
    use kappa, only : opalz

    implicit double precision(a-h,o-z)

    real*8, dimension(4,5) :: k, y, q, dd
    real*8, dimension(5) :: ax, bx, cx
    real*8, dimension(7) :: xavec
    real*8, parameter :: pecx_factor = 0.40 !scaling factor for convection
    integer :: newtime
    integer, pointer :: j
    
    data ax/0.e0,5.e-1,2.928932190e-1,1.7071067810e0,1.6666666667e-1/ 
    data bx/0.e0,2.e0,1.e0,1.e0,2.e0/
    data cx/0.e0,5.e-1,2.928932190e-1,1.7071067810e0,5.e-1/

    j => jj

    dept=deptx
    if ( ko .le. 1 ) then
       q(1,1)=0.
       q(2,1)=0.
       q(3,1)=0.
       q(4,1)=0.
    else
       q(1,1)=q(1,5)
       q(2,1)=q(2,5)
       q(3,1)=q(3,5)
       q(4,1)=q(4,5)
    endif

    h=dr
    n=ko
    y(1,1)=r(n)
    y(2,1)=p(n)
    y(3,1)=t(n)
    y(4,1)=mm(n)
    dd(1,2)=1.
    dd(2,2) =0.4342945*dpdr*10.**(-y(2,1))
    dd(3,2)=0.4342945*dtdr*10.**(-y(3,1))
    dd(4,2)=0.4342945*dmdr/(1.989e+33*10.**y(4,1))
    do jx=2,4
       mx=jx-1
       kx=jx+1
       do ix=1,4
          k(ix,jx)=h*dd(ix,jx)
          y(ix,jx)=y(ix,mx)+ax(jx)*(k(ix,jx)-bx(jx)*q(ix,mx))
          q(ix,jx)=q(ix,mx)+3.*ax(jx)*(k(ix,jx)-bx(jx)*q(ix,mx))-cx(jx)*k(ix,jx)
       enddo
       pl=y(2,jx) 
       tl=y(3,jx) 
       xmass=10.**y(4,jx)
       call comp_env 
       call estatem
       call check_time(newtime)
       if ((newtime - t0) .gt. xmaxtime) goto 806
       if (.not. kap_new ) then
          ol1x = opac3(rl1x,tl)
       else
          xavec(:) = 0.0
          xavec(1) = x(1)
          xavec(2) = x(2)
          xavec(3) = x(3)
          call opalz(10**tl,10**rl1x, xavec, fkappa, opr, opt)
          ol2x = dlog10(fkappa)
! write(*,'(10f8.4,6es13.5,4f8.3)') rl1x, tl, ol1x, ol2x, ol1x-ol2x, xavec(1:3)
          ol1x = ol2x
       endif

       rtgx = 7.732749e+9*10.**(ol1x+pl-4.*tl)*sl(j)/xmass
       dd(1,kx)=1. 
       dd(2,kx)=-0.4342945*27398. * 10.**(y(4,jx)-y(2,jx)) * rho * &
            (6.96e+10/y(1,jx))**2 

!     the mixing length is evaluated locally (shell n, halfway, shell n+1)

       if ( rtgx .ge. atgx ) then
          phsx=dabs(0.4342945/dd(2,kx)) 
          dept=dept+dabs(y(1,mx)-y(1,jx))
          if ( phsx .gt. dept ) phsx = dept
          pecx = 1.601743e+07 * cpx * 10.**(ol1x+2.5*rl1x+y(4,jx)-3.*y(3,jx) &
               -0.5*y(2,jx)) * dsqrt(dabs(drdtpx))*phsx**2*(6.96e+10/y(1,jx))**2 
          pecx=9./8./pecx
          pecx=pecx/aml/aml 

! pick between ml1 and ml3 versions of mixing-length theory 

          if(mls.eq.0)then
             !*** ml1 ***
             call gradt
          else
             !*** ml3 ***
! Some scaling to get the right convection zone depth
             pecx=(8./9.)*3./4.242440687*pecx  ! original, mysterious scaling
             pecx = pecx_factor * pecx ! new scaling
             call gradt2
          endif
       else
          ttgx=rtgx
       endif
       dd(3,kx)=ttgx*dd(2,kx) 
       dd(4,kx)=0.4342945*4.*pi*rho*y(1,jx)**2/(1.989e+33*10.**y(4,jx))
    enddo
    do  ix=1,4
       k(ix,5)=h*dd(ix,5)
       y(ix,5)=y(ix,4)+ax(5)*(k(ix,5)-bx(5)*q(ix,4))
       q(ix,5)=q(ix,4)+3.*ax(5)*(k(ix,5)-bx(5)*q(ix,4))-cx(5)*k(ix,5)
    enddo

    r(n+1)=y(1,5) 
    p(n+1)=y(2,5) 
    t(n+1)=y(3,5) 
    mm(n+1)=y(4,5) 
806 return

  end subroutine intg

  !************************************************************************

  subroutine write4

    !     print the structure of envelopes and convective transfer properties
    !     both scaled down,  also interpolate between shells to find exact
    !     physical properties when eta=0,20, at the bottom of the convection
    !     zone and when the mass of the envelope is=1e-6* stellar mass

    ! Subroutines
    use utils_subroutines, only : isnan

    ! Common blocks
    use misc_const
    use shells2
    use tenv
    use jee
    use comp
    use xcompp
    use slr
    use cc
    use d 
    use dirac
    use grad
    use surf
    use terp
    use wprop
    use writeitt
    use oldval
    use mixl
    use crash
    use flags, only: verbose, evoloutput

    implicit double precision(a-h,o-z)

    real*8, dimension(600) :: dd,amfrac, s
    integer, dimension(600) :: nbn2
    real*8, pointer :: sm
    integer, pointer :: j

    sm => smm
    j => jj

12  format(/13h log(m/msun)=,f6.3,2x,13h log(r/rsun)=,f6.3,2x, &
         11h log(grav)=,f6.3,2x,11h log(teff)=,f6.3,2x,13h log(l/lsun)=, &
         f6.3,2x,12h log(l*/m*)=,f6.3// 12h iterations=,i2, &
         59h  nc=-1 if convv + nondeg, 0 if convv + deg, +1 if nonconvv/)
37  format(4x,6hradius,5x,11hlog(1-r/rs),2x,11hlog(1-m/ms),2x,4hlogp, &
         4x,4hlogt,4x,6hlogrho,4x,7hopacity,5x,3heta,5x,3httg,8x,3hatg,8x, &
         3hrtg,6x,2hnc,3x,4hne/z)
38  format(2x,1pe11.5,3x,0pf8.4,5x,f8.4,1x,2(1x,f7.3),2x,f7.3,3x, &
         1pe9.3,2x,0pf6.2,3(2x,e9.3),2x,i2,3x,f5.3) 
101 format(/1x,8henvelope,2x,9hstructure,16x,6hscale=,i2,5x, &
         ' total shells: ',i5)
103 format('log(1-r/rs)=',f8.4,2x,'log(1-m/ms)=',f8.4,2x,'logp=', &
         f7.3,2x,'logt=',f6.3,2x,'logrho=',f6.3,2x,'logkappa=',f6.3,2x, &
         'eta=',f6.2)
104 format(//57h physical conditions at the bottom of the convection zone/)
106 format(///////31h convective transfer properties,16x,6hscale=i2)
107 format(/4x,6hz=rs-r,5x,11hlog(1-r/rs),2x,11hlog(1-m/ms),3x, &
         7hlocal g,6x,        6hpeclet,6x,3httg,6x,5hsuper,5x, &
         10hmix length,2x,5hfc/ft,3x,9hconv velo,3x,5hv/vth)
108 format(2x,1pe11.5,3x,0pf8.4,5x,f8.4,4x,1pe9.3,4x,e9.3,3x,0pf5.3, &
         3x,1pe9.3,3x,e9.3,3x,0pf5.3,3x,1pe9.3,3x,0pf5.3)
109 format(//' physical conditions when m(env)/m(star)=',1pd14.6/)
900 format(/,4x,6hradius,5x,11hlog(1-r/rs),2x,11hlog(1-m/ms),10x,3hbn2 &
         ,11x,3hac2,8x,2hje,4x,4hnbn2)
902 format(2x,1pe11.5,3x,0pf8.4,5x,f8.4,2(2x,1pe16.9),0p,3f12.9,i5)

    ! envelope structure
    rads=r(1)
    nmod=1
    if ( evoloutput .and. writeit )  write(26,101) nmod,nlim
    alc4pig = -6.07666 - 2. * dlog10(gs)
    amass = 10.**sm * amsun
    do n=1,nlim 
       if( 1.-r(n)/rads .ge. 1.d-18 ) then
          dd(n) = dlog10( 1.-r(n)/rads )
       else
          dd(n) = -10000.
       endif
       !*******************************************
       !*** check here for NaN's to avoid crash ***
       !*******************************************
       if(isnan(mm(n))) then
          ihotspot = 1
          goto 906
       endif
       !*******************************************
       amfrac(n)=10.**mm(n)/sm
       ame1=1.d0-amfrac(n)
       if (ame1 .ge. 1.d-18) then
          s(n)=dlog10(ame1)
       else
          s(n) = -10000.
       endif
    enddo

    xs=dlog10(sl(j))
    xr=dlog10(rs) 
    xg=dlog10(gs) 
    xt=dlog10(te) 
    xm=dlog10(sm) 
    ratio=xs-xm

    if ( evoloutput .and. writeit ) then
       write(26,12) xm,xr,xg,xt,xs,ratio,iter
    end if

33  format(' too few shells')
    if ( nsuff .le. 0 ) then
       if ( evoloutput ) then
          write(26,33)
          write(20,33)
       end if
       if ( verbose ) write(*,33)
    endif
    if ( evoloutput .and. writeit ) write(26,37)

    ion=-1
    do n=1,nlim,nmod
       if ( atg(n) .gt. rtg(n) ) then 
          nc=1
       else
          if ( el1(n) .lt. 0. ) then
             nc=-1
          else
             nc=0 
          endif
       endif
       opy=10.**ol1(n)

       !     compute degree of ionization xion 

       if ( ion .gt. 0. ) then
          xion=1.0
       else
          zzz = xcomp(n,1)*z(1) + xcomp(n,2)*z(2) + xcomp(n,3)*z(3) 
          !***************************************************
          !*** it will mess up here if something goes awry ***
          !***************************************************
          if (el1(n).gt.1e6) then
             ihotspot = 1
             goto 906
          endif

          if ( xion .gt. 1.0) then
             xion=1.0
             ion=1
          endif
       endif
       if ( evoloutput .and. writeit) then
          write(26,38) r(n),dd(n),s(n),p(n),t(n),rl1(n), &
               opy,el1(n),ttg(n),atg(n),rtg(n),nc,xion
       end if
    end do

    !* physical properties at places of interest

    it=1
    stop=dlog10(1.-10.**stpms)
    do n=1,nlim
       l=n
       if ( s(n) .gt. stop) goto 91
    end do
91  slope=(s(l-1)-stop)/(s(l-1)-s(l)) 
    ! returns to 22 from below

22  dfin=dd(l-1)+(dd(l)-dd(l-1))*slope
    sfin=s(l-1)+(s(l)-s(l-1))*slope
    pfin=p(l-1)+(p(l)-p(l-1))*slope
    tfin=t(l-1)+(t(l)-t(l-1))*slope
    rfin=rl1(l-1)+(rl1(l)-rl1(l-1))*slope
    ofin=ol1(l-1)+(ol1(l)-ol1(l-1))*slope
    efin=el1(l-1)+(el1(l)-el1(l-1))*slope
    if ( it .le. 1 ) then
       stp = 10.**stop
       if ( evoloutput .and. writeit ) then
          write(26,109) stp
          write(26,103) dfin,sfin,pfin,tfin,rfin,ofin,efin
       end if
       slum=sl(1) 
       u(1,is)=pfin
       u(2,is)=tfin
       v(1,is)=dlog10(slum) 
       v(2,is)=dlog10(1-10.**dfin)+ww(2,is)
       it=it+1

       !  skip over the next section with the following statement. 
       !  the next section computes the convective properties.
       !  it's caused some problems in the hot models

       if ( nlim1 .le. 0 ) go to 30
       do n=kind,nlim
          l=n
          if ( atg(n) .gt. rtg(n) ) goto 26
       end do
26     slope=(atg(l-1)-rtg(l-1))/(atg(l-1)-rtg(l-1)-atg(l)+rtg(l))
       go to 22
    endif

    if ( evoloutput .and. writeit ) then
       write(26,104)
       write(26,103) dfin,sfin,pfin,tfin,rfin,ofin,efin
    end if

    ! convection zone

    nmod1=((nlim1-kind)/100)+1

    if ( evoloutput .and. writeit ) then
       write(26,106) nmod1 
       write(26,107)
    end if

    ! compute geometrical depth zc, local gravity glocl, superadiabaticity
    ! super, ratio of convective to total flux fcf, convective velocity
    ! conv, ratio of convective to thermal velocities vcvt and peclet
    ! number

    do n=kind,nlim1,nmod1
       zc=rads-r(n)
       glocl=27398.*10.**mm(n)*(6.96e+10/r(n))**2
       super=ttg(n)-atg(n)
       fcf=1.-ttg(n)/rtg(n) 
       if ( fcf .le. 0.) then
          fcf=0.
          conv=0. 
       else
          conv=27398.*(te/5800.)**4*10.**(mm(n)-t(n)-rl1(n))*3.90e+33*phs(n) &
               *drdtp(n)*fcf/(16.*pi*cp(n)*r(n)**2)
          conv=conv*aml

          ! choose mixing-length version here
          ! (mls) =0 Bohm-Vitense or (mls=1) Bohm & Cassinelli
          ! default is Bohm-Vitense version

          if(mls.eq.1)then
             !*** ml2 or ml3 version ***
             conv=2.*dsqrt(2.d0)*conv
          endif
          if ( conv .gt. 0.0 ) conv = conv**0.333333
       endif
       azave = xcomp(n,1)*az(1)+xcomp(n,2)*az(2)+xcomp(n,3)*az(3) 
       vcvt=conv/(3.*bk_const*an0/azave*10.**t(n))**0.5

       !*** also need to adjust the peclet #. default is Bohm-Vitense. ***
       peclet=cp(n)*phs(n)*conv/4.535688e-04*10.**(2.*rl1(n)+ol1(n)-3.*t(n))
       peclet=peclet*aml*aml

       !*** Bohm & Cassinelli version  ***
       if(mls.eq.1)then
          peclet=peclet/dsqrt(2.d0)
       endif
       if ( evoloutput .and. writeit ) then
          write(26,108) zc,dd(n),s(n),glocl,peclet,ttg(n),super,phs(n),fcf, &
               conv,vcvt 
       end if
    end do

30  continue
    if ( evoloutput .and. writeit ) write(26,900)
    do n=1,nlim,nmod
       nbn2(n)=bn2(n)/dabs(bn2(n))
       if (bn2(n) .ge. 0.) then
          bn2(n)=dlog10(bn2(n))
       else
          bn2(n) = -7.
       endif
       if (ac2(n) .ge. 0.) then
          ac2(n)=dlog10(ac2(n))
       else
          ac2(n) = -7.
       endif

       if ( evoloutput .and. writeit ) then
          write(26,902) r(n),dd(n),s(n),bn2(n),ac2(n),xcomp(n,1), &
               xcomp(n,2),xcomp(n,3),nbn2(n)
       end if
    end do
906 return

  end subroutine write4

  !************************************************************************

  function opac3(rhol,tl)

    !  this acts like opac, but interpolates between opacities.
    !  the vector x(1-3) must be set for this to work!
    !  13 Nov 89 MAW
    !  Added in calls to opacder under certain conditions. 4-91 PAB

    use opacfin
    use jee
    use comp

    implicit double precision(a-h,o-z)

    real*8, dimension(3) :: z,az

    ! initialize constants

    az(1)=1.0
    az(2)=4.0
    az(3)=12.0
    z(1)=1.0
    z(2)=2.0
    z(3)=6.0
    opl = 0.
    jesav = je

    ! pure composition (comp. > 0.999) 

    do i = 1,3
       xx = x(i)
       if(xx .gt. 0.999)then
          je = i
          if(ifinal.eq.0)then
             opac3 = opac(rhol,tl)
          else
             call opacder(rhol,tl,opcs,dlklr,dlklt)
             opac3 = opcs
          endif
          je = jesav
          goto 999
       end if
    end do

    ! mixtures

    do i = 1,3
       xje = x(i)
       if ( xje.gt. 0.001) then
          je = i
          opl = opl + xje * 10.**opac(rhol,tl)
       end if
    end do
    je = jesav
    opac3   = dlog10(opl)

999 return
  end function opac3

  !************************************************************************

  function opac(d,t)

    !  compute the total (rad. + cond.) opacity of a specific rho-t point. 
    !  This function assumes a pure composition only. (See opac3 for mixtures.)
    !  NOTE:  in this version, the Hubbard-Lampe opacities have been replaced
    !  by the conductive opacities of itoh et al (see top of code for
    !  references.)

    ! Subroutines
    use ocon_functions, only : ocon

    ! Common blocks
    use opacs
    use jee
    use comp

    implicit double precision(a-h,o-z)

    real*8, dimension(2) :: two
    integer, parameter :: jrow = 8, krow = 29

    ! trap point on staggered rho,T grid
    j=1 
    k=krow
6   i=(j+k)/2
    if ( t - cst(i) ) 7,7,8 
7   k=i 
    go to 9
8   j=i 
9   if(iabs(j-k)-1) 10,10,6 
10  do 14 i=1,2
       l=k+i-2
       delta=d-csr(l)+1
       j=delta+1
       if(j-1) 11,11,12
11     j=2 
       go to 14
12     if(j-jrow) 14,14,13
13     j=jrow
14     two(i)=cso(je,l,j-1)+(cso(je,l,j)-cso(je,l,j-1))*(delta-j+1)

       orad=two(1)+(two(2)-two(1))*(t-cst(k-1))/(cst(k)-cst(k-1))
       if(dabs(orad).gt.40.) orad=dsign(40.d0,orad)
       ocond = ocon(d,t)
       opac=dlog10(1./(10.**(-orad)+10.**(-ocond)))
       return

     end function opac

     !************************************************************************

     subroutine opacder(d,t,opcs,dlklr,dlklt)

       ! opacities, but from Lagrange re-interpolator
       ! Call for each chemical species. If pure, use derivatives from
       ! within the routine. Otherwise, use 5 point differencing.
       ! See opac3 for further details.

       ! Subroutines
       use utils_subroutines, only : fh4
       use ocon_functions, only : ocon

       use jee
       use comp
       use opcswch
       use flags

       implicit double precision(a-h,o-z)

       real*8, dimension(4) :: op,dkrx,ytb,xtb
       real*8, dimension(30) :: ta,tb
       real*8, dimension(30,30,3) :: ty

       ! read in the opacities if we haven't already

       if ( first ) then
          if (verbose) write(*,*) 'entering opac 1st time'
          first = .false.
          read(25,*) iaend,jbend
          read(25,705) (ta(li),li=1,iaend)
          read(25,705) (tb(li),li=1,jbend)
705       format(8f10.5)
710       format(10f8.4)
          do lk=1,3
             do li=1,iaend
                read(25,710) (ty(li,lj,lk),lj=1,jbend) 
             enddo
          enddo
          if (verbose) write(*,*) 'opacs read'
       endif

       ! locate our place on the grid.

       ider=1
       if(km.lt.1)then
          i = 9
          j = 7
          km = 0
       endif
65     continue
       km=km+1

       ! set temp = xa and dens = xb

       xa = t
       xb = d

       if( xa .lt. ta(1) ) goto 50
       if( xa .gt. ta(iaend)) goto 50
1      if ( xa - ta(i)) 4, 6, 2
2      if (xa - ta(i+1)) 6, 5, 3
3      i = i + 2
       goto 1
4      i = i - 1
       goto  1
5      i = i + 1
6      continue

       if( xb .lt. tb(1) ) goto 50
       if( xb .gt. tb(jbend)) goto 50
11     if (xb - tb(j)) 14, 16, 12
12     if(xb-tb(j+1)) 16, 15, 13
13     j = j + 2
       goto 11
14     j = j - 1
       goto 11
15     j = j + 1
16     continue
       itop=0
       it10=i-1
       if(it10 .ge. 1) goto 150
       it10=it10+1
       itop=1
150    continue
       if(ty(it10,j,je).gt.-7.99 .or. ty(it10,j+1,je).gt.-7.99) goto 160
       it10=it10+1
       itop=1
       if(it10 .gt.i) goto 400
       goto 150
160    it1f=it10+3
       if(it1f.le.iaend) goto 170
       if(itop .ne. 0) goto 400
       it10=it10-1
       goto 150
170    continue
       if(ty(it1f,j,je).gt.-7.99 .or. ty(it1f,j+1,je).gt.-7.99) goto 180
       if(itop .ne. 0) goto 400
       it10 = it10-1 
       if(it10.le. 0) goto 400
       if(it1f .le. i+1) goto 400
       goto 150
180    continue

       int = 0
       do 
250       int = int + 1
          if (int .gt. 4) exit
          it1=it10+int-1
          if(ty(it1,j,je).le.-7.99 .or. ty(it1,j+1,je).le.-7.99) goto 200
          js1=j-1
          if(js1 .le. 0) js1=js1+1
          if(ty(it1,js1,je) .le. -7.99) js1=js1+1
          if(js1+3 .gt. jbend) js1=js1-1
          if(ty(it1,js1+3,je) .le. -7.99) js1=js1-1
          do k=1,4
             k1=js1+k-1
             !*** add conductive opac. to rad. opac. ***
             ocond = ocon(tb(k1),ta(it1))
             ytb(k)=dlog10(1./(10.**(-ty(it1,k1,je))+10.**(-ocond)))
             xtb(k)=tb(k1) 
          enddo
          call fh4(xb,ytb,xtb,op(int),dkrx(int),ider) 
          goto 250

200       if(ty(it1,j,je) .le. -7.99) js1=j+1
          if(ty(it1,j+1,je) .le. -7.99) js1=j-1
          !*** add conductive opac. to rad. opac. ***
          ocond1 = ocon(tb(js1+1),ta(it1))
          ocond2 = ocon(tb(js1),ta(it1))
          op1=dlog10(1./(10.**(-ty(it1,js1+1,je))+10.**(-ocond1)))
          op2=dlog10(1./(10.**(-ty(it1,js1,je))+10.**(-ocond2)))
          dkrx(int)=(op1-op2)/(tb(js1+1)-tb(js1))
          op(int)=op1+dkrx(int)*(xb-tb(js1+1))
       enddo

       do k=1,4
          xtb(k)=ta(it10+k-1)
       enddo
       call fh4(xa,op,xtb,opc,dlklt,ider)
       opcs=opc
       if(ider .ne. 1)then
          goto 1009
       endif
       call fh4(xa,dkrx,xtb,dlklr,dumfh,0)

       goto 1009
400    continue

       do k=1,2
          it1=i+k-1
          if(ty(it1,j,je).le.-7.99 .and. ty(it1,j+1,je).le.-7.99) goto 450
          js1=j
          if(ty(it1,j,je).le.-7.99) js1=j+1
          if(ty(it1,j+1,je).le.-7.99) js1=j-1
          frxs=(xb-tb(js1+1))/(tb(js1+1)-tb(js1))
          !*** add conductive opac. to rad. opac. ***
          ocond1 = ocon(tb(js1+1),ta(it1))
          ocond2 = ocon(tb(js1),ta(it1))
          op1=dlog10(1./(10.**(-ty(it1,js1+1,je))+10.**(-ocond1)))
          op2=dlog10(1./(10.**(-ty(it1,js1,je))+10.**(-ocond2)))
          op(k)=op1+(op1-op2)*frxs
          if(ider .ne. 1) goto 500
          !*** add conductive opac. to rad. opac. ***
          dop1=(op1-op2)/(tb(js1+1)-tb(js1))
          dop2=dop1
          if(js1-1 .le. 0) goto 420
          if(ty(it1,js1-1,je) .gt. -7.99) then
             !*** add conductive opac. to rad. opac. ***
             ocond1 = ocon(tb(js1-1),ta(it1))
             ocond2 = ocon(tb(js1),ta(it1))
             op1=dlog10(1./(10.**(-ty(it1,js1-1,je))+10.**(-ocond1)))
             op2=dlog10(1./(10.**(-ty(it1,js1,je))+10.**(-ocond2)))
             dop1=0.5*(dop1+(op2-op1)/(tb(js1)-tb(js1-1)))
          endif
420       if(js1+2 .gt. jbend) goto 430
          if(ty(it1,js1+2,je).gt.-7.99) then
             !*** add conductive opac. to rad. opac. ***
             ocond1 = ocon(tb(js1+1),ta(it1))
             ocond2 = ocon(tb(js1+2),ta(it1))
             op1=dlog10(1./(10.**(-ty(it1,js1+1,je))+10.**(-ocond1)))
             op2=dlog10(1./(10.**(-ty(it1,js1+2,je))+10.**(-ocond2)))
             dop2=0.5*(dop2+(op3-op1)/(tb(js1+2)-tb(js1+1)))
          endif
430       dkrx(k)=dop2+(dop2-dop1)/(tb(js1+1)-tb(js1))*(xb-tb(js1+1))
          goto 500

450       continue
          jc=1
          jco=0
          js1=j
          if(j.gt. 11) jc=-1
455       do ljs=1,jbend
             js1=js1+jc
             if(js1.ge.jbend .or. js1.le.1) goto 465
             if(ty(it1,js1,je) .gt. -7.99) goto 470
          enddo
465       if(jco .ne. 0) goto 50 
          jco=jco+1
          jc=-jc
          goto 455
470       if(jc .eq. -1) js1=js1-1
          !*** add conductive opac. to rad. opac. ***
          ocond1 = ocon(tb(js1+1),ta(it1))
          ocond2 = ocon(tb(js1),ta(it1))
          op1=dlog10(1./(10.**(-ty(it1,js1+1,je))+10.**(-ocond1)))
          op2=dlog10(1./(10.**(-ty(it1,js1,je))+10.**(-ocond2)))
          dkrx(k)=(op1-op2)/(tb(js1+1)-tb(js1))
          op(k)=op1+dkrx(k)*(xb-tb(js1+1))

       enddo
500    continue
       dlklt=(op(2)-op(1))/(ta(i+1)-ta(i))
       opcs=op(2)+dlklt*(xa-ta(i+1))
       if(ider .eq. 0) goto 1009
       dlklr=dkrx(2)+(dkrx(2)-dkrx(1))*(xa-ta(i+1))/(ta(i+1)-ta(i))
       goto 1009

50     continue

620    format(1h ,52hargument of opacity table is range out for log(t) =&
            ,f8.5, 16h, and log(rho) = , f8.4 )
       if (verbose) write(*,620) xa, xb
       stop

1009   return
     end subroutine opacder

     !**************************************************************************
     ! compute modified Ledoux B term to include composition gradient
     ! effects in the Brunt-Vaisala frequency
     ! Pressure is the dependent variable throught and I sum partial pressures.
     ! In this version I change only the He abundance at the H/He interface.
     ! while I keep the sum of the abundances unity at the He/C interface.
     ! This duplicates the numerical differencing results, but the physical 
     ! justification eludes me at present.
     ! PAB Sept 10, 1991
     ! Changed the way the partial pressure additions are done (MHM, Jan 1999).
     ! These changes can be found by searching for "MHM" in comments.
     !**************************************************************************
     subroutine ledoux

       ! Subroutines
       use utils_subroutines, only: hunt, spline, splint

       ! Common blocks
       use misc_const
       use cc
       use tenv
       use jee
       use comp
       use d 
       use alpha
       use bldx
       use crash

       implicit double precision(a-h,o-z)

       integer, parameter :: nd=38
       integer, pointer :: j
       real*8, dimension(nd) :: rhot1,rhot2,rhot3,rhot4,ptmp1,ptmp2,ptmp3,ptmp4
       real*8, dimension(nd) :: y21,y22,y23,y24
       real*8, dimension(60) :: temp
       real*8 :: yp1

       j => jj

       !.................................
       !  initialize summation variables
       !.................................
       bled    = 0.0
       ! step size increments in composition
       del     = 0.01

       !.............................
       ! This routine is complicated because we have to consider the
       ! H/He interface and He/C interfaces.
       ! We also have to make sure to avoid non-physical compositions
       ! such as X<0 or X>1.
       !.............................
       ! H/He interface
       if(x(1).gt.0.0 .and. x(2).gt.0.0)then
          if(x(3) .gt. x(1))then
             ! Very thin He buffer layer. Call H/He or He/C transition zone
             ! depending on case involved.
             go to 45
          endif
          ! first, locate T point on grid for Helium
          do i=1,lmax
             temp(i) = tk(2,i)
          enddo
          call hunt(temp,lmax,tl,ktlo)
          kthi = ktlo + 1
          ! set size of temporary arrays
          kmax1 = ml(2,ktlo)
          kmax2 = ml(2,kthi)
          ! fill up temporary 1-D arrays for rho and T He arrays
          do k=1,kmax1
             rhot1(k) = rhok(2,ktlo,k)
             ptmp1(k) =   pp(2,ktlo,k)
          enddo
          do k=1,kmax2
             rhot2(k) = rhok(2,kthi,k)
             ptmp2(k) =   pp(2,kthi,k)
          enddo
          ! first, locate T point on grid for Hydrogen
          do i=1,lmax
             temp(i) = tk(1,i)
          enddo
          call hunt(temp,lmax,tl,ktlo)
          kthi = ktlo + 1
          ! set size of temporary arrays
          kmax3 = ml(1,ktlo)
          kmax4 = ml(1,kthi)
          ! H arrays
          do k=1,kmax3
             rhot3(k) = rhok(1,ktlo,k)
             ptmp3(k) =   pp(1,ktlo,k)
          enddo
          do k=1,kmax4
             rhot4(k) = rhok(1,kthi,k)
             ptmp4(k) =   pp(1,kthi,k)
          enddo
          ! set BC's so second der is zero there
          yp1 = 1.0e+40
          ypn = 1.0e+40
          ! call for spline coefficients
          call spline(rhot1,ptmp1,kmax1,yp1,ypn,y21)
          call spline(rhot2,ptmp2,kmax2,yp1,ypn,y22)
          call spline(rhot3,ptmp3,kmax3,yp1,ypn,y23)
          call spline(rhot4,ptmp4,kmax4,yp1,ypn,y24)
          ! now compute P for each isotherm at a given rho,comp point
          call splint(rhot1,ptmp1,y21,kmax1,rho,p11)
          call splint(rhot2,ptmp2,y22,kmax2,rho,p12)
          call splint(rhot3,ptmp3,y23,kmax3,rho,p13)
          call splint(rhot4,ptmp4,y24,kmax4,rho,p14)
          ! linear interpolation for P at the actual temp.
          frac=(tl-tk(je,ktlo))/(tk(je,kthi)-tk(je,ktlo))
          pone=p11+frac*(p12-p11)
          ptwo=p13+frac*(p14-p13)
          pone=10.**(pone)
          ptwo=10.**(ptwo)

          !.................................................................
          ! now to consider specific composition mixtures!
          ! CASE 1: x(1) or x(2) is less than del (use a 3 point formula)
          ! CASE 1.1: Nearly pure helium: Treat first point as pure He.
          !           treat second point as H / He mixture 
          !.................................................................
          if(x(1).lt.del)then
             diff = 1.0-x(2)
             ! changed to take care of diff=0.0 case (MHM Jan. 1999)
             if (diff.lt.1.d-05) diff=1.d-05
             !             xone = x(2) - diff
             ! since x(1)<del, we are "close enough" to pure x(2) that
             ! I'm setting xtwo=1.0 and xone=1.0-diff. Now the total change
             ! in abundance is just diff instead of 2*diff (MHM Jan 1999)
             xone = 1.0 - diff
             xtwo = 1.0
             ! compute P for mixture by partial pressures.
             ! pmix1 is for x(2)-diff; pmix2 is for x(2)+diff
             !             pmix1 = xone*pone + (diff)*ptwo
             ! changed above line to the following so that respective mass
             ! fractions for pone and ptwo add up to one (MHM Jan. 1999)
             pmix1 = xone*pone + (1.-xone)*ptwo
             pmix1 = dlog(pmix1)
             ! same change as for pmix1 (MHM)
             pmix2 = xtwo*pone + (1.-xtwo)*ptwo
             pmix2 = dlog(pmix2)
             ! compute d ln P/d ln Y and return
             !             bled=(x(2)/(2.*diff))*(pmix2-pmix1)
             ! changed previous line to the following (MHM)
             bled=(1./(diff))*(pmix1-pmix2)
             goto 1109
             !.................................................................
             ! CASE 1.2: Nearly pure hydrogen: Treat first point as pure H.
             !           treat second point as H / He mixture 
             !.................................................................
          else if(x(2).lt.del)then
             diff=x(2)
             ! changed to take care of diff=0.0 case (MHM Jan. 1999)
             if (diff.lt.1.d-05) diff=1.d-05
             ! same comments as for Case 1.1, except now we set xone=0.0, 
             ! and xtwo=diff (MHM, Jan 1999)
             xone = 0.0
             xtwo = diff
             ! compute P for mixture by summing partial pressures.
             ! pmix1 is for x(2)-diff; pmix2 is for x(2)+diff
             pmix1 = xone*pone + (1.-xone)*ptwo
             pmix1 = dlog(pmix1)
             pmix2 = xtwo*pone + (1.-xtwo)*ptwo
             pmix2 = dlog(pmix2)
             bled=(x(2)/(diff))*(pmix1-pmix2)
             goto 1109
          endif

          !.................................................................
          ! CASE 2: We are not near 0 or 1 in composition
          ! we can use 3-point differencing formula with impunity.
          !.................................................................
          xone = x(2) - del
          xtwo = x(2) + del
          pmix1 = xone*pone + (1.-xone)*ptwo
          pmix2 = xtwo*pone + (1.-xtwo)*ptwo
          pmix1 = dlog(pmix1)
          pmix2 = dlog(pmix2)
          ! compute d ln P/d ln Y and return
          bled=(x(2)/(2.*del))*(pmix1-pmix2)
          goto 1109
       endif
       ! finished with H/He mixtures

       !.................................................................
       ! He/C interface
       if(x(2).gt.0.0 .and. x(3).gt.0.0)then
          ! first, locate T point on grid for Helium
45        do i=1,lmax
             temp(i) = tk(2,i)
          enddo
          call hunt(temp,lmax,tl,ktlo)
          kthi = ktlo + 1
          ! set size of temporary arrays
          kmax1 = ml(2,ktlo)
          kmax2 = ml(2,kthi)
          ! fill up temporary 1-D arrays for rho and T
          ! He arrays
          do k=1,kmax1
             rhot1(k) = rhok(2,ktlo,k)
             ptmp1(k) =   pp(2,ktlo,k)
          enddo
          do k=1,kmax2
             rhot2(k) = rhok(2,kthi,k)
             ptmp2(k) =   pp(2,kthi,k)
          enddo
          ! first, locate T point on grid for Carbon
          do i=1,lmax
             temp(i) = tk(3,i)
          enddo
          call hunt(temp,lmax,tl,ktlo)
          kthi = ktlo + 1
          ! set size of temporary arrays
          kmax3 = ml(3,ktlo)
          kmax4 = ml(3,kthi)
          ! C arrays
          do k=1,kmax3
             rhot3(k) = rhok(3,ktlo,k)
             ptmp3(k) =   pp(3,ktlo,k)
          enddo
          do k=1,kmax4
             rhot4(k) = rhok(3,kthi,k)
             ptmp4(k) =   pp(3,kthi,k)
          enddo
          ! set BC's so second der is zero there
          yp1 = 1.0e+40
          ypn = 1.0e+40
          ! call for spline coefficients
          call spline(rhot1,ptmp1,kmax1,yp1,ypn,y21)
          call spline(rhot2,ptmp2,kmax2,yp1,ypn,y22)
          call spline(rhot3,ptmp3,kmax3,yp1,ypn,y23)
          call spline(rhot4,ptmp4,kmax4,yp1,ypn,y24)
          ! the following line is a kludge to fix an unknown array problem
          ! in the previous 4 lines
          ! now compute P for each isotherm at a given rho,comp point
          call splint(rhot1,ptmp1,y21,kmax1,rho,p11)
          call splint(rhot2,ptmp2,y22,kmax2,rho,p12)
          call splint(rhot3,ptmp3,y23,kmax3,rho,p13)
          call splint(rhot4,ptmp4,y24,kmax4,rho,p14)
          ! linear interpolation for P at the actual temp.
          frac=(tl-tk(je,ktlo))/(tk(je,kthi)-tk(je,ktlo))
          pone=p11+frac*(p12-p11)
          ptwo=p13+frac*(p14-p13)
          pone=10.**(pone)
          ptwo=10.**(ptwo)

          !.................................................................
          ! now to consider specific composition mixtures!
          ! CASE 3: x(3) or x(2) is less than del (use a 3 point formula)
          ! CASE 3.1: Nearly pure helium: Treat first point as pure He.
          !           treat second point as He / C mixture 
          !.................................................................
          if(x(3).lt.del)then
             diff = 1.0-x(2)
             ! changed to take care of diff=0.0 case (MHM Jan. 1999)
             if (diff.lt.1.d-05) diff=1.d-05
             !           xone = x(2) - diff
             ! since x(3)<del, we are again "close enough" to pure x(2) that
             ! I'm setting xtwo=1.0 and xone=1.0-diff. Now the total change
             ! in abundance is just diff instead of 2*diff (MHM Jan 1999)
             xone = 1. - diff
             xtwo = 1.0
             ! compute P for mixture by summing partial pressures.
             ! pmix1 is for x(2)-diff; pmix2 is for x(2)+diff
             ! note that I have not changed the following lines. They were
             ! already "correct". (MHM, Jan. 1999)
             pmix1 = xone*pone + (1.0-xone)*ptwo
             pmix2 = xtwo*pone + (1.0-xtwo)*ptwo
             pmix1 = dlog(pmix1)
             pmix2 = dlog(pmix2)
             bled=(1./(diff))*(pmix1-pmix2)
             goto 1109

             !.................................................................
             ! CASE 3.2: Nearly pure Carbon: Treat first point as pure C.
             !           treat second point as He / C mixture 
             !.................................................................
          else if(x(2).lt.del)then
             diff=x(2)
             ! changed to take care of diff=0.0 case (MHM Jan. 1999)
             if (diff.lt.1.d-05) diff=1.d-05
             xone = 0.0
             xtwo = diff
             ! compute P for mixture by additive volume method.
             ! pmix1 is for x(2)-diff; pmix2 is for x(2)+diff
             pmix1 = xone*pone + (1.0-xone)*ptwo
             pmix2 = xtwo*pone + (1.0-xtwo)*ptwo
             pmix1 = dlog(pmix1)
             pmix2 = dlog(pmix2)
             ! compute d ln P/d ln Y and return.
             bled=(x(2)/(diff))*(pmix1-pmix2)
             goto 1109
          endif

          !.................................................................
          ! CASE 2: We are not near 0 or 1 in composition
          ! we can use 3-point differencing formula with impunity.
          !.................................................................
          xone = x(2) - del
          xtwo = x(2) + del
          pmix1 = xone*pone + (1.0-xone)*ptwo
          pmix2 = xtwo*pone + (1.0-xtwo)*ptwo
          pmix1 = dlog(pmix1)
          pmix2 = dlog(pmix2)
          ! compute d ln P/d ln Y and return.
          bled=(x(2)/(2.*del))*(pmix1-pmix2)
          goto 1109
       endif
       ! finished with He/C mixtures.

       !.................................................................
       ! If no transition zone, then set bled to 0.0
       bled = 0.0
1109   return

     end subroutine ledoux

     !************************************************************************

     subroutine smooth(nl,mver)

       ! this subroutine is a quick and dirty way to fix the jump in
       ! in derivative quants. at the core-envelope boundary.
       ! If the jump is greater than a pre-set limit (see eprep),
       ! we come here and replace the jump by a linear interpolation
       ! of delad between the last core point and some predetermined
       ! point out in the envelope.
       ! This is not a terribly physical subroutine, but adjusting the
       ! stop mass (where it's possible) show this assumption isn't
       ! too awful.

       use dprep

       implicit double precision(a-h,o-z)

       !print *, 'smooth'
       !*** set upper limit to smoothing function at 10 points (for now) ***
1      nup = nl + 20

       ! determine quant. to be smoothed and mass 
       ! at top and bottom of smoothing interval

       if(mver.eq.1)then
          delbot = aa(16,nl)
          deltop = aa(16,nup)
       endif
       if(mver.eq.2)then
          delbot = aa(15,nl)
          deltop = aa(15,nup)
       endif
       if(mver.eq.3)then
          delbot = aa(10,nl)
          deltop = aa(10,nup)
       endif
       if(mver.eq.4)then
          delbot = aa(9,nl)
          deltop = aa(9,nup)
       endif
       botmass = aa(2,nl)
       topmass = aa(2,nup)

       !*** determine size of mass and delad difference ***
       rmass = topmass - botmass
       deladr = dabs(delbot - deltop)
       do n=nl+1,nup

          !*** special case, deltop = delbot ***
          if(deltop.eq.delbot)then
             if(mver.eq.1)then
                aa(16,n) = delbot
             endif
             if(mver.eq.2)then
                aa(15,n) = delbot
             endif
             if(mver.eq.3)then
                aa(10,n) = delbot
             endif
             if(mver.eq.4)then
                aa(9,n) = delbot
             endif
          endif

          !*** determine step size in mass ***
          frac = (aa(2,n) - botmass) / rmass

          !*** determine stepsize in delad ***
          delfrac = deladr * frac

          ! compute new value of delad
          ! delad(top) < delad(bot)

          if(deltop.lt.delbot)then
             if(mver.eq.1)then
                aa(16,n) = delbot - delfrac
             endif
             if(mver.eq.2)then
                aa(15,n) = delbot - delfrac
             endif
             if(mver.eq.3)then
                aa(10,n) = delbot - delfrac
             endif
             if(mver.eq.4)then
                aa(9,n) = delbot - delfrac
             endif
          endif

          !*** delad(top) > delad(bot) ***
          if(deltop.gt.delbot)then
             if(mver.eq.1)then
                aa(16,n) = delbot + delfrac
             endif
             if(mver.eq.2)then
                aa(15,n) = delbot + delfrac
             endif
             if(mver.eq.3)then
                aa(10,n) = delbot + delfrac
             endif
             if(mver.eq.4)then
                aa(9,n) = delbot + delfrac
             endif
          endif
       enddo

       return

     end subroutine smooth

     !************************************************************************

   end module eprep_subroutines
