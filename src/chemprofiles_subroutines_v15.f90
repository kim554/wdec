module chemprofiles_subroutines

contains

!************************************************************************

subroutine comp_core(i,amr)

! Returns the oxygen, carbon, and helium mass fraction in the core.
! For now, we don't have hydrogen in the core, though that could be coded in
! by someone clever. There are no limitations as far as EOS tables and opacities
! are concerned. ABK, June 2016.

!Subroutines

!Common blocks
  use comp, only : x
  use flags, only : idif, chemprofmode
  use coeff, only : cfac1, cfac2
  use misc_const, only : amsun
  use shells, only: s
  use startmod, only : M_env
  use terp, only : stpms
  use xcompp

  implicit none
  
  real*8 :: amr, mmr, ame, mr_hyhe, mr_heca, logqi, xhe, xo2
  real*8 :: stpms10
  real*8 :: menv_mr,boundary
  integer :: i
  
! We have three regions to consider:
! I  . The tail end of the oxygen profile, mixed with Carbon
! II . The region where we have O/C/He all mixed
! Hydrogen does not reach into the core and is zero everywhere in this routine
! If it looks like there should be hydrogen in the core, reduce stop mass stpms
! log10(.99) is a good number
  stpms10 = 10.0d0**stpms
!  amr = amr*stpms10

  ame = 1.d0 - amr
  menv_mr = 1.03d0 - 10**(-M_env)
!  menv_mr = 10.0d0**stpms - 10.0d0**(-M_env)
!  print *, M_env, menv_mr
!  stop
  boundary = menv_mr-buffer_inner
!  boundary = 1.d0 - 10.d0**(-1.d0*boundary)

  if (menv_mr-buffer_inner .lt. 0.) then
     print *, "buffer_inner too large. Make it less than ", menv_mr 
     stop
  end if
! On the flip side, if the helium abundance takes off sharply from zero
! (appears chopped off on its inner edge, and/or there is a weird spike in the
! carbon abundance profile), try increasing buffer_inner in inputprof.

! 1 = Hydrogen
! 2 = Helium
! 3 = Carbon
! 4 = Oxygen
  
! Make up Oxygen profile
  call profsm1(amr,ams_o,corat_o,ndimo,ao,xo2)
  xcomp(i,4) = xo2

! Calculate core Helium abundance
! In region I, it's zero. In region II, calculate it
  if (amr .ge. boundary) then ! Region II
     ame = 1.d0 - amr
     call profche(ame,xhe)
  else                        ! Region I
     xhe = 0.d0
  end if

! Conclude
     xcomp(i,2) = xhe
! Carbon abundance is leftover, but no more than 1
     xcomp(i,3) = 1.d0 - xo2 - xhe
     if (xcomp(i,3) .gt. 1.d0) xcomp(i,3) = 1.d0
     xcomp(i,1) = 0.d0 ! Hydrogen abundance set to zero

  x(:) = xcomp(i,:)
!  print *, mmr, xcomp(i,:)
!  print *, -log10(1-mmr), xcomp(i,:)

  return

end subroutine comp_core

!************************************************************************
! Similar routine to above, but gets called in the envelope part of the code.
! At the base of the helium layer, we use profsm again to calculate the 
! helium abundance, similarly to how we did for the oxygen abundance (user
! defined in chemprinput_points.dat). At the base of the hydrogen layer, 
! use diffusive equilibrium.
! ABK, June 2016

subroutine comp_env

!Common blocks
  use cc, only : xmass, smm 
!xmass is Mr in solar masses, smm is stellar mass in solar masses
  use xcompp
  use comp
  use comptmp, only : hecdexp_tmp
  use flags, only : verbose
  use startmod, only: M_env
  use terp, only: stpms, stpms_orig

  implicit none

  real*8 :: logqi, mmr, ame, xh, xhe, xc, xo, menv_mr
  real*8 :: stpms10, stpms10_orig
  real*8 :: alpha3, alpha4, qhyhe

  alpha3 = hecdexp_tmp
  alpha4 = -hecdexp_tmp
  qhyhe = abs(log10(amr_hyhe))
  menv_mr = 1.d0 - 10.d0**(-M_env)

  mmr = xmass/smm !Mr now in stellar mass (1.0 being the surface)
!  mr = mr*0.873596
  ame = 1.d0 - mmr
  if ( ame .ne. 0.d0) then
     logqi = log10(ame)
  else
     logqi = -18.d0 !Somewhat arbitrary value that puts us close to the surface
  end if

! Below, we have three regions to consider:
! I  . The tail end of the oxygen profile, mixed with Carbon
! II . The region where we have O/C/He all mixed
! III. The region where C and O have gone to zero, and we have He/H mix

  if (menv_mr-buffer_inner .lt. 0.) then
     print *, "buffer_inner too large. Choose a value smaller than M_env"
     stop
  end if

  if (qhyhe-buffer_outer .lt. M_env) then
     print *, "buffer_outer too large. Choose a value smaller &
          than", qhyhe-M_env
     stop
  end if
  
! Region I is below M_env, staying away from M_env by buffer_inner
  if ( mmr .lt. (menv_mr-buffer_inner) ) then
     call profsm4(mmr,ams_o,corat_o,ndimo,ao,xo)
     x(1) = 0.d0      ! No Hydrogen allowed
     x(2) = 0.d0      ! No Helium yet either
     x(3) = 1.d0 - xo ! Carbon
     x(4) = xo        ! Oxygen
  end if

! Region II starts buffer_inner below M_env and ends buffer_outer below qhyhe
  if ( (mmr .ge. (menv_mr-buffer_inner)) .and. &
       (abs(logqi) .lt. (qhyhe - buffer_outer)) ) then
! Use profche for base of He layer
     call profche(ame,xhe)
! And we may still have some oxygen in that region
     call profsm1(mmr,ams_o,corat_o,ndimo,ao,xo)
     x(1) = 0.d0 ! H set to zero
     x(2) = xhe
     x(3) = 1.d0 - xo - xhe ! Carbon is 1 - oxygen - helium
     x(4) = xo
  end if

! Region III starts buffer_outer below qhyhe
! Use diffusive equilbrium at the He/H interface. We assume there is no
! oxygen or carbon in that region
  if ( abs(logqi) .ge. (qhyhe - buffer_outer) ) then
     call DIFUSE(amhyhe,amheca,alpha3,alpha4,ame,xhe)
     x(1) = 1.d0 - xhe
     x(2) = xhe
     x(3) = 0.d0 ! No more carbon
     x(4) = 0.d0 ! No more oxygen
  end if

!  print *, mmr, x(1), x(2), x(3), x(4)
!  print *, abs(logqi), x(1), x(2), x(3), x(4)

end subroutine comp_env

!***************************************************************
!                    subroutine profsm
! This routine uses the profile specified
! in the header in ams() and corat() and computes a smoothed version
! of it so that the derivatives of this profile are also smooth.
! The parameter "a" controls the sharpness of the changes in 
! slope, where a larger "a" means a more rapid change in slope
! (MHM March 2006)
!***************************************************************

  subroutine profsm1(amr,ams,corat,ndim,a,x2)

    implicit none
    
    integer :: ndim,i,npoints
    real*8, dimension(ndim) ::  ams,corat
    real*8 :: amr,x2,efac,a
    real*8, dimension(ndim-1) :: am
    real*8, dimension(ndim) :: dm !ndim-2
    real*8, parameter :: emax=25.d0

    do i=1,ndim
       if (ams(i).ge.0.99999999) goto 10
    enddo
10  continue
    npoints = i

    do i=1,npoints-1
       am(i) = (corat(i+1)-corat(i))/(ams(i+1)-ams(i))
    enddo
    do i=1,npoints-2
       dm(i) = am(i+1)-am(i)
    enddo

    x2=corat(1)
    do i=1,npoints-2
       efac=a*(amr-ams(i+1))
       if (efac.lt.-emax) then
          x2 = x2 + 0.0
       elseif (efac.gt.emax) then
          x2 = x2 + efac*dm(i)/a
       else
          x2 = x2 + log(1.d0+exp(efac))*dm(i)/a
       endif
    enddo
    return

  end subroutine profsm1

!*****************************************************************************
! Cubic Spline interpolation
! Things get weird at very small oxygen abundances
! Would need to implement the kind of kludge I did in the linear interpolation
! routine in order for this routine to be OK to use.
  subroutine profsm3(amr,xvals,yvals,ndim,a,x2)

    use utils_subroutines, only : spline, splint, hunt

    implicit none
    
    integer :: ndim,i
    real*8, dimension(ndim) ::  xvals, yvals,y2
    real*8 :: amr,x2,a,yp1,ypn

    yp1 = -40.d0
    ypn = 0.d0

    call spline(xvals,yvals,ndim,yp1,ypn,y2)
    call splint(xvals,yvals,y2,ndim,amr,x2)
 
    return

  end subroutine profsm3

!*****************************************************************************
! Linear interpolation
! Not smooth enough  

  subroutine profsm4(amr,xvals,yvals,ndim,a,x2)

    implicit none
    
    integer :: ndim,i,i0,i1
    real*8, dimension(ndim) ::  xvals, yvals
    real*8 :: amr,x2,a,slope,b

! Find points from input profile on either side of point of interest
    do i=1,ndim
       if ( (xvals(i) .gt. amr) .or. (i .eq. ndim) ) then
          i1 = i
          exit 
       end if
    end do

    if ( i1 .lt. ndim) then
       i0 = i1 + 1
! Get equation of the line defined by (xvals(i0),yvals(i0)) and
! (xvals(i1),yvals(i1))
       slope = (yvals(i1) - yvals(i0))/(xvals(i1)-xvals(i0))
       b = yvals(i1) - slope*xvals(i1)
       x2 = slope * amr + b
    else
       x2 = 0.d0
    end if
 
    return

  end subroutine profsm4

!*****************************************************************************
  subroutine profsm2(amr,xvals,yvals,ndim,a,x2)

    implicit none
    
    integer :: ndim,i
    real*8, dimension(ndim) ::  xvals, yvals
    real*8 :: amr,x2,xsum,a,efac
    real*8, dimension(ndim-1) :: mvals
    real*8, dimension(ndim-2) :: dm
    real*8, parameter :: bigfac = 10.d0

    do i=1,ndim-1
       mvals(i) = (yvals(i+1)-yvals(i))/(xvals(i+1)-xvals(i))
    enddo
   
    x2 = yvals(1) + mvals(1)*amr
    xsum = 0.d0

    do i=1,ndim-2
       efac = a*(amr-xvals(i+1))
! Exploit behavior of ln(1+exp(x)) for "large" x to prevent overflow
       if (efac .lt. bigfac) then
          dm(i) = (mvals(i+1)-mvals(i))*log(1.d0+exp(efac))
       else
          dm(i) = (mvals(i+1)-mvals(i))*(efac)
       end if
       xsum = xsum + dm(i)
    end do

    x2 = x2 + xsum/a

    return

  end subroutine profsm2

!*****************************************************************************
! Borrowing code from subroutine ccomp of makedb (ff_dk.f)

  subroutine profche(ame,xhe)

    use comp
    use startmod

    implicit none

    real*8 :: ame, xame, xamhyhe, xamheca, xc, xhe

!*** kludge to calculate smooth Dehner and Kawaler-type profiles ***
    if (ame.gt.0.0) then
       xame=-log10(ame)
    else
       xame=18.0
    endif

    xhe = 1.0-xhebar*fmike(alph1,M_env,xame)- &
         (1.0-xhebar)*fmike(alph2,M_he,xame)
   
    return

  contains

!************************************************************************

  function fmike(a0,x0,x)

!  simple/stupid analytical formula for calculating
!  diffusion profiles as in Dehner & Kawaler 1995, ApJ, L141

    implicit none

    real*8 :: fmike, a0, x0, x

    fmike = 1.d+00/(1.d+00 + exp(a0*(x-x0)))

    return

  end function fmike

end subroutine profche

!*****************************************************************************
! Revised profsm, rewritten in Phython by MHM to work more generally (no longer
! requires a flat start to profiles). Translated by ABK.
! Clever smoothing algorithm detailed in Agnes' research notebook 05/26/16

  subroutine profsm2_old(amr,xvals,yvals,ndim,a,xc2)

    implicit none
    
    integer :: ndim,i
    real*8, dimension(ndim) ::  xvals, yvals
    real*8 :: amr,xc2,xcsum,a,efac
    real*8, dimension(ndim-1) :: mvals
    real*8, dimension(ndim-2) :: dm
    real*8, parameter :: bigfac = 10.d0

    do i=1,ndim-1
       mvals(i) = (yvals(i+1)-yvals(i))/(xvals(i+1)-xvals(i))
    enddo
   
    xc2 = yvals(1) + mvals(1)*amr
    xcsum = 0.d0

    do i=1,ndim-2
       efac = a*(amr-xvals(i+1))
! Exploit behavior of ln(1+exp(x)) for "large" x to prevent overflow
       if (efac .lt. bigfac) then
          dm(i) = (mvals(i+1)-mvals(i))*log(1.d0+exp(efac))
       else
          dm(i) = (mvals(i+1)-mvals(i))*(efac)
       end if
       xcsum = xcsum + dm(i)
    end do

    xc2 = xc2 + xcsum/a

    return

  end subroutine profsm2_old

!*****************************************************************************

  subroutine profsm_orig(amr,ams,corat,ndim,xc2)

    use xcompp, only : xcomp_orig

    implicit none
    
    integer :: ndim,i,npoints
    real*8, dimension(ndim) ::  ams,corat
    real*8 :: amr,xc2,efac
    real*8, dimension(ndim-1) :: am
    real*8, dimension(ndim-2) :: dm
    real*8, parameter :: a=150.d0, emax=25.d0

    do i=1,ndim
       if (ams(i).ge.0.99999999) goto 10
    enddo
10  continue
    npoints = i

    do i=1,npoints-1
       am(i) = (corat(i+1)-corat(i))/(ams(i+1)-ams(i))
    enddo
    do i=1,npoints-2
       dm(i) = am(i+1)-am(i)
    enddo

    xc2=corat(1)
    do i=1,npoints-2
       efac=a*(amr-ams(i+1))
       if (efac.lt.-emax) then
          xc2 = xc2 + 0.0
       elseif (efac.gt.emax) then
          xc2 = xc2 + efac*dm(i)/a
       else
          xc2 = xc2 + log(1.d0+exp(efac))*dm(i)/a
       endif
    enddo
    return

  end subroutine profsm_orig

!****************************************************************************

  subroutine DIFUSE(qx1,qy1,alpha3,alpha4,ame,y)

! Computes abundances assuming equilibrium diffusion profiles, using
! exact solutions of equation (A5) in Arcoragi and Fontaine. The
! relationship between qm, the SURFACE mass fraction at the midpoint of
! the profile, and qh and qhe, the total mass fraction of hydrogen
! and/or helium, is from numeric integration over the abundance profiles.

! This routine assumes that the He/H profile is in diffusive equilibrium
! (i.e., alpha(1) and alpha(2) are not used). If alpha(4) > -1, then it
! warns about significant He in the core, which may not be a good thing.
! alpha(3) and alpha(4) are used as pseudo-diffusion exponents to compute
! the upper and lower parts of the C/He profile, respectively.

! Note that qx1 and qy1 are mass fractions, NOT surface mass fractions
! MHM

! This version stripped to only compute helium abundance, and to only
! calculate abundances at the He/H interface. The C/He interface is handled by 
! profsm.
! ABK 05/31/2016

! Subroutines
    use utils_subroutines, only : zbrent, beta, hypgeo

! Common blocks
    use flags
    use coeff
    use crash

    implicit none

    real*8, parameter :: smallnegative = -1.0d-4
    real*8 ::qxtot,qytot,qx1,qy1,ame,y
    real*8 :: gam,tmp,tol,g1,g2,alpha3,alpha4,tmp2,r
!    real*8, external :: fheh, fche
    complex*16 :: h1,h2,h3,h4,cunit,h5,h6,h7,h8
  
    !print *, 'difuse'
    if (idif.eq.0) then
       cfac1=1.76351076
       if (alpha4.ge.-1.0) then
          aa=1./alpha3
          bb=-(1./alpha3 + 1./alpha4)
          cfac2=0.0
          qs= 1.0
          qmy=1.0
          cfac2=fche(1.d0)
       else
          aa=1./alpha3
          bb=-(1./alpha3 + 1./alpha4)
          r=1./3.0           ! A_He/A_C = 4/12 = 1/3
          cunit=(1.,0.)
          h1=-bb*cunit
          h2=aa*cunit
          h3=(1.-bb)*cunit
          h4=(1.-r)*cunit
          h5=(1.-bb)*cunit
          h6=(aa+1.)*cunit
          h7=(2.-bb)*cunit
          h8=h4
          cfac2 = aa*r**(aa)*beta(aa,1.-bb-aa)*hypgeo(h1,h2,h3,h4) + &
               bb*r**(aa+1.)*beta(aa+1.,1.-bb-aa)*hypgeo(h5,h6,h7,h8)
      endif
       idif=1
    endif
    
    qxtot=qx1
    qytot=qy1
    
    qmx=qx1
    qmy=qy1
    
    qs=ame
    
    if (qs.le.0.d00) then
       if (qxtot.le.0.d0) then
          y=1.0
          goto 206
       else
          y=0.0
       endif
    endif
    
    tol=1.d-06

    if (qxtot.gt.0.d0) then
       g1=0.0
       g2=100.*(cfac1 * (qs/qmx))**(1.25)
       gam=zbrent(fheh,g1,g2,tol)
       tmp=4.*gam/(1.+4.*gam)
       y=tmp
    else
       y=1.d0
    endif
          
    if (xheneg.eq.1) then
       if (y.lt.0.0 .and. y.gt.smallnegative) y = 0.0
    endif
    
206 return  

    contains

!**********************************************************************

      function fheh(x)
        
        use coeff
        
        implicit none
        real*8 :: x,fheh
        
        !print *, 'fheh'
        fheh=x**0.2*(1.+x)**0.72*(1.+2.*x)**(-0.1)*(1.+6.*x)**(-0.02)- &
             cfac1*(qs/qmx)
        
        return
        
      end function fheh
      
!****************************************************************************
      
      function fche(x)
        
        use coeff
        
        implicit none
        
        real*8 :: x,fche
    
        !print *, 'fche'
        fche=x**aa * (1.+x)**bb - cfac2 * (qs/qmy)
    
        return
  
      end function fche

!*****************************************************************************

    end subroutine DIFUSE

end module chemprofiles_subroutines
