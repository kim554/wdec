module chemprofiles_subroutines

contains

!************************************************************************

subroutine comp_core(i,amr)

! This being a helium core white dwarf, simply returns pure helium.
  
!Common blocks
  use xcompp
  use comp, only : x
  use comptmp, only : sm_tmp
  use flags, only : idif, chemprofmode
  use coeff, only : cfac1, cfac2
  use misc_const, only : amsun

  implicit none
  
  real*8 :: amr, ame, xh
  integer :: i
  
  xcomp(i,1) = 0.d0 ! hydrogen
  xcomp(i,2) = 1.d0 ! pure helium
  xcomp(i,3) = 0.d0 ! no carbon
  xcomp(i,4) = 0.d0 ! no oxygen

  x(:) = xcomp(i,:)
!  print *, amr, xcomp(i,:)
  
  return

end subroutine comp_core

!************************************************************************
! Similar routine to above, but gets called in the envelope part of the code.
! ABK, June 2016

subroutine comp_env

!Common blocks
  use cc, only : xmass, smm
  use xcompp
  use comp
  use comptmp, only : sm_tmp

  implicit none

  real*8 :: logqi, mr, ame, xh

  mr = xmass/smm  
  ame = 1.d0 - mr   

  call profche(ame,xh)

  x(1) = xh        !hydrogen
  x(2) = 1.d0 - xh !helium
  x(3) = 0.d0      ! No carbon
  x(4) = 0.d0      ! No oxygen

!  write(8000,*) mr, x(1), x(2), x(3), x(4)

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
! Revised profsm, rewritten in Phython by MHM to work more generally (no longer
! requires a flat start to profiles). Translated by ABK.
! Clever smoothing algorithm detailed in Agnes' research notebook 05/26/16

  subroutine profsm2(amr,xvals,yvals,ndim,a,xc2)

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

  end subroutine profsm2

!*****************************************************************************
! Borrowing code from subroutine ccomp of makedb (ff_dk.f)
! Names for helium core WD may be confusing. Essentially we are calculating
! the hydrogen abundance in the He --> H transtion zone instead of the helium
! abundance in the C --> He transition zone. Fortran doesn't care.

  subroutine profche(ame,xhe)

    use comp
    use startmod

    implicit none

    real*8 :: ame, xame, xamhyhe, xamheca, xc, xhe

! ame = 1 - Mr/M*
! xame = log(1-Mr/M*)
! xhebar = abundance of hydrogen in mixed He/H region
! M_env = location of the base of the Hydrogen/Helium envelope    
! M_he = location of the base of the pure Hydrogen layer
    
    if (ame.gt.0.0) then
       xame=-log10(ame)
    else
       xame=18.0
    endif

!*** kludge to calculate smooth Dehner and Kawaler-type profiles ***
    xhe = 1.0-xhebar*fmike(alph1,M_env,xame)- &
         (1.0-xhebar)*fmike(alph2,M_he,xame)
!    write(8000,*) xame, 1.0-xhebar*fmike(alph1,M_env,xame), (1.0-xhebar)*fmike(alph2,M_he,xame)
    
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

!*****************************************************************************

  subroutine profsm_orig(amr,ams,corat,ndim,xc2)

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

end module chemprofiles_subroutines
