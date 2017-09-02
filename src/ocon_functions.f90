module ocon_functions

contains

!************************************************************************

  function ocon(d,t)

!  compute the conductive opacity of a specific rho-t point. 
!  note:  in this version, the hubbard lampe opacities have been replaced
!  by the conductive opacities of itoh et al (see top of code for
!  references.)

    use jee
    use comp

    implicit double precision(a-h,o-z)
    
    real*8 :: d, t

!  the itoh opacities are only good down to log(rho) = 1.6
!  therefore, maw put in this kludge where if the density is
!  less than 1.8 and greater than 1.5, the interpolate between
!  opacity tables.  if d < 1.5, then use hubbard-lampe conductive
!  opacities.
    if (d .ge. 1.8) then
       ocon = oconit(d,t)
    elseif ( d .ge. 1.5 ) then
       dfrac = (d - 1.5) / 0.3
       ocon = dfrac * oconit(d,t) + (1.-dfrac) * oconhl(d,t)
    else
       ocon = oconhl(d,t)
    endif
    if(dabs(ocon).gt.40.)then
       ocon=dsign(40.d0,ocon)
    endif
    return
  end function ocon

!************************************************************************

  function oconit(d,t)

!  compute the conductive opacities using the method of itoh et al.
!  see top of program for references.

    use jee
    use comp
    use temp

    implicit double precision(a-h,o-z)

    real*8 :: d, t
    real*8, parameter :: alac43 = -3.519188
!  gamma constants
    real*8, parameter :: gamc3c = 3.5772e+06, gamc3o = 5.7782e+06, &
         gamc2 = 573264., gamc1 = 2.275e+05

! compute gamma, and call the appropriate routine (liquid metal or
! crystalline).  if gamma > 5000, then set gamma = 5000.
! This is the suggested proceedure from itoh et al (1984), since
! the results barely change when gamma > 5000.
! NOTE: there is no fit to xtal hydrogen.
    if ( je .eq. 3 ) then
       gammac = 10.**(d/3. - t) * gamc3c
       gammao = gammac / gamc3c * gamc3o
       gamma = xc2 * gammac + xo2 * gammao

! Carbon contribution

       if ( xc2 .gt. 0. ) then
          if ( gamma .lt. 171. ) then 
             condc = cndlc(d,t,gammac)
          elseif ( gamma .ge. 171 ) then
             if ( gammac .gt. 5000.) gammac = 5000.
             condc = cndsc(d,t,gammac)
          endif
       else
          condc = 0.
       endif

! Oxygen contribution

       if ( xo2 .gt. 0. ) then
          if ( gamma .lt. 171. ) then 
             condo = cndlo(d,t,gammao)
          elseif ( gamma .ge. 171 ) then
             if ( gammao .gt. 5000.) gammao = 5000.
             condo = cndso(d,t,gammao)
          endif
       else
          condo = 0.
       endif
       oconit = alac43 + (t*3. - d - (xc2*condc + xo2*condo) )

! Helium contribution

    elseif ( je .eq. 2 ) then
       gamma = gamc2 * 10.**(d/3. - t)
       if ( gamma .lt. 171. ) then
          conduct = cndlhe(d,t,gamma)
       elseif ( gamma .ge. 171 ) then 
          if ( gamma .gt. 5000.) gamma = 5000.
          conduct = cndshe(d,t,gamma)
       endif
       oconit = alac43 + (t*3. - d - conduct)

! Hydrogen contribution

    elseif ( je .eq. 1 ) then
       gamma = gamc1 * 10.**(d/3. - t)
       conduct = cndlhy(d,t,gamma)
       oconit = alac43 + (t*3. - d - conduct)
    endif

    return
  end function oconit

!****************************************************************************

  function oconhl(d,t)

! get Hubbard-Lampe opacities for the appropriate composition

    use jee
    use comp

    implicit double precision(a-h,o-z)

    real*8 :: d, t

    if ( je .eq. 3 ) then
       polytmp = hlcnc(d,t)
    elseif ( je .eq. 2 ) then
       polytmp = hlcnhe(d,t)
    elseif ( je .eq. 1 ) then
       polytmp = hlcnh(d,t)
    endif
    
    oconhl = polytmp
    
    return
  end function oconhl

!************************************************************************

  function cndlc(dens,t,gamma)

! this function returns the (carbon) thermal conductivity computed for a
! specific dens, temp, and gamma using the results if itoh et al (1983)
! for the liquid metal regime.

    use cndlc_data

    implicit double precision (a-z)
 
    real*8 :: dens, t, gamma
    real*8, parameter :: at = 12.d0

    t8 = 10.**(t -8)
    d6 = 10.**(dens - 6 )
    d6mu = d6 / 2.
    d6mu23 = d6mu**(2./3.)
 
    rs    = 1.388e-02 / d6mu**(1./3.) 
    
    x = 0.45641 * dlog( gamma ) - 1.31636
    x2 = x*x
    x3 = x2 * x
    
    asum = a(1) + a(2) * x + a(3) * x2 + a(4) * x3
    bsum = b(1) + b(2) * x + b(3) * x2
    csum = c(1) + c(2) * x + c(3) * x2
    dsum = d(1) + d(2) * x + d(3) * x2 + d(4) * x3
    esum = e(1) + e(2) * x + e(3) * x2
    fsum = f(1) + f(2) * x + f(3) * x2
    
    sminus1 = asum * ( 1. + bsum * rs + csum * rs*rs )
    splus1  = dsum * ( 1. + esum * rs + fsum * rs*rs )
    
    cstd6  = 1.018 * d6mu23 
    tmpsum = 1. + cstd6
    s = sminus1 - splus1 * cstd6 / tmpsum
    
    opacity = 2.363e+17 * d6 * t8 / ( at * tmpsum * s )
    cndlc  = dlog10( opacity )
    
    return
  end function cndlc

!************************************************************************

  function cndsc(dens,t,gamma)

! this function returns the thermal conductivity computed for a specific
! dens, temp, and gamma using the results if itoh et al (1984) for carbon
! in the xtal regime.
    
    use cndsc_data
    use misc_const, only: pi

    implicit double precision (a-z)
    
    real*8 :: dens, t, gamma

!  declare variables (find the values in the block data piece of code).

    real*8, parameter :: at = 12.d0, z = 6.d0

    twopi = 2.*pi
    pisqr = pi**2
    x13 = -1./3.
    x23 = -2./3.
    t8 = 10.**t / 1.e+08
    d6 = 10.**dens / 1.e+06 
    d6mu = d6 * z / at
    d6mu23 = d6mu**(2./3.)

! compute the interpolation coefficients v and w

    u = twopi * ( dens - 3. ) / 9.
    v = alpha(1) + alpha(2)*gamma**x13 + alpha(3)*gamma**x23 + alpha(4)/gamma
    w = beta(1) + beta(2) * gamma**x13 + beta(3) * gamma**x23 + beta(4) / gamma

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

!  here's the loops, unwrapped

    asum = a(1) + a(2) * cu    + a(3) * c2u + a(4) * c3u + a(5) * c4u
    bsum = b(1) * su + b(2) * s2u + b(3) * s3u + c*u + d
    esum = e(1) + e(2) * cu    + e(3) * c2u + e(4) * c3u + e(5) * c4u
    fsum = f(1) * su + f(2) * s2u + f(3) * s3u + ggg*u + h
    isum = i(1) + i(2) * cu    + i(3) * c2u + i(4) * c3u + i(5) * c4u
    jsum = j(1) * su + j(2) * s2u + j(3) * s3u + k*u + l
    psum = p(1) + p(2) * cu    + p(3) * c2u + p(4) * c3u + p(5) * c4u
    qsum = q(1) * su + q(2) * s2u + q(3) * s3u + r*u + s
 
    isig171 = asum + bsum
    isig5k  = esum + fsum
    ikap171 = isum + jsum
    ikap5k  = psum + qsum
 
    isig    = (1.-v) * isig171 + v * isig5k
    ikap    = (1.-w) * ikap171 + w * ikap5k
    
    smgamma = .07832 * (z/at) * dsqrt(d6) / t8
    sg2     = smgamma * smgamma
 
    g0      = 13.00 / dsqrt( 1. + 0.0174 * sg2 ) 
    g2   = (sg2 / pisqr) * (1. + 0.0118 * sg2)**(-3./2.)
 
    fkap    = isig * g0 + ikap * g2
    cstd6   = 1.018 * d6mu23
    nukap   = 9.554e+16 * t8 * dsqrt( 1. + 1./cstd6 ) * fkap
    
    opacity = 4.146d+15 * d6mu / dsqrt( 1. + cstd6 ) * t8 * 1.e+18 / nukap
 
    cndsc  = dlog10( opacity )
 
    return
  end function cndsc

!************************************************************************

  function cndlo(dens,t,gamma)

! this function returns the (oxygen) thermal conductivity computed for a
! specific dens, temp, and gamma using the results if itoh et al (1983)
! for the liquid metal regime.

    use cndlo_data

    implicit double precision (a-z)
 
    real*8 :: dens, t, gamma
    real*8, parameter :: at = 16.d0

!  the following data from itoh et al. are for oxygen (z=8) 

    t8 = 10.**(t -8)
    d6 = 10.**(dens - 6 )
    d6mu = d6 / 2.
    d6mu23 = d6mu**(2./3.)
    
    rs    = 1.388e-02 / d6mu**(1./3.) 
    
    x = 0.45641 * dlog( gamma ) - 1.31636
    x2 = x*x
    x3 = x2 * x
 
    asum = a(1) + a(2) * x + a(3) * x2 + a(4) * x3
    bsum = b(1) + b(2) * x + b(3) * x2
    csum = c(1) + c(2) * x + c(3) * x2
    dsum = d(1) + d(2) * x + d(3) * x2 + d(4) * x3
    esum = e(1) + e(2) * x + e(3) * x2
    fsum = f(1) + f(2) * x + f(3) * x2
    
    sminus1 = asum * ( 1. + bsum * rs + csum * rs*rs )
    splus1  = dsum * ( 1. + esum * rs + fsum * rs*rs )
    
    cstd6  = 1.018 * d6mu23 
    tmpsum = 1. + cstd6
    s = sminus1 - splus1 * cstd6 / tmpsum
    
    opacity = 2.363e+17 * d6 * t8 / ( at * tmpsum * s )
    cndlo  = dlog10( opacity )
    
    return
  end function cndlo

!************************************************************************

  function cndso(dens,t,gamma)

! this function returns the thermal conductivity computed for a specific
! dens, temp, and gamma using the results if itoh et al (1984) for oxygen
! in the xtal regime.

    use cndso_data
    use misc_const, only: pi

    implicit double precision (a-z)

!  declare variables (find the values in the block data piece of code).

    real*8 :: dens, t, gamma
    real*8, parameter :: at = 16.d0, z = 8.d0

    twopi = 2.*pi
    pisqr = pi**2
    x13 = -1./3.
    x23 = -2./3.
    t8 = 10.**t / 1.e+08
    d6 = 10.**dens / 1.e+06 
    d6mu = d6 / 2.
    d6mu23 = d6mu**(2./3.)

! compute the interpolation coefficients v and w

    u = twopi * ( dens - 3. ) / 9.
    v = alpha(1) + alpha(2)*gamma**x13 + alpha(3)*gamma**x23 + alpha(4)/gamma
    w = beta(1) + beta(2)*gamma**x13 + beta(3)*gamma**x23 + beta(4)/gamma

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

!  here's the loops, unwrapped

    asum = a(1) + a(2) * cu    + a(3) * c2u + a(4) * c3u + a(5) * c4u
    bsum = b(1) * su + b(2) * s2u + b(3) * s3u + c*u + d
    esum = e(1) + e(2) * cu    + e(3) * c2u + e(4) * c3u + e(5) * c4u
    fsum = f(1) * su + f(2) * s2u + f(3) * s3u + ggg*u + h
    isum = i(1) + i(2) * cu    + i(3) * c2u + i(4) * c3u + i(5) * c4u
    jsum = j(1) * su + j(2) * s2u + j(3) * s3u + k*u + l
    psum = p(1) + p(2) * cu    + p(3) * c2u + p(4) * c3u + p(5) * c4u
    qsum = q(1) * su + q(2) * s2u + q(3) * s3u + r*u + s
 
    isig171 = asum + bsum
    isig5k  = esum + fsum
    ikap171 = isum + jsum
    ikap5k  = psum + qsum
    
    isig    = (1.-v) * isig171 + v * isig5k
    ikap    = (1.-w) * ikap171 + w * ikap5k
    
    smgamma = .07832 * (z/at) * dsqrt(d6) / t8
    sg2     = smgamma * smgamma
    
    g0      = 13.00 / dsqrt( 1. + 0.0174 * sg2 ) 
    g2   = (sg2 / pisqr) * (1. + 0.0118 * sg2)**(-3./2.)
    
    fkap    = isig * g0 + ikap * g2
    cstd6   = 1.018 * d6mu23
    nukap   = 9.554e+16 * t8 * dsqrt( 1. + 1./cstd6 ) * fkap
    
    opacity = 4.146d+15 * d6mu / dsqrt( 1. + cstd6 ) * t8 * 1.e+18 / nukap
    
    cndso  = dlog10( opacity )
    
    return
  end function cndso

!************************************************************************

  function cndlhe(dens,t,gamma)

! this function returns the (helium) thermal conductivity computed for a
! specific dens, temp, and gamma using the results if itoh et al (1983)
! for the liquid metal regime.

    use cndlhe_data

    implicit double precision (a-z)
 
    real*8 :: dens, t, gamma
    real*8, parameter :: at = 4.d0
    
    t8 = 10.**(t - 8.)
    d6 = 10.**(dens - 6. )
    d6mu = d6 / 2.
    d6mu23 = d6mu**(2./3.)
    
    rs    = 1.388e-02 / d6mu**(1./3.) 
 
    x = 0.45641 * dlog( gamma ) - 1.31636
    x2 = x*x
    x3 = x2 * x
    
    asum = a(1) + a(2) * x + a(3) * x2 + a(4) * x3
    bsum = b(1) + b(2) * x + b(3) * x2
    csum = c(1) + c(2) * x + c(3) * x2
    dsum = d(1) + d(2) * x + d(3) * x2 + d(4) * x3
    esum = e(1) + e(2) * x + e(3) * x2
    fsum = f(1) + f(2) * x + f(3) * x2
    
    sminus1 = asum * ( 1. + bsum * rs + csum * rs*rs )
    splus1  = dsum * ( 1. + esum * rs + fsum * rs*rs )
    
    cstd6  = 1.018 * d6mu23 
    tmpsum = 1. + cstd6
    s = sminus1 - splus1 * cstd6 / tmpsum
    
    opacity = 2.363e+17 * d6 * t8 / ( at * tmpsum * s )
    cndlhe = dlog10( opacity )
    
    return
  end function cndlhe

!************************************************************************

  function cndshe(dens,t,gamma)

! this function returns the thermal conductivity computed for a specific
! dens, temp, and gamma using the results if itoh et al (1984) for helium
! in the xtal regime.
    
    use cndshe_data
    use misc_const, only: pi

    implicit double precision (a-z)

!  declare variables (find the values in the block data piece of code).

    real*8 :: dens, t, gamma
    real*8, parameter :: at = 4.d0, z = 2.d0

    twopi = 2.*pi
    pisqr = pi**2
    x13 = -1./3.
    x23 = -2./3.
    t8 = 10.**t / 1.e+08
    d6 = 10.**dens / 1.e+06 
    d6mu = d6 * z / at
    d6mu23 = d6mu**(2./3.)

! compute the interpolation coefficients v and w

    u = twopi * ( dens - 3. ) / 9.
    v = alpha(1) + alpha(2)*gamma**x13 + alpha(3)*gamma**x23 + alpha(4)/gamma
    w = beta(1) + beta(2) * gamma**x13 + beta(3) * gamma**x23 + beta(4) / gamma

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

!  here's the loops, unwrapped

    asum = a(1) + a(2) * cu    + a(3) * c2u + a(4) * c3u + a(5) * c4u
    bsum = b(1) * su + b(2) * s2u + b(3) * s3u + c*u + d
    esum = e(1) + e(2) * cu    + e(3) * c2u + e(4) * c3u + e(5) * c4u
    fsum = f(1) * su + f(2) * s2u + f(3) * s3u + ggg*u + h
    isum = i(1) + i(2) * cu    + i(3) * c2u + i(4) * c3u + i(5) * c4u
    jsum = j(1) * su + j(2) * s2u + j(3) * s3u + k*u + l
    psum = p(1) + p(2) * cu    + p(3) * c2u + p(4) * c3u + p(5) * c4u
    qsum = q(1) * su + q(2) * s2u + q(3) * s3u + r*u + s
 
    isig171 = asum + bsum
    isig5k  = esum + fsum
    ikap171 = isum + jsum
    ikap5k  = psum + qsum
 
    isig    = (1.-v) * isig171 + v * isig5k
    ikap    = (1.-w) * ikap171 + w * ikap5k
 
    smgamma = .07832 * (z/at) * dsqrt(d6) / t8
    sg2     = smgamma * smgamma
    
    g0      = 13.00 / dsqrt( 1. + 0.0174 * sg2 ) 
    g2   = (sg2 / pisqr) * (1. + 0.0118 * sg2)**(-3./2.)
 
    fkap    = isig * g0 + ikap * g2
    cstd6   = 1.018 * d6mu23
    nukap   = 9.554e+16 * t8 * dsqrt( 1. + 1./cstd6 ) * fkap
    
    opacity = 4.146d+15 * d6mu / dsqrt( 1. + cstd6 ) * t8 * 1.e+18 / nukap
 
    cndshe = dlog10( opacity )
 
    return
  end function cndshe

!************************************************************************

  function cndlhy(dens,t,gamma)

! this function returns the (hydrogen) thermal conductivity computed for a
! specific dens, temp, and gamma using the results if itoh et al (1983)
! for the liquid metal regime.
    
    use cndlhy_data

    implicit double precision (a-z)

    real*8 :: dens, t, gamma
    real*8, parameter :: at = 1.d0

    t8 = 10.**(t -8.)
    d6 = 10.**(dens - 6. )
    d6mu = d6
    d6mu23 = d6mu**(2./3.)
 
    rs    = 1.388e-02 / d6mu**(1./3.) 
    
    x = 0.45641 * dlog( gamma ) - 1.31636
    x2 = x*x
    x3 = x2 * x
    
    asum = a(1) + a(2) * x + a(3) * x2 + a(4) * x3
    bsum = b(1) + b(2) * x + b(3) * x2
    csum = c(1) + c(2) * x + c(3) * x2
    dsum = d(1) + d(2) * x + d(3) * x2 + d(4) * x3
    esum = e(1) + e(2) * x + e(3) * x2
    fsum = f(1) + f(2) * x + f(3) * x2
    
    sminus1 = asum * ( 1. + bsum * rs + csum * rs*rs )
    splus1  = dsum * ( 1. + esum * rs + fsum * rs*rs )
    
    cstd6  = 1.018 * d6mu23 
    tmpsum = 1. + cstd6
    s = sminus1 - splus1 * cstd6 / tmpsum
    
    opacity = 2.363e+17 * d6 * t8 / ( at * tmpsum * s )
    cndlhy = dlog10( opacity )
    
    return
  end function cndlhy
  
!************************************************************************

  function hlcnc(d,t)

    use hlcnc_data

! Don Lamb's fit to the Hubbard-Lampe carbon conductive opacities

    implicit double precision(a-h,o-z)
    
    real*8 :: d, t
    real*8, dimension(25) :: w

    w(1)=1.
    m=0 
    n=1 
    do jmax=1,4 
       n=n+1
       do j=1,jmax 
          m=m+1
          w(n)=d*w(m)
          n=n+1
       end do
       w(n)=t*w(m)
    end do
    do i=16,25
       w(i)=d*w(i-5) 
    end do
    hlcnc=0
    do i=1,25
       hlcnc=hlcnc+a(i)*w(i) 
    end do
    return
  end function hlcnc

!************************************************************************

  function hlcnhe(d,t)

! Gilles Fontaine's fit to the Hubbard-Lampe Helium conductive opacities

    use hlcnhe_data

    implicit double precision(a-h,o-z)

    real*8 :: d, t
    real*8, dimension(5) :: temp, rho
  
    sum = 0.
    temp(1) = 1.
    rho(1)  = 1.
    do i = 2,5 
       temp(i) = temp(i-1) * t
       rho(i)  = rho(i-1) * d
    end do
    do i=1,5
       do j=1,5
          sum = sum + a(i,j)*temp(j)*rho(i)
       end do
    end do
    hlcnhe = sum 
    return
  end function hlcnhe

!************************************************************************

  function hlcnh(d,t)

! Gilles Fontaine's fit to the Hubbard-Lampe Hydrogen conductive opacities

    use hlcnh_data

    implicit double precision(a-h,o-z)

    real*8 :: d, t
    real*8, dimension(5) :: temp, rho

    sum = 0.
    temp(1) = 1.
    rho(1)  = 1.
    do i = 2,5 
       temp(i) = temp(i-1) * t
       rho(i)  = rho(i-1) * d
    end do
    do i=1,5
       do j=1,5
          sum = sum + a(i,j)*temp(j)*rho(i)
       end do
    end do
    hlcnh = sum
    return
  end function hlcnh

!************************************************************************

end module ocon_functions
