module utils_subroutines

contains

!************************************************************************

  subroutine arline(a,b,c,ja,jb,kc,ip8,p)

!  array linear extrapolation.  use a and b to extrapolate to c, but
!  first interpolate a to the same mass mesh points as in b. 

    use flags, only: evoloutput

    implicit double precision(a-h,o-z)
 
    real*8, dimension(600,1) ::  a,b,c
    real*8, dimension(7) :: aj

9     format (21h a does not contain b//)
    
    ip8=1
    check=(b(1,1)-a(1,1))/a(1,1)-1.d-07
 
    if( check .le. 0. )then
       l=2
       do  j=1,jb
7         if ( l .le. ja ) then
             if ( b(j,1) .le. a(l,1) ) then
                sfrac=(b(j,1)-a(l,1))/(a(l,1)-a(l-1,1))
                do k=1,kc
                   aj(k)=a(l,k)+(a(l,k)-a(l-1,k))*sfrac
                   c(j,k)=b(j,k)+(b(j,k)-aj(k))*p
                enddo
             else 
                l=l+1
                goto 7
             endif
          else
             
!*** if a does not contain b, then simply set c = b (no extrapolation) ***
             if ( evoloutput ) write(17,9)
             do m=j,jb
                do k=1,kc
                   c(m,k)=b(m,k)
                enddo
             enddo
             goto 116
          endif
       enddo
       goto 116
    else
       write(17,9)
       ip8=0
       goto 116
    endif
    
116 return
  end subroutine arline

!************************************************************************

  subroutine armove(x,y,j,k)

!  copy array y into array x
  
    implicit double precision(a-h,o-z)
  
    real*8, dimension(600,k):: x,y
  
    do  l=1,k
       do  i=1,j
          x(i,l)=y(i,l)
       enddo
    enddo
    
    return
  end subroutine armove

!****************************************************************************
!Some kind of integrator from numerical recipes

  SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)

    INTEGER :: nv
    real*8 ::  eps,hdid,hnext,htry,x,eps1,errmax,fact,h,red,scale,work, &
         wrkmin,xest,xnew,epsold
    real*8, dimension(nv) :: dydx,y,yscal
    integer , parameter :: NMAX=50, KMAXX=8, IMAX=KMAXX+1
    real*8, parameter :: SAFE1=.25,SAFE2=.7,REDMAX=1.e-5,REDMIN=.7, &
         TINY=1.e-30,SCALMX=.1
    INTEGER :: i,iq,k,kk,km,kmax,kopt
    integer, dimension(IMAX) :: nseq
    real*8, dimension(IMAX) :: a
    real*8, dimension(KMAXX,KMAXX) :: alf
    real*8, dimension(KMAXX) :: err
    real*8, dimension(NMAX) :: yerr,ysav,yseq
    logical :: reduct, first

    SAVE a,alf,kmax,kopt,nseq,xnew
    EXTERNAL derivs

    DATA nseq /2,4,6,8,10,12,14,16,18/

    first = .true.
    epsold = -1.
    if(eps.ne.epsold)then
       hnext=-1.e29
       xnew=-1.e29
       eps1=SAFE1*eps
       a(1)=nseq(1)+1
       do k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
       enddo
       do iq=2,KMAXX
          do k=1,iq-1
             alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
          enddo
       enddo
       epsold=eps
       do kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
       enddo
1      kmax=kopt
    endif
    h=htry
    do i=1,nv
       ysav(i)=y(i)
    enddo
    if(h.ne.hnext.or.x.ne.xnew)then
       first=.true.
       kopt=kmax
    endif
    reduct=.false.
2   do k=1,kmax
       xnew=x+h
       !        if(xnew.eq.x)pause 'step size underflow in bsstep'
       call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
       xest=(h/nseq(k))**2
       call pzextr(k,xest,yseq,y,yerr,nv)
       if(k.ne.1)then
          errmax=TINY
          do i=1,nv
             errmax=max(errmax,abs(yerr(i)/yscal(i)))
          enddo
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
       endif
       if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
             red=SAFE2/err(km)
             goto 3
          else if(k.eq.kopt)then
             if(alf(kopt-1,kopt).lt.err(km))then
                red=1./err(km)
                goto 3
             endif
          else if(kopt.eq.kmax)then
             if(alf(km,kmax-1).lt.err(km))then
                red=alf(km,kmax-1)*SAFE2/err(km)
                goto 3
             endif
          else if(alf(km,kopt).lt.err(km))then
             red=alf(km,kopt-1)/err(km)
             goto 3
          endif
       endif
    enddo
3   red=min(red,REDMIN)
    red=max(red,REDMAX)
    h=h*red
    reduct=.true.
    goto 2
4   x=xnew
    hdid=h
    first=.false.
    wrkmin=1.e35
    do kk=1,km
       fact=max(err(kk),SCALMX)
       work=fact*a(kk+1)
       if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
       endif
    enddo
    hnext=h/scale
    if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
       fact=max(scale/alf(kopt-1,kopt),SCALMX)
       if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
       endif
    endif
    return
  END SUBROUTINE bsstep

!************************************************************************
!Some kind of numerical subroutine (Fixed point iteration)
  SUBROUTINE xlocate(xx,n,x,j)

    INTEGER :: j,n
    real*8 :: x
    real*8, dimension(n) :: xx
    INTEGER :: jl,jm,ju

    jl=0
    ju=n+1
10  if(ju-jl.gt.1)then
       jm=(ju+jl)/2
       if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
       else
          ju=jm
       endif
       goto 10
    endif
    if(x.eq.xx(1))then
       j=1
    else if(x.eq.xx(n))then
       j=n-1
    else
       j=jl
    endif
    return
  END SUBROUTINE xlocate

!************************************************************************

  subroutine hunt(xx,n,x,jlo)

! Given an array XX of length N, and given a value X, returns a 
! value jlo such that x is between xx(jlo) and xx(jlo+1).  XX must
! be monotonic, either increasing or decreasing.

    implicit double precision(a-h,o-z)
    
    real*8, dimension(n) :: xx
    real*8 :: x
    logical :: ascnd
    
!*** Is table ascending (.TRUE.) or decending (.FALSE.)? ***
    ascnd = xx(n) .gt. xx(1)
!*** Input guess not useful, go straight to bisection ***
    if ( jlo .le. 0 .or. jlo .gt. n) then
       jlo = 0
       jhi = n + 1
       goto 3
    endif
!*** set hunting increment ***
    inc = 1
!*** hunt up. ***
    if ( x .ge. xx(jlo) .eqv. ascnd) then
1      jhi = jlo + inc
       if (jhi .gt. n) then
          jhi = n+1
       else if ( x .ge. xx(jhi) .eqv. ascnd) then
          jlo = jhi
          inc = inc + inc
          goto 1
       endif
    else
!*** hunt down ***
       jhi = jlo
2      jlo = jhi - inc
       if ( jlo .lt. 1 ) then
          jlo = 0
       else if ( x.lt.xx(jlo) .eqv.ascnd) then
          jhi = jlo
          inc= inc+inc
          goto 2
       endif
    endif
3   if ( jhi-jlo .eq. 1) goto 216
    jm = (jhi+jlo) / 2
    if ( x .gt. xx(jm) .eqv. ascnd) then
       jlo=jm
    else
       jhi = jm
    endif
    goto 3
216 return
  end subroutine hunt

!************************************************************************

  subroutine spline(x,y,n,yp1,ypn,y2)

! Given arrays X & Y (Y=f(X)) of length n containing a tabulated
! function, and given the values yp1 and ypn of the first 
! derivatives at the endpoints, this routine returns an array Y2
! of length n containing the second derivatives ofthe interpolating
! function at the tabulated points of X.
! If yp1 and/or ypn are 10^{30} or larger, the routine sets the
! BC's for a natural spline, with zero second derivatives at that
! boundary.
! Routine is from Numerical Recipies

    implicit double precision(a-h,o-z)
    
    integer, parameter :: nmax=10000
    real*8, dimension(n) :: x,y,y2
    real*8, dimension(nmax) :: u
    real*8 :: yp1, ypn
    
    if(yp1.gt.0.99e30) then
       y2(1)=0.
       u(1)=0.
    else
       y2(1)=-0.5
       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.
       y2(i)=(sig-1.)/p
       u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
            /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if(ypn.gt.0.99e30) then
       qn=0.
       un=0.
    else
       qn=0.5
       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do k=n-1,1,-1
       y2(k)=y2(k)*y2(k+1)+u(k)
    enddo
    
    return
  end subroutine spline
  
!************************************************************************

  subroutine splint(xa,ya,y2a,n,x,y)
    
! given the input arrays xa & ya (y=f(x)) of length n, tabulating
! a function, and given the spline coefficients y2a (from spline)
! and given an evaluation point x, this routine returns a cubic-spline
! interpolated value y.
! Routine is from Numerical Recipes

    use crash
    use flags, only: verbose

    implicit double precision(a-h,o-z)
    
    real*8, dimension(n) :: xa,ya,y2a
    real*8 :: x,y
    
!*** find place in table ***
    call hunt(xa,n,x,klo)
    khi = klo + 1
    h=xa(khi)-xa(klo)
    if (h.eq.0.)then
!*** xa points are not distinct ***
       if (verbose) write(*,*)  'bad xa input.'
       ihotspot = 1
       goto 126
    endif
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
!*** evaluate spline ***
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

126 return
  end subroutine splint

!************************************************************************

  subroutine fh4(x,yt,xt,y,dydx,indh)

! Lagrange interpolating polynomial

    implicit double precision(a-h,o-z)

    real*8, dimension(4) ::  yt,xt,xmt
    real*8 :: x,y,dydx

    a1=yt(1)/((xt(1)-xt(2))*(xt(1)-xt(3))*(xt(1)-xt(4)))
    a2=yt(2)/((xt(2)-xt(1))*(xt(2)-xt(3))*(xt(2)-xt(4)))
    a3=yt(3)/((xt(3)-xt(1))*(xt(3)-xt(2))*(xt(3)-xt(4)))
    a4=yt(4)/((xt(4)-xt(1))*(xt(4)-xt(2))*(xt(4)-xt(3)))
  
    do k=1,4
       xmt(k)=x-xt(k)
    end do
    y=a1*xmt(2)*xmt(3)*xmt(4)+a2*xmt(1)*xmt(3)*xmt(4) + &
         a3*xmt(1)*xmt(2)*xmt(4)+a4*xmt(1)*xmt(2)*xmt(3)
    if( indh .eq. 0) goto 136
    fr=xmt(2)/(xt(3)-xt(2)) 
    if(xmt(4)*xmt(3) .lt. 0.) goto 150
    dydx1=(xt(1)-xt(4))*a1*(xmt(2)+xmt(3)) &
         +(xt(2)-xt(4))*a2*(xmt(1)+xmt(3)) &
         +(xt(3)-xt(4))*a3*(xmt(2)+xmt(1))
    goto 200
150 dydx1=(yt(4)-yt(3))/(xt(4)-xt(3)) 
    fr=-xmt(4)/(xt(4)-xt(3))
200 if(xmt(1)*xmt(2) .lt. 0.) goto 250
    dydx2 = (xt(2)-xt(1))*a2*(xmt(3)+xmt(4)) + &
         (xt(3)-xt(1))*a3*(xmt(2)+xmt(4)) + &
         (xt(4)-xt(1))*a4*(xmt(3)+xmt(2))
    goto 300
250 dydx2=(yt(2)-yt(1))/(xt(2)-xt(1)) 
    fr=-xmt(2)/(xt(2)-xt(1))
300 dydx=dydx2*fr+dydx1*(1.-fr)
136 return 

  end subroutine fh4

!************************************************************************
! Finds the zero of a function
  
  FUNCTION zbrent(func,x1,x2,tol)
  
    INTEGER :: iter
    real*8 :: zbrent,tol,x1,x2,func
    real*8 :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    EXTERNAL func
    integer, parameter :: ITMAX=100
    real*8, parameter :: EPS=3.d-8

    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    c=b
    fc=fb
    do iter=1,ITMAX
       if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
       endif
       if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2.*EPS*abs(b)+0.5*tol
       xm=.5*(c-b)
       if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          goto 146
       endif
       if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
             p=2.*xm*s
             q=1.-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
             q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
          else
             d=xm
             e=d
          endif
       else
          d=xm
          e=d
       endif
       a=b
       fa=fb
       if(abs(d) .gt. tol1) then
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif
       fb=func(b)
    enddo
    zbrent=b
146 return
  END FUNCTION zbrent

!****************************************************************************

  FUNCTION beta(z,w)
      
    real*8 :: beta,w,z

    beta=exp(gammln(z)+gammln(w)-gammln(z+w))
    
    return
  END FUNCTION beta

!****************************************************************************

  FUNCTION gammln(xx)
  
    real*8 :: gammln,xx
    INTEGER :: j
    real*8 :: ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
      
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
         24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         -.5395239384953d-5,2.5066282746310005d0/

    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
    end do
    gammln=tmp+log(stp*ser/x)
    
    return
  END FUNCTION gammln
!****************************************************************************
!Something related to hypergeometric functions

  FUNCTION hypgeo(a,b,c,z)
 
! Common blocks 
    use hypg
    use path

    complex*16 ::  hypgeo,a,b,c,z
    real*8, parameter :: EPS=1.d-6
    INTEGER :: nbad,nok
    complex*16 :: y(2)

    kmax=0
    if (real(z)**2+aimag(z)**2.le.0.25) then
       call hypser(a,b,c,z,hypgeo,y(2))
       goto 156
    else if (real(z).lt.0.) then
       z0=cmplx(-0.5,0.)
    else if (real(z).le.1.0) then
       z0=cmplx(0.5,0.)
    else
       z0=cmplx(0.,sign(0.5,aimag(z)))
    endif
    aa=a
    bb=b
    cc=c
    dz=z-z0
    call hypser(aa,bb,cc,z0,y(1),y(2))
    call odeint(y,4,0.,1.,EPS,.1,.0001,nok,nbad,hypdrv,bsstep)
    hypgeo=y(1)
156 return

  END FUNCTION hypgeo

!****************************************************************************
!Hypergeometric series obviusly related to hypergeometric funtions

  SUBROUTINE hypser(a,b,c,z,series,deriv)

    INTEGER :: n
    complex*16 :: a,b,c,z,series,deriv,aa,bb,cc,fac,temp
 
    deriv=cmplx(0.,0.)
    fac=cmplx(1.,0.)
    temp=fac
    aa=a
    bb=b
    cc=c
    do n=1,1000
       fac=((aa*bb)/cc)*fac
       deriv=deriv+fac
       fac=fac*z/n
       series=temp+fac
       if (series.eq.temp) goto 166
       temp=series
       aa=aa+1.
       bb=bb+1.
       cc=cc+1.
    enddo
 !      pause 'convergence failure in hypser'
166 return
  END SUBROUTINE hypser

!****************************************************************************

  SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)

    use path

    INTEGER :: nbad,nok,nvar,i,nstp
    complex*16 :: ystart(nvar)
    real*8 :: eps,h1,hmin,x1,x2,h,hdid,hnext,x,xsav
    real*8, dimension(NMAX) :: dydx,y,yscal
    real*8, parameter :: TINY=1.d-30
    
    EXTERNAL derivs,rkqs
    
    x=x1
    h=sign(h1,x2-x1)
    nok=0
    nbad=0
    kount=0
    do  i=1,nvar
       y(i)=ystart(i)
    enddo
    if (kmax.gt.0) xsav=x-2.*dxsav
    do  nstp=1,MAXSTP
       call derivs(x,y,dydx)
       do i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
       enddo
       if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
             if(kount.lt.kmax-1)then
                kount=kount+1
                xp(kount)=x
                do i=1,nvar
                   yp(i,kount)=y(i)
                enddo
                xsav=x
             endif
          endif
       endif
       if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
       call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
       if(hdid.eq.h)then
          nok=nok+1
       else
          nbad=nbad+1
       endif
       if((x-x2)*(x2-x1).ge.0.)then
          do i=1,nvar
             ystart(i)=y(i)
          enddo
          if(kmax.ne.0)then
             kount=kount+1
             xp(kount)=x
             do i=1,nvar
                yp(i,kount)=y(i)
             enddo
          endif
          goto 176
       endif
!        if(abs(hnext).lt.hmin) pause
!     *'stepsize smaller than minimum in odeint'
       h=hnext
    enddo
!      pause 'too many steps in odeint'
176 return

  END SUBROUTINE odeint

!****************************************************************************

  SUBROUTINE hypdrv(s,y,dyds)

    use hypg

    real*8 :: s
    complex*16 :: y(2),dyds(2),z

    z=z0+s*dz
    dyds(1)=y(2)*dz
    dyds(2)=((aa*bb)*y(1)-(cc-((aa+bb)+1.)*z)*y(2))*dz/(z*(1.-z))
    return

  END SUBROUTINE hypdrv

!****************************************************************************

  SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)

    INTEGER :: nstep,nvar, i, n
    real*8 :: htot, xs, h, h2, swap, x
    real*8, dimension(50) :: ym, yn
    real*8, dimension(nvar) :: dydx, y, yout

    EXTERNAL derivs
     
    h=htot/nstep
    do i=1,nvar
       ym(i)=y(i)
       yn(i)=y(i)+h*dydx(i)
    enddo
    x=xs+h
    call derivs(x,yn,yout)
    h2=2.*h
    do n=2,nstep
       do i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
       enddo
       x=x+h
       call derivs(x,yn,yout)
    enddo
    do i=1,nvar
       yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
    enddo
    return

  END SUBROUTINE mmid

!***********************************************************************

  SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)

    INTEGER :: iest,nv,j,k1
    integer, parameter :: IMAX=13,NMAX=50
    real*8 :: xest,delta,f1,f2,q
    real*8, dimension(nv) :: dy,yest,yz
    real*8, dimension(NMAX) :: d
    real*8, dimension(NMAX,IMAX) :: qcol
    real*8, dimension(IMAX) :: x
  
    SAVE qcol,x

    x(iest)=xest
    do j=1,nv
       dy(j)=yest(j)
       yz(j)=yest(j)
    enddo
    if(iest.eq.1) then
       do j=1,nv
          qcol(j,1)=yest(j)
       enddo
    else
       do j=1,nv
          d(j)=yest(j)
       enddo
       do k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do j=1,nv
             q=qcol(j,k1)
             qcol(j,k1)=dy(j)
             delta=d(j)-q
             dy(j)=f1*delta
             d(j)=f2*delta
             yz(j)=yz(j)+dy(j)
          enddo
       enddo
       do j=1,nv
          qcol(j,iest)=dy(j)
       enddo
    endif
    return

  END SUBROUTINE pzextr

!************************************************************************

  function isnan(f)

    implicit none

    real*8 :: f
    logical :: isnan

    if (f.eq.f) then
       isnan=.FALSE.
    else
       isnan=.TRUE.
    endif
    
    return
    
  end function isnan

!************************************************************************

  subroutine dali(x,arg,val,y)

! four-point lagrange interpolation via interpolating polynomial.

    implicit double precision(a-h,o-z)
    
    real*8 :: x
    real*8, dimension(4) :: arg,val 
 
    do j=2,4
       iend=j-1
       do i=1,iend
          h=arg(i)-arg(j)
          val(j)=(val(i)*(x-arg(j))-val(j)*(x-arg(i)))/h
       enddo
    enddo
    y=val(4)
 
    return
  
  end subroutine dali

!****************************************************************************
! Home made Runge-Kutta routines
!****************************************************************************

  subroutine rkfcow(rt,y,yp)

    use misc
    use dmisc
    use splinq

    implicit double precision(a-h,o-z)

    real*8, dimension(2) :: y, yp

    call splintl(rt,vq)
    yp(1)=y(1)*(vq(2)-3.+lindex)+y(2)*(lhat*vq(1)/eigt-vq(2))
    yp(2)=y(1)*(eigt/vq(1)-vq(3))+y(2)*(1.-vq(4)+vq(3)+lindex)
    do i=1,2
       yp(i)=yp(i)*vq(5)
    end do
    
    return
    
  end subroutine rkfcow

!************************************************************************

  subroutine rkfliq(rt,y,yp)

    use misc
    use dmisc
    use splinq
    use ray

    implicit double precision(a-h,o-z)

    real*8, dimension(8) :: y, yp
    real*8, dimension(4) :: z

    call splintl(rt,vq)
    yp(1)=y(1)*(vq(2)-3.+lindex)+y(2)*(lhat*vq(1)/eigt-vq(2))+y(3)*vq(2)
    yp(2)=y(1)*(eigt/vq(1)-vq(3))+y(2)*(1.-vq(4)+vq(3)+lindex)-y(3)*vq(3)
    yp(3)=y(3)*(1.-vq(4)+lindex)+y(4) 
    yp(4)=y(1)*vq(4)*vq(3)+y(2)*vq(4)*vq(2)+y(3)*(lhat-vq(4)*vq(2)) &
         +y(4)*(lindex-vq(4))
    do i=1,4
       yp(i)=vq(5)*yp(i)
    end do
    if (iray.eq.0) goto 186
    rrt=1./vq(6)**lindex
    do i=1,4
       z(i)=y(i)*rrt 
    end do
    rhot=vq(4)*vq(1)/pi4/grav
    yp(5)=rhot*vq(6)**5*(z(1)**2+lhat*z(2)**2*(vq(1)/eigt)**2)
    rrt=rhot*vq(1)*vq(6)**5 
    yp(6)=rrt*vq(2)*(z(2)-z(3))**2
    yp(7)=rrt*vq(3)*z(1)**2 
    yp(8)=-rrt*(z(4)+(l+1.)*z(3))**2/vq(4)
    do i=5,8
       yp(i)=yp(i)*vq(5)
    end do

186 return
     
  end subroutine rkfliq

!************************************************************************

  subroutine rkf(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)

    implicit double precision(a-h,o-z)

    real*8, dimension(8) :: y, yp, f1, f2, f3, f4, f5
    real*8, dimension(51) :: work
    integer, dimension(5) :: iwork
    external f

    k1m=neqn+1
    k1=k1m+1
    k2=k1+neqn
    k3=k2+neqn
    k4=k3+neqn
    k5=k4+neqn
    k6=k5+neqn
    do i=1,neqn
       yp(i)=work(i)
    end do
    h=work(k1m)
    do j=1,neqn
       f1(j)=work(k1m+j)
    end do
    do k=1,neqn
       f2(k)=work(k2-1+k)
    end do
    do j=1,neqn
       f3(j)=work(k3-1+j)
    end do
    do k=1,neqn
       f4(k)=work(k4-1+k)
    end do
    do k=1,neqn
       f5(k)=work(k5-1+k)
    end do
    savre=work(k6)
    savae=work(k6+1)
    call rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,yp,h,f1,f2,f3,f4,f5,savre, &
         savae,iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))
    work(k1m)=h
    work(k6)=savre
    work(k6+1)=savae
    do k=1,neqn
       work(k)=yp(k)
    end do
    do j=1,neqn
       work(k1m+j)=f1(j)
    end do
    do k=1,neqn
       work(k2-1+k)=f2(k)
    end do
    do k=1,neqn
       work(k3-1+k)=f3(k)
    end do
    do k=1,neqn
       work(k4-1+k)=f4(k)
    end do
    do k=1,neqn
       work(k5-1+k)=f5(k)
    end do
    return
     
  end subroutine rkf

!************************************************************************

  subroutine rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,yp,h,f1,f2,f3,f4,f5, &
       savre,savae,nfe,kop,init,jflag,kflag)

    implicit double precision(a-h,o-z)
    
    logical :: hfaild,output
    real*8, dimension(8) :: y,yp,f1,f2,f3,f4,f5
    external f

    u26=2.0d-13
    remin=1.0d-12
    maxnfe=3000
    if (neqn .lt. 1) goto 10
    if ((relerr .lt. 0.)  .or.  (abserr .lt. 0.)) goto 10
    mflag=iabs(iflag)
    if ((mflag .ge. 1) .and. (mflag .le. 7)) goto 20
10  iflag=7
    goto 196
20  if (mflag .eq. 1) goto 50
    if (t .eq. tout) goto 10
    if(mflag .ne. 2) goto 25
    if (init .eq. 0) goto 45
    if (kflag .eq. 3) goto 40
    if ((kflag .eq. 4) .and.  (abserr .eq. 0.)) goto 30
    if ((kflag .eq. 5)  .and. (relerr .le. savre) .and. &
         (abserr .le. savae)) goto 30
    goto 50
25  if (iflag .eq. 3) goto 40
    if ((iflag .eq. 4) .and. (abserr .gt. 0.)) goto 45
30  continue
    print *, "mysterious stop in subroutine rkfs"
    stop
40  nfe=0
    if (mflag .eq. 2) goto 50
45  iflag=jflag
50  jflag=iflag
    kflag=0
    savre=relerr
    savae=abserr
    rer=max(relerr,remin) 
    dt=tout-t
    if (mflag .eq. 1) goto 60
    if (init .eq. 0) goto 65
    goto 80
60  init=0
    kop=0
    a=t 
    call f(a,y,yp)
    nfe=1
    if (t .ne. tout) goto 65
    iflag=2
    goto 196
65  init=1
    ymax=0.
    ypn=0.
    do k=1,neqn
       ypn=max(abs(yp(k)),ypn)
       ymax=max(abs(y(k)),ymax)
    end do
    etn=rer*ymax+abserr
    h=dabs(dt)
    if(etn.ge.ypn*h**5) goto 80
    h=max((etn/ypn)**0.2,u26*max(dabs(t),h)) 
80  h=sign(h,dt)
    if (dabs(h).ge.dabs(dt)) kop=kop+1
    if (kop.ne.100) goto 85
    iflag=6
    goto 196
85  if (dabs(dt).gt.u26*dabs(t))goto 95
    do k=1,neqn
       y(k)=y(k)+dt*yp(k)
    end do
    a=tout
    call f(a,y,yp)
    nfe=nfe+1
    goto 300
95  output=.false.
    scale=2./rer
    ae=scale*abserr
100 hfaild=.false.
    hmin=u26*dabs(t)
    dt=tout-t
    if (dabs(dt).ge.2.*dabs(h)) goto 200
    if (dabs(dt).gt.dabs(h)/0.9) goto 150
    output=.true. 
    h=dt
    goto 200
150 h=0.5*dt
200 if (nfe.le.maxnfe) goto 220
    iflag=3
    kflag=3
    goto 196
220 continue
    call fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,f1)
    nfe=nfe+5
    eeoet=0.
    do k=1,neqn
       et=dabs(dble(y(k)))+dabs(dble(f1(k)))+ae
       if (et.gt.0.) goto 240
       iflag=4
       kflag=4
       goto 196
240    ee=dabs((-2090.d0*yp(k)+(21970.*f3(k)-15048.*f4(k)))+ &
            (22528.*f2(k)-27360.*f5(k)))
       eeoet=max(eeoet,ee/et)
    end do
    esttol=dabs(h)*eeoet*scale/752400. 
    if (esttol.le.1.) goto 260
    hfaild=.true. 
    output=.false.
    s=0.1
    if (esttol.lt.59049.)s=0.9/esttol**0.2
    h=s*h
    if (dabs(h).gt.hmin) goto 200
    iflag=5
    kflag=5
    goto 196
260 t=t+h
    do  k=1,neqn
       y(k)=f1(k)
    end do
    a=t 
    call f(a,y,yp)
    nfe=nfe+1
    if (hfaild) goto 290
    s=5.
    if (esttol.gt.1.889568e-04) s=0.9/esttol**0.2
    h=sign(max(s*dabs(h),hmin),h)
290 if (output) goto 300
    if (iflag.gt.0) goto 100
    iflag=-2
    goto 196
300 t=tout
    iflag=2

196 return
    
  end subroutine rkfs

!************************************************************************

  subroutine splintl(xint,yout)
    
    use puls
    use savingl
    use setup
    use mfac
    use yfacs
    use rs

!  interpolates for liquid 

    implicit double precision(a-h,o-z)
    
    real*8, dimension(6) :: yout
    real*8, dimension(650) :: y6
    real*8, dimension(4,650) :: c1, c2, c3, c4, c5, c6

! compute spline coefficients for a model if we haven't done so

    if (isetup.eq.1)then
       do i=1,m
          y6(i)=r(i)
       end do
       call consl(x,y1,m,c1)
       call consl(x,y2,m,c2)
       call consl(x,y3,m,c3)
       call consl(x,y4,m,c4)
       call consl(x,y5,m,c5)
       call consl(x,y6,m,c6)
       isetup=0
    endif
    mm=m-1
    if (xint.ge.x(1).and.xint.lt.x(m)) goto 40 
    if (xint .lt. (1.+1.e-8)*x(1)) goto 110
    if (xint .ge. x(m)) goto 20
    k=1 
    goto 70
20  if (xint .gt. (1.+1.e-6)*x(m)) goto 110
    k=mm
    goto 70
40  il=1
    ir=m
50  k=il+((ir-il)/2)
    if (xint.ge.x(k)) goto 60
    ir=k
    goto 50
60  if (xint.lt.x(k+1)) goto 70
    il=k
    goto 50
70  continue
    x1=x(k+1)-xint
    xx=xint-x(k)
    x12=x1*x1
    xx2=xx*xx
    yout(1)=x1*(c1(1,k)*x12 +c1(3,k))+xx*(c1(2,k)*xx2 +c1(4,k))
    yout(2)=x1*(c2(1,k)*x12 +c2(3,k))+xx*(c2(2,k)*xx2 +c2(4,k))
    yout(3)=x1*(c3(1,k)*x12 +c3(3,k))+xx*(c3(2,k)*xx2 +c3(4,k))
    yout(4)=x1*(c4(1,k)*x12 +c4(3,k))+xx*(c4(2,k)*xx2 +c4(4,k))
    yout(5)=x1*(c5(1,k)*x12 +c5(3,k))+xx*(c5(2,k)*xx2 +c5(4,k))
    yout(6)=x1*(c6(1,k)*x12 +c6(3,k))+xx*(c6(2,k)*xx2 +c6(4,k))
    goto 216
110 continue
    goto 216
700 continue

    print *, "mysterious stop in subroutine spintl"
    stop 
216 return
  end subroutine splintl

!************************************************************************

  subroutine consl(x,y,m,c)

!  interpolates for liquid 

    implicit double precision(a-h,o-z)

    real*8, dimension(4,m) :: c
    real*8, dimension(m) :: x, y
    real*8, dimension(650,3) :: a
    real*8, dimension(650) :: d, b, z, p

    mm=m-1
    do k=1,mm
       d(k)=x(k+1)-x(k)
       p(k)=d(k)/6.
       z(k)=(y(k+1)-y(k))/d(k)
    end do
    do k=2,mm
       b(k)=z(k)-z(k-1)
    end do
    a(1,2)=-1.-d(1)/d(2)
    a(1,3)=d(1)/d(2)
    a(2,3)=p(2)-p(1)*a(1,3) 
    a(2,2)=2.*(p(1)+p(2))-p(1)*a(1,2) 
    a(2,3)=a(2,3)/a(2,2)
    b(2)=b(2)/a(2,2)
    do k=3,mm
       a(k,2)=2.*(p(k-1)+p(k))-p(k-1)*a(k-1,3)
       b(k)=b(k)-p(k-1)*b(k-1)
       a(k,3)=p(k)/a(k,2) 
       b(k)=b(k)/a(k,2)
    end do
    q=d(m-2)/d(m-1)
    a(m,1)=1.+q+a(m-2,3)
    a(m,2)=-q-a(m,1)*a(m-1,3)
    b(m)=b(m-2)-a(m,1)*b(m-1)
    z(m)=b(m)/a(m,2)
    mn=m-2
    do i=1,mn
       k=m-i
       z(k)=b(k)-a(k,3)*z(k+1)
    end do
    z(1)=-a(1,2)*z(2)-a(1,3)*z(3)
    do k=1,mm
       q=1./(6.*d(k))
       c(1,k)=z(k)*q
       c(2,k)=z(k+1)*q
       c(3,k)=y(k)/d(k)-z(k)*p(k)
       c(4,k)=y(k+1)/d(k)-z(k+1)*p(k)
    end do

    return
     
  end subroutine consl

!************************************************************************

  subroutine fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,s)

    implicit double precision(a-h,o-z)
      
    real*8, dimension(8) :: y, yp, f1, f2, f3, f4, f5, s

    ch=0.25*h
    do k=1,neqn
       f5(k)=y(k)+ch*yp(k)
    end do
    call f(t+0.25*h,f5,f1)
    ch=0.09375*h
    do  k=1,neqn
       f5(k)=y(k)+ch*(yp(k)+3.*f1(k))
    end do
    call f(t+0.375*h,f5,f2) 
    ch=h/2197.
    do k=1,neqn
       f5(k)=y(k)+ch*(1932.*yp(k)+(7296.*f2(k)-7200.*f1(k)))
    end do
    call f(t+12./13.*h,f5,f3)
    ch=h/4104.
    do k=1,neqn
       f5(k)=y(k)+ch*((8341.*yp(k)-845.*f3(k))+(29440.*f2(k)-32832.*f1(k)))
    end do
    call f(t+h,f5,f4)
    ch=h/20520.
    do k=1,neqn
       f1(k)=y(k)+ch*((-6080.*yp(k)+(9295.*f3(k)-5643.*f4(k)))+ &
            (41040.*f1(k)-28352.*f2(k)))
    end do
    call f(t+0.5*h,f1,f5)
    ch=h/7618050. 
    do k=1,neqn
       s(k)=y(k)+ch*((902880.*yp(k)+(3855735.*f3(k)-1371249.*f4(k)))+ &
            (3953664.*f2(k)+277020.*f5(k)))
    end do

    return

  end subroutine fehl

!****************************************************************************

  subroutine check_time(newtime)
    
    use crash

    character(10) :: date, time, zone
    integer, dimension(8) :: values

    call date_and_time(date,time,zone,values) 
!values(3) = day, values(5) = hours, values(6) = minutes, values(7) = seconds
    newtime = 24*3600*values(3) + 3600*values(5) + 60*values(6) + values(7)
    return
   
  end subroutine check_time
    
!****************************************************************************
! Numerical integration by Trapezoid rule. (Called by QSIMP)
! From Numerical Recipies
!****************************************************************************

  subroutine trapzd(func,a,b,s,n) 
        
    use flags, only: verbose

    implicit double precision(a-h,o-z)
    
    real*8 :: func
    logical :: debug

    debug = .false.

    if (n.eq.1) then
       s = 0.5*(b-a)*(func(a)+func(b))
       it = 1
    else
       tnm = it 
       del = (b-a)/tnm
       x = a + 0.5 * del
       sum = 0. 
       do j = 1,it
          dum = func(x)
          sum = sum + dum 
          x = x + del
31        format(1x,2f10.5)
       end do
       s = 0.5 * (s + (b-a) * sum/tnm )
       it = 2 * it
       if (debug .and. verbose) write(*,*) 'trapzd: best guess = ',s 
    endif
    return
    
  end subroutine trapzd

!************************************************************************
! Numerical integration via Simpson's rule. (From Numerical Recipies)
!************************************************************************

  subroutine qsimp(func,a,b,s)
    
    use flags, only: verbose
    implicit double precision(a-h,o-z)
    
    external func
    real*8, parameter :: eps = 1.d-6
    integer, parameter :: jmax = 25

    ost = -1.e-30
    os =  -1.e-30

    do j = 1,jmax
       call trapzd(func,a,b,st,j)
       s = (4.*st-ost)/3. 
       if ( dabs(s-os) .lt. eps*dabs(os)) goto 226
       os = s
       ost = st 
    end do

    if (verbose) then
       write(*,*) 'qsimp: too many steps'
       write(*,*) 'exiting.'
       stop
    end if
226 return

  end subroutine qsimp

!************************************************************************
  subroutine dosim(result,f,h,n)
        
    use lager
    use delta
    use resvalue
    use flags

    implicit double precision(a-h,o-z)

    real*8, dimension(n) :: f,h
    real*8, dimension(650) :: fpp,res1,res2
    
    do ii=1,650
       fpp(ii) = 0.0d0
    end do
    
    sum=0.
    nm=n-1
    sumc=0.
    res1(1)=0.0
    res2(1)=0.0
    res2(2)=0.0
    do j=2,n
       sum=sum+h(j)*(f(j)+f(j-1))
       res1(j)=h(j)*(f(j)+f(j-1))/2.
    end do
    fpp(2)=2.*((f(3)-f(2))/h(3) - (f(2)-f(1))/h(2))/(h(2)+h(3))
    trapez=sum/2. 
    do j=3,nm
       fpp(j)=2.*((f(j+1)-f(j))/h(j+1)-(f(j)-f(j-1))/h(j))/(h(j)+h(j+1))
       dp=h(j)+h(j+1)
       dm=h(j)+h(j-1)
       res2(j)=(dm*fpp(j)+dp*fpp(j-1))*h(j)**3/(dm+dp)
       sumc=sumc+(dm*fpp(j)+dp*fpp(j-1))*h(j)**3/(dm+dp)
    end do
    endcor=(fpp(2)+(h(2)+h(3))*(fpp(2)-fpp(3))/(h(2)+2.*h(3)+h(4)))*h(2)**3
    endcor=endcor+(fpp(n-1)+(h(n)+h(n-1))*(fpp(n-1)-fpp(n-2))/ &
         (h(n)+2.*h(n-1)+h(n-2)))*h(n)**3
    sumc=sumc+endcor
    result=trapez-sumc/12.
    do j=2,n
       res(j)=res(j-1)+res1(j)-res2(j)/12.
    enddo
1002  format(20x,'endcor',e12.5,5x,'ekin',e12.5)
    if ( pulsoutput ) write(7,1002) endcor, result
    
    return
  end subroutine dosim

!************************************************************************

  SUBROUTINE TRIP (N,A,B,C,Y,Z,INT)
    IMPLICIT REAL*8 (A-H,O-Z)
    real*8, dimension(1) :: A,B,C,Y
    real*8, dimension(650) :: Z


! ARITHMETIC STATEMENT FUNCTION USED TO LOCATE ENTRIES IN ARRAY Z.

    II(INDEX)=(INDEX-1)*INT+1

! GAUSSIAN ELIMINATION

    BN = B(N)
    YN = Y(N)
    V = C(N)
    Y(1) = Y(1)/B(1)
    A(1) = A(1)/B(1)
    B(1) = C(1)/B(1)
    NM2 = N-2
    DO J=2,NM2
       DEN = B(J)-A(J)*B(J-1)
       B(J) = C(J)/DEN
       Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
       A(J) = -A(J)*A(J-1)/DEN
       BN = BN-V*A(J-1)
       YN = YN-V*Y(J-1)
       V = -V*B(J-1)
    end DO
    DEN = B(N-1)-A(N-1)*B(N-2)
    B(N-1) = (C(N-1)-A(N-1)*A(N-2))/DEN
    Y(N-1) = (Y(N-1)-A(N-1)*Y(N-2))/DEN
    BN = BN-V*A(N-2)
    YN = YN-V*Y(N-2)
    V = A(N)-V*B(N-2)
! BACK SUBSTITUTION
    IIN = II(N)
    Z(IIN) = (YN-V*Y(N-1))/(BN-V*B(N-1))
    IIN2 = II(N-1)
    Z(IIN2) = Y(N-1)-B(N-1)*Z(IIN)
    NM1 = N-1
    IN = II(N)
    DO J=2,NM1
       K = N-J
       IK = II(K)
       IKT = IK+INT
    Z(IK) = Y(K)-B(K)*Z(IKT)-A(K)*Z(IN)

    end DO
    RETURN

  END SUBROUTINE TRIP

!************************************************************************

  SUBROUTINE TERP1 (N,X,F,W,Y,INT,TAB,ITAB)
    IMPLICIT REAL*8(A-H,O-Z)
    real*8, dimension(1) :: X
    real*8, dimension(650) :: F,W
    real*8, dimension(3) :: TAB
    integer, dimension(3) :: ITAB

    CALL SEARCH (Y,X,N,I)
    CALL INTERPS (N,X,F,W,Y,I,INT,TAB,ITAB)
    RETURN
  END SUBROUTINE TERP1

!************************************************************************
  SUBROUTINE SEARCH (XBAR,X,N,I)
    IMPLICIT REAL*8 (A-H,O-Z)
    real*8, dimension(N) :: X
    DATA B/.69314718D0/

! IF XBAR IS OUTSIDE RANGE OF X TABLE EXTRAPOLATE

    IF (XBAR .GT. X(2)) GO TO 101
    I = 1
    RETURN
101 CONTINUE
    IF (XBAR .LT. X(N-1)) GO TO 102
    I = N-1
    RETURN
102 CONTINUE

! FIND MAXIMUM POWER OF TWO LESS THAN N

    M = INT((DLOG(DFLOAT(N)))/B)
    I = 2**M
    IF (I .GE. N) I = I/2
    K = I
    NM1 = N-1

! CONDUCT BINARY SEARCH.

103 CONTINUE
    K = K/2
    IF (XBAR .GE. X(I)) GO TO 104
    I = I-K
    GO TO 103
104 CONTINUE
    IF (XBAR .LE. X(I+1)) RETURN
    I = MIN0(I+K,NM1)
    GO TO 103
  END SUBROUTINE SEARCH

!************************************************************************

  SUBROUTINE INTERPS (N,X,F,W,Y,I,INT,TAB,ITAB)
    IMPLICIT REAL*8 (A-H,O-Z)
    real*8, dimension(1) :: X
    real*8, dimension(650) :: F, W
    real*8, dimension(3) :: TAB
    integer, dimension(3) :: ITAB

! ARITHMETIC STATEMENT FUNCTION USED TO LOCATE ENTRIES IN F AND W ARRAYS

    II(INDEX)=(INDEX-1)*INT+1

! PERFORM INTERPOLATION OR EXTRAPOLATION

    FLK = X(I+1)-X(I)
    FLP = X(I+1)-Y
    FL0 = Y-X(I)
    I0 = II(I)
    IP = I0+INT
    IF (ITAB(1) .NE. 0) GO TO 101
    GO TO 102
101 CONTINUE

! CALCULATE F(Y)

    A = (W(I0)*FLP**3+W(IP)*FL0**3)/(6.D0*FLK)
    B = (F(IP)/FLK-W(IP)*FLK/6.D0)*FL0
    C = (F(I0)/FLK-W(I0)*FLK/6.D0)*FLP
    TAB(1) = A+B+C

102 CONTINUE
    IF (ITAB(2) .NE. 0) GO TO 103
    GO TO 104
103 CONTINUE

! CALCULATE FIRST DERIVATIVE AT Y
    A = (W(IP)*FL0**2-W(I0)*FLP**2)/(2.D0*FLK)
    B = (F(IP)-F(I0))/FLK
    C = (W(I0)-W(IP))*FLK/6.D0
    TAB(2) = A+B+C
104 CONTINUE
    IF (ITAB(3) .NE. 0) GO TO 105
    GO TO 106
105 CONTINUE

! CALCULATE SECOND DERIVATIVE AT Y

    TAB(3) = (W(I0)*FLP+W(IP)*FL0)/FLK
106 CONTINUE
    RETURN
  END SUBROUTINE INTERPS

!************************************************************************
!computes 2-point logarithmic dervatives.

  subroutine dert(y,h,der,np)
  
    implicit double precision(a-h,o-z)
    real*8, dimension(650) :: y,h,der

! program allows y or h = 0 for n = 1 only
    if (y(1) .eq. 0.) go to 5
    if (h(1) .eq. 0.) go to 5
    ym=dlog(y(2)/y(1))
    if (dabs(ym) .lt. 1.d-6) ym=0.
    hm=dlog(h(2)/h(1))
5   npm2=np-2
    do n=1,npm2
       if(y(n+2).eq.0.0 .or. y(n+1).eq.0.0)then
          der(n+1) = 0.0
          go to 10
       endif
       yone=y(n+2)/y(n+1)
       if(yone.le.0.0)then
          yone=-yone
          nfac=-1
       endif
       yp=dlog(yone)
       if (dabs(yp) .lt. 1.e-6) yp=0.
       hp=dlog(h(n+2)/h(n+1))
       if(y(n) .eq. 0.) go to 9
       if(h(n) .eq. 0.) go to 9
       htot=dlog(h(n+2)/h(n))
       der(n+1)=(hm*yp/hp+hp*ym/hm)/htot 
       if(nfac.eq.-1)then
          der(n+1)=-der(n+1)
       endif
       ym=yp
       hm=hp 
      go to 10
9     der(2)=yp/hp
      ym=yp
      hm=hp 
10 end do
   if(y(1).eq.0.0 .and. y(2).eq.0.0)then
      der(2)=0.0
   endif
   der(np)=ym/hm 
   return
 end subroutine dert

!************************************************************************


end module utils_subroutines
