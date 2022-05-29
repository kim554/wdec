module pulse_subroutines

contains

!****************************************************************************

  subroutine openfi

    use outfile
    use flags, only: makeplots

    implicit none

    open(29,file='tape29.dat', status='unknown')
    open(35,file=fname,status='old')

    rewind(29)
  
    if ( makeplots) then
       open(8,file='output.dat',status='unknown')        !done
       open(18,file='tape18.dat',status='unknown')       !done but not right
       open(19,file='tape19.dat',status='unknown')       !done
       open(28,file='tape28.dat',status='unknown')       !done
       open(50,file='modelp.dat',status='unknown')       !done
       open(60,file='prop.dat',status='unknown')         !done
       open(61,file='deld.dat',status='unknown')         !done
       open(62,file='kaprt.dat',status='unknown')        !done
       open(63,file='xhe.dat',status='unknown')          !done
       open(64,file='lrad.dat',status='unknown')         !done
       open(65,file='temprp.dat',status='unknown')       !done
       open(66,file='chirt.dat',status='unknown')        !done
       open(67,file='epsrt.dat',status='unknown')        !done
       open(68,file='cpvtga1.dat',status='unknown')      !done
       open(71,file='struc.dat',status='unknown')        !done
       open(72,file='pform.dat',status='unknown')        !done
!       open(73,file='gam_xtal.dat',status='unknown')    !not used
       open(74,file='corsico.dat',status='unknown')      !done
!       open(75,file='reflection.dat',status='unknown')  !not reintroducing
       open(77,file='thorne.dat',status='unknown')       !done

       rewind(8)
       rewind(18)
       rewind(19)
       rewind(25)
       rewind(26)
       rewind(28)
       rewind(50)
       rewind(60)
       rewind(61)
       rewind(62)
       rewind(63)
       rewind(64)
       rewind(65)
       rewind(66)
       rewind(67)
       rewind(68)

    end if

    return

  end subroutine openfi

!****************************************************************************

  subroutine pulse

    ! pulsate the prepped model

    use utils_subroutines, only: dosim

    use misc_const
    use puls
    use misc
    use dmisc
    use ray
    use rs
    use ekint
    use setup
    use modect
    use eig1
    use ms
    use rot1
    use rot2
    use perds
    use fracxtal
    use modelp1
    use modelp2
    use cperiods
    use flags
    use flagtmp, only: mlver

    implicit double precision(a-h,o-z)

    real*8 :: ekintt
    integer, parameter :: sequen = 600502 !This is a filler and means
                                          !absolutely nothing.
    integer, dimension(3) :: modes
    integer, dimension(2) :: num

    if ( pulsoutput ) then
       !Open some files
       open(unit=1,file='discr.dat',status='unknown')
       open(unit=7,file='check.dat',status='unknown')
       open(unit=51,file='results.dat',status='unknown')
       open(unit=83,file='xir.dat',status='unknown')
       !       if ( kernels ) then
       !           open(unit=52,file='reflection2.dat',status='unknown')
       !          open(unit=77,file=fname2,status='unknown')
       !          open(unit=78,file='kernels/modelist.dat',status='unknown')
       !       end if
    end if

    !Re-assign some former common block variables
    np = np_tmp
    num = ncalc

    call init(model)
    
9006  format(10x,'white dwarf sequence = Tc',i6,1x,'ML',i1/,10x, &
           'temperature of model =',i7,' K')
9900  format(5x,'model number =',i5,3x,'age (yrs) = ',1pe12.5,3x, &
           'l/l sun = ',e12.5/,5x,'number of points in model =',0p,i4,3x, &
           'r/r sun = ',1pe12.5/,5x,' DE = +3.0, -3.0; Mod. Ledoux of &
           the BV freq, Reint opacs in WDXD')
9901 format(10x,'Nu_o',10x,'Pi_o',10x,'t dyn')
9902 format(2(7x,f8.4),6x,f8.4)
9903 format(5x,'l',2x,'k',2x,'period',3x, &
          'int. per',3x,'log K.E.',6x,'y3',3x,'nodes',1x,'= y1', &
          3x,'= y2',3x,'=plmode',3x,'rint 1'/)

    if ( pulsoutput ) then
       write(51,9006) sequen,mlver,idnint(teff)
       write(51,9900) model,age,llsun,np,rrsun
       write(51,9901)
       write(51,9902) gnu0,per0,tdyn
       write(51,9903)
    end if

    ! First, calculate surface discriminant only
    ! in cowling approximation to obtain period guesses

    delper=(permax-permin)/nper

    do lind=1,lmax
       l=float(lind)
       lindex=2.-l
       lhat=l*(l+1.)

       do idisc=1,nper+1
          periods(idisc)=permin+delper*(idisc-1)
          eig=(2.*pi/periods(idisc))**2
          eigt=eig
          iray=0
          call bump(discr(idisc))
          discr(idisc)=discr(idisc)/(r(nsurf)**lindex)
          if ( pulsoutput ) then
             write(1,290) periods(idisc),discr(idisc)
          end if
       end do  ! end search for periods
290    format(1pe14.6,1pe14.6) 

       call perg(lind,modes,nper,permax)

       if ( pulsoutput ) then
          write(1,*) '    period  yguess'
          do i=1,modes(lind)
             write(1,290) pguess(lind,i),yguess(i)
          end do
       end if

    end do   ! end search for modes

    ! begin loop over l 

    do li=1,lmax

       ! initialize index for calc_per array

       num(li) = 0

       limit=modes(li)
       l=float(li)
       lindex=2.-l
       lhat=l*(l+1.)

       ! begin loop over periods

       do nbnum=1,limit
100       nsearch=10
          numb=nbnum

          ! stop execution if period = 0.0

          if (pguess(li,numb).eq.0.0) then
             if(li.eq.lup)then
                if ( pulsoutput ) then
                   ! this replaces write(7,1001) in old version of code
                   write(7,*) '     calculation completed'
                   stop
                endif
             end if
             goto 7500
          endif

          ! continue with period guess
          eig=(2.*pi/pguess(li,numb))**2
          y3i=yguess(numb)
          if ( pulsoutput ) then
             write(7,1002) pguess(li,numb)
          end if
1002      format ('guessed period=',1pe14.6,4h sec)
          nconv=0
          iray=0

          do ntry=1,nsearch
             eigt=eig
             y3t=y3i
             b1 = 0.0
             b2 = 0.0
             call grind(b1,b2)
             eigt=(1.+eps)*eig
             dum1 = 0.0
             dum2 = 0.0
             call grind(dum1,dum2)
             deigb1=dum1-b1
             deigb2=dum2-b2
             eigt=eig
             y3t=(1.+eps)*y3i
             call grind (dum1,dum2)
             dy3b1=dum1-b1 
             dy3b2=dum2-b2 
             d2=(deigb1*dy3b2-deigb2*dy3b1)
             dume=b2*dy3b1-b1*dy3b2
             dumy=b1*deigb2-b2*deigb1
             deig=eps*eig*dume/d2
             dy3=eps*y3i*dumy/d2
             if ( pulsoutput ) then
                write(7,1003) ntry,eig,y3i,deig,dy3
                write(7,1020) b1,b2,yliq(3,nsurf),yliq(1,nsurf)
             end if
             adeig=dabs(deig/eig)
             ady3=dabs(dy3/y3i)
1003         format(6h ntry=,i3,5h eig=,1pe20.12,4h y3=,1pe14.6,6h deig=, &
                  1pe14.6,5h dy3=,1pe14.6)
1020         format(6x,8h and b1=,1pe12.4,4h b2=,1pe12.4,4h y3=,1pe12.4, &
                  4h y1= ,1pe12.4)
             if (adeig.lt.verg.and.ady3.lt.verg) then
                nconv=1 
             endif
             epseig=deig/eig
             epsy3=dy3/y3i
             if (adeig.gt.0.1) then
                epseig=0.1*deig/dabs(deig) 
             endif
             if (ady3.gt.2.0) then
                epsy3=2.*dy3*dabs(y3i)/dabs(dy3)/y3i-1.
             endif
             eig=eig*(1.+epseig)
             y3i=y3i*(1.+epsy3)

! if we converged to a solution, skip down a few lines to continue

             if (nconv.eq.1) goto 201
             
          end do
          if ( pulsoutput ) then
             ! replaces write(7,1004) in old code
             write(7,*) '  not converged'
          end if
          goto 7600

!  continue with converged period guess

201       eigt=eig
          y3t=y3i
          iray=1
          call grind (b1,b2)
          period=(2.*pi)/dsqrt(eig)
          if ( pulsoutput ) then
                ! replaces write(7,9879) in old code
             write(7,*) ' ***************************************************'
             write(7,1005) period,eig,y3i
1005         format (25h final values are,period=,1pe17.9,5h eig=,1pe17.9, &
                  4h y3=,1pe14.6)
             write(7,1007) b1,b2
1007         format (28h boundary conditions are b1=,1pe11.3,4h b2=,1pe11.3) 
             ! count nodes in y1 
          end if

          do i=1,nsurf
             y1m(i)=yliq(1,i)
             y2m(i)=yliq(2,i)
          end do
          nfine=nsurf
          call modeid
          if ( pulsoutput ) then
             write(7,1010) nodes1,nodes2,modep
          end if
1010      format(i3,' nodes in y1   ',i3," nodes in y2 "/ &
               '     phase diagram mode ',i3)

             ! copy integrated t,c,n,g to array rint(4,600)
          do  iint=1,nsurf
             do  jint=1,4
                rint(jint,iint)=rayy(jint,iint)
             end do
          end do
!**** integrated frequency and period ******************************
          eigtry=(rayy(2,nsurf)+rayy(3,nsurf)+rayy(4,nsurf))/rayy(1,nsurf)
          pertry=(2.*pi)/dsqrt(eigtry)
          if ( pulsoutput ) then
             write(7,1011) eigtry,pertry
          end if
1011      format('integrated sig**2=',e12.4,2x,'integrated period=',f14.10)
             
          call eigenf
          facke=r(nsurf)**(2.d0*lindex)/yliq(1,nsurf)**2
          flke = dlog10(rint(1,nsurf)*facke*2.d0*pi*eig)
          
!!$      if ( kernels ) then
!!$      call kern_write(modep,li,nsurf,x,period)
!!$   end if

!*** compute kinetic energy (ekin)
          do i=1,nsurf
             if(i.ne.1)then 
                h(i)=r(i)-rmi
             endif
14           rmi=r(i)
             aa=r(i)*yliq(1,i)
1012  format(I6,2F24.6,F24.15)
             if ( pulsoutput ) write(83,1012) i, rmi, aa, yliq(1,i)
             bb=(gg(i)/eigt)*yliq(2,i)
             front=rho(i)*r(i)**2
             f(i)=front*((aa*aa)+(bb*bb)*lhat)
          end do
          call dosim(ekin,f,h,nsurf)
          ! normalize kinetic energy to \xi_h/R=1 at the surface
          ekintt=2.*pi*eigt*ekin
          ekl=dlog10(ekintt)
          !**** compute rotational splitting coefficients
          call rotcoef

!**** print out results ************************************************
1006      format(20x,'kinetic energy',e12.5,' log kinetic energy',f10.5, &
               2x,f10.5)
1025      format(4x,i2,1x,i2,1x,f8.3,1x,f8.3,3x,f8.4,2x,e12.5,4x,i3,4x, &
               i3,4x,i3,1x,2f8.3)
1030      format('n, r/r* ,  y1  ,  y2  ,   tint , wint')
1021      format(i4,f10.6,2(1pe12.4),3x,2(1pe12.4))
1022      format ('y3  , y4  , t(i) , c(i) , n(i) , g(i) ')
1023      format (i4,6(1pe12.4))
          if ( pulsoutput ) then
             write(7,1006) ekintt,ekl,flke
             nl=li
             write(51,1025) nl,modep,period,pertry,ekl,y3i,nodes1,nodes2, &
                  modep, flke, clk
             write(7,1030)
             do i=1,nsurf
                wint=rint(2,i)+rint(3,i)+rint(4,i)
                write(7,1021) i,r(i)/r(nsurf),yliq(1,i), &
                     yliq(2,i),rint(1,i),wint
             end do
             write(7,1021)
             write(7,1022)
             do i=1,nsurf
                write(7,1023) i,yliq(3,i),yliq(4,i),rayy(1,i),rayy(2,i), &
                     rayy(3,i),rayy(4,i)
             end do
          end if

7600      continue
7500      continue
          
          num(li) = num(li) + 1
          calc_per(li,num(li)) = period
          
          write(102,*) li, period
          
       end do   ! end loop over periods
       
    end do   ! end loop over l
       
    ! Rename some common block variables before returning
    np = np_tmp
    num = ncalc
    
    return
       
  end subroutine pulse

!************************************************************************

  subroutine init(model)
    
    use flags
    use puls
    use misc_const
    use eig1
    use ms
    use misc
    use dmisc
    use rs
    use setup
    use mfac
    use yfacs
    use tape28a
    use tape28b
    use tape29
    use fracxtal
    
    implicit double precision(a-h,o-z)
    
    real*8 :: mr_tmp

! Reassign some common block variables
    x = xi2

    nsurf = nnsurf
    isetup = 1
    verg=1.00d-07 
    eps=1.00d-03
    p43=4.*pi/3.
    pi4=4.*pi
    grav = g
      
1003  format(6h verg=,1pe10.2,5h eps=,1pe10.2)
    if ( pulsoutput ) write(7,1003) verg,eps

    do i=1,nsurf
       gg(i) = ggrav(i+1)
       rho(i) = rrho(i+1)
       mr(i)= mr2(i+1)
       r(i)=rn(i+1)
    end do
    
    m = nsurf
      
!  if there are too many shells for the arrays, terminate execution.

    if(m.gt.650)then
       stop
    endif

! transfer contents of COMMON block

    do i=1,m
       y1(i) = yy1(i+1)
       y2(i) = voga1(i+1)
       y3(i) = yy3(i+1)
       y4(i) = u(i+1)
       y5(i) = yy5(i+1)
    end do

!***** the following lines added to mimic core crystallization
    nxtal=0
    do n=1,nsurf
       frac1=mr(n)
       frac2=mr(n+1)
       if ((frac2 .ge. fracm) .and. (frac1 .lt. fracm)) nxtal=n+1
    enddo
    
    if (nxtal .gt. 0) then
       fracm_new=mr(nxtal)
    else
       fracm_new = 0.0
    endif
    
    do n=1,nsurf-nxtal
       x(n)=x(nxtal+n)
       r(n)=r(nxtal+n)
       gg(n)=gg(nxtal+n)
       rho(n)=rho(nxtal+n)
       mr(n)=mr(nxtal+n)
       y1(n)=y1(nxtal+n)
       y2(n)=y2(nxtal+n)
       y3(n)=y3(nxtal+n)
       y4(n)=y4(nxtal+n)
       y5(n)=y5(nxtal+n)
    enddo
      
!    nsurf is now the number of uncrystallized shells
    nsurf=nsurf-nxtal
    m=m-nxtal
!***** the previous lines added to mimic core crystallization

    return
  end subroutine init

!************************************************************************

  subroutine perg(li,modes,nper,permax)

! this subroutine computes trial periods for the pulsation programs 
! from the discriminant values. It uses a simple linear interpolation
! between data points where the sign changes. 

    use perds
    use flags, only: pulsoutput

    implicit double precision(a-h,o-z)

    integer, dimension(3) :: modes

! load up Y3 guess array

    yguess(1) = -0.5
    yguess(2) = -0.02
    yguess(3) = -0.02
    yguess(4) = -0.02
    yguess(5) = -0.02
    yguess(101) = 0.0
    do i=6,18
       yguess(i) = -0.001
    end do
    do i=19,100
       yguess(i) = -0.0001
    end do

! start computing period guesses
    num = 0
    do i=1,nper-1

! if we have more periods than the array allows, 
! truncate and use what we have already.

       if(num.eq.100)then
          pguess(li,num+1) = 0.0
          modes(li) = num+1
          if ( pulsoutput ) then
             write(51,*) '  NOTE: There are more than 100 period guesses in &
                  this scan'
          end if
          goto 99
       endif

! pass through zero going down 
       if(discr(i).gt.0.0.and.discr(i+1).lt.0.0)then
          num = num + 1
          if(periods(i).lt.permax)then
             frac=discr(i)/(discr(i)+dabs(discr(i+1)))
             pguess(li,num)=periods(i)+(periods(i+1)-periods(i))*frac
          endif
       endif

! pass through zero going up

       if(discr(i).lt.0.0.and.discr(i+1).gt.0.0)then
          num = num + 1
          if(periods(i).lt.permax)then
             frac=discr(i)/(discr(i)-discr(i+1))
             pguess(li,num)=periods(i)+(periods(i+1)-periods(i))*frac
          endif
       endif
    end do
    pguess(li,num+1)=0.0
    modes(li)=num+1

99  return 
    
  end subroutine perg

!************************************************************************

  subroutine bump(b2)
 
! Subroutines
    use utils_subroutines, only : rkfcow, rkf, splintl

! Common blocks
    use misc_const
    use puls
    use misc
    use dmisc
    use rs
    use splinq
    use fracxtal
    
    implicit double precision(a-h,o-z)

    real*8, dimension(8) :: y
    real*8, dimension(51) :: work
    real*8, parameter :: small = 1.0d-02
    integer, dimension(5) :: iwork

! put in if statement to treat crystallized and uncrystallized state
    if (fracm_new .lt. small) then
       yliq(1,1)=1.0
    else
       yliq(1,1)=0.0
    endif
    yliq(2,1)=eigt/(l*grav*p43*rho(1))
    y(1)=yliq(1,1)
    y(2)=yliq(2,1)
    iflag=1
    neqn=2
4   nsurf1=nsurf-1
    do j=1,nsurf1
       k=j+1
       xstart=x(j)
       xfin=x(k)
       relerr=1.d-8
       abserr=dabs(y(1))
       t=dabs(y(2))
2      abserr=dmin1(abserr,t)
       abserr=1.d-8*abserr
       abserr=max(abserr,1.d-20)
       call rkf(rkfcow,neqn,y,xstart,xfin,relerr,abserr,iflag,work,iwork)
       iflag=1
       do i=1,2
          yliq(i,k)=y(i)
       end do
    end do
    rt=x(nsurf)
    call splintl(rt,vq)
    temp=1./vq(5)-1.
    b2=yliq(1,nsurf)*((4.+r(nsurf)*eigt/gg(nsurf))/temp-1.) &
         +(1.-lhat*gg(nsurf)/r(nsurf)/eigt/temp)*yliq(2,nsurf) 
    return

  end subroutine bump

!************************************************************************

  subroutine grind(b1,b2)

! Subroutines
    use utils_subroutines, only : rkfliq, rkf, splintl

! Common blocks  
    use misc
    use puls
    use misc_const
    use dmisc
    use ray
    use rs
    use splinq
    use fracxtal
    
    implicit double precision(a-h,o-z)

    real*8, dimension(8) :: y
    real*8, dimension(51) :: work
    real*8, parameter :: small = 1.0d-02
    integer, dimension(5) :: iwork

!    print *, 'grind'
! put in if statement to treat crystallized and uncrystallized state
    if (fracm_new .lt. small) then
       yliq(1,1)=1.0
    else
       yliq(1,1)=0.0
    endif
    
    yliq(2,1)=eigt/(l*grav*p43*rho(1))
    y(1)=yliq(1,1)
    y(2)=yliq(2,1)
    yliq(3,1)=y3t 
    yliq(4,1)=l*yliq(3,1)
    y(3)=yliq(3,1)
    y(4)=yliq(4,1)
    iflag=1
    neqn=4
    if(iray.eq.0) goto 4
    do j=1,4
       y(j+4)=0.
       rayy(j,1)=0.
    end do
    neqn=8
4   nsurf1=nsurf-1
    do j=1,nsurf1
       k=j+1
       xstart=x(j)
       xfin=x(k)
       relerr=1.d-8
       abserr=dabs(y(1))
       do i=2,neqn 
          t=dabs(y(i))
          abserr=dmin1(abserr,t)
       end do
       abserr=1.d-8*abserr
       abserr=max(abserr,1.d-20)
       call rkf(rkfliq,neqn,y,xstart,xfin,relerr,abserr,iflag,work,iwork)
       iflag=1
       do i=1,4
          yliq(i,k)=y(i)
       end do
       if (iray.eq.0) goto 1
       do i=1,4
          rayy(i,k)=y(i+4)
       end do
1      continue
    end do

    rt=x(nsurf)
    call splintl(rt,vq)
    temp=1./vq(5)-1.
    b1=yliq(1,nsurf)*vq(4)+yliq(3,nsurf)*(l+1.)+yliq(4,nsurf)
    
! these are saio and cox bc's, changed from osaki and hansen's
    
    b2=yliq(1,nsurf)*((4.+r(nsurf)*eigt/gg(nsurf))/temp-1.) &
         +(1.-lhat*gg(nsurf)/r(nsurf)/eigt/temp)*yliq(2,nsurf) &
         + ((l+1.)/temp-1.)*yliq(3,nsurf) 

    return

  end subroutine grind

!************************************************************************

  SUBROUTINE COEFF1 (N,X,F,W,IOP,INT,WK)

    use utils_subroutines, only : trip

    IMPLICIT REAL*8 (A-H,O-Z)
    real*8, dimension(N) :: X
    real*8, dimension(650) :: F, W
    real*8, dimension(N,4) :: wk
    integer, dimension(2) :: IOP

! ARITHMETIC STATENENT FUNCTION USED TO LOCATE ENTRIES IN F AND W ARRAYS

    II(INDEX)=(INDEX-1)*INT+1

! START TO SET UP TRIDIAGONAL SYSTEM

    J0 = 1
    DO I=2,N
       JM = J0
       J0 = J0+INT
       WK(I,1) = X(I)-X(I-1)
       WK(I,2) = (F(J0)-F(JM))/WK(I,1)
       WK(I,3) = WK(I,1)/6.D0
       WK(I,1) = WK(I,1)/3.D0
    end DO
    NN = N
    MK = IOP(1)
    ML = IOP(2)

! APPLY BOUNDARY CONDITIONS AT BOUNDARY 1
    select case (MK)
    case (1) 
       goto 102
    case (2) 
       goto 103
    case (3) 
       goto 104
    case (4) 
       goto 105
    end select

! SECOND DERIVATIVE GIVEN AT BOUNDARY 1

102 CONTINUE
    WK(2,2) = WK(3,2)-WK(2,2)-WK(2,3)*W(1)
    WK(2,3) = 0.D0
    WK(2,1) = WK(2,1)+WK(3,1)
    I1 = 2
    NN = NN-1
    GO TO 106

! FIRST DERIVATIVE GIVEN AT BOUNDARY 1

103 CONTINUE
    WK(1,2) = WK(2,2)-W(1)
    WK(2,2) = WK(3,2)-WK(2,2)
    WK(1,3) = 0.D0
    WK(1,1) = WK(2,1)
    WK(2,1) = WK(2,1)+WK(3,1)
    I1 = 1
    GO TO 106

! PERIODIC BOUNDARY CONDITION

104 CONTINUE
    Y2 = WK(2,2)
    B2 = WK(2,1)
    WK(2,2) = WK(3,2)-WK(2,2)
    WK(2,1) = WK(3,1)+WK(2,1)
    I1 = 2
    NN = NN-1
    GO TO 106

! FIRST DERIVATIVE AT BOUNDARY 1 FROM 4 POINT INTERPOLATION.

105 CONTINUE
    A12 = X(1)-X(2)
    A13 = X(1)-X(3)
    A14 = X(1)-X(4)
    A23 = X(2)-X(3)
    A24 = X(2)-X(4)
    A34 = X(3)-X(4)
    J1 = 1
    J2 = J1+INT
    J3 = J2+INT
    J4 = J3+INT
    W(1)    = (1.D0/A12+1.D0/A13+1.D0/A14)*F(J1)- &
         A13*A14/(A12*A23*A24)*F(J2)+A12*A14/(A13*A23*A34)*F(J3)- &
         A12*A13/(A14*A24*A34)*F(J4)
    GO TO 103

! COMPUTE TRIDIAGONAL ARRAYS

106 CONTINUE
    I2 = N-2
    DO I=3,I2
       WK(I,2) = WK(I+1,2)-WK(I,2)
       WK(I,1) = WK(I+1,1)+WK(I,1)
    end DO

! APPLY BOUNDARY CONDITIONS AT BOUNDARY 2.

    IN = II(N)
    
    select case (ML)
    case (1) 
       goto 108
    case (2) 
       goto 109
    case (3) 
       goto 110
    case (4) 
       goto 111
    end select

!SECOND DERIVATIVE GIVEN AT BOUNDARY 2.

108 CONTINUE
    WK(N-1,2) = WK(N,2)-WK(N-1,2)-WK(N,3)*W(IN)
    WK(N,3) = 0.D0
    WK(N-1,1) = WK(N-1,1)+WK(N,1)
    NN = NN-1
    GO TO 112

! FIRST DERIVATIVE GIVEN AT BOUNDARY 2.

109 CONTINUE
    WK(N-1,2) = WK(N,2)-WK(N-1,2)
    WK(N,2) = -WK(N,2)+W(IN)
    WK(N-1,1) = WK(N-1,1)+WK(N,1)
    WK(1,4) = 0.D0
    GO TO 112

! PERIODIC BOUNDARY CONDITION

110 CONTINUE
    WK(N-1,2) = WK(N,2)-WK(N-1,2)
    WK(N,2) = Y2-WK(N,2)
    WK(N-1,1) = WK(N-1,1)+WK(N,1)
    WK(N,1) = WK(N,1)+B2
    WK(1,4) = WK(2,3)
    GO TO 112

! FIRST DERIVATIVE AT BOUNDARY 2 FROM 4 POINT INTERPOLATION.

111 CONTINUE
    A12 = X(N)-X(N-1)
    A13 = X(N)-X(N-2)
    A14 = X(N)-X(N-3)
    A23 = X(N-1)-X(N-2)
    A24 = X(N-1)-X(N-3)
    A34 = X(N-2)-X(N-3)
    J1 = IN
    J2 = J1-INT
    J3 = J2-INT
    J4 = J3-INT
    W(IN)   = (1.D0/A12+1.D0/A13+1.D0/A14)*F(J1)- &
         A13*A14/(A12*A23*A24)*F(J2)+A12*A14/(A13*A23*A34)*F(J3)- &
         A12*A13/(A14*A24*A34)*F(J4)
    GO TO 109
112 CONTINUE
    IW1 = II(I1)
    CALL TRIP (NN,WK(I1,3),WK(I1,1),WK(I1+1,3),WK(I1,2),W(IW1),INT)
    select case (MK)
    case (1) 
       goto 114
    case (2) 
       goto 114
    case (3) 
       goto 113
    case (4) 
       goto 114
    end select

113 CONTINUE
    W(1) = W(IN)
114 CONTINUE
    RETURN
  END SUBROUTINE COEFF1

!************************************************************************
  subroutine plots(np,totm)

    use values !See comments for this module in commonblocks.f90 for a 
               !description of what each row of xx(:,:) represents.
    use flags
    use freq
    use stuff
    use heat
    use lnrp
    use bvorig
    use xcompp
    
    implicit double precision(a-h,o-z)
    
    real*8, dimension(650) :: fr

    do l=2,np
       fr(l) = dlog10(xx(1,l))! / xx(1,np)
       rem = 1. - mr(l)
       sumup = 0.
       do i=l,np
          sumup = piece(i)+sumup
       end do
       rem = sumup/totm
       if(rem.gt.1.d-18)then 
          partt(l)=-dlog10(rem)
       else
          partt(l)=18.0
       endif
       if (partt(l).lt.17.999) then
          write(60,2108) partt(l),acous(l),bvfreq(l),bvfrq(l),fr(l),xi(l), &
               tfreq(l)
       endif
       write(61,2103) partt(l),xx(16,l),xx(17,l)
       gam3m1=xx(17,l)*ga1(l)
       drivkap = xx(15,l) + xx(14,l)/gam3m1
       write(62,2104) partt(l),xx(14,l),xx(15,l),dlog10(xx(8,l)),fr(l),xi(l), &
            drivkap,tthl(l),gam3m1
       tsix = xx(4,l)/1.0d+6
       tsix3 = tsix**3
       rlog = dlog10(xx(5,l)/tsix3)
       if(xx(3,l).gt.0.0)then
          radl = dlog10(xx(3,l))
       else
          radl = 0.0
       endif
       write(64,2103) partt(l),radl,tthl(l),xi(l)
       write(65,2109) partt(l),dlog10(xx(4,l)),dlog10(xx(5,l)), &
            dlog10(xx(6,l)),fr(l),xi(l),radl
       write(66,2103) partt(l),xx(10,l),xx(11,l)
       epsgen = dabs(xx(7,l))
       write(67,2104) mr(l),xx(12,l),xx(13,l),epsgen,fr(l)
       write(68,2104) partt(l),dlog10(xx(9,l)),dlog10(cp(l)),ga1(l),fr(l)
    end do
    call calc_crystalization(np,partt,xtal1,xtal2)
!    call profgen(np)
    write(74,2115) 
2115 format('# -log(1-Mr/M*)    X_O        X_C        X_He       X_H', &
           '         B          N^2')
    do l=2,np     
       write(63,2107) partt(l), mr(l), xheli(l),xx(20,l),xoxyg(l), &
            xcarb(l),fr(l),xhydr(l),xi(l)
       write(74,2116) (partt(l)-partt(2)),xoxyg(l),xcarb(l),xheli(l), &
            xhydr(l),xx(20,l),bvfreq3(l)
    enddo
    
2103 format(4(f10.5,1x))
2104 format(11(1x,f10.6))
2109 format(7(f12.6))
2105 format(2f10.5)
2106 format(2(f10.5,1x)1pe12.5,1x,0pf10.5,2(1x,f10.5))
2107 format(9(1x,f10.6))
2116 format(2x,6(1x,f10.6),2x,e12.5)
2108 format(7(1x,f10.6))
2114 format(1x,f10.6,1x,'/',1x,f10.6)
    return

  end subroutine plots
  
!************************************************************************
subroutine profgen(ndim)

!Variables
    use misc_const
    use shells, only : nshint
    use values
    use xcompp, only : xcomp_orig, xhydr, xheli, xcarb, xoxyg
    
    implicit none
    
    integer :: i,j,ndim,nsub,imin,imax
    integer, parameter :: nsm=5
    real*8 :: xox,xca,xhe,xh,dhe,yp1,ypn,interp_value
    real*8, dimension(2*nsm) :: xsub, ysub, yder

! Fill the composition arrays used by pulsation routines
    do i=2,ndim
       xhydr(i) = xcomp_orig(i,1)
       xheli(i) = xcomp_orig(i,2)
       xcarb(i) = xcomp_orig(i,3)
       xoxyg(i) = xcomp_orig(i,4)
       xh  = xhydr(i)
       xhe = xheli(i)
       xca = xcarb(i)
       xox = xoxyg(i)
    enddo

  end subroutine profgen
!************************************************************************

  subroutine calc_crystalization(ndim,partt,xtalm1,xtalm2)

!Variables
    use misc_const
    use values
    use xcompp, only : xhydr, xheli, xcarb, xoxyg
    
    implicit none
    
    integer :: i,j,iflag,ndim,l
    integer, parameter :: n=650
    real*8 :: aavg,zavg,t,rho,gammac,gammao, gammahe,gammah,u2,t2, &
         gam1,gam2,xtalm1,xtalm2
    real*8, parameter :: gammin=150.0d0,gammax=200.d0,gamc=3577200.d0, &
         gamo=5778200.d0,gamhe=573264.d0,gamh=227500.d0,eps=1.d-6
    real*8, dimension(n) :: gammatot,xtal,partt

    iflag=0

! Calculate the Coulomb gamma, gamma_xtal (about 160--180) for crystallization)
    do l=2,ndim
       t=xx(4,l)
       rho=xx(5,l)
! the following is lifted from the evolution code
       u2=dlog10(rho)
       t2=dlog10(t)
       gammah  =  gamh * 10.**(u2/3. - t2)
       gammahe = gamhe * 10.**(u2/3. - t2)
       gammac  =  gamc * 10.**(u2/3. - t2)
       gammao  =  gamo * 10.**(u2/3. - t2)
! now kludge up a total gamma using additive volume technique
       zavg=xcarb(l)*6.0 +xoxyg(l)*8.0 +xheli(l)*2.0+xhydr(l)*1.0
       aavg=xcarb(l)*12.0+xoxyg(l)*16.0+xheli(l)*4.0+xhydr(l)*1.0
       gammatot(l)=xheli(l)*gammahe+xcarb(l)*gammac+xoxyg(l)*gammao &
            + xhydr(l)*gammah
       xtal(l)=1.d+0 - 10**(-(partt(l)-partt(2)))
    enddo

! loop through shells to find massfrac corresponding to gammin & gammax

    xtalm1=0.0
    xtalm2=1.0
    do i=2,ndim-1
       gam1=gammatot(i)
       gam2=gammatot(i+1)
       if (gam1.gt.gammax .and. gam2.lt.gammax) xtalm1=xtal(i)
       if (gam1.gt.gammin .and. gam2.lt.gammin) xtalm2=xtal(i+1)
    enddo

  end subroutine calc_crystalization
!************************************************************************
! this subroutine computes the asymptotic frequency and period
! spacings for a given equilibrium model.
! Formulae taken from 5 horsemen.
! Added in computation for Pi_H and Pi_He for trapping cycles
!************************************************************************

  subroutine asymp(np,r,sound,totr,totm,akrdsig2,mu,mr,rho)

! Subroutines
    use utils_subroutines, only: dosim

    use freq
    use resvalue
    use misc_const
    use element
    use modelp1, only: gnu0, per0, tdyn

    implicit double precision(a-h,o-z)

    real*8, dimension(650) ::  bvn,r,sound,h,bvnh,bvnhe,mu,akrdsig2, &
         akrdsig,avgne,rho,mr,dpxtal,mxtal
    real*8 :: perH,perHe,alpha2,ane,torfreq0
    integer :: l

    do ii=1,650
       bvn(ii)=0.0d0
       h(ii)=0.0d0
       bvnh(ii)=0.0d0
       bvnhe(ii)=0.0d0
    end do
       
! get sound speed and brunt-vaisala frequency
      nmax=0
      tdyn=dsqrt(g*totm/totr**3)
      do i=2,np
! integrand for torsional frequency spacing (constant, like p-modes)
         if(akrdsig2(i).gt.0.0 .and. mr(i)/mr(np).lt..99)then
            akrdsig(i)=dsqrt(akrdsig2(i))
         else
            akrdsig(i)=0.0
         endif
! integrand for Pi_o
         if(bvfreq(i).gt.0.0)then
            bvn(i)=dsqrt(bvfreq(i))
            bvn(i)=bvn(i)/totr/r(i)
! integrand for Pi_H
            if(xhe(i).lt.0.97 .and. xhe(i-1).gt.xhe(i))then
               bvnh(i)=dsqrt(bvfreq(i))
               bvnh(i)=bvnh(i)/totr/r(i)
               if (ihflag.eq.0) ihyzone=i
               if (i.lt.np) then
                  avgne(i)=(r(i+1)-1.)*log(rho(i+1)/rho(i))/(r(i+1)-r(i))
               else
                  avgne(i)=0.0
               endif
               ihflag = 1
            elseif(ihflag.eq.1)then
               bvnh(i)=dsqrt(bvfreq(i))
               bvnh(i)=bvnh(i)/totr/r(i)
               if (i.lt.np) then
                  avgne(i)=(r(i+1)-1.)*log(rho(i+1)/rho(i))/(r(i+1)-r(i))
               else
                  avgne(i)=0.0
               endif
            else
               bvnh(i)=0.0
               avgne(i)=0.0
            endif
! integrand for Pi_He
            if(xhe(i).gt.0.03)then
               bvnhe(i)=dsqrt(bvfreq(i))
               bvnhe(i)=bvnhe(i)/totr/r(i)
            elseif(xhe(i).lt.0.03 .and. xhe(i-1).gt.xhe(i))then
               bvnhe(i)=dsqrt(bvfreq(i))
               bvnhe(i)=bvnhe(i)/totr/r(i)
               iheflag = 1
            elseif(iheflag.eq.1)then
               bvnhe(i)=dsqrt(bvfreq(i))
               bvnhe(i)=bvnhe(i)/totr/r(i)
            else
               bvnhe(i)=0.0
            endif
! add zero to integrand in convection zone
         else
            bvn(i)=0.0
            bvnh(i)=0.0
            bvnhe(i)=0.0
            avgne(i)=0.0
         endif
      end do
      do n=2,np
         sound(n)=1./sound(n)
         rs=r(n)*totr
         if(n.ne.1)then
            h(n)=rs-rm1
         else
            h(1)=r(1)*totr
         endif
         rm1=rs
      end do
      npm1=np-1
      call dosim(rslt,akrdsig,h,npm1)
      torint=rslt
      torfreq0=(1./2.)/torint
      pertor=1./torfreq0
!      write(*,*) 'torfreq0 pertor=',torfreq0,pertor
      call dosim(rslt,sound,h,npm1)
      cint=rslt
      call dosim(rslt,bvn,h,npm1)
      gint=rslt
      gnu0=1./(2.*cint)
      pgnu0=1./gnu0
      per0=2.*pi*pi/gint

      do i=3,npm1
         mxtal(i)=mr(i)/mr(np)
         if (mxtal(i).lt.0.98) then
            dpxtal(i)=2.*pi*pi/(res(npm1)-res(i))
            write(9,1004) mxtal(i),dpxtal(i)
         endif
1004     format(2(e12.5,2x))
      enddo

      if(ihflag .eq. 1)then
         call dosim(rslt,bvnh,h,npm1)
         ghint=rslt
         perH=2.*pi*pi/ghint
         call dosim(ane,avgne,h,npm1)
         dr=totr*(1.-r(ihyzone))
         ane=ane/dr
      else
         perH=0.0
      endif
      if(iheflag .eq. 1)then
         call dosim(rslt,bvnhe,h,npm1)
         gheint=rslt
         perHe=2.*pi*pi/gheint
      else
         perHe=0.0
      endif
!      write(*,*) per0,perH,ane,dr
      write(8,1005) gnu0,per0,perH,perHe,tdyn
      write(50,1003) gnu0,per0,tdyn
1003  format(3(2x,f8.4))
1005  format(2x,'asymptotic frequency spacing =',1pe12.5,/ &
           2x,'asymptotic period spacing =',3x,1pe12.5,/ &
           2x,'Characteristic period spacing for Hydrogen =',3x,1pe12.5,/ &
           2x,'Characteristic period spacing for Helium =',3x,1pe12.5,/ &
           2x,'dynamical time scale =',8x,1pe12.5)
      return
    end subroutine asymp

!************************************************************************
    subroutine modeid

      use modect

      implicit double precision(a-h,o-z)

      integer :: quad2,quad1

      modep=0

! first count the nodes in y1 and y2

      nodes1=0
      nodes2=0
      do i=2,nfine
         fit1=y1m(i)*y1m(i-1)
         fit2=y2m(i)*y2m(i-1)
         if (fit1.le.0) nodes1=nodes1+1
         if (fit2.le.0) nodes2=nodes2+1
      end do

! count crossings in the phase diagram

      quad1=quad(y1m(1),y2m(1))
      do i=2,nfine
         quad2=quad(y1m(i),y2m(i))
         idq=quad2-quad1
         if (idq.eq.0) go to 30
!  see if the quadrant change is a crossing in y1 
         if (abs(idq).eq.3) go to 100
         if (quad1.eq.3.and.quad2.eq.2) go to 100
         if (quad1.eq.2.and.quad2.eq.3) go to 100
!  not a y1 crossing
         quad1=quad2
         go to 30

! if crossing is clockwise, increment mode count; 
! if crossing is negative, decrement mode count

100      if (idq.eq.-3.or.idq.eq.1) kdmod=-1
         if (idq.eq.3.or.idq.eq.-1) kdmod=1
         modep=modep+kdmod
         quad1=quad2
! get next point
30       continue         
      end do
      return

    contains

      integer function quad(y1,y2)
        
        implicit none
        
        real*8, intent(in) :: y1, y2
        
        if (y1.gt.0.and.y2.ge.0) quad=1
        if (y1.le.0.and.y2.gt.0) quad=2
        if (y1.lt.0.and.y2.le.0) quad=3
        if (y1.ge.0.and.y2.lt.0) quad=4
        
        return
      end function quad
      
    end subroutine modeid

!************************************************************************
    subroutine eigenf

      use utils_subroutines, only: splintl, spline, splint

      use misc_const
      use puls
      use eig1
      use misc
      use dmisc
      use rs
      use splinq
      use gkern
      use if

      implicit double precision(a-h,o-z)

 
      real*8, dimension(650) :: dpartphi(650)

! save integration of mode inertia
      facke=r(nsurf)**(2.d0*lindex)/yliq(1,nsurf)**2
      xinertia=rayy(1,nsurf)*facke
      
      do i=1,nsurf
         t1=1./r(i)**lindex
         do j=1,4
            yliq(j,i)=yliq(j,i)*t1
         end do
      end do
      t1=yliq(1,nsurf)
      do i=1,nsurf
         do j=1,4
            yliq(j,i)=yliq(j,i)/t1
         end do
      end do
      do i=1,nsurf
         rt=x(i)
         call splintl(rt,vq)
         rayy(1,i)=rho(i)*r(i)**2*(yliq(1,i)**2+lhat*(vq(1)/eigt)**2* &
              yliq(2,i)**2)
         temp=rho(i)*gg(i)*r(i)
         rayy(2,i)=temp*vq(2)*(yliq(2,i)-yliq(3,i))**2
         rayy(3,i)=temp*vq(3)*yliq(1,i)**2
         rayy(4,i)=-temp*(yliq(4,i)+(l+1.)*yliq(3,i))**2/vq(4)
! *** the g kernel is per unit x=ln(r/p), which will be the independent
! *** variable for the inversions
! *** the g-mode kernels **for period shifts (in sec)** must be 
! *** properly normalized
         gkernel(i)=(r(i)*vq(5)*rayy(3,i)*r(i)**2)* &
              (-period**3/(8.*pi**2*xinertia))
         temp2= sqrt(abs(vq(1)*vq(3)))*vq(5)
         gkernel2(i)=gkernel(i)/temp2
      end do

! compute normalized mode phase as a function of radius (shell number)
      if (ifirst.eq.0) then
         phi(1)=0.0
         small=1.d-06
         do i=2,nsurf
            rt=x(i)
            call splintl(rt,vq)
            dphi=sqrt(abs(vq(1)*vq(3)))*vq(5)*(x(i)-x(i-1))
            if (dphi.lt.small) dphi=small
            phi(i)=phi(i-1) + dphi
         enddo
         do i=1,nsurf
            phi(i)=phi(i)/phi(nsurf)
         enddo
! call spline routine to set up spline coefficients
         ybig=1.d+32
         call spline(x,phi,nsurf,ybig,ybig,dphix)
         call spline(phi,x,nsurf,ybig,ybig,dxphi)
         call spline(phi,part,nsurf,ybig,ybig,dpartphi)
         open(unit=52,file='reflection2.dat',status='unknown')
         do i=1,nsurf
            x0=x(i)
            phirev=1.d0 - phi(i)
            call splint(phi,x,dxphi,nsurf,phirev,xrev)
            call splint(phi,part,dpartphi,nsurf,phirev,prev)
            partrev(i)=prev
            write(52,10) x0,xrev,phi(i),phirev,part(i),partrev(i), &
                 1.d0-10**(-part(i))
         enddo
10       format(7(e12.5,2x))

         ifirst=1
      endif

      sum=0.0

      return
    end subroutine eigenf

!************************************************************************
!** compute integrands for rotational splitting coeffs
!************************************************************************
    subroutine rotcoef

      use utils_subroutines, only: dosim

      use misc_const
      use puls
      use misc
      use dmisc
      use rs
      use ekint
      use eig1
      use rot1
      use rot2

      implicit double precision(a-h,o-z)

      rstar = r(nsurf)
      do  i=1,nsurf
         if(i.eq.1)goto 14
         h(i)=r(i)-rmi
14       rmi=r(i)
         rorstr = rmi / rstar
         aa=r(i)*yliq(1,i)
         bb=(gg(i)/eigt)*yliq(2,i)
         front = rho(i)*r(i)*r(i)
         rone(i)=front*(2.*aa*bb+(bb*bb))
         rpone(i)=front*(2.-rorstr)*(aa**2-(2.*aa*bb)-((1.-lhat)*bb**2))
         rptwo(i)=rpone(i)*(5.-rorstr)
         rpthr(i)=rpone(i)/dsqrt(2.-rorstr)
         rpfour(i)=rpone(i)/(2.-rorstr)**0.75
         rpfive(i)=rpone(i)/(2.-rorstr)**0.90
         angfac(i)=rho(i)*r(i)**4
         f(i)=rho(i)*r(i)*r(i)*((aa*aa)+(bb*bb)*lhat)
      end do
!***compute moment of inertia ********
      call dosim(angmom,angfac,h,nsurf)
      angmom=(8.*pi/3.)*angmom
!*** compute uniform rotational splitting factor clk
      call dosim(clk,rone,h,nsurf)
      clk = clk / ekin
!*** compute nonuniform splitting portion based on (2-r) power law
      call dosim(crone,rpone,h,nsurf)
      crone = crone / ekin
      betone = rpone(1)/rpone(nsurf)
!*** compute nonuniform splitting portion based on (2-r)(5-r) power law
      call dosim(crtwo,rptwo,h,nsurf)
      crtwo = crtwo / ekin
      bettwo = rptwo(1)/rptwo(nsurf)
!**** compute splitting based on (2-r)1/2 power law
      call dosim(crthr,rpthr,h,nsurf)
      crthr=crthr/ekin
      betthr = rpthr(1)/rpthr(nsurf)
!**** now rpfour in splitting based on a (2-r)1/4 power law 
      call dosim(crfour,rpfour,h,nsurf)
      crfour=crfour/ekin
      betfour = rpfour(1)/rpfour(nsurf) 
!**** compute rotational splitting based on (2-r)0.1 power law
      call dosim(crfive,rpfive,h,nsurf)
      crfive = crfive/ekin
      betfive = rpfive(1)/rpfive(nsurf) 
!*** print out results ******************************************
!      write(*,1940) clk,crone,crtwo,crthr
!      write(*,1942) angmom,ekin,crfour,crfive
!      write(*,1943) betone,bettwo,betthr,betfour,betfive
!     write(51,1941) idint(l),modep,period,clk,crone,crtwo,crthr
!     write(51,1942) angmom,ekin,crfour,crfive
!     write(51,1943) betone,bettwo,betthr,betfour,betfive
1940  format(10x,'clk =',e12.5,5x,'c (r+1)=',e12.5,3x,'c (2-r)(5-r)', & 
           e12.5,5x,'c (r+0.5)',e12.5)
1941  format(4x,i2,1x,i2,2x,f10.3,1x,'c lk',e12.5,2x,'c(r-1)',e12.5, &
           x,'c(r+1)',e12.5,2x,'c(r-0.5)',e12.5)
1942  format(10x,'mom. of inertia',e12.5,'bottom',e12.5,'c (r+1/4)', &
           e12.5,'c (+0.1)',e12.5)
1943  format(5x,'omega(cent)/omega(surf)',5(2x,e12.5))
      
      return
      
    end subroutine rotcoef
    
!************************************************************************


end module pulse_subroutines
