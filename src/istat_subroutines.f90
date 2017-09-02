module istat_subroutines

contains

  !************************************************************************

  subroutine istatco(pp2,tt2,iw,converg)

    ! compute the thermodynamic variables for a specific P, T
    ! point.  This version searches both the c and o eos tables,
    ! and then interpolates to the desired composition using the
    ! additive volume technique

    ! Subroutines
    use phase_func
    use utils_subroutines, only : dali
    use nuax_subroutines, only: gnutrino, axions
    use eprep_subroutines, only: opac
    use wd_eos_mod, only : wd_eos
    use kappa, only : opalz


    ! Common blocks
    use flags, only: verbose, evoloutput, ios_new, kap_new
    use misc_const
    use temp
    use tablecc
    use tableoo
    use sizec
    use sizeo
    use thermo 
    use thermoc
    use thermoo
    use neut
    use iiii
    use ixswchh
    use jee
    use ax

    implicit double precision(a-h,o-z)

    real*8, dimension(4) :: porderc,torderc,xorderc,yorderc
    real*8, dimension(4) :: pordero,tordero,xordero,yordero
    real*8, dimension(7) :: thermoc_array, thermoo_array
    real*8, pointer :: d2,dp,dt,e,ep,et,psi,op,ot,fcc2,fa2
    real*8 :: d3, dt3, dp3, e3, et3, ep3, psi3, o3, ot3, op3
    real*8, dimension(7) :: xavec
    real*8, dimension(11) :: res

    integer, dimension(4) :: iorderc,jorderc,iordero,jordero
    logical :: converg

    real*8, parameter :: delta = 2d-3

    d2   => u2
    dp   => up2
    dt   => ut2
    e    => e2 
    ep   => ep2
    et   => et2
    psi  => psi2
    op   => op2
    ot   => ot2
    fcc2 => fcc
    fa2  => fa

    ! initialize some variables

    xc = xc2
    xo = xo2
    xhe = xhe2
    je=3
    iii=0
    ci=1/(xc/12.+xo/16.)
    p=pp2
    t=tt2

    if (.not. ios_new) then

       ! find eos box and interpolate.
       ! note that we're doing both tables in parallel.

       call wdselect(t,tminc,deltc,torderc,iorderc,1,itc,converg)
       call wdselect(t,tmino,delto,tordero,iordero,1,ito,converg)

       kmax=6
       if ( iw .ge. 1 )  kmax=1

       do k=1,kmax 
          do i=1,4 
             indexc=iorderc(i) 
             indexo=iordero(i) 
             if ( ixswch .eq. 1 ) then
                pfrzc = phasec(t)
                pfrzo = pfrzc
             elseif ( ixswch .eq. 2 ) then
                pfrzo = phaseo(t)
                pfrzc = pfrzo
             elseif ( ixswch .eq. 3 ) then
                pfrzc = xc * phasec(t) + xo * phaseo(t)
                pfrzo = pfrzc
             elseif ( ixswch .eq. 4 ) then
                pfrzc = phasec(t)
                pfrzo = phaseo(t)
             else
                if ( evoloutput ) then
                   write(17,*) 'ixswch not set. exiting.'
                   write(20,*) 'ixswch not set. exiting.'
                end if
                if ( verbose ) print *, 'ixswch not set. exiting.'
                stop 
             endif

             if ( p .lt. pfrzc ) then
                jminc=iphasec(indexc)+1
                jmaxc=ipc(indexc)
                pstartc=pmaxc-ioverc*delpc
             else
                jminc=1
                jmaxc=iphasec(indexc)
                pstartc=pmaxc
             endif

             if ( p .lt. pfrzo ) then
                jmino=iphaseo(indexo)+1
                jmaxo=ipo(indexo)
                pstarto=pmaxo-iovero*delpo
             else
                jmino=1
                jmaxo=iphaseo(indexo)
                pstarto=pmaxo
             endif

             call wdselect(p,pstartc,delpc,porderc,jorderc,jminc,jmaxc,converg)
             call wdselect(p,pstarto,delpo,pordero,jordero,jmino,jmaxo,converg)

             do j=1,4
                xorderc(j)=tablec(k,indexc,jorderc(j))
                xordero(j)=tableo(k,indexo,jordero(j))
             end do

             call dali(p,porderc,xorderc,yorderc(i))
             call dali(p,pordero,xordero,yordero(i))

          end do

          l=k
          if(k.gt.4) l=k+1

          call dali(t,torderc,yorderc,thermoc_array(l))
          call dali(t,tordero,yordero,thermoo_array(l))

       end do

       dc = thermoc_array(1)
       dpc = thermoc_array(2)
       dtc = thermoc_array(3)
       ec = thermoc_array(4)
       epc = thermoc_array(5)
       etc = thermoc_array(6)
       psic = thermoc_array(7)

       do = thermoo_array(1)
       dpo = thermoo_array(2)
       dto = thermoo_array(3)
       eo = thermoo_array(4)
       epo = thermoo_array(5)
       eto = thermoo_array(6)
       psio = thermoo_array(7)

       ! c/o interpolations finished.  at this time we've interpolated in the c and o 
       ! tables to get d, dp, dt, e, et, and psi for each composition. Use additive 
       ! volume technique to interpolate the therm quantities for the mixture.
       ! see Fontaine, Graboske, & Van Horn 1977, ApJS, 35, 293
       ! see my log 4 for derivation of these interpolations.

       d2 = -dlog10( xc/10.**dc + xo/10.**do )

       d  = d2

       dt  = 1./( d * ( xc/(dc*dtc) + xo/(do*dto) ) )
       dp  = 1./( d * ( xc/(dc*dpc) + xo/(do*dpo) ) )
       e  = xc * 10.**ec + xo * 10.**eo
       et = xc * etc * 10.**ec + xo * eto * 10.**eo
       ep = 10.**(p-t)/fnat * (xc*dtc*10.**(-dc) + xo*dto*10.**(-do))
       psi = xc * psic + xo * psio

    else

       xavec(:) = 0.d0
       xavec(1) = xh
       xavec(2) = xhe
       xavec(3) = xc
       xavec(5) = xo
       call wd_eos(10**p,10**t,xavec,'p',res)

       d3 = dlog10(res(1))
       ep3 = -dlog(10.)*res(2)*res(4)
       e3  = res(10)
       dt3 = -res(6)/res(5)
       dp3 = 1./res(5)
       et3 = dlog(10.)*res(2)
       psi3 = res(9)

       ! write(*,'(8f8.4,10e13.5)') p,t,d2,d3,dt,dt3,dp,dp3,e,e3,ep,ep3,et,et3,psi,psi3
       ! write(*,'(8f8.4,10e13.5)') p,t,d,d/d3,dt/dt3,dp/dp3,e/e3,ep/ep3,et/et3

       d2 = d3
       dt = dt3
       dp = dp3
       e = e3
       et = et3
       ep = ep3
       psi = psi3

    endif

    d  = d2

    ! Now that we've interpolated the main quantities, we need to calculate the 
    ! opacity, carbon burning and neutrino luminosities and their derivitives wrt 
    ! rho and t. Note that the carbon burning luminosity is reduced by a factor
    ! of 5 arbitrarily (otherwise, runaway).

    n  = -1
20  x1=49.4+d-2*t/3-36550*(1+7e-1*10.**(t-10))**(1./3)*10.**(-t/3)

    if ( x1 .le. -30. ) then
       ffcc=0
    else
       ffcc=4.368*xc2*xc2*10.**x1
       ffcc=ffcc/5. 
    endif

    if (iw .ge. 1) then
       fcc2=ffcc
       goto 66
    endif

    if (nu .le. 0) then
       fn=0
    elseif(nu .gt. 0)then
       call gnutrino(d,t)
    endif

    if ( maxion .gt. 0.0d0 ) then
       call axions(xc,d,t,epsaxion)
       fax=epsaxion
    else
       fax=0.0d0
    end if

    f=fn-ffcc*w + fax
    if (.not. kap_new ) then
       o=10.**opac(d,t)
    endif
    if (n .lt. 0) goto 14
    if (n .eq. 0) goto 15
    if (n .gt. 0) goto 16

14  o2=o
    ff2=f
    fn2=fn
    fcc2=ffcc
    fa2=fax
    do i=1,5
       en2(i)=en(i)
    enddo


    pg=dlog10(10.**p-10.**(4*t-14.59836))
    d=d2+delta
    n=0 
    go to 20

15  od=(o-o2)/delta
    fd=(f-ff2)/delta
    d=d2
    t=tt2+delta
    n=1 
    go to 20

16  ot=(o-o2)/delta
    ft=(f-ff2)/delta
    op=od*dp
    ot=ot+od*dt
    fp=fd*dp
    ft=ft+fd*dt

    if (kap_new) then
       call opalz(10**tt2,10**d2, xavec, fkappa, opr, opt)
       o3 = fkappa
       ot3 = opt*fkappa*dlog(10.)
       op3 = opr*fkappa*dlog(10.)*dp
       ! write(*,'(2f8.4,6es13.5,4f8.3)') d2, tt2, o2, fkappa, ot, ot3, op, op3
       ! write(*,'(10f8.4,6es13.5,4f8.3)') d2, tt2, o2/fkappa, ot/ot3, op/op3
       
       op = op3
       ot = ot3
       o2 = o3

    endif

66  return

  end subroutine istatco

  !************************************************************************

  subroutine istat1(pp2,tt2,iw,iswch,converg)

    ! compute the thermodynamic variables for a specific P, T
    ! point.  This version searches either the C or O eos table, depending
    ! on the composition and returns the thermodynamic properties.

    ! Subroutines
    use utils_subroutines, only : dali
    use nuax_subroutines, only : gnutrino, axions
    use phase_func
    use eprep_subroutines, only: opac
    use wd_eos_mod, only : wd_eos
    use kappa, only : opalz

    ! Common blocks
    use misc_const
    use temp
    use tablecc
    use tableoo
    use sizec
    use sizeo
    use thermo
    use neut
    use iiii
    use ixswchh
    use jee
    use ax
    use flags

    implicit double precision(a-h,o-z)

    logical :: converg
    real*8, dimension(4) :: porder,torder,xorder, yorder
    real*8, pointer :: d2,dp,dt,e,ep,et,psi,op,ot,fcc2,fa2

    real*8 :: d3, dt3, dp3, e3, et3, ep3, psi3
    real*8, dimension(7) :: xavec
    real*8, dimension(11) :: res

    real*8, dimension(24) :: thermo_array
    integer, dimension(4) :: iorder,jorder
    real*8, parameter :: delta = 2.d-3

    d2   => u2
    dp   => up2
    dt   => ut2
    e    => e2 
    ep   => ep2
    et   => et2
    psi  => psi2
    op   => op2
    ot   => ot2
    fcc2 => fcc
    fa2  => fa

    ! initialize some variables

    xc = xc2
    xo = xo2
    je=3
    iii=0
    p=pp2
    t=tt2

    if (.not. ios_new) then

       if ( iswch .eq. 12 ) then
          ci=1./(xc/12.+xo/16.) 
       else if ( iswch .eq. 16 ) then
          ci=1./(xo/16.+(1-xo)/21.45)
       else
          if (verbose) write(*,*) 'istat1: iswch not set.  exiting.'
          stop
       endif
       thermo_array(23) = ci

       ! find eos box and interpolate.
       ! We search carbon if ISWCH is 12 and oxygen if ISWCH is 16.

       if( iswch .eq. 12 ) then
          call wdselect(t,tminc,deltc,torder,iorder,1,itc,converg)
          kmax=6
          if( iw .ge. 1 )  kmax=1
          do k=1,kmax
             do i=1,4
                indexc=iorder(i)
                pfrzc = phasec(t)
                if( p .lt. pfrzc ) then 
                   jminc=iphasec(indexc)+1
                   jmaxc=ipc(indexc)
                   pstartc=pmaxc-ioverc*delpc
                else 
                   jminc=1
                   jmaxc=iphasec(indexc) 
                   pstartc=pmaxc
                endif

                call wdselect(p,pstartc,delpc,porder,jorder,jminc,jmaxc, converg)

                do j=1,4
                   xorder(j)=tablec(k,indexc,jorder(j))
                end do
                call dali(p,porder,xorder,yorder(i))
             enddo
             l=k
             if(k.gt.4) l=k+1
             call dali(t,torder,yorder,thermo_array(l))
             ! Transfer thermo_array terms to thermo common block
             d2     = thermo_array( 1)
             dp     = thermo_array( 2)
             dt     = thermo_array( 3)
             e      = thermo_array( 4)
             ep     = thermo_array( 5)
             et     = thermo_array( 6)
             psi    = thermo_array( 7)
             !Potential source of bug
             !Not sure whether call to dali always changes only first 7 elements of array
!!$          pg     = thermo_array( 8)
!!$          o2     = thermo_array( 9)
!!$          op     = thermo_array(10)
!!$          ot     = thermo_array(11)
!!$          fp     = thermo_array(12)
!!$          ft     = thermo_array(13)
!!$          en2(1) = thermo_array(14)
!!$          en2(2) = thermo_array(15)
!!$          en2(3) = thermo_array(16)
!!$          en2(4) = thermo_array(17)
!!$          en2(5) = thermo_array(18)
             fn2    = thermo_array(19)
!!$          fcc2   = thermo_array(20)
!!$          fa2    = thermo_array(21)
!!$          ce     = thermo_array(22)
!!$          ci     = thermo_array(23)
!!$          cif    = thermo_array(24)
          enddo

          !  or if oxygen, then come here

       else if ( iswch .eq. 16 ) then
          call wdselect(t,tmino,delto,torder,iorder,1,ito,converg)
          kmax=6
          if ( iw .ge. 1 )  kmax=1
          do k=1,kmax
             do i=1,4
                indexo=iorder(i)
                pfrzo = phaseo(t)
                pfrzc = pfrzo
                if ( p .lt. pfrzo ) then 
                   jmino=iphaseo(indexo)+1
                   jmaxo=ipo(indexo)
                   pstarto=pmaxo-iovero*delpo
                else 
                   jmino=1
                   jmaxo=iphaseo(indexo) 
                   pstarto=pmaxo
                endif

                call wdselect(p,pstarto,delpo,porder,jorder,jmino,jmaxo,converg)

                do j=1,4
                   xorder(j)=tableo(k,indexo,jorder(j))
                enddo
                call dali(p,porder,xorder,yorder(i))
             end do
             l=k
             if(k.gt.4) l=k+1
             call dali(t,torder,yorder,thermo_array(l))
             !Potential source of bug
             !Not sure whether call to dali always changes only first 7 elements of array
             ! Transfer thermo_array terms to thermo common block
             d2     = thermo_array( 1)
             dp     = thermo_array( 2)
             dt     = thermo_array( 3)
             e      = thermo_array( 4)
             ep     = thermo_array( 5)
             et     = thermo_array( 6)
             psi    = thermo_array( 7)
!!$          pg     = thermo_array( 8)
!!$          o2     = thermo_array( 9)
!!$          op     = thermo_array(10)
!!$          ot     = thermo_array(11)
!!$          fp     = thermo_array(12)
!!$          ft     = thermo_array(13)
!!$          en2(1) = thermo_array(14)
!!$          en2(2) = thermo_array(15)
!!$          en2(3) = thermo_array(16)
!!$          en2(4) = thermo_array(17)
!!$          en2(5) = thermo_array(18)
             fn2    = thermo_array(19)
!!$          fcc2   = thermo_array(20)
!!$          fa2    = thermo_array(21)
!!$          ce     = thermo_array(22)
!!$          ci     = thermo_array(23)
!!$          cif    = thermo_array(24)
          end do

       endif

       d  = d2

       e  = 10.**e
       et = e * et
       ep = (dt/fnat)*10.**(p-d-t)

    else

       xavec(:) = 0.d0
       if ( iswch .eq. 12 ) then
          xavec(3) = 1.d0
       else if ( iswch .eq. 16 ) then
          xavec(5) = 1.d0
       else
          write(*,*) 'istat1: iswch not set.  exiting.'
          stop
       endif
       call wd_eos(10**p,10**t,xavec,'p',res)

       d3 = dlog10(res(1))
       ep3 = -dlog(10.)*res(2)*res(4)
       e3  = res(10)
       dt3 = -res(6)/res(5)
       dp3 = 1./res(5)
       et3 = dlog(10.)*res(2)
       psi3 = res(9)

!debug
 ! write(*,'(8f8.4,10e13.5)') p,t,d,d3,dt,dt3,dp,dp3,e,e3,ep,ep3,et,et3,psi,psi3
 ! write(*,'(8f8.4,10e13.5)') p,t,d,d/d3,dt/dt3,dp/dp3,e/e3,ep/ep3,et/et3

       d2 = d3
       dt = dt3
       dp = dp3
       e = e3
       et = et3
       ep = ep3
       psi = psi3

    endif

    d  = d2

    ! Now that we've interpolated the main quantities, we need to calculate the 
    ! opacity, carbon burning and neutrino luminosities and their derivitives wrt 
    ! rho and t. Note that the carbon burning luminosity is reduced by a factor
    ! of 5 arbitrarily (otherwise, runaway).

    n  = -1
20  x1=49.4+d-2*t/3-36550*(1+7e-1*10.**(t-10))**(1./3)*10.**(-t/3)

    if ( x1 .le. -30. ) then
       ffcc=0
    else
       xtmp = .999999
       ffcc=4.368*xtmp*xtmp*10.**x1
       ffcc=ffcc/5. 
    endif

    if (iw .ge. 1) then
       fcc2=ffcc
       thermo_array(20) = fcc2
       goto 76
    endif

    if (nu .le. 0) then
       fn=0
    elseif(nu .gt. 0)then
       call gnutrino(d,t)
    endif

    if ( maxion .gt. 0.0d0 ) then
       call axions(xc,d,t,epsaxion)
       fax=epsaxion
    else
       fax=0.0d0
    end if

    f=fn-ffcc*w + fax
    if (.not. kap_new ) then
       o=10.**opac(d,t)
    endif

    if (n.lt.0) goto 14
    if (n.eq.0) goto 15
    if (n.gt.0) goto 16

14  o2=o
    thermo_array(9) = o2
    ff2=f
    fn2=fn
    fcc2=ffcc
    fa2=fax
    do i=1,5
       en2(i)=en(i)
    enddo
!    e  = 10.**e
!    et = e * et
!    ep = (dt/fnat)*10.**(p-d-t)
    pg=dlog10(10.**p-10.**(4*t-14.59836))
    thermo_array(8) = pg
    d=d2+delta
    n=0 
    thermo_array(4) = e
    thermo_array(5) = ep
    thermo_array(6) = et
    thermo_array(14) = en2(1)
    thermo_array(15) = en2(2)
    thermo_array(16) = en2(3)
    thermo_array(17) = en2(4)
    thermo_array(18) = en2(5)
    thermo_array(19) = fn2
    thermo_array(20) = fcc2
    thermo_array(21) = fa2

    goto 20

15  od=(o-o2)/delta
    fd=(f-ff2)/delta
    d=d2
    t=tt2+delta
    n=1 
    goto 20

16  ot=(o-o2)/delta
    ft=(f-ff2)/delta
    op=od*dp
    ot=ot+od*dt
    fp=fd*dp
    ft=ft+fd*dt

    thermo_array(10) = op
    thermo_array(11) = ot
    thermo_array(12) = fp
    thermo_array(13) = ft


    if (kap_new) then
       call opalz(10**tt2,10**d2, xavec, fkappa, opr, opt)
       o3 = fkappa
       ot3 = opt*fkappa*dlog(10.)
       op3 = opr*fkappa*dlog(10.)*dp
       ! write(*,'(2f8.4,6es13.5,4f8.3)') d2, tt2, o2, fkappa, ot, ot3, op, op3
       ! write(*,'(10f8.4,6es13.5,4f8.3)') d2, tt2, o2/fkappa, ot/ot3, op/op3
       
       op = op3
       ot = ot3
       o2 = o3

    endif

76  return

  end subroutine istat1

  !************************************************************************

  subroutine wdselect(x,z,dz,xorder,iorder,jmin,jmax,converg)

    !  find the 'coordinates' in the P, T grid of four points that enclose 
    !  a specified pressure and temperature.
    !  NOTE: if you try to go off the grid, failsafe is visited.
    !  do it too many times, and the program quits.

    use iiii
    use ifaill

    implicit double precision(a-h,o-z)

    logical :: converg
    real*8, dimension(4) :: xorder
    integer, dimension(4) :: iorder

    !*** changed itest from 50 to 100 (TSM) ***
    itest=100
    j=(x-z)/dz+1.5e0
    zmax = jmax * dz + z
    if ( j .le. jmin ) then 
       j=jmin
       iii=iii+1
       if ( converg .and. iii .ge. itest ) then 
          stop
       endif
       ifail = ifail + 1
    elseif ( j .gt. jmax ) then
       j=jmax
       iii=iii+1
       if ( iii .ge. itest ) then
          stop
       endif
       ifail = ifail + 1
    endif

    ii=j
    jl=0
    jr=0
    do i=1,4
       xorder(i)=z+(ii-1)*dz
       iorder(i)=ii
       if(j+jr-jmax) 12,15,12
12     if(j-jl-jmin) 13,14,13
13     if((xorder(i)-x)*dz) 14,15,15
14     jr=jr+1
       ii=j+jr
       goto 16
15     jl=jl+1
       ii=j-jl
16  end do

    return

  end subroutine wdselect

  !************************************************************************

end module istat_subroutines







