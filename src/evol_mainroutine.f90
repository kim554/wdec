subroutine evolve (tnot)
!************************************************************************
!*  This is the white dwarf evolution code as written by Martin         *
!*  Schwarzschild. Subsequent modifications made by:                    *
!*  Version         Author(s)                                           *
!*  WDEC 1.0        Kutter & Savedoff   (1969, ApJ, 157, 1021)          *
!*  WDEC 2.0        D.Q. Lamb           (1975, ApJ, 200, 306)           *
!*  WDEC 3.0        D.E. Winget                                         *
!*  WDEC 4.0        S.D. Kawaler                                        *
!*  WDEC 5.0        M.A. Wood           (aka WDXD)                      *
!*  WDEC 6.0        P.A. Bradley        (aka WDXDL, WDEC6)              *
!*  WDEC 7.0        M.H. Montgomery     (aka WDXDFIT)                   *
!*  Specialized for Metacomputer use by Travis Metcalfe, 1998/1999      *
!************************************************************************

!Subroutines
  use evol_subroutines, only : openem, read1, read2, write1, begin, calc, &
       gshell, write2, interp, end, setprof
  use utils_subroutines, only : check_time

!Common blocks
  use comptmp
  use contrl
  use crash
  use deldgg
  use dprep
  use flags
  use flagtmp
  use idiffuss
  use idone
  use ifaill
  use ixswchh
  use kill
  use mixl
  use oldval
  use shells2
  use startmod
  use temp
  use tfitt
  use threshh
  use writeitt
  use wprep
  use xcompp, only: xcomp !debugging

  implicit double precision(a-h,o-z)

  integer :: newtime, tnot
  logical :: callcshell
     
  icomp = 1
  ihotspot = 0
  irestart = 0
  t0 = tnot
  
666 killit = .false.
  ikill = 0
  itconverged = 0
  itcrashed = 0

!*** fix temperature *** 
  tfit(1) = Teff

!*** fix H and He boundaries ***
  amhyhe_tmp = 10**(-1.*M_h)
  amheca_tmp = 10**(-1.*M_he)
  amheca2_tmp = 10**(-1.*M_env)

  sm_tmp = SM_star

  itfit=0
  ifit=1

  call openem
  if ( evoloutput ) then
     write(17,321) 
     write(20,321) 
  end if
  if (verbose) then
     print 321
321  format(/,20x,'white dwarf evolution code wdxd'/ &
          //22x'last modified August 94 (PAB)'//)
  end if

!*** read from input file ***
  read(7,*) thresh, deldg, ixswch, iprep, idiffus
  read(7,*) fmls, fmls !We ignore these

!Toggle MLT convection version
mls = 1 !Bohm & Cassinelli 
!mls = 0 !Bohm-Vitense 

  nxmods = 5

  if ( verbose ) then
     print *,    '    nxmods entered:',nxmods
     if(iprep .eq. 1)then
        print *,    'Prepped models? YES!' 
     else
        print *,    'Prepped models? NO!' 
     endif
     if(idiffus .eq. 1)then
        print *,    'Diffusionesque transition zones'
     else
        print *,    'Discontinuous transition zones'
     endif
     if (icomp.eq.1) then
        write(*,337) 
337     format(/' Layering of envelope is C/He4/H -- normal (DA) case'/)
     elseif (icomp.eq.2) then
        write(*,333) 
333     format(/' Layering of envelope is C/He4/He3 -- a DB with a', &
             ' surface Helium3 layer'/)
     elseif (icomp.eq.3) then
        write(*,334) 
334     format(/' Layering of envelope is C/He4/He4. This is for test', &
             ' purposes--it',/' should be equivalent to a normal DB model', &
             ' with just C/He layering'/)
     else
        write(*,*) 'Flag icomp must be set...Exiting'
	stop
     end if

     if(mls .eq. 0)then      
        print *, 'Bohm-Vitense ML Version (inefficient)'
        print *, 'mixing length/pressure scale height ratio = ',aml
     else
        print *, 'Bohm & Cassinelli ML Version (efficient)'
        print *, 'mixing length/pressure scale height ratio = ',aml
     endif
  end if

  if ( evoloutput ) then
     write(17,*) ' threshold entered:',thresh
     write(20,*) ' threshold entered:',thresh
     write(17,*) '     deldg entered:',deldg
     write(20,*) '     deldg entered:',deldg
     write(17,*) '    ixswch entered:',ixswch
     write(20,*) '    ixswch entered:',ixswch
     write(17,*) '    nxmods entered:',nxmods
     write(20,*) '    nxmods entered:',nxmods
     if(iprep .eq. 1)then
        write(17,*) 'Prepped models? YES!' 
        write(20,*) 'Prepped models? YES!' 
     else
        write(17,*) 'Prepped models? NO!' 
        write(20,*) 'Prepped models? NO!'
     endif
     if(idiffus .eq. 1)then
        write(17,*) 'Diffusionesque transition zones'
        write(20,*) 'Diffusionesque transition zones'
     else
        write(17,*) 'Discontinuous transition zones'
        write(20,*) 'Discontinuous transition zones'
     endif
     if (icomp.eq.1) then
        write(17,337) 
!337     format(/' Layering of envelope is C/He4/H -- normal (DA) case'/)
     elseif (icomp.eq.2) then
        write(17,333) 
!333     format(/' Layering of envelope is C/He4/He3 -- a DB with a', &
!             ' surface Helium3 layer'/)
     elseif (icomp.eq.3) then
        write(17,334) 
!334     format(/' Layering of envelope is C/He4/He4. This is for test', &
!             ' purposes--it',/' should be equivalent to a normal DB model', &
!             ' with just C/He layering'/)
     else
	stop
     end if

     if(mls .eq. 0)then
        write(17,*) 'Bohm-Vitense ML Version (inefficient)'
        write(17,*) 'mixing length/pressure scale height ratio = ',aml
        write(20,*) 'Bohm-Vitense ML Version (inefficient)'
        write(20,*) 'mixing length/pressure scale height ratio = ',aml
     else
        write(17,*) 'Bohm & Cassinelli ML Version (efficient)'
        write(17,*) 'mixing length/pressure scale height ratio = ',aml
        write(20,*) 'Bohm & Cassinelli ML Version (efficient)'
        write(20,*) 'mixing length/pressure scale height ratio = ',aml
     endif
  end if

!*** initialize variables ***
  ifail = 0
  bold = 0.
  teold = 0.

!*** main program ***
  call read1
  call read2

  if ( idiffus .eq. 1 ) then
     call setprof
  endif

1 call write1
2 call begin
3 call calc(line)

!*** look for crashes ***
!Keep trying until we took too long
!t0 is the time in seconds that we called ff_da from the main program
  if (itcrashed.eq.1) then
     call check_time(newtime)
     if ((newtime - t0) .le. xmaxtime) then
        irestart = irestart + 1
     else
        ff=0.0
        goto 10
     endif
!*** look for hotspots ***
     if (irestart .gt. 10) then
        ihotspot = 1
!*** assign fitness of 0.0 ***
        goto 10
     endif
!*** resubmit with new timestep ***
     goto 666
  endif
  if ( line .eq. 2 ) then 
     callcshell = .true.
     call gshell(callcshell)
  elseif ( line .eq. 1 ) then
     callcshell = .false.
     call gshell(callcshell)
  endif
  if ( nite.eq.1) then
     call write2
  end if
  if ( it .lt. 0 ) then
     call interp
     goto 3
  else
     call end(line)
     call check_time(newtime)
     if ((newtime - t0) .gt. xmaxtime) then
        ittimedout=1
        goto 10
     endif
!*** avert hotspot crash if possible ***
     if (ihotspot.eq.1) then
        irestart = irestart + 1
        if (irestart .gt. 10) then
           goto 10
        endif
        !print *, 'evol 3'
        call check_time(newtime)
        if ((newtime - t0) .gt. xmaxtime) then
           goto 10
        endif
        goto 666
     endif
!***************************************
  endif
  if ( itconverged .eq. 1 ) goto 100
  if ( killit ) then
     if (verbose) then
        write(*,*) 'End of the line.'
        write(*,*) 'Terminus.'
     end if
     stop
  endif

  if (line .eq. 1) goto 1
  if (line .eq. 2) goto 2
  if (line .eq. 10) goto 10

10 continue
  
100  return

end subroutine evolve

