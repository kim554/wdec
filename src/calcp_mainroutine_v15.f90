subroutine calculate_periods(data)

  !Common blocks
  use comp
  use outfile 
  use cperiods
  use mixl
  use startmod
  use ax
  use flagtmp
  use crash
  use fracxtal
  use flags, only: pulsemodel,makeplots,evoloutput,alphachoice
  use terp, only: stpms

  implicit none

  real*8, dimension(36), intent(in) :: data
  real*8 :: M_x, xtal, sigma, ff, best
  real*8 :: best_sigma
  integer :: nans,li,ell,type,i, stpmsid, M_tot, tnot
  character(15) :: naname
  character(80) :: cmd_string  
  
  ittimedout=0

! Set input parameters

  mode=0 
  fracm = 0.

  Teff = data(1)
  SM_star = data(2)/1000.
  M_env = data(3)/100.
  M_he = data(4)/100.
  M_h = data(5)/100.
!  M_h = MAX(data(5)/100., M_he+2.0)
  xhebar = data(6)/100.
  alph1 = data(7)
  alph2 = data(8)
  
!alpha for convection. It's complicated
!alphachoice = 1 - alpha is a free parameter, set in gridparameters
!alphachoice = 2 - For DAs, alpha is taken from Tremblay et al. 2015.
!                  For DBs, alpha is taken from non-linear light curve fitting
!alphachoice = 3 - alpha is taken from non-linear light curve fitting both for
!                  DBs and for DAs.

  select case (alphachoice)
  case(1)
     aml = data(9)
  case(2)
     if (M_h .lt. 20.) then
        call tremblay(Teff, SM_star)
     else
        aml = 0.96
     end if
  case(3)
     if (M_h .lt. 20.) then
        call nonlinearfitDA(Teff)
     else
        aml = 0.96
     end if
  end select

!data(10) not used
!data(11) not used
  permin = data(12)
  permax = data(13)
  lmax = INT(data(14))
  nper = INT(data(15))
  neuts = INT(data(16))
!data(17) not used
  maxion = data(18)
  stpms = data(19)
  tnot = INT(data(24))
      
  maxion = maxion*1.e-3 !meV -> eV
      
! Open proper starter model file

  M_tot = 10*INT(100.*SM_star)

11 format('masses/0',I3.3,'.',I4,'.dat')
  stpmsid = 9900   ! up to 99% stop mass
!  stpmsid = 9995   ! up to 99.94% stop mass
  
  if (SM_star .lt. 0.550) then
     write(start_file,11) 500, stpmsid
  else
     if (SM_star .lt. 0.650) then
        write(start_file,11) 600, stpmsid
     else
        if (SM_star .lt. 0.750) then
           write(start_file,11) 700, stpmsid
        else
           if (SM_star .lt. 0.850) then
              write(start_file,11) 800, stpmsid
           else
              if (SM_star .le. 1.440) then
                 write(start_file,11) 900, stpmsid
              else
                 print *, "White dwarf mass is super Chandrasekhar. Exiting."
                 stop
              end if
           end if
        end if
     end if
  end if

!*** evolve and pulsate the model ***

!!$write(*,*) "fname ", fname
!!$write(*,*) "start_file ", start_file
!!$write(*,*) "maxion ", maxion
!!$write(*,*) "MLT alpha ", aml
!!$write(*,*) "neuts ", neuts
!!$write(*,*) "itcrashed,ittimedout,irestart,ihotspot,istpms,mode,xheneg"
!!$write(*,*) itcrashed,ittimedout,irestart,ihotspot,istpms,mode,xheneg
!!$write(*,*) "xmaxtime, t0", xmaxtime, t0
!!$write(*,*) "permin, permax ", permin, permax
!!$write(*,*) "nper, lmax ", nper, lmax
!!$stop

!!$write(*,*) "Teff, SM_star ", Teff, SM_star
!!$write(*,*) "M_h, M_he, XO_c, XO_fm ", M_h, M_he, XO_c, XO_fm
!!$write(*,*) "hecdexp, w1, w3 ", hecdexp, w1, w3
!!$write(*,*) "h1, h2, tnot ", h1frac, h2frac,tnot
!!$stop

call evolve (tnot)

if (ittimedout.eq.1) then
   print *,  "Timed out, exiting"
   goto 100
endif
!*** avoid crashes by checking for NANs ***
if ( evoloutput )   then   
32 format("NANs.in",A8)
   write(naname,32) fname(8:15)
33 format("grep -ic NAN ",A15," > ",A15)
   write(cmd_string,33) fname,naname
   call system(cmd_string)
end if
!open(98,file=naname,status='old')
!read(98,*) nans
!close(98,status='delete')
!if (nans .gt. 0) ihotspot=1
if (ihotspot .eq. 1) then
   print *, "Model not converging, exiting"
   ff=0.0
   goto 100
endif
if (ihotspot .eq. 2) then
   print *, "Model not converging at evolution stage, exiting"
   ff=0.0
   goto 100
endif

!*** loop over crystallized mass fraction ***

best=10000.
best_sigma=10000.
M_x = 0.0

if ( pulsemodel .or. makeplots) then
   call pulsate       
!   write(*,*) calc_per, ncalc
end if

ff = 1.

100 return

end subroutine calculate_periods

!******************************************************************************
!Computes a value for alpha based on Tremblay et al. 2015, ApJ, 799, 142
!What's coded below is their fitting function given in equations 7, 8, and 9
!One major difficulty is that the formula requires logg, not known at this point
!The routine calls another subroutine that makes a guess for logg based on
!an analytic fit using DA models previously ran.

subroutine tremblay(teff,mass)

  use mixl, only : aml !The coveted alpha
  use startmod, only : logg_init
  
  implicit none

  real*8, intent(in) :: teff
  real*8, intent(in) :: mass
  real*8 :: go, to
  real*8, dimension(13) :: a

  data a / 1.1989083d+00, -1.8659403d+00, 1.4425660d+00, 6.4742170d-02, &
       -2.9996192d-02, 6.0750771d-02, -5.2572772d-02, 5.4690218d+00, &
       -1.6330177d-01, 2.8348941d-01, 1.7353691d+01, 4.3545950d-01, &
       -2.1739157d-01/ 

  call getloggda(teff,mass)

  go = logg_init - 8.0d0
  to = (teff - 12000.d0)/1000.d0 - 1.6d0 * go

  aml = (a(1) + (a(2) + a(3)*exp(a(4)*to + a(5)*go)) * &
       exp((a(6) + a(7)*exp(a(8)*to))*to + a(9)*go)) + &
       a(10) * exp(-a(11)*((to-a(12))**2 + (go-a(13))**2))
    
end subroutine tremblay

!******************************************************************************
!Provides an inital guess for logg based on the mass and effective temperature
!I ran a grid of DA models and for models of a fixed teff but varying masses
!graphed logg vs logmass (cgs). I made graphs for a range of temperatures.
!I empirically found that at around and below 20000, the best fit curve was
!a parabola, while things were best described by a line at higher temperatures
!Those two cases are coded up below.

!I determined the best fit parameters for all the curves generated. I got
!values for the curves at each effective temperature. I then took say all
!the a's of the quadratic fits and graphed a vs teff. Lastly I perforemd a
!linear fit to that.

!The result is a 2 variable fitting function. It fits a parabola for models
!20000K and below, and a line for hotter models. The parameters of the fitting
!themselves are calculated as a function of effective temperature.

!Since this was calibrated off DA models, it works best with DA models. It's
!off for DB models. A similar method could be used for DBs, if required.

subroutine getloggda(teff,mass)

  use misc_const, only : amsun
  use startmod, only : logg_init
  
  implicit none

  real*8, intent(in) :: teff,mass
  real*8 :: logm
  real*8 :: a, b, c

  logm = dlog10(mass*amsun)
  
  if (teff .le. 20000.) then
     a = -2.21581504e-05 * teff + 1.77725016e+00
     b =  1.49175069e-03 * teff - 1.15340771e+02
     c = -2.51058836e-02 * teff + 1.87871260e+03
     logg_init = a * logm**2 + b * logm + c
  else
     a =  2.87099333e-05 * teff + 2.31671805e+00
     b = -9.57665911e-04 * teff - 6.85177342e+01
     logg_init = a * logm + b
  end if

  return
  
end subroutine getloggda

!******************************************************************************

subroutine nonlinearfitDA(teff)

  use mixl, only : aml !The coveted alpha

  implicit none

  real*8, intent(in) :: teff
  real*8, dimension(3) :: a
  real*8, parameter :: small = 5.d-1

  data a / -3.37364797d-08,   1.00707231d-03,  -6.46188378d+00 /
  
  aml = a(1)*teff**2 + a(2)*teff + a(3)
  if ((teff .lt. 15000.d0) .and. (aml .lt. 0.d0)) aml = small
  if ((teff .ge. 15000.d0) .and. (aml .lt. 1.d0)) aml = 1.0d0
  
end subroutine nonlinearfitDA

!******************************************************************************
