program getpar

  use utils_subroutines
  use flags

  implicit none

  real*8, dimension(36) :: datain
  real*8, dimension(7) :: minval, maxval, step
  real*8 :: fitness,ff
  integer :: starttime
  character(1) :: response

  open(unit=100,file='controlparams')
  open(unit=101,file='gridparameters')
  open(unit=102,file='calcperiods')

!Read control parameters (where user makes decision on how much output he/she
!wants
  read(100,*)
  read(100,*) response
!If yes, then use the new EOS tables from MESA
  if ( (response .eq. 'Y') .or. (response .eq. 'y') ) then
     ios_new = .TRUE.
  else
     ios_new = .FALSE.
  end if

  read(100,*)
  read(100,*) response
!If yes, then use the new opacities from MESA
  if ( (response .eq. 'Y') .or. (response .eq. 'y') ) then
     kap_new = .TRUE.
  else
     kap_new = .FALSE.
  end if

  read(100,*)
  read(100,*) alphachoice
!=1 - alpha is a free parameter, set in gridparameters
!=2 - For DAs, alpha is taken from Tremblay et al. 2015.
!    For DBs, alpha is taken from non-linear light curve fitting
!=3 - alpha is taken from non-linear light curve fitting both for DBs and DAs
  
  read(100,*)
  read(100,*) response
!If yes, then screen output from evolution code
  if ( (response .eq. 'Y') .or. (response .eq. 'y') ) then
     verbose = .TRUE.
  else
     verbose = .FALSE.
  end if

  read(100,*)
  read(100,*) response
!If yes, then write tape files from evolution code
  if ( (response .eq. 'Y') .or. (response .eq. 'y') ) then
     evoloutput = .TRUE.
  else
     evoloutput = .FALSE.
  end if

  read(100,*)
  read(100,*) response
!If yes, then write output files of prep code
  if ( (response .eq. 'Y') .or. (response .eq. 'y') ) then
     makeplots = .TRUE.
  else
     makeplots = .FALSE.
  end if
  read(100,*)
  read(100,*) response
!If yes, then we are prepping and pulsing the model
  if ( (response .eq. 'Y') .or. (response .eq. 'y') ) then
     pulsemodel = .TRUE.
  else
     pulsemodel = .FALSE.
  end if
  read(100,*)
  read(100,*) response
!If yes, then write output files of pulsation code
  if ( (response .eq. 'Y') .or. (response .eq. 'y') ) then
     pulsoutput = .TRUE.
  else
     pulsoutput = .FALSE.
  end if
!Kernels, under construction
  read(100,*)
  read(100,*) !kernels !If yes, then write out kernel files

!default values
!!$datain(1)  = 11500.!Effective temperature
!!$datain(2)  = 600. !Mass
!!$datain(3)  = 150  !Menv
!!$datain(4)  = 200. !Mhe
!!$datain(5)  = 400. !Mh
!!$datain(6)  = 0.70 !Helium abundance in mixed C/He/H region
!!$datain(7)  = 2.0  !Diffusion coefficient for He at the base of the envelope
!!$datain(8)  = 2.0  !Diffusion coefficient for He at the base of pure He
!!$datain(9)  = 0.65 !Alpha of MLT convection. Default prescription is
                     !ml2 (Bohm & Cassinelli). To toggle to Bohm-Vitense
                     !edit source file evol_mainroutine.f90 and search 
                     !for the variable mls and its assignment.
!!$datain(10) = Not used at the moment
!!$datain(11) = Not used at the moment
  datain(12) = 100.         !Minimum period to compute
  datain(13) = 1500.        !Maximum period to compute
  datain(14) = 2.           !Maximum ell to compute
  datain(15) = 1000.        !Leave that at 1000, no less than 500.
  datain(16) = 3.           !Neutrinos. Anything other than 0 turns them on.
! datain(17) = Not used at the moment        

  datain(18) = 0.           !Axion mass.
  datain(19) = -4.365d-3    !(10**datain(19) = 99.0%) Stop mass we want to use
                            !Furthest out allowed = 99.94%. Best to stay around
                            !99%. If wanting to use more than 99, you need
                            !to edit the source file calcp_mainroutine.f90 and
                            !change the value of the variable stpmsid to 9995.

10 format(8F9.1)

firstmod = .true.
!read(101,*)  !First line contains header for the user
  do
     open(unit=8000,file='testoutput',position='append')
     read(101,*,end=1) datain(1:8)
     write(*,10) datain(1:8)
     write(102,10) datain(1:8)
     call check_time(starttime)
     datain(24) = real(starttime)
     call calculate_periods(datain)  
     write(102,*) 0.0
     write(102,*) 100000
     close(5)
     close(8000)
     firstmod = .false.
  end do

1 continue                

end program getpar
