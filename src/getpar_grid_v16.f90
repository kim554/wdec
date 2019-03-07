program getpar

  use datain
  use flags
  use utils_subroutines
  
  
  implicit none

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
!!$inputdata(1)  = 11500.!Effective temperature
!!$inputdata(2)  = 600. !Mass
!!$inputdata(3)  = 150  !Menv
!!$inputdata(4)  = 200. !Mhe
!!$inputdata(5)  = 400. !Mh
!!$inputdata(6)  = 0.70 !Helium abundance in mixed C/He/H region
!!$inputdata(7)  = 2.0  !Diffusion coefficient for He at the base of the envelope
!!$inputdata(8)  = 2.0  !Diffusion coefficient for He at the base of pure He
!!$inputdata(9)  = 0.65 !Alpha of MLT convection. Default prescription is
                     !ml2 (Bohm & Cassinelli). To toggle to Bohm-Vitense
                     !edit source file evol_mainroutine.f90 and search 
                     !for the variable mls and its assignment.
!!$inputdata(10) = Not used at the moment
  inputdata(11) = 99         !Edge of oxygen core (x100)
  inputdata(12) = 100.         !Minimum period to compute
  inputdata(13) = 1500.        !Maximum period to compute
  inputdata(14) = 2.           !Maximum ell to compute
  inputdata(15) = 1000.        !Leave that at 1000, no less than 500.
  inputdata(16) = 3.           !Neutrinos. Anything other than 0 turns them on.
! inputdata(17) = Not used at the moment        

  inputdata(18) = 0.           !Axion mass.
  inputdata(19) = -4.365d-3    !(10**inputdata(19) = 99.0%) Stop mass we want to use
                            !Furthest out allowed = 99.94%. Best to stay around
                            !99%. If wanting to use more than 99, you need
                            !to edit the source file calcp_mainroutine.f90 and
                            !change the value of the variable stpmsid to 9995.
!Core parameters (see Fig. 1 of https://zenodo.org/record/1715917#.XH8E6XXwYrg)
!w4 defined in terms of w1, w2, w3 and edge of oxygen core (stop mass), which is fixed by inputdata(19)  
!!$  inputdata(20) = 70.  !h1 x 100
!!$  inputdata(21) = 50.  !h2, expressed as a percent of h1
!!$  inputdata(22) = 40.  !h3, expressed as a percent of h2
!!$  inputdata(23) = 50.  !w1 x 100
!!$  inputdata(24) = compiler bug prevents this one from getting passed properly. Don't use.
!!$  inputdata(25) = 10.  !w2 x 100
!!$  inputdata(26) = 20.  !w3 x 100
  
10 format(15F9.1)

firstmod = .true.
!read(101,*)  !First line contains header for the user
  do
     open(unit=8000,file='testoutput',position='append')
     read(101,*,end=1) inputdata(1:9), inputdata(20:23), inputdata(25:26)
     write(*,10) inputdata(1:9), inputdata(20:23), inputdata(25:26)
     write(102,10) inputdata(1:9), inputdata(20:23), inputdata(25:26)
     call check_time(starttime)
     inputdata(24) = real(starttime)
     call calculate_periods(inputdata)  
     write(102,*) 0.0
     write(102,*) 100000
     close(5)
     close(8000)
     firstmod = .false.
  end do

1 continue                

end program getpar
