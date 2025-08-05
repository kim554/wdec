program getpar

  use datain
  use flags
  use outfile
  use pdotoutput
  use utils_subroutines
  
  
  implicit none

  real*8 :: fitness,ff
  integer :: starttime, i, wconfig
  character(1), parameter :: hconfig = "1"
  character(1) :: response, wconfig_char
  character(2) :: h3
  character(3) :: menv, mhe

  open(unit=90,file='radii.dat', position='append')
  open(unit=100,file='controlparams')
  open(unit=101,file='gridparameters')
  open(unit=102,file='calcperiods', position='append')

!Read control parameters (where user makes decision on how much output they
!wants
  read(100,*)
  read(100,*) response
!If yes, show parameters running on the screen
  if ( (response .eq. 'Y') .or. (response .eq. 'y') ) then
     screenoutput = .TRUE.
  else
     screenoutput = .FALSE.
  end if

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

!Special output for calculating pdots
  read(100,*)
  do i=1,2
     read(100,*) litosave(i), period1(i), period2(i)
  end do

!default values
!!$inputdata(1)  = 11500.!Effective temperature
!!$inputdata(2)  = 600. !Mass
!!$inputdata(3)  = 161  !Menv
!!$inputdata(4)  = 212. !Mhe
!!$inputdata(5)  = 412. !Mh
!!$inputdata(6)  = 0.69 !Helium abundance in mixed C/He/H region
!!$inputdata(7)  = 13.0  !Diffusion coefficient for He at the base of the envelope
!!$inputdata(8)  = 10.0  !Diffusion coefficient for He at the base of pure He
!!$inputdata(9)  = 0.96 !Alpha of MLT convection. Default prescription is
                     !ml2 (Bohm & Cassinelli). To toggle to Bohm-Vitense
                     !edit source file evol_mainroutine.f90 and search 
                     !for the variable mls and its assignment.
!!$inputdata(10) = Not used at the moment
  inputdata(11) = 99         !Edge of oxygen core (x100)
  inputdata(12) = 50.         !Minimum period to compute
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
!!$  inputdata(20) = 67.  !h1 x 100
!!$  inputdata(21) = 56.  !h2 x 100
!!$  inputdata(22) = 38.  !h3 x 100
!!$  inputdata(23) = 52.  !w1 x 100
!!$  inputdata(24) = compiler bug prevents this one from getting passed properly. Don't use.
!!$  inputdata(25) = 2.  !w2 x 100
!!$  inputdata(26) = 41.  !w3 x 100
!!$  inputdata(27) = 2. !w4 x 100
  
10 format(16F9.1)

  firstmod = .true.
  wconfig = 1

!read(101,*)  !First line contains header for the user
  open(unit=8000,file='testoutput',position='append')
  read(101,*) inputdata(1:9), inputdata(20:23), inputdata(25:27)
  write(*,10) inputdata(1:9), inputdata(20:23), inputdata(25:27)
  write(102,10) inputdata(1:9), inputdata(20:23), inputdata(25:27)
!build core description string
  write(h3,'(I2)') int(inputdata(22))
  write(wconfig_char, '(I1)') wconfig
  core_description = 'hset' // trim(hconfig) // '.h3.' // trim(h3) // '.mass' // trim(wconfig_char)
  write(menv,'(I3)') int(inputdata(3))
  write(mhe,'(I3)') int(inputdata(4))
  env_description = 'menv.' // trim(menv) // '.mhe.' // trim(mhe)
!  print *, 'xhe.' // env_description // '.dat'
  call check_time(starttime)
  inputdata(24) = real(starttime)
  call calculate_periods(inputdata)  
  write(102,*) 0.0
  write(102,*) 100000
  close(5)
  close(8000)
  firstmod = .false.
  wconfig = wconfig + 1

end program getpar
