module kappa
  use kap_lib
  use kap_def
  use chem_lib
  use chem_def
  use const_def
  use const_lib
  use math_lib
  use utils_lib, only: mesa_error

  implicit none
  !this program demonstrates how to use mesa/kap in a stellar structure code
  !it reads in a mesa/star model of an AGB star so that it makes full use
  !of mesa/kap's capabilities. modules kap_lib and kap_def contain access to
  !all of the pieces required to set up and use mesa/kap.

  !  character(len=256) :: model_file, output_file
  integer :: handle2, ifirst
  integer, parameter :: NSpec = 7
  logical :: use_cache
  type (Kap_General_Info), pointer :: rq2
  integer :: ierr

  integer, parameter :: h1 = 1
  integer, parameter :: h2 = 2
  integer, parameter :: he3 = 3
  integer, parameter :: he4 = 4
  integer, parameter :: li7 = 5
  integer, parameter :: be7 = 6
  integer, parameter :: b8 = 7
  integer, parameter :: c12 = 8
  integer, parameter :: c13 = 9
  integer, parameter :: n13 = 10
  integer, parameter :: n14 = 11
  integer, parameter :: n15 = 12
  integer, parameter :: o16 = 13
  integer, parameter :: o17 = 14
  integer, parameter :: o18 = 15
  integer, parameter :: f19 = 16
  integer, parameter :: ne20 = 17
  integer, parameter :: ne21 = 18
  integer, parameter :: ne22 = 19
  integer, parameter :: na22 = 20
  integer, parameter :: na23 = 21
  integer, parameter :: mg24 = 22
  integer, parameter :: mg25 = 23
  integer, parameter :: mg26 = 24
  integer, parameter :: al26 = 25
  integer, parameter :: al27 = 26
  integer, parameter :: si28 = 27
  integer, parameter :: si29 = 28
  integer, parameter :: si30 = 29
  integer, parameter :: p31 = 30
  integer, parameter :: s32 = 31

!  real(dp) :: Mstar, Xc, Xn, Xo, Xne, xc_base, xn_base, xo_base, xne_base, &
!       lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
!       eta, d_eta_dlnRho, d_eta_dlnT
  real(dp) :: Z_init
!  real(dp) :: lnRho(maxpts), lnT(maxpts), logRho(maxpts), logT(maxpts), X(maxspec,maxpts)
!  real(dp) :: X(maxspec,maxpts)
!  real(dp) :: lnR(maxpts), L, dq(maxpts)
!  real(dp) :: kappa(maxpts), dlnkap_dlnRho(maxpts), dlnkap_dlnT(maxpts)
!  real(dp) :: kappa_fracs(num_kap_fracs), dlnkap_dxa(maxspec)
!  real(dp) :: kappaCO(maxpts), dlnkapCO_dlnRho(maxpts), dlnkapCO_dlnT(maxpts)
!  real(dp) :: kappaCO_fracs(num_kap_fracs), dlnkapCO_dxa(maxspec)
  character (len=32) :: my_mesa_dir

  integer, dimension(:), pointer :: chem_id, net_iso

contains

  !*****************************************************************************************************
  ! This subroutine is a wrapper for the MESA opacity/kappa subroutines. It is designed to be used with 
  ! the white dwarf evolution code (WDEC) so it only considers the elements, H1, He4, C12, N14, O16, Ne20,
  ! and Mg24 . It would be easy to generalize it to larger element lists. 
  ! The input parameters are
  ! t = T (in Kelvin)
  ! ro = Rho (in cgs)
  ! xavec(:) contains the mass fractions of different elements, as given below:
  !    xavec(1)  = H1   mass fraction
  !    xavec(2)  = He4  mass fraction
  !    xavec(3)  = C12  mass fraction
  !    xavec(4)  = N14  mass fraction
  !    xavec(5)  = O16  mass fraction
  !    xavec(6)  = Ne20 mass fraction
  !    xavec(7)  = Mg24 mass fraction
  !
  ! The output values are the following:
  !  fkappa = opacity (cgs)
  !  opr = d_ln_opacity/d_ln_rho 
  !  opt = d_ln_opacity/d_ln_T 
  !*****************************************************************************************************
  subroutine opalz(t, ro, xa, fkappa, opr, opt)
    real(dp), intent(in) :: t, ro
    real(dp), dimension(NSpec), intent(in) :: xa
    real(dp), intent(out) :: fkappa, opr, opt
    real(dp) :: logt, logro
    real(dp) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, eta, d_eta_dlnRho, d_eta_dlnT
    real(dp) :: kappa_fracs(num_kap_fracs), dlnkap_dxa(NSpec)
    integer :: i

    ! initialization and setup
    if (ifirst .ne. 1) then
       write(*,*) 'loading MESA opacity tables...'
       
       ierr = 0
       use_cache = .false.
       Z_init = 0.0

       my_mesa_dir = ''
       call const_init(my_mesa_dir,ierr)     
       if (ierr /= 0) then
          write(*,*) 'const_init failed'
          call mesa_error(__FILE__,__LINE__)
       end if

       call math_init()

       call chem_init('isotopes.data', ierr)
       call kap_init(use_cache, '', ierr) 
       if(ierr/=0) call mesa_error(__FILE__,__LINE__,'problem in kap_init')

       !next it is necessary to create a 'handle' for the general kap structure
       !using handles, it is possible to simultaneously access more than one
       !working copy of the opacity subroutines with different settings
!       handle1 = alloc_kap_handle_using_inlist('inlist_sample', ierr)
!       call kap_ptr(handle1, rq1, ierr)
!       rq1% use_Type2_opacities = .false.
!       rq1% Zbase = Z_init

       handle2 = alloc_kap_handle_using_inlist('inlist_sample', ierr)
       call kap_ptr(handle2, rq2, ierr)
       rq2% use_Type2_opacities = .true.
       rq2% Zbase = Z_init

       if(ierr/=0) call mesa_error(__FILE__,__LINE__,'problem in alloc_kap_handle')

       allocate(chem_id(NSpec), net_iso(num_chem_isos))
       !chem_id(:) = (/ ih1, ih2, ihe3, ihe4, ili7, ibe7, ib8, ic12, &
       !     ic13, in13, in14, in15, io16, io17, io18, if19, ine20, ine21, &
       !     ine22, ina22, ina23, img24, img25, img26, ial26, ial27, &
       !     isi28, isi29, isi30, ip31, is32 /)
       chem_id(:) = (/ ih1, ihe4, ic12, in14, io16, ine20, img24 /)
       net_iso(:) = 0
       do i = 1, NSpec
          net_iso(chem_id(i)) = i
       end do

       ifirst = 1

    endif

    ! these should come from an eos call
    lnfree_e = 0d0            ! needed for Compton at high T
    d_lnfree_e_dlnRho = 0d0
    d_lnfree_e_dlnT = 0d0
    eta = 0d0            ! needed for Compton at high T
    d_eta_dlnRho = 0d0
    d_eta_dlnT = 0d0

    logt  = log10(t)
    logro = log10(ro)

    call kap_get( &
         handle2, Nspec, chem_id, net_iso, xa, logro, logt, lnfree_e, d_lnfree_e_dlnRho, &
         d_lnfree_e_dlnT, eta, d_eta_dlnRho, d_eta_dlnT, kappa_fracs, fkappa, &
         opr, opt, dlnkap_dxa, ierr)
    if(ierr/=0) write(*,*) 'kap_get (Type 2) failed at i=', i
    ! write(*,*) 'kappa_fracs:',kappa_fracs

  end subroutine opalz

end module kappa
