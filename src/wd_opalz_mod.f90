module kappa

  use kap_lib
  use kap_def
  use chem_lib
  use chem_def
  use const_def
  use const_lib
  use crlibm_lib

  implicit none

  character(len=256) :: kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix
  logical :: use_cache, cubic_X_interpolation, &
       cubic_Z_interpolation, include_electron_conduction, &
       use_Zbase_for_Type1
  integer :: handle2, ikfirst
  !  integer :: handle, i, ii
  integer :: ierr

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
  subroutine opalz(t,ro, xa, fkappa, opr, opt)

    real(dp) :: Xc, Xn, Xo, Xne, xc_base, xn_base, xo_base, xne_base, &
         zbar, frac_Type2, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
    real(dp), intent(in) :: t, ro
    real(dp), intent(out) :: fkappa, opr, opt
    real(dp) :: logRho, logT
    real(dp) :: Z_init, Z, type2_logT_lower_bdy, Zbase
    real(dp) :: kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy
    real(dp) :: kap_Type2_full_off_X = 0.71d0, kap_Type2_full_on_X = 0.70d0
    real(dp) :: kap_Type2_full_off_dZ = 0.001d0, kap_Type2_full_on_dZ = 0.01d0
    real(dp), dimension(7), intent(in) :: xa
    real(dp), dimension(7) :: amass, nion
    real(dp), parameter, dimension(7) :: Aion = [1., 4., 12., 14., 16., 20.18, 24.31]
    real(dp), parameter, dimension(7) :: Zion = [1., 2., 6., 7., 8., 10., 12.]
    character (len=32) :: my_mesa_dir

    kappa_file_prefix = 'gn93'
    !kappa_lowT_prefix = 'lowT_Freedman11'
    kappa_lowT_prefix = 'lowT_fa05_gn93'

    kappa_CO_prefix = '' ! use default
    kappa_lowT_prefix = '' ! use default
    kappa_blend_logT_upper_bdy = 0 ! use default
    kappa_blend_logT_lower_bdy = 0 ! use default
    type2_logT_lower_bdy = 0 ! use default
    use_cache = .false.
    ierr = 0
    use_Zbase_for_Type1 = .false.

    ! initialization and setup

    if (ikfirst .ne. 1) then

       write(*,*) 'loading MESA opacity tables...'
       ! my_mesa_dir = '../..'         
       my_mesa_dir = ''         
       call const_init(my_mesa_dir,ierr)     
       if (ierr /= 0) then
          write(*,*) 'const_init failed'
          stop 1
       end if

       call crlibm_init

       call chem_init('isotopes.data', ierr)
       call kap_init( &
            kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
            kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
            type2_logT_lower_bdy, use_cache, '', ierr) 
       if(ierr/=0) stop 'problem in kap_init'

       !next it is necessary to create a 'handle' for the general kap structure
       !using handles, it is possible to simultaneously access more than one
       !working copy of the opacity subroutines with different settings
       handle2 = alloc_kap_handle(ierr)
       if(ierr/=0) stop 'problem in alloc_kap_handle'


       !the final step is to set the type of X and Z interpolation, either linear or cubic spline
       !this is achieved by two logical variables, in this example both X and Z use cubic splines
       cubic_X_interpolation = .true.
       cubic_Z_interpolation = .true.
       include_electron_conduction = .true.
       call kap_set_choices( &
            handle2, cubic_X_interpolation, cubic_X_interpolation, &
            include_electron_conduction, &
            kap_Type2_full_off_X, kap_Type2_full_on_X, &
            kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
            ierr)
       if(ierr/=0) stop 'problem in kap_set_interpolation_choices'
       !end of initialization and setup

       ikfirst=1
    end if

!    zbar = 1.5d0              ! needed for electron conduction at high rho
    lnfree_e = 0d0            ! needed for Compton at high T
    d_lnfree_e_dlnRho = 0d0
    d_lnfree_e_dlnT = 0d0

    XC  = xa(3)
    XN  = xa(4)
    XO  = xa(5)
    XNe = xa(6)
    Z = 1d0 - (xa(1)+xa(2))

    nion(:) = xa(:)/Aion(:)
    nion(:) = nion(:)/sum(nion(:))
    zbar = sum(nion(:)*Zion(:))
    ! write(*,'(10f9.4)') xa(:)
    ! write(*,'(10f9.4)') nion(:)
    ! write(*,'(10f9.4)') Zion(:)
    ! write(*,'(10f9.4)') Aion(:)
    ! write(*,*) 'zbar = ',zbar

    logT=log10(t)
    logRho=log10(ro)
    Zbase = 0.d0

    !  subroutine opalz(t,ro, xa, fkappa, opr, opt)
    call kap_get_Type2( &
         handle2, zbar, xa(1), Z, Zbase, XC, XN, XO, XNe, logRho, logT, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, use_Zbase_for_Type1, &
         frac_Type2, fkappa, opr, opt, ierr)
    if(ierr/=0) write(*,*) 'kap_get_Type2 failed'



    !all finished? then deallocate the handle and unload the opacity tables
!    call free_kap_handle(handle2)
!    call kap_shutdown

  end subroutine opalz

end module kappa
