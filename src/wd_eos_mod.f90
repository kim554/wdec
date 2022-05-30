
module wd_eos_mod
  use eos_def
  use eos_lib
  use chem_def
  use chem_lib
  use const_lib
  use crlibm_lib

  implicit none

  real(dp) :: X, Z, Y, abar, zbar, z2bar, ye
  integer, parameter :: species = 7
  integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
  integer, pointer, dimension(:) :: net_iso, chem_id
  real(dp) :: xa(species)
  character (len=256) :: my_mesa_dir

contains

  !*****************************************************************************************************
  ! This subroutine is a wrapper for the MESA EOS subroutines. It is designed to be used with the 
  ! white dwarf evolution code (WDEC) so it only considers the elements, H1, He4, C12, and O16. It would
  ! be easy to generalize it to larger element lists. 
  ! The input parameters are
  ! in1  = Rho or P, depending on eflag
  ! temp = T (in Kelvin)
  ! xavec(:) contains the mass fractions of different elements, as given below:
  !    xavec(1)  = H1   mass fraction
  !    xavec(2)  = He4  mass fraction
  !    xavec(3)  = C12  mass fraction
  !    xavec(4)  = N14  mass fraction
  !    xavec(5)  = O16  mass fraction
  !    xavec(6)  = Ne20 mass fraction
  !    xavec(7)  = Mg24 mass fraction
  ! eflag  = 'd' or 'P'
  !    This parameter "eflag" controls whether the tables are called with (rho,T) or (P,T):
  ! eflag = 'd'  --- (rho,T) is assumed, so "in1" is taken to be the density (in cgs).
  !                  In this case the returned value of res2(1) = Pgas.
  ! eflag = 'P'  --- (P,T) is assumed, so "in1" is taken to be the density (in cgs).
  !                  In this case the returned value of res2(1) = Rho.
  !
  ! The output values are stored in the array res2(), and are defined as 
  !  res2(1)  = Pgas or Rho (depending on eflag) 
  !  res2(2)  = Cp 
  !  res2(3)  = Cv 
  !  res2(4)  = grad_ad
  !  res2(5)  = chiRho
  !  res2(6)  = chiT
  !  res2(7)  = cs
  !  res2(8)  = E  (internal energy)
  !  res2(9)  = eta (degeneracy parameter)
  !  res2(10) = S (entropy) 
  !  res2(11) = gamma1
  !*****************************************************************************************************

  subroutine wd_eos(in1,temp,xavec,eflag,res2)

    integer, save :: handle, ifirst
    real(dp) :: Rho, T, Pgas, log10Rho, log10T
    real(dp), intent(in) :: in1, temp, xavec(species)
    real(dp), dimension(11), intent(out) :: res2
    real(dp) :: dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, d_dlnRho_const_T, d_dlnT_const_Rho
    real(dp) :: p, dpdt, dpdro, gamma1, cs, out
    real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT, d_dabar, d_dzbar
    integer :: ierr
    character (len=1) :: eflag

    ierr = 0

    if (ifirst.ne.1) then
       my_mesa_dir = '' ! if empty string, uses environment variable MESA_DIR
       call const_init(my_mesa_dir,ierr)     
       if (ierr /= 0) then
          write(*,*) 'const_init failed'
          stop 1
       end if

       call crlibm_init

       call chem_init('isotopes.data', ierr)
       if (ierr /= 0) then
          write(*,*) 'failed in chem_init'
          stop 1
       end if

       ! allocate and initialize the eos tables
       call Setup_eos(handle)

       allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
       if (ierr /= 0) stop 'allocate failed'

       net_iso(:) = 0

       chem_id(h1) = ih1; net_iso(ih1) = h1
       chem_id(he4) = ihe4; net_iso(ihe4) = he4
       chem_id(c12) = ic12; net_iso(ic12) = c12
       chem_id(n14) = in14; net_iso(in14) = n14
       chem_id(o16) = io16; net_iso(io16) = o16
       chem_id(ne20) = ine20; net_iso(ine20) = ne20
       chem_id(mg24) = img24; net_iso(img24) = mg24

       ifirst=1
    endif

    xa(:) = xavec(:)

    if (sum(xa(:)) > 1.0d0) then
       !write(*,*) 'Warning: Sum(xa(:)) > 1',sum(xa(:))
       !write(*,*) 'xa(:):',xa(:)
       xa(:) = xa(:)/sum(xa(:))
    endif

    X = xa(1)
    Z = 1.d0 - X - xa(2)
    
    call Init_Composition


    T = temp
    log10T = log10_cr(T)

    if (eflag == 'd') then
       Rho = in1

       call eosDT_get( &
            handle, Z, X, abar, zbar,  &
            species, chem_id, net_iso, xa, &
            Rho, log10_cr(Rho), T, log10_cr(T),  &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)

       Pgas = exp_cr(res(i_lnPgas))
       out = Pgas

    else if (eflag == 'p') then
       Pgas = in1

       call eosPT_get(  &
            handle, Z, X, abar, zbar,   &
            species, chem_id, net_iso, xa, &
            Pgas, log10_cr(Pgas), T, log10_cr(T),  &
            Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,  &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)

       out = Rho
    else
       write(*,*) 'Error: unrecognized value for eflag: ',eflag
       write(*,*)
       stop

    endif

    if (ierr /= 0) then
       write(*,*) 'bad result from eos_get'
       write(*,*) T, in1
       write(*,*) Rho, Pgas
       write(*,*) xa(:)
       stop 1
    end if

    p = Pgas
    dpdro=(p/rho)*res(i_chiRho)
    dpdt =  (p/t)*res(i_chiT)
    gamma1=res(i_gamma1)
    cs=sqrt(gamma1*p/rho)
         
    res2(1)=out
    res2(2)=res(i_Cp)
    res2(3)=res(i_Cv)
    res2(4)=res(i_grad_ad)
    res2(5)=res(i_chiRho)
    res2(6)=res(i_chiT)
    res2(7)=cs
    res2(8)=exp(res(i_lnE))
    res2(9)=res(i_eta)
    ! res2(9)=res(i_chiT)/res(i_chiRho)
    res2(10)=exp(res(i_lnS))
    res2(11)=gamma1
    

  end subroutine wd_eos


  subroutine Setup_eos(handle)
    ! allocate and load the eos tables
    use eos_def
    use eos_lib
    integer, intent(out) :: handle

    character (len=256) :: eos_file_prefix
    integer :: ierr
    logical, parameter :: use_cache = .true.

    eos_file_prefix = 'mesa'

    call eos_init(eos_file_prefix, '', '', '', use_cache, ierr)
    if (ierr /= 0) then
       write(*,*) 'eos_init failed in Setup_eos'
       stop 1
    end if

    write(*,*) 'loading MESA eos tables...'

    handle = alloc_eos_handle(ierr)
    if (ierr /= 0) then
       write(*,*) 'failed trying to allocate eos handle'
       stop 1
    end if

  end subroutine Setup_eos


  subroutine Shutdown_eos(handle)
    use eos_def
    use eos_lib
    integer, intent(in) :: handle
    call free_eos_handle(handle)
    call eos_shutdown
  end subroutine Shutdown_eos


  subroutine Init_Composition
    use chem_lib

!    real(dp), parameter :: Zfrac_C = 0.173312d0
!    real(dp), parameter :: Zfrac_N = 0.053177d0
!    real(dp), parameter :: Zfrac_O = 0.482398d0
!    real(dp), parameter :: Zfrac_Ne = 0.098675d0

    real(dp) :: frac, dabar_dx(species), dzbar_dx(species),  &
         sumx, xh, xhe, xz, mass_correction, dmc_dx(species)

!    net_iso(:) = 0

!    chem_id(h1) = ih1; net_iso(ih1) = h1
!    chem_id(he4) = ihe4; net_iso(ihe4) = he4
!    chem_id(c12) = ic12; net_iso(ic12) = c12
!    chem_id(n14) = in14; net_iso(in14) = n14
!    chem_id(o16) = io16; net_iso(io16) = o16
!    chem_id(ne20) = ine20; net_iso(ine20) = ne20
!    chem_id(mg24) = img24; net_iso(img24) = mg24

!    Y = 1 - (X + Z)

!    xa(h1) = X
!    xa(he4) = Y
!    xa(c12) = Z * Zfrac_C
!    xa(n14) = Z * Zfrac_N
!    xa(o16) = Z * Zfrac_O
!    xa(ne20) = Z * Zfrac_Ne
!    xa(species) = 1 - sum(xa(1:species-1))

    call composition_info( &
         species, chem_id, xa, xh, xhe, xz, abar, zbar, z2bar, ye,  &
         mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)

  end subroutine Init_Composition


end module wd_eos_mod

