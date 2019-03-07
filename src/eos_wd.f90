
      module wd_eos
      use eos_def
      use eos_lib
      use chem_def
      use chem_lib
      use const_lib
      use crlibm_lib

      implicit none

      double precision fkappa,opr,opt
      double precision :: X, Z, Y, abar, zbar, z2bar, ye, Rho, T
      integer, parameter :: species = 7
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
      integer, pointer, dimension(:) :: net_iso, chem_id
      double precision :: xa(species)
      character (len=256) :: data_dir


      contains
      
      subroutine wd_eos_DT(Rho,T,xh,xhe,xc,xo,res2)
         
        double precision, intent(in) :: Rho, T, xh, xhe, xc, xo
        double precision :: Pgas, log10Rho, log10T
        double precision :: dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, d_dlnRho_const_T, d_dlnT_const_Rho
        double precision, dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT, d_dabar, d_dzbar
        double precision, dimension(11) :: res2
        double precision :: p,dpdt,dpdro,gamma1,cs
        double precision :: dabar_dx(species), dzbar_dx(species), sumx, &
             xz, mass_correction, dmc_dx(species), xxh, xxhe
        integer, save :: handle, ifirst
        integer :: ierr

        if (ifirst.ne.1) then
           !data_dir = ''
           data_dir = '/Users/mikemon/mesa'
           call const_init(data_dir,ierr)
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
            
!     allocate and initialize the eos tables
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

        log10T = log10(T)

        xa(h1) = xh
        xa(he4) = xhe
        xa(c12) = xc
        xa(n14) = 0.0
        xa(o16) = xo
        xa(ne20) = 0.0
        xa(species) = 0.0

        if (sum(xa(:)) > 1.0) then
           write(*,*) 'Error: Sum(xa(:)) > 1',sum(xa(:))
        endif
         
        call composition_info( &
             species, chem_id, xa, xxh, xxhe, xz, abar, zbar, z2bar, ye,  &
             mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)

        X = xh
        Z = xz
         
        ! get a set of results for given temperature and density
        call eosDT_get( &
             handle, Z, X, abar, zbar, &
             species, chem_id, net_iso, xa, &
             Rho, log10_cr(Rho), T, log10T,  &
             res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         
        
1       format(a20,3x,e20.12)

        Pgas = exp_cr(res(i_lnPgas))

        ! deallocate the eos tables
        ! call Shutdown_eos(handle)

        p=exp_cr(res(i_lnPgas))

        dpdro=(p/rho)*res(i_chiRho)
        dpdt =  (p/t)*res(i_chiT)
        gamma1=res(i_gamma1)
        cs=sqrt(gamma1*p/rho)
         
        res2(1)=p
        res2(2)=res(i_Cp)
        res2(3)=res(i_Cv)
        res2(4)=res(i_grad_ad)
        res2(5)=dpdro
        res2(6)=dpdt
        res2(7)=cs
        res2(8)=exp(res(i_lnE))
        res2(9)=res(i_chiT)/res(i_chiRho)
        res2(10)=exp(res(i_lnS))
        res2(11)=gamma1
         
        if (ierr /= 0) then
           write(*,*) 'bad result from eos_get'
          ! stop 1
        end if

      end subroutine wd_eos_DT
      

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
         
         write(*,*) 'loading eos tables'
         
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

         real(dp), parameter :: Zfrac_C = 0.173312d0
         real(dp), parameter :: Zfrac_N = 0.053177d0
         real(dp), parameter :: Zfrac_O = 0.482398d0
         real(dp), parameter :: Zfrac_Ne = 0.098675d0
         
         real(dp) :: frac, dabar_dx(species), dzbar_dx(species),  &
               sumx, xh, xhe, xz, mass_correction, dmc_dx(species)
         
         net_iso(:) = 0
         
         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         
         Y = 1 - (X + Z)
               
         xa(h1) = X
         xa(he4) = Y
         xa(c12) = Z * Zfrac_C
         xa(n14) = Z * Zfrac_N
         xa(o16) = Z * Zfrac_O
         xa(ne20) = Z * Zfrac_Ne
         xa(species) = 1 - sum(xa(1:species-1))
         
         call composition_info( &
               species, chem_id, xa, xh, xhe, xz, abar, zbar, z2bar, ye,  &
               mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
         
      end subroutine Init_Composition

    end module wd_eos

