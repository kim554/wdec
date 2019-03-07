module evol_subroutines

contains

!************************************************************************

  subroutine openem

    use startmod
    use names
    use flags, only : evoloutput, firstmod

    implicit double precision(a-h,o-z)
    
    character(30) :: t40,out,log
  
10  format(a30)

    open(unit=7,file=start_file,status='old')

!*** read from input file ***
    read(7,10) t1
    read(7,10) t6
    read(7,10) t9
    read(7,10) out
    read(7,10) log
    read(7,10) t40
    read(7,10) t50
    
    open(unit=2,file='EEOSH',status='old')
    open(unit=3,file='EEOSHE',status='old')
    open(unit=4,file='EEOSC',status='old')
    open(unit=111,status='scratch',form='unformatted')
    open(unit=13,file='IEOSC',status='old',form='unformatted')
    open(unit=14,file='IEOSO',status='old',form='unformatted')
    open(unit=15,file='AUXIN5',status='old') 
    open(unit=25,file='SQOPAC',status='old') 

    if (evoloutput) then
       open(unit=11,file=t1)
       open(unit=26,file=t6)
       open(unit=9,file=t9)
       open(unit=17,file=out)
       open(unit=20,file=log)
       open(unit=40,file=t40)
       open(unit=51,file='neut.out',status='unknown') 
    end if
    
    open(unit=5,file='inputprof')
!    open(unit=10,file='initprof.dat')

    rewind(1)
    rewind(9)
    rewind(11)
    rewind(12)
    rewind(17)
    rewind(20)
    rewind(26)
    rewind(15)
    rewind(25)

    return
  end subroutine openem

!************************************************************************

  subroutine read1

!  read in the Heubner opacities and the envelope eos.
!  then compute the optical depths.

! Subroutines
    use eprep_subroutines, only: opac !a function

! Common blocks
    use misc_const
    use mixl
    use cc
    use tenv
    use jee
    use comp
    use comptmp
    use xcompp, only : amr_hyhe
    use flagtmp
    use d 
    use opacs
    use dirac
    use acc
    use aoo
    use alpha
    use opacfin
    use opcswch
    use flags
    use corats, only: irdold
    use startmod

    implicit double precision(a-h,o-z)
    
    real*8, dimension(7) :: dummy(7)
    integer, pointer :: j

    j => jj
  
! format statements 

2    format(i5,i3,f7.4,2i3)
31   format(1x,f5.1,2(1x,f10.6),40x)
1004 format(f16.6,8f8.2)
1005 format(8f10.4)
1006 format(1p,2e12.2,0p,4f8.4,' < amhyhe, amheca, alph(1-4)')

!*** transfer H and He boundaries from COMMON block ***
    amhyhe = amhyhe_tmp
    amheca = amheca_tmp
    amr_hyhe = 1.d0-10**-amhyhe_tmp

! read the mass boundaries and alpha's for diffusion calculation
! (changed, now read from input (first line))

    read(7,*) amhyhe_dum,amheca_dum,alph(1),alph(2),alph(3),alph(4)
! Using our own diffusion coefficients
    alph(1) = alph1
    alph(2) = alph2
    if ( evoloutput ) then
       write(11,1006) amhyhe,amheca, alph(1),alph(2),alph(3),alph(4)
       write(9,1006) amhyhe,amheca, alph(1),alph(2),alph(3),alph(4)
       write(17,1006) amhyhe,amheca, alph(1),alph(2),alph(3),alph(4)
       write(20,1006) amhyhe,amheca, alph(1),alph(2),alph(3),alph(4)
       write(26,1006) amhyhe,amheca, alph(1),alph(2),alph(3),alph(4)
    end if

!***If Salaris, set He/C diffusion exponents ***

    if (irdold .eq. 3) then
       alph(3) = hecdexp
       alph(4) = -hecdexp
    end if

!     read control parameters and the atomic number z and mass az
!     for hydrogen, helium, and carbon. 

    read(15,2) nmax,lmax,err,lmam,nlum 
    first = .true.
    ifinal = 0
    idif = 0
    az(1)=1.0
    az(2)=4.0
    az(3)=12.0
    z(1)=1.0
    z(2)=2.0
    z(3)=6.0

!     define the eos grid

    do nel = 1,3
       lunit = nel + 1
       read ( lunit, * ) dum
       do l = 1, lmax
          read (lunit,*) ml(nel,l)
          kmax = ml(nel,l)
          do k = 1,kmax 
             read (lunit,*) tk(nel,l),rhok(nel,l,k),pp(nel,l,k),uu(nel,l,k), &
                  xtt(nel,l,k),xrt(nel,l,k),datg(nel,l,k),eta(nel,l,k) 
          enddo
       enddo
       rewind lunit
    enddo

! read the Heubner radiative opacities cso on the density-
! temperature grid csr-cst
! Used to be Cox-Stewart opacities.

    read(15,1004) (cst(l),csr(l),dummy,l=1,29)
    read(15,1005) (((cso(jm,l,k),k=1,8),l=1,29),jm=1,3)

! compute optical depths(log) od for first lmam isotherms

    gs=1.e+08

    do jm=1,3
       je=jm
       do l=1,lmam
          kmax=ml(jm,l)
          tau=10.**(opac(rhok(jm,l,1),tk(jm,l))+pp(jm,l,1))/gs*0.5
          od(jm,l,1)=dlog10(tau)
          do k=2,kmax
             tau=tau+0.5/gs*(10.**pp(jm,l,k)-10.**pp(jm,l,k-1)) &
                  *(10.**opac(rhok(jm,l,k),tk(jm,l)) &
                  +10.**opac(rhok(jm,l,k-1),tk(jm,l)))
             od(jm,l,k)=dlog10(tau)
          enddo
       enddo
    enddo

!     read fermi-dirac f12,f32 functions in the range of degeneracy
!     parameter -4 < psi < 10 

    do k=1,141 
       read(15,31) psi(k),f32(k),f12(k)
    enddo
    
    if (verbose) print *, ' eos tables and tape5 read'
    
    return
  end subroutine read1

!************************************************************************

  subroutine read2

!  read model and some more control parameters

!Subroutines
    use istat_subroutines, only : istatco, istat1
    use utils_subroutines, only : xlocate, armove
    use chemprofiles_subroutines, only : comp_core, profsm_orig

!Common blocks
    use cc, only : smm   ! Stellar mass in solar masses
    use comp, only : x
    use comptmp
    use contrl
    use corats
    use crash
    use dfin
    use dvca
    use fin
    use flags, only : verbose, evoloutput, firstmod
    use flagtmp
    use gne
    use mixl
    use outfile
    use read2input
    use rw2
    use shells
    use temp
    use terp
    use tfitt
    use thermo 
    use times
    use surf
    use vca
    use xbnd
    use xcompp
    use xxnew

    implicit double precision(a-h,o-z)

    real*8, dimension(600,7) :: read2x, y
    real*8, pointer :: fn

! y has 7 columns. The column headers are:
! y(:,1) = log(m/M*)   = s(:)
! y(:,2) = log(r)      = r(:)
! y(:,3) = l/Lsun      = b(:)
! y(:,4) = log(P)      = p(:)
! y(:,5) = log(T)      = t(:)
! y(:,6) = s (entropy) = e(:)
! y(:,7) = xc          = xc(:)

! read2x and z are copies of y. read2x <-> sa..., y <-> s..., z <-> sk...
! z not used in this subroutine.

    fn   => fn2
    fname = "evolved"
    gnew=0.0

! Format Statements for READ2

!1 format(12i4)
2   format(f6.3,f10.6,e14.7,e10.3,f5.2,2e9.2,i4,f9.6)
3   format(2e16.8,4f8.5,i4,f5.2)
4   format(7f7.5) 
7   format(0pE15.8,F9.6,ES11.4,F10.6,F9.6,0pE14.7,F8.6)
8   format(4i10/2f10.6/(6d13.7))
20  format(e15.8,f9.6,e15.8,f9.6,2e15.8)
21  format(2f9.3,f10.6)

! subroutine TAPE reads in carbon and oxygen interior eos table

    call tape

! set main points of the core oxygen profile

    if (firstmod) then
       allocate(ams_o(ndimo), corat_o(ndimo))
    else
       ams_o(:) = 0.0d0
       corat_o(:) = 0.0d0
    end if

    ams_o(:)   = core_points(:,1)
    corat_o(:) = core_points(:,2)
    
! read some more parameters that shape the chemical profiles
    read(5,*)
    read(5,*)
    read(5,*) ao  !smoothing parameter for oxygen profile
    read(5,*)
    read(5,*) buffer_inner  !have to do with stitching of Helium profile
    read(5,*) buffer_outer  !see inputprof and documentation for more details

    rewind(5)
    
    if (evoloutput) then
       do i=1,ndimo
          write(10,*) -log10(1.d0-ams_o(i)), corat_o(i)
       end do
    end if

!  read in the model, and set some variables
!  note that the model is now tagged on the end of the input file

    !*** read from input file ***
    read(7,*)ip5,ip6,ip7,nmod,nite1,m,md,nu,ip8,ip1,irdold,ip40
    if (neuts .EQ. 0) nu = 0
    nite1 = 80
    !overide what is in the starter model file to force Salaris profiles
    irdold = 3

    if (verbose) then
       print *, 'nite ', nite1
       print *,    '             nmod:',nmod
       print *,    '             ip40:',ip40
       if (nu .eq. 0) then
          print *,    'neutrino-less models'
       else
          print *,    'Itoh etal neutrino rates'
       end if
    end if
    if ( evoloutput ) then
       write(20,*)    '             nmod:',nmod
       write(20,*)    '             ip40:',ip40
    end if

    iii = 0
888 iii = iii + 1
    
    

    if (iii .le. 2) read(7,*) ams(iii), corat(iii)

    if ( iii .gt. 1 ) then
       dcorat(iii) = corat(iii) - corat(iii - 1)
       delmass(iii) = ams(iii) - ams(iii-1)
    endif

    if ( ams(iii) .lt. 1.0 ) goto 888

    !*** hardwire npoints in C/O profile ***
    if (mode .NE. 0) iii = 4
    !***************************************
    ncore = iii

    !*** read from input file ***
    read(7,3) sm_dum,rat1,c,xh,yh,yyh,khomo,time
    sm = int(1.e6*dlog10(sm_tmp*1.989e33)+0.5)/1.e6
    read(7,4) ce,cif,sin,sout,smid,grid1,grid2
    read(7,2) dg,gx,ssg,stpms,f,tmax1,tmax2,modnr,g3

    ! set stop mass to whatever you want here
      if(istpms.eq.2) stpms = -0.0043648012 !(99.0%)
    ! if(istpms.eq.2) stpms = -0.007784 ! (98.2%)
    ! if(istpms.eq.2) stpms = -0.0132   ! (97%)
    !      if(istmps.eq.2) stpms = -0.0223   ! (95%)
    !      if(istpms.eq.2) stpms = -0.0269   ! (94%)

    ! change time step for crashed jobs

!    dg = dg - 0.1*float(irestart)
    
    tmax3=10.*tmax2
    gold=gx

    !*** read from input file ***
    read(7,20) speak,rpeak1,vpeak1,rlast,vexp1,egrav1
    read(7,21) bgrav,w,g1
    read(7,8) ks,ls,ms,kstart,rm,bm,u,v,ww

    !  here is where we read in the model.  if we are reading in
    !  an old model, then we also have to specify the run of the
    !  desired c/o profile.  (corat = 1 ==> pure c).

    read(7,*) jb
    do j = 1,jb
! special input format needed if working with an old input file
! use if you get the error "Bad value during floating point read" for unit 7
!       read(7,7) y(j,1), y(j,2), y(j,3), y(j,4), y(j,5), y(j,6), y(j,7)
       read(7,*) y(j,1), y(j,2), y(j,3), y(j,4), y(j,5), y(j,6), y(j,7)
    end do

! Save starter model values to common block "shells"
       
    s(:) = y(:,1)
    r(:) = y(:,2)
    b(:) = y(:,3)
    p(:) = y(:,4)
    t(:) = y(:,5)
    e(:) = y(:,6)
    xc(:)= y(:,7)

!  if irdold = 1, then do linear interpolation between the specified
!      mass points
!  if irdold = 2, then use the relation given in Barrat, Hansen,
!      & Mochkovitch 1988, A+A, 199, L15
!  if irdold = 3, for flexible chemical profiles (see comments under subroutine
!  ccomp). Even though with MESA equations of state and opacities we are no 
!  longer contrived by the core-envelope boundary, it still exists numericaly.
!  Here, we are only taking care of the core (but it is OK to put Helium in the
!  core). ABK (June 2016)
    
    if ( irdold .eq. 3 ) then

       do i = 1,jb
          amr  = 10.**y(i,1)      !!! mass coordinate in stellar masses
          call comp_core(i,amr)
          xcomp_orig = xcomp
          xc2 = xcomp(i,3)
          xo2 = xcomp(i,4)
          xc(i) = xcomp(i,3)
          xhe(i) = xcomp(i,2)
          y(i,7) = xcomp(i,3)

          if ( (xc2 .gt. .00001) .and. (xc2 .lt. .99999))then
             call istatco(p(i),t(i),0,.false.)
          elseif ( xc2 .ge. .999999) then
             call istat1(p(i),t(i),0,12,.false.)
          elseif ( xc2 .le. .000001) then
             call istat1(p(i),t(i),0,16,.false.)
          endif
          e(i) = e2
          y(i,6) = e2
       enddo
       close(7)

!  write header line to tape1 and tape9 
    elseif ( irdold .eq. 1 ) then
       
       do i = 1,jb
          amr  = 10.**y(i,1) 
          call profsm_orig(arm,ams(iii),corat(iii),ndim,xc2)
          iii = 2 
999       if ( amr .gt. ams(iii) ) then
             iii = iii + 1
             goto 999
          endif
          xc2 = dcorat(iii)/delmass(iii) * (amr - ams(iii-1)) + corat(iii-1)
          xo2 = 1 - xc2
          xc(i) = xc2
          y(i,7) = xc2
          if ( xc2 .gt. .00001 .and. xc2 .lt. .99999)then
             call istatco(p(i),t(i),0,.false.)
          elseif ( xc2 .ge. .999999) then
             call istat1(p(i),t(i),0,12,.false.)
          elseif ( xc2 .le. .000001) then
             call istat1(p(i),t(i),0,16,.false.)
          endif
          e(i) = e2
          y(i,6) = e2

       end do

       if ( verbose ) write(*,*) 'model read in and rescaled to c/o model'
       if ( evoloutput ) write(20,*) 'model read in and rescaled to c/o model'
       
    elseif ( irdold .eq. 2 ) then
 
       ndim = 10
       do  i = 1,jb
          amr  = 10.**y(i,1) 
          call profsm_orig(arm,ams(iii),corat(iii),ndim,xc2)
          if ( i .ne. jb ) then
             xo2 = .5 * ( 1. + alph ) * ( 1. - amr )**alph
          else
             xo2 = 0.
          endif
          xc2 = 1 - xo2
          xc(i) = xc2
          y(i,7) = xc2
          if ( xc2 .gt. .000001 .and. xc2 .lt. .999999) then
             call istatco(p(i),t(i),0,.false.)
          elseif ( xc2 .ge. .999999) then
             call istat1(p(i),t(i),0,12,.false.)
          elseif ( xc2 .le. .000001) then
             call istat1(p(i),t(i),0,16,.false.)
          endif
          e(i) = e2
          y(i,6) = e2
       end do
       if ( verbose ) write(*,*) 'model read in and rescaled to bhm c/o model'
       if (evoloutput) write(20,*) 'model read in and rescaled to bhm c/o model'
       corat(1) = xc(1)
       corat(2) = xc(jb)
       ncore = 2

    else

       do j=1,jb
          xc(j)=y(j,7)
       enddo
       
    endif
          
    smass=10.**(sm-33.298635)
6   format('*** h/he/c-o ',f4.2,'msun wd: dg=',f6.3,' stpms='e9.3, &
         ' co1=',f3.1,' co2=',f3.1,' ***')
    if ( evoloutput ) then
       write(9,6) smass, dg, stpms, corat(1), corat(ncore)
       write(26,6) smass, dg, stpms, corat(1), corat(ncore)
       write(11,6) smass, dg, stpms, corat(1), corat(ncore)
    end if

! Copy y into read2x
    call armove(read2x,y,jb,7)

!Copy read2x to common block
    sa(:) = read2x(:,1)
    ra(:) = read2x(:,2)
    ba(:) = read2x(:,3)
    pa(:) = read2x(:,4)
    ta(:) = read2x(:,5)
    ea(:) = read2x(:,6)
!    xca(:)= read2x(:,7)
    xca(:) = xcomp(:,3)  ! Carbon abundance
    xxhe(:) = xcomp(:,2) ! Helium abundance

! perform a homology transform if requested, store new model in name.out
    if(khomo.gt.0) then
       call homo(xh,yh,yyh,jb)
    endif
    ja=jb
    time1=time
    dg1=dg
    if ( w.lt.1. ) ip8=-1
    npass=0

    if ( evoloutput ) then
       write(40,*) '10**s2     r2     t2      p2    e2    b2      fn'
       write(40,*) '   cp(k)   cv(k)  o2  wc   w2   theta    psi2  rho'
    end if
    
    return

  end subroutine read2

!****************************************************************************

  subroutine tape

! read in the specification for the interior eos, and the eos itself. 
! new version here reads c eos from tape13 and o eos from tape14,
! while the control parameters are still read from tape5 (MAW)

    use tablecc
    use tableoo
    use sizec
    use sizeo
    use flags, only: verbose, evoloutput

    implicit double precision(a-h,o-z)

1000 format(8f10.2)
1001 format(8i10)

! read c interior eos

    read(15,1001) itc,ioverc 
    read(15,1000) (pmc(i),i=1,itc)
    read(15,1001) (ipc(i),i=1,itc)
    read(15,1001) (iphasec(i),i=1,itc) 
    read(15,1000) pminc,pmaxc,delpc
    read(15,1000) tminc,tmaxc,deltc
    read(13) tablec
    rewind 13
    close(13)
    if (verbose) write(*,100)
    if (evoloutput) write(20,100)
100 format('  carbon interior eos table read')
    
! read o interior eos
    
    read(15,1001) ito,iovero 
    read(15,1000) (pmo(i),i=1,ito)
    read(15,1001) (ipo(i),i=1,ito)
    read(15,1001) (iphaseo(i),i=1,ito) 
    read(15,1000) pmino,pmaxo,delpo
    read(15,1000) tmino,tmaxo,delto
    read(14) tableo
    rewind 14
    close(14)
    if (verbose) write(*,200) 
    if (evoloutput) write(20,200)
200 format('  oxygen interior eos table read')

    return

  end subroutine tape

!****************************************************************************

  subroutine homo(xh,yh,yyh,nshell) 
    
! homology transformation of wdec model 

! first group of relations are those that came with the
! code.  The others are empirically derived from
! pre-wd sequences of .601 and .7795 mo 

! Subroutines
    use istat_subroutines, only: istatco, istat1

! Common blocks
    use shells 
    use contrl
    use thermo 
    use temp
    use flags, only: verbose

    implicit double precision(a-h,o-z)

    real*8, pointer :: fn

    fn => fn2
    stt=2.*xh+5.*yh-yyh
    str=-xh-4.*yh+yyh
    stp=6.*xh+16.*yh-4.*yyh 
    stb=10.**(3.*xh+4.*yh)
  
    open(18,file='homo.dat',status='unknown') 

    do j=1,nshell
       t(j)=t(j)+stt
       r(j)=r(j)+str
       p(j)=p(j)+stp
       b(j)=b(j)*stb
       xc2 = xc(j)
       xhe2 = xhe(j)
       xo2 = 1.d0 - xc2 - xhe2
       if ( xc2 .gt. .000001 .and. xc2 .lt. .999999) then 
          call istatco(p(j),t(j),0,.true.)
       elseif ( xc2 .ge. .999999) then
          call istat1(p(j),t(j),0,12,.true.)
       elseif ( xc2 .le. .000001) then
          call istat1(p(j),t(j),0,16,.true.)
       endif
       e(j)=e2
       write(18,11)s(j),r(j),b(j),p(j),t(j),e(j),xc2
    enddo

11  format(e15.8,f9.6,1pe11.4,0pf10.6,f9.6,e14.7,f8.6)
    close(18)
    if (verbose) write(*,*) 'homology transform stored in homo.dat'
    stop
  end subroutine homo

!************************************************************************

  subroutine write1 

! write a (header) summary of a model, then set some
! variables to start a new sequence.

! Subroutines
    use istat_subroutines, only: istatco, istat1

! Common blocks
    use flags, only: verbose, evoloutput
    use rw2
    use rw
    use shells 
    use thermo 
    use surf
    use temp
    use contrl
    use xxnew
    use times
    use comptmp
    use flagtmp
    use terp
    use vca
    use dvca
    use xbnd
    use corats
    use fin
    use dfin
    use tfitt
    use gne
    use crash
    use outfile
    use read2input
    use xcompp

    implicit double precision(a-h,o-z)
    
    real*8, pointer, dimension(:) :: zp
    real*8, pointer :: fn
    real*8, dimension(600,7) :: z

103 format(///20h    **** model nr.  ,i4,7h  ***  ,f6.3, &
         43h solar masses (h - he - c/o ))  ***  age = ,1pd14.8, &
         11h years ****/)
104 format(//6h  ip5=i4,6h  ip6=i4,6h  ip7=i4,7h  nmod=i4,8h  nite1= &
         i4,4h  m=i4,5h  md=i4,5h  nn=i4,6h  ip8=i4,6h    w=f8.3)
114 format(60x,4hsin=f7.4,6h sout=f7.4,7h  smid=f7.4,7h grid1=f7.4, &
         7h grid2=f7.4/)
105 format(4h dg=f10.6,4h  g=f10.6,4h sg=e16.8,7h stpms=e11.4, &
         4h  f=f10.6,8h  tmax1=e10.4,8h  tmax2=e10.4//)
106 format  (5h  sm=e16.8,6h rat1=e16.8,5h   c=f6.3,6h   xh=f9.5, &
         6h   yh=f9.5,6h  yyh=f9.5,6h   ja=i4//)
    
    zp => sk
    fn => fn2
    z(:,1) = zp

    nite=nite1
    modnr=modnr+1 
    smass=10.**(sm-33.298635)
    if ( evoloutput ) then
       write(17, 103)modnr,smass,ssg
       write(17, 104) ip5,ip6,ip7,nmod,nite1,m,md,nu,ip8,w
       write(17, 114) sin,sout,smid,grid1,grid2
       write(17, 105) dg,gx,ssg,stpms,f,tmax1,tmax2
       write(17, 106) sm,rat1,c,xh,yh,yyh,ja
    end if
    fcc1=0.
    xdel=0.
    sc1=0.
    fn1=0.
    sn1=0.
    k1=0
    kon=0
    kom=0
    total=0.
    tmass=0.
    egrav=0.
    bcc=1e-10
    do i=1,5
       bn(i)=1e-10
    enddo
    bnt=1.e-10
    bax=1.e-30
    cste=c
    g2=dg-g3  
    do j=1,jb 
!       xc2 = xc(j)
!       xo2 = 1 - xc2
       xc2 = xcomp(j,3)
       xo2 = xcomp(j,4)
       xhe2 = xcomp(j,2)
       xh2 = xcomp(j,1)
       if ( xc2 .gt. .000001 .and. xc2 .lt. .999999) then 
          call istatco(p(j),t(j),0,.false.)
       elseif ( xc2 .ge. .999999) then
          call istat1(p(j),t(j),0,12,.false.)
       elseif ( xc2 .le. .000001) then
          call istat1(p(j),t(j),0,16,.false.)
       endif
       fca(j)=fcc/((1.-2.2892*fcc*10.**(gx-18.)/0.999999)**2)
    enddo
    
! fca is fractional C abundance at time t + delta t
! See eqn 3 of Kutter & Savedoff 1969, ApJ, 157, 1021.
    zp = z(:,1)

  end subroutine write1

!************************************************************************

  subroutine begin

!  initialize variables for central boundary condition

    use shells 
    use contrl
    use temp
    use fin

    implicit double precision(a-h,o-z)

    j=2 
    k=0 
    l=2 
    s1=-10.
    r1=-10.
    s2=s(1)
    r2=r(1)
    b2=b(1)
    p2=p(1)
    t2=t(1)
    ifin = 0
    return
  end subroutine begin
  
!************************************************************************

  subroutine calc(line)

!  set up variables for the step to the next shell.

! Subroutines
    use istat_subroutines, only : istatco, istat1
    use eprep_subroutines, only : opac

! Common blocks
    use misc_const, only: an0, bk_const, ck, pi, fnat, amsun, alsun, g
    use shells 
    use thermo 
    use temp
    use cal
    use xxnew
    use contrl
    use jee
    use comp
    use corats
    use wprep
    use dprep
    use bldxc
    use crash
    use xcompp
    use flags, only: verbose, evoloutput

    implicit double precision(a-h,o-z)
    
    real*8, pointer :: fn, gg
    logical :: converg


    fn => fn2
    gg => gx

    line = 2
    je = 3
    k=k+1
    if ( nite .eq. 1 ) then 
       converg = .true.
    else
       converg = .false.
    endif

! if memory full, then make program dump
402 format(' quit memory full',/'k = ',i4)
    if ( k .ge. 400) then
       if ( evoloutput ) then
          write(20,402) k
          write(17,402) k
       end if
       if ( verbose ) print 402, k
       itcrashed = 1
       goto 16
    endif

! otherwise, go on about our business.

    sk(k)  = s2
    rk(k)  = r2
    bk(k)  = b2
    pk(k)  = p2
    tk(k)  = t2
    kk=k

420 if ( l .le. ja ) then
       if ( s2 .le. sa(l) ) then
          sfrac = (s2-sa(l))/(sa(l)-sa(l-1))
          ea2   = ea(l)+(ea(l)-ea(l-1))*sfrac
          xc2   = xca(l)+(xca(l)-xca(l-1))*sfrac
          xhe2  = xxhe(l)+(xxhe(l)-xxhe(l-1))*sfrac
          xo2   = 1.d0-xc2-xhe2
          fca2  = fca(l)+(fca(l)-fca(l-1))*sfrac
       else
          l=l+1
          goto 420
       endif
    elseif ( w1 .lt. 0.4 ) then
       ea2=e1-1.914535e+08*dlog10(1.+ds/s1)*(1.-2.5*w1)/cei
    endif

600 c=cste      
    xc(k) = xc2
    xhe(k) = xhe2

    if ( s2 .gt. -sin .and. s2 .lt. -sout ) then
       c=c*(dabs(smid+s2)*grid2+grid1) 
    endif
    
    cei=1./(.5+1./ci)

    if ( xc2 .gt. .000001 .and. xc2 .lt. .999999) then
       call istatco(p2,t2,0,converg)
    elseif ( xc2 .ge. .999999) then
       call istat1(p2,t2,0,12,converg)
    elseif ( xc2 .le. .000001) then
       call istat1(p2,t2,0,16,converg)
    endif

    e(k)=e2
    fcc=(1.-w)*fca2+w*fcc
    if(k.le.1) ucent=u2
    q2=10.**(-3.*r2+s2-u2-1.099210+sm)
    f2=10.**(-4.*r2-p2+2.*s2-8.275214+2.*sm)
    w2=o2*b2*10.**(p2-4.*t2-s2+43.18734-sm)
    wc=-ep2/et2
    if (w2 .ge. wc) then
       w2 = wc
       if (kon .le. 0) then 
          c = c/ck
          kon = k 
       endif
    else
       if ( kon .gt. 0 ) then
          kom = kon
          kon = 0 
          c = c * ck
       endif
    endif
    qx=10.**(s2-33.22885+sm)
    qn2=(fcc-fn-fa)*qx
    qnp2=-fp*qx
    qnt2=-qx*ft
    qe2=qx*10.**(t2-gg)
    if ( k .gt. 1 ) line = 1
!  Now load up the variables in aa() if we're calculating 
!  the converged model

    if ( converg .and. iprep .eq. 1 ) then
       aa( 1,k) = 10.**r2
       star=10.0**(sm)
       aa( 2,k) = 10.**s2 * star
       aa( 3,k) = b2 * alsun
       aa( 4,k) = 10.**t2
       aa( 5,k) = 10.**u2
       aa( 6,k) = 10.**p2
       aa( 7,k) = fn
       cpx  = ut2/wc*10.**(-u2+p2-t2)
       xt   = -ut2/up2
       cv   = cpx*(1.e0-(xt*wc))
       aa( 8,k) = -cv
       aa( 9,k) = 1./up2
       aa(10,k) = xt
       aa(11,k) = ft
       aa(12,k) = fp/up2

! Paul's opac business (Will have to go to opac3 when comp is mixed)
! only calling opac for carbon here, because we lack O table. 
! Rad opac contribution is zilch here anyway.
! We use 5pt differencing in the core.

       cent = 0.01
       dent = 0.02
       oplo4=opac(u2+dent,t2)
       oplo1=opac(u2+cent,t2)
       oplo2=opac(u2-cent,t2)
       oplo3=opac(u2-dent,t2)
!*** kapr ***
       aa(13,k)=(oplo3-(8.*oplo2)+(8.*oplo1)-oplo4)/(12.*0.01)
       oplo2=opac(u2,t2+dent)
       oplo3=opac(u2,t2+cent)
       oplo4=opac(u2,t2-cent)
       oplo1=opac(u2,t2-dent)
       !*** kapt ***
       aa(14,k)=(oplo1-(8.*oplo4)+(8.*oplo3)-oplo2)/(12.*0.01)
       aa(15,k) = w2
       aa(16,k) = wc
       aa(17,k) = xxhe(k) ! Helium abundance
       aa(18,k) = o2
       aa(20,k) = xo2
!*** this is core ledoux term. (set to zero, then reset if need be) ***
       aa(19,k) = 0.0
!*** if at center, set Ledoux term to 0 by definition. ***
       if(k.eq.1)then
          goto 16
       endif
       if(xc2.gt.0.000001 .and. xc2.lt.0.999999)then
!*** if O fraction is constant, set Ledoux to zero. ***
          if(aa(20,k).eq.aa(20,k-1))then
             goto 16
          else
!*** otherwise, compute a Ledoux term. ***
             call ledouxc
             if(aa(20,k).le.0.0 .or. aa(20,k-1).le.0.0)then
                goto 16
             endif
             dlo = ( dlog(aa(20,k)) - dlog(aa(20,k-1)) )
             dlp = ( dlog(aa(6,k)) - dlog(aa(6,k-1)) )
             chtor = -aa(9,k)/aa(10,k)
             aa(19,k) = chtor*(dlo/dlp)*bledc

! Note: there is a minus sign difference between this ledoux term and
! the one in Brassard et al. (1991, ApJ, 367, 601).
! This was "corrected" by the absolute value below, but I want to be
! able to see if certain configurations are convectively unstable in the
! core, so I am going to remove the dabs() (MHM April 1998).
             aa(19,k) = -aa(19,k)
          endif
       endif
    endif
16  return
  end subroutine calc

!************************************************************************

  subroutine ledouxc

! compute modified Ledoux B term to include composition gradient
! effects in the Brunt-Vaisala frequency
! This is the version that computes the B term in the core.
! I vary both the C and O abundance at the C/O interface.
! PAB 1 October 1991

! Subroutines
    use istat_subroutines, only: istatco

! Common blocks
    use temp
    use thermo 
    use bldxc
    use flags, only: verbose

    implicit double precision(a-h,o-z)

    real*8, dimension(4) ::  x
    real*8, pointer :: fn



! DEL is step size increment in composition
    fn => fn2
    bledc = 0.0
    del  = 0.01
    x(1) = 0.0
    x(2) = 0.0
    x(3) = xc2
    x(4) = xo2

! C/O interface
! I am computing partial ln P/partial ln O here.
! now to consider specific composition mixtures!
! CASE 3: x(4) or x(3) is less than del (use a 3 point formula)
! CASE 3.1: Nearly pure oxygen: Treat first point as pure O.
!           treat second point as C/O mixture 

    if(x(3).lt.del)then
       diff = 1.0-x(4)
       xone = x(4) - diff
       xc2 = 1.0 - xone
       xo2 = xone
       call istatco(p2,t2,1,.false.)
       dens1 = 10.**(u2)
       xtwo = 1.0
       xc2 = 0.0
       xo2 = xtwo
       call istatco(p2,t2,1,.false.)
       dens2 = 10.**(u2)
       xc2 = x(3)
       xo2 = x(4)
       bledc=(xo2/(2.*diff))*(dlog(dens2)-dlog(dens1))
       goto 26

! CASE 3.2: Nearly pure Carbon: Treat first point as pure C.
!           treat second point as C/O mixture 

    else if(x(4).lt.del)then
       diff = x(4)
       xone = 0.0
       xc2 = 1.0 - xone
       xo2 = xone
       call istatco(p2,t2,1,.false.)
       dens1 = 10.**(u2)
       xtwo = x(4) + diff
       xc2 = 1.0 - xtwo
       xo2 = xtwo
       call istatco(p2,t2,1,.false.)
       dens2 = 10.**(u2)
       xc2 = x(3)
       xo2 = x(4)
       bledc=(xo2/(2.*diff))*(dlog(dens2)-dlog(dens1))
       goto 26
    endif

! CASE 2: We are not near 0 or 1 in composition
! we can use 3-point differencing formula with impunity.

    xone = x(4) - del
    xo2 = xone
    xc2 = 1.0 - xone
    call istatco(p2,t2,1,.false.)
    dens1 = 10.**(u2)
    xtwo = x(4) + del
    xo2 = xtwo
    xc2 = 1.0 - xtwo
    call istatco(p2,t2,1,.false.)
    dens2 = 10.**(u2)
    xc2 = x(3)
    xo2 = x(4)
    bledc0=(xo2/(2.*del))*(dlog(dens2)-dlog(dens1))
    call istatco(p2,t2,1,.false.)
    dens3 = 10.**(u2)
    bledc1=((xo2+del/2.)/(del))*(dlog(dens2)-dlog(dens3))
    bledc2=((xo2-del/2.)/(del))*(dlog(dens3)-dlog(dens1))

! the following lines are a kludge to correct for a discontinuity
! in the density-composition relationship across the crystallization
! boudary. dlog(rho)/dlnX is artificially high when X=xone and X=xtwo
! are on different sides of the crystallization boundary. We check for
! this case with the following line and compute on just one or the other
! side of the boundary (bledc1 *or* bledc2)

    if ( abs(bledc1-bledc2)/bledc0 .gt. 0.5 ) then
       if (bledc1 .lt. bledc2 ) then
          bledc=bledc1
       else
          bledc=bledc2
       endif
       if (verbose) write(*,*) bledc,u2,t2,abs(bledc1-bledc2)/bledc0
    else
       bledc=bledc0
    endif

26  return
  end subroutine ledouxc

!*************************************************************************

!  subroutine cshell

!  cshell steps from the center shell to the next shell.
!  gshell steps from shell n to shell n+1
!  See Schwarzschild & Harm 1965, ApJ, 142, 855
!  for a partial explanation of variable meanings

!  see subroutine gshell, as cshell = gshell + intro bit  

!*************************************************************************

  subroutine gshell(docshell)
!  formerly "entry gshell"

!  gshell steps from shell n to shell n+1
!  See Schwarzschild & Harm 1965, ApJ, 142, 855
!  for a partial explanation of variable meanings

    use misc_const
    use coef 
    use cal
    use temp
    use thermo 
    use contrl

    implicit double precision(a-h,o-z)

    real*8, pointer :: fn
    logical :: docshell

1244 Format(5F16.4)
    fn => fn2
 
    if ( docshell ) then
       hrp(1)=-(up2)*.333333
       hrt(1)=-(ut2)*.333333
       hra(1)=-.333333*(u2+3.*r2-s2+.622089-sm)
       hpp(1)=0
       hpt(1)=0
       hpa(1)=0
       hbp(1)=(-qe2*ep2+qnp2)*fnat
       hbt(1)=-qe2*(e2-ea2+et2*fnat)+qnt2*fnat
       hba(1)=-b2+(-qe2*(e2-ea2)+qn2)*fnat
       htp(1)=0
       htt(1)=0
       hta(1)=0
       it=-1
       return
    end if

    gs=2.*fnat/(s2-s1)
    ar=gs*(r1-r2)+fnat*(q1+q2)
    ap=gs*(p1-p2)-fnat*(f1+f2)
    gp=2.*fnat/(p1-p2)
    gt=gp*(t1-t2)/(p1-p2)
    at=gp*(t1-t2)-fnat*(w1+w2)
    ab=gs*(b2-b1)+fnat*(qe1*(e1-ea1)+qe2*(e2-ea2)-qn1-qn2)
    c1=q1*up1+hrp(k-1)*(-gs+3.*q1)
    c2=q1*ut1+hrt(k-1)*(-gs+3.*q1)
    c3=gs+3.*q2
    c4=q2*up2
    c5=q2*ut2
    car=ar+hra(k-1)*(gs-3.*q1)
    c6=-gs-f1-4.*f1*hrp(k-1)
    c7=-4.*f1*hrt(k-1)
    c8=-4.*f2
    c9=gs-f2
    cap=ap+4.*f1*hra(k-1)
    c10=gt+w1+fnat*w1*(op1/o1+hbp(k-1)/b1)
    c11=-gp-4.*w1+fnat*w1*(ot1/o1+hbt(k-1)/b1)
    c12=fnat*w2/b2
    c13=-gt+w2+w2*op2*fnat/o2
    c14=gp-4.*w2+w2*ot2*fnat/o2
    cat=at-fnat*w1*hba(k-1)/b1
    c15=gs*hbp(k-1)+fnat*(qnp1-qe1*ep1)
    c16=-qe1*(e1-ea1)+fnat*(qnt1-qe1*et1)+gs*hbt(k-1)
    c17=-gs
    c18=fnat*(qnp2-qe2*ep2) 
    c19=-qe2*(e2-ea2)+fnat*(qnt2-qe2*et2)
    cab=ab-gs*hba(k-1)
    c20=c1*c8-c3*c6
    c21=c2*c8-c3*c7
    c22=c4*c8-c3*c9
    c23=c5*c8
    carp=c8*car-c3*cap
    c24=c12*c15-c10*c17
    c25=c12*c16-c11*c17
    c26=c12*c18-c13*c17
    c27=c12*c19-c14*c17
    cabt=c12*cab-c17*cat
    c28=c20*c25-c21*c24
    c29=c22*c25-c21*c26
    c30=c23*c25-c21*c27
    cp=c25*carp-c21*cabt
    c32=c20*c26-c22*c24
    c33=c20*c27-c23*c24
    ct=c20*cabt-c24*carp
    cd1=-1./c28
    cd2=-1./c17
    cd3=-1./c3
    hpp(k)=c29*cd1
    hpt(k)=c30*cd1
    hpa(k)=-cp*cd1
    htp(k)=c32*cd1
    htt(k)=c33*cd1
    hta(k)=-ct*cd1
    hbp(k)=( c15*hpp(k)+c16*htp(k)+c18)*cd2
    hbt(k)=( c15*hpt(k)+c16*htt(k)+c19)*cd2
    hba(k)=( c15*hpa(k)+c16*hta(k)-cab)*cd2
    hrp(k)=( c1*hpp(k)+c2*htp(k)+c4)*cd3
    hrt(k)=( c1*hpt(k)+c2*htt(k)+c5)*cd3
    hra(k)=( c1*hpa(k)+c2*hta(k)-car)*cd3

    return
    
  end subroutine gshell

!************************************************************************

  subroutine write2

    use flags, only: evoloutput
    use rw2
    use shells 
    use thermo 
    use surf
    use temp
    use contrl
    use xxnew
    use times
    use comptmp
    use flagtmp
    use terp
    use vca
    use dvca
    use xbnd
    use corats
    use fin
    use dfin
    use tfitt
    use gne
    use crash
    use outfile
    use read2input

    implicit double precision(a-h,o-z)

!    real*8, dimension(5) :: bn
    real*8, dimension(600) :: cp, cv
    real*8, pointer, dimension(:) :: x1, y1, z1
    real*8, dimension(600,7) :: x, y, z
    real*8, pointer :: fn
   
    fn => fn2
    x1 => sa
    y1 => s
    z1 => sk
    
309 format(42h         **** convection in mass-shells j=i3, &
         7h, s(j)=, 1e16.8,12h, through j=i3,7h, s(j)=e16.8)
310 format(26h              containing  f11.5, &
         23h solar masses with x12=f8.6,'  x16=',f8.6) 
319 format(/' s2=',1pe16.8,'  r2=',0pf12.7,'  b2=', 1pe12.4, &
         '  p2=', 0pf11.6,'  t2=',f11.6,'  e2=',1pe14.7, &
         '  xc2=',0pf8.6,'  xo2=',f8.6)
321 format(' q2=', 1pe12.4, '  f2=', e12.4, '  w2=', e12.4,'  ea2=', &
         e12.4,'  psi=',e12.4,'  cp= ',e12.4,'  cv = ',e12.4) 
323 format(4h u2= e12.4,5h  o2= e12.4,4h  b=f11.6,5h  wc=f11.6, &
         5hvesc=e12.4,9h  eratio=e12.4)
325 format(5h fte=e12.4,5h fcc=e12.4,5h  fn=e12.4,23x,2hr=e12.6, &
         5h   m=f9.6)
327 format(5h epl=1pe12.4,5h eph=e12.4,5h epa=e12.4,5h erc=e12.4, &
         5h ebr=e12.4,6h etot=e12.4)
    
    x(:,1) = x1
    y(:,1) = y1
    z(:,1) = z1

!  called by end, write2 writes a shell-by-shell summary of 
!  the computed model

!  compute specific heats cp and cv, and store in 
!  cp() and cv().  write out both here and along with
!  the model (to tape9).

    cp(k) = et2
    cv(k) = et2 - ep2 * ut2 / up2
    xmass=10.**s2-10.**s1
    smr=10.**(sm-33.591065)
    xr=xmass*smr
    do i=1,5
       bn(i)=bn(i)+en2(i)*xr
    enddo
    bnt=bnt+fn*xr 
    bax=bax+fa*xr
    egrav=egrav-xmass*10.**(2.*sm-7.176004)*(10.**s2+10.**s1)/(10.**r2+10.**r1)
    
! compute radius, velocity, acceleration, and shock at location of max.
! expansion. (At least that's what it used to do.)
! Now it computes this for m/M* = 0.1 (fixed) 

    if ( speak .le. s2 .and. s1 .le. speak ) then
       rpeak=r1+(r2-r1)*(speak-s1)/(s2-s1)
       vpeak=(10.**rpeak-10.**rpeak1)*10.**(-gx) 
       apeak=2.*(vpeak-vpeak1)/(10.**gx+10.**g1) 
       shock=apeak*10.**(rpeak*2.-speak-sm+7.176004)
    endif

! mess with core convection zone, if one is present.

    if(xc2.ne.0.) xctmp=xc2/(1.+2.2892*fcc*10.**(g2-18.)/xc2) 
    if ( xc(k) .lt. 0.000001 ) xc(k)=0.000001
    if ( kon .gt. 0 ) then
       k1=k
       tmass=tmass+xmass
       total=total+xctmp*xmass
    elseif ( k1 .gt. 0 ) then
       do k3=kom,k1
          xctmp=total/tmass
       enddo
       xotmp = xo2
       tmass=tmass*10.**(sm-33.2986)
       if ( evoloutput ) then
          write(17, 309)kom,s(kom),k1,s(k1)
          write(17, 310) tmass,xctmp,xotmp  
          write(20,309) kom, s(kom), k1, s(k1)
          write(20,310) tmass, xctmp, xotmp
       end if
       k1=0
       kom=0
       total=0.
       tmass=0.
    endif

! carbon burning stuff

    if ( fcc .gt. fcc1 ) then
       fcc1=fcc
       xdel=xc(k)-xc2
       sc1=s2
    endif
 
    if ( fn .gt. fn1 ) then 
       fn1=fn
       sn1=s2
    endif
    if(it.le.0 .and. mod(modnr,ip1).ne.0 .and. mod(k,10).ne.1)then
       if ( ip40 .eq. 1 ) then
          goto 1980
       else
          goto 36
       endif
    endif
    
    pz=10.**(pg-p2)
    eratio=10.**(u2+s2+sm-p2-r2-7.352095)/(2.-pz)
    vesc=10.**(.5*(s2+sm-r2-6.874974))
    fte=(e2-ea2)*10.**(t2-gx)
    r3=10.**r2
    s3=10.**s2
    if ( evoloutput ) then
       write(17,319) s2,r2,b2,p2,t2,e2,xc2,xo2
       write(17,321) q2,f2,w2,ea2,psi2,cp(k),cv(k)
       write(17,323) u2,o2,pz,wc,vesc,eratio 
       write(17,325) fte,fcc,fn,r3,s3
       write(17,327) (en2(i),i=1,5),fn
    end if
    theta = 7832.3 * dsqrt( 10.**u2 ) / 10.**t2
    gammac = 10.**(u2/3. -t2) * 3.5772e6
    gammao = gammac / 3.5772e6 * 5.7782e6
    gamma = xc2 * gammac + xo2 * gammao
    if ( evoloutput ) then
       write(17, 332) theta, gamma
    end if
332 format(' theta= ',1pe12.5,' gamma= ',e12.5/)
    
! before exiting, write to tape40 if ip40=1

1980 if ( ip40 .eq. 0 ) then
       goto 36
    endif

    bl=dlog10(b(jb))
    tel=bl/4.-rm/2.+9.18458 

    if ( evoloutput ) then
       if ( k .eq. 1 ) then 
          write(40,*) '---------------------------------------------'
          write(40,1999) modnr,ssg,10.**rm,10.**tel,bl
          write(40,700) jb
       end if
       write(40,2001) 10.**s2,r2,t2,p2,e2,b2,fn,cp(k),cv(k),o2,wc,w2, &
            theta,psi2,u2
       if ( k .eq. jb ) then
          write(40,1999) modnr,ssg,10.**rm,10.**tel,bl,dlog10(bnt)
       endif
    end if
700 format(i10)
1999 format(i4,1pe13.5,1pe13.5,0p,f9.0,f9.4,f9.4)
2001 format(f9.7,3f9.5,1p,3e13.5,/,3x,3e11.3,0p,2f9.5,1p,2e11.3,0pf9.4)

    x1 = x(:,1)
    y1 = y(:,1)
    z1 = z(:,1)
    
36  return
    
  end subroutine write2

!************************************************************************

  subroutine interp
    
! interpolate model from old to new meshes

    use shells 
    use cal
    use thermo 
    use temp
    use terp
    use contrl
    use corats
    
    implicit double precision(a-h,o-z)

    real*8, pointer :: fn, g

    fn => fn2
    g => gx

    if (nite.ge.20) niteo = 30
    s1=s2
    r1=r2
    b1=b2
    p1=p2
    t1=t2
    up1=up2
    ut1=ut2
    o1=o2
    ot1=ot2
    op1=op2
    e1=e2
    ep1=ep2
    et1=et2
    ea1=ea2
    q1=q2
    f1=f2
    w1=w2
    qe1=qe2
    qn1=qn2
    qnp1=qnp2
    qnt1=qnt2

! now that we've saved all the variables for the previous shell,
! let's go about stepping to the next shell.

    fac = 1.
60  ds=c/max(fac,q1,f1)

! 0.008 was the value I used before 2/28/91
! Matt uses .005 here. Gets lots more shells.
! Signed either Travis or Mikemon

! Agnes Kim, 2017:
! The most obscure bug ever. Took me a week. If core appears chopped off
! before the stop mass, try reducing the number below (e.g. 0.008) to a
! slightly smaller value. If other funny business occurs at the C/He transition
! try increasing it. Careful with this, it's sensitive. Save last working 
! version before modifying.

    dsx=dlog10(10.**s1+.008)-s1

! Matt's new value as of 2/28/91 
! want -1.0 instead of -4.0 to  prevent undesirable jumps in core
    if(s1.gt.-1.0.and.s1.lt.-.02)then
       ds=dsx
    endif
    if( nite .le. 20 .and. nite .ge. 10 &
         .or.(nite.eq.0.and.niteo.le.20.and. niteo .ge. 10)) then
       ds = .90*ds
       niteo = nite
    elseif ( nite .lt. 10 .and. nite .ge. 2  &
         .or.(nite.eq.0.and.niteo.lt.10.and. niteo .ge. 2)) then
       ds = .80*ds
       niteo = nite
    endif

    s2=s1+ds

! make sure that corat boundaries are mesh points 

200 if ( s2 .ge. stpms ) then
       it=1
       s2=stpms
       ds=stpms-s1
    endif

    do
       if ( j .gt. jb ) exit
       if ( s2 .le. s(j) ) then
          sfrac=(s2-s(j))/(s(j)-s(j-1))
          p2=p(j)+(p(j)-p(j-1))*sfrac 
          b2=b(j)+(b(j)-b(j-1))*sfrac 
          r2=r(j)+(r(j)-r(j-1))*sfrac 
          t2=t(j)+(t(j)-t(j-1))*sfrac 
          goto 46
       endif
       j=j+1
    end do

    cs=dlog10(1.+ds/s1)
    p2=p1+cs
    t2=t1+w1*cs
    r2=r1+ds*q1
    b2=b1

46  return
  end subroutine interp

!************************************************************************

  subroutine end(line)

!  writes out the model according to the switches (ipn's) that have
!  been set. [PROFILE: approximately 10% runtime]

! Subroutines
    use utils_subroutines, only : check_time, arline, armove
    use eprep_subroutines, only : eprep

    use shells 
    use coef 
    use contrl
    use times
    use vca
    use dvca
    use deldgg
    use fin
    use tfitt
    use idone
    use surf
    use stuck
    use gne
    use crash
    use flags, only: verbose, evoloutput

    implicit double precision(a-h,o-z)
    
    real*8, dimension(8000) :: stash2
    real*8, dimension(600,8) :: x, y
    real*8, dimension(600,5) :: z
    integer :: newtime
    
603 format(' sk= ',1pe16.8,'  dr2=',0pf12.7,'  db2=',1pe12.4, &
         '  dp2= ',0pf11.6,'  dt2= ',f11.6)
701 format('  nite=',i12,'  jmax=',i12,'  smax=',1pe16.8, &
         '  dtmax=',e12.4)
704 format(' s=',1pe16.8,'  r=',0pf12.7,'  b=',1pe12.4,'  p=',0pf11.6, &
         '  t=',0pf11.6,'  e= ',1pe14.7,'  xc=',0pf8.6,'  xo=',f8.6)
718 format (//14h dtmax too big//)
720 format(' COARSE approximation, time= ',f8.3) 
727 format(' VERY COARSE approximation, time= ',f8.3)

    numtries = 0
    line = 3
    call stitch(hbp(k),hbt(k),hba(k),hrp(k),hrt(k),hra(k),dp2,dt2)
    call check_time(newtime)
    if ((newtime - t0) .gt. xmaxtime) goto 806
    if (ihotspot.eq.1) goto 806
    jb=k
    dtmax=0.0

    do j=1,jb 
       k=jb-j+1
       dr2=hrp(k)*dp2+hrt(k)*dt2+hra(k)
       db2=hbp(k)*dp2+hbt(k)*dt2+hba(k)
       if ( evoloutput ) then
          if ( ip5.gt.0 .or. (ip5.eq.0 .and. nite .le. 1)) then
             write(17, 603)sk(k),dr2,db2,dp2,dt2
          endif
       end if
       s(k)=sk(k) 
       r(k)=rk(k)+f*dr2
       b(k)=bk(k)+f*db2
       p(k)=pk(k)+f*dp2
       t(k)=tk(k)+f*dt2
       if(dabs(dtmax).lt.dabs(dt2)) then
          dtmax=dt2
          smax=sk(k)
       endif
       dp1=hpp(k)*dp2+hpt(k)*dt2+hpa(k)
       dt1=htp(k)*dp2+htt(k)*dt2+hta(k)
       dp2=dp1
       dt2=dt1
    enddo

    if ( evoloutput ) write(17,701) nite,jb,smax,dtmax
    nite=nite-1

! Tried to modernize the following. Hopefully execution branching is preserved.
! if ( ip6 ) 706,702,703
! 702 if ( nite ) 703,703,714 
! 703 continue
!
! means
! 
!    if ( ip6 < 0 ) goto 706
!    if ( ip6 = 0 ) goto 702 
!    if ( ip6 > 0 ) goto 703
!702 test nite:
!         if ( nite < 0 ) goto 703
!         if ( nite = 0 ) goto 703
!         if ( nite > 0 ) goto 714
!703 continue

    if ( ip6 < 0 ) goto 706
    
    if ( ip6 == 0 ) then
       if ( nite > 0 ) then
          goto 714
       end if
    end if

    if ( evoloutput ) then
       do j = 1, jb
          write(17,704) y(j,1), y(j,2), y(j,3), y(j,4), y(j,5), y(j,6), &
               y(j,7), y(j,8)
       end do
    end if

    if ( ip5 .gt. 0 .and. ip6 .gt. 0 ) then
       print *, "mysterious stop in subroutine end"
       stop
    end if
!  Stash the interior shells
!  Call eprep to calculate converged envelope, store it in aa,
!  and write it to tape50!

706 if(nite .le. 0 ) then

!***  Read the interior shells back in! ***
       nshint = jb
       ifin = 1
! Save common block variables to stash2. I think just sa, ra, ...
! Add Helium abundance (xxhe) to the stash
       stash2(   1: 600) = sa(:)
       stash2( 601:1200) = ra(:)
       stash2(1201:1800) = ba(:)
       stash2(1801:2400) = pa(:)
       stash2(2401:3000) = ta(:)
       stash2(3001:3600) = ea(:)
       stash2(3601:4200) = xca(:)
       stash2(4201:4800) = fca(:)
       stash2(4801:5400) = xxhe(:)
! Debugging
!!$       open(unit=1631,file='stash')
!!$       do i=1,8000
!!$          write(1631,*) stash2(i)
!!$       enddo
!!$       print *, 'reached debugging stop statement in subroutine end'
!!$       stop

       write(111) stash2
       rewind 111
!*** write out model to tape9 ***
       call write3
       if (ihotspot .eq. 1) goto 806

!*** write out complete model to tape 50 ***
       call eprep
       call check_time(newtime)
       if ((newtime - t0) .gt. xmaxtime) goto 806
       if (ihotspot .eq. 1) goto 806
       ifin = 0
       read(111) stash2
       rewind 111

! Retrieve variables from stash
       sa(:)  = stash2(   1: 600)
       ra(:)  = stash2( 601:1200)
       ba(:)  = stash2(1201:1800)
       pa(:)  = stash2(1801:2400)
       ta(:)  = stash2(2401:3000)
       ea(:)  = stash2(3001:3600)
       xca(:) = stash2(3601:4200)
       fca(:) = stash2(4201:4800)
       xxhe(:)= stash2(4801:5400)
! Save sa, ra, ... to x and s, r, ... to y.
       x(:,1) = sa(:)
       x(:,2) = ra(:)
       x(:,3) = ba(:)
       x(:,4) = pa(:)
       x(:,5) = ta(:)
       x(:,6) = ea(:)
       x(:,7) = xca(:)
       x(:,8) = xxhe(:)

       y(:,1) = s(:)
       y(:,2) = r(:)
       y(:,3) = b(:)
       y(:,4) = p(:)
       y(:,5) = t(:)
       y(:,6) = e(:)
       y(:,7) = xc(:)
       y(:,8) = xxhe(:)
      
       if ( ip8 .ge. 0 ) then
          pext=10.**(gx-g1)
          call arline(x,y,z,ja,jb,5,ip8,pext)

! Save z to common block
          sk(:) = z(:,1)
          rk(:) = z(:,2)
          bk(:) = z(:,3)
          pk(:) = z(:,4)
          tk(:) = z(:,5)
       endif
       call armove(x,y,jb,8)

! Save x to common block
       sa(:)   = x(:,1)
       ra(:)   = x(:,2)
       ba(:)   = x(:,3)
       pa(:)   = x(:,4)
       ta(:)   = x(:,5)
       ea(:)   = x(:,6)
       xca(:)  = x(:,7) 
       xxhe(:) = x(:,8)
      
       ja=jb
       if(ip8.gt.0) then
          call armove(y,z,jb,5)

! Save y to common block
          s(:)   = y(:,1)
          r(:)   = y(:,2)
          b(:)   = y(:,3)
          p(:)   = y(:,4)
          t(:)   = y(:,5)
          e(:)   = y(:,6)
          xc(:)  = y(:,7)
          xxhe(:)= y(:,8)
       end if

       line=1
       if (tfit(ifit).lt.1.d-03) then
!***************************************************************
!*** THIS IS WHERE IT HAS CONVERGED TO THE RIGHT TEMPERATURE ***
!***************************************************************
          itconverged = 1
          if (verbose) then
             write(*,*) 'Converged to the',ifit-1,' specified temperatures'
             write(*,*) 'Exiting'
          end if
!Debugging
!Changed a stop statement in old code to a return statement. Not sure why
!or whether that's the right thing to do.
          goto 806
       end if
       goto 806
    endif
714 if ( nite .gt. 1 ) then 
       if ( dabs(dtmax) .ge. tmax1 ) then
          line=2
          goto 806
       else
          dg=dmin1(dg1,dg+.5*time)
          time=dmin1(time1,time+.1)
!            nite = 1
!            line = 2
!*** put convergence criteria for specific temperatures here ***

          bl=dlog10(b(jb))
          tel=bl/4.-rm/2.+9.18458

          if ( (itfit.eq.0 .and. tel.gt.log10(tfit(ifit)) ) .or. &
               abs(10**tel-tfit(ifit)).lt.0.1 ) then
             if (verbose) write(*,*) 'Done: tel tfit ifit itfit',&
                  10**tel,tfit(ifit),ifit,itfit
             told2= told
             told = tel
             if ( abs(10**tel-tfit(ifit)).lt.0.1 ) itfit=1
             if (itfit.eq.1) then
                ifit=ifit+1
                dg = dg + 0.1
             endif
             itfit=0
             istuck=0
             nite = 1
             line = 2
             gnew=gx
             if (itfit.eq.0 .and. abs(10**tel-tfit(ifit)) .lt. 50.0 ) then
                if (verbose) then
                   write(*,*) 'lengthening timestep to include next t(i)'
                end if
                nite = nite1
                line = 2
                itfit=1
                told=told2
             endif
          else
             nite = nite1
             line = 2
             gcor=(10**told-tfit(ifit))/(10**told-10**tel)
             if (istuck.gt.25 .and. abs(gcor).gt.0.05) then
                fac=float(istuck-25)/100.
                gcorarg=1.+(1.0+fac)*(gcor - 1.)
                if ( gcorarg .gt.0 ) then
                   gcor=log10(gcorarg)
                else
                   ihotspot=2
                   goto 806
                endif
             else
                gcor=log10(1.+0.40*(gcor - 1.))
             endif
             gx = gx + gcor
707          format('Last iteration: tel tfit gcor ',2f9.2,x,f7.4)
             if (verbose) write(*,707) 10**tel,tfit(ifit),gcor
             itfit=1
             istuck=istuck+1
             goto 806
          endif
          goto 806
       endif
    else
       if ( dabs(dtmax) .lt. tmax2 ) then
          time=time*.8
          dg=dg-time
          if ( evoloutput ) then
             write(17,720) time
             write(20,720) time
          end if
          if (verbose) write(*,720) time
          modnr = modnr - 1 
          sg=sg-10.**(gx-7.499095)
          gx = gx + deldg
          nite = 1
          line = 1
          goto 806
       elseif ( dabs(dtmax) .lt. tmax3 ) then
          time=time*.8
          dg=dg-time
          if ( evoloutput ) then
             write(17,727) time
             write(20,727) time
          end if
          if ( verbose ) write(*,727) time
          modnr = modnr - 1 
          sg=sg-10.**(gx-7.499095)
          gx = gx + deldg
          nite = 1
          line = 1
          goto 806
       else
          write(*,*) ' dtmax too big'
          modnr = modnr - 1 
          sg=sg-10.**(g-7.499095)
          gx = gx + deldg
          ihotspot = 1
          nite = 1
          line = 1
          goto 806
       endif
    endif
806 return

  end subroutine end

!************************************************************************

  subroutine stitch(hbp,hbt,hba,hrp,hrt,hra,dp2,dt2)

!  stitch together the envelope and the interior. 

! Subroutines
    use eprep_subroutines, only: env
    use utils_subroutines, only: check_time

! Common blocks
    use misc_const
    use tablecc
    use tableoo
    use shells3
    use temp
    use surf
    use crash
!Debugging
    use shells2 

    implicit double precision(a-h,o-z)

    real*8 :: g1, g2, g3
    real*8, dimension(2,2) :: w1,ww2,w3,w4,h
    real*8, dimension(2) :: del(2)
    integer, dimension(3) :: kk(3)
    integer, pointer :: i, k, l, m
    integer :: newtime

    i => is
    k => ks
    l => ls
    m => ms   

! initialize variables
    bm=dlog10(b2)
    do i=1,3
       kk(i)=0
    end do
! start loop where we try to trap luminosity and radius 2 
! is the top of the loop.
2   g1=k*((ww(1,2)-ww(1,3))*(rm-ww(2,2))+(ww(2,3)-ww(2,2))*(bm-ww(1,2)))
    if ( g1 .lt. 0. ) then
       do i=1,2 
          ww(i,1)=ww(i,3)+ww(i,2)-ww(i,1) 
       end do
       k=-k
       ms=ms-3*ls
       ls=-ls
       kk(1)=1
       goto 2
    endif
    g2=k*((ww(1,3)-ww(1,1))*(rm-ww(2,3))+(ww(2,1)-ww(2,3))*(bm-ww(1,3)))
    if ( g2 .lt. 0. ) then
       do i=1,2 
          ww(i,2)=2*ww(i,1)-ww(i,2)
       end do
       k=-k
       kk(2)=1
       goto 2
    endif
    g3=real(k)*((ww(1,1)-ww(1,2))*(rm-ww(2,1))+(ww(2,2)-ww(2,1))*(bm-ww(1,1)))
    if ( g3 .lt. 0. ) then
       do i=1,2
          ww(i,3)=2*ww(i,1)-ww(i,3)
       end do
       k=-k
       kk(3)=1
       goto 2
    endif

! this is the bottom of the loop
! now that we've trapped the values, let's get on with it
    if ( kstart .eq. 1 ) then
       kk(1)=1
       kk(2)=1
       kk(3)=1
       kstart=0
    endif
    if ( kk(1)+kk(2)+kk(3) .gt. 0 ) then
       write(111) stash2
       rewind 111
       do i=1,3
          if (kk(i) .gt. 0) then
             call env
             call check_time(newtime)
             if ((newtime - t0) .gt. xmaxtime) goto 56
             if (ihotspot.eq.1) goto 56
          endif
       end do
       read(111) stash2
       rewind 111
    elseif ( npass .gt. 0 ) then
       goto 17
    endif

! goes through here just after starting up the sequence. 
    do i=1,2
       do j=1,2
          n=ms-ls*j 
          cd=ww(j,n)-ww(j,1)
          w1(i,j)=(u(i,n)-u(i,1))/cd
          ww2(i,j)=(v(i,n)-v(i,1))/cd 
       end do
    end do

    call solve(w3(1,1),w3(2,1),w1,-1.0d0,0.0d0)
    call solve(w3(1,2),w3(2,2),w1,0.0d0,-1.0d0)

    do i=1,2
       do j=1,2
          w4(i,j)=ww2(i,1)*w3(1,j)+ww2(i,2)*w3(2,j)
       end do
    end do
17  del(1)=p2-u(1,1)
    del(2)=t2-u(2,1)
    rm=ww(2,1)
    h1=v(1,1)
    h2=v(2,1)
    do i=1,2
       rm=rm+w3(2,i)*del(i)
       h1=h1+w4(1,i)*del(i) 
       h2=h2+w4(2,i)*del(i) 
    end do
    h1=hba+b2-10.**h1
    h2=hra+r2-h2
    h(1,1)=hbp-w4(1,1)*b2/fnat
    h(1,2)=hbt-w4(1,2)*b2/fnat
    h(2,1)=hrp-w4(2,1)
    h(2,2)=hrt-w4(2,2)
    call solve(dp2,dt2,h,h1,h2)
    npass=npass+1 
56  return
    
  end subroutine stitch

!************************************************************************

  subroutine solve(u1,u2,t,c1,c2)

    implicit double precision(a-h,o-z)
    real*8, dimension(2,2) :: t
 
    d=t(1,1)*t(2,2)-t(2,1)*t(1,2)
    u1=(t(1,2)*c2-t(2,2)*c1)/d
    u2=(t(2,1)*c1-t(1,1)*c2)/d

    return
  end subroutine solve

!************************************************************************

  subroutine write3
  
!Subroutines
    use eprep_subroutines, only: eprep
    use phase_func

!Common blocks
    use comptmp
    use contrl
    use corats
    use crash
    use dfin
    use dvca
    use fin
    use flags, only : verbose, evoloutput
    use flagtmp
    use gne
    use ixswchh
    use misc_const, only : amsun
    use outfile
    use read2input
    use rw2
    use shells 
    use surf
    use tape50
    use temp
    use terp
    use tfitt
    use thermo 
    use times
    use vca
    use xbnd
    use xxnew

    implicit double precision(a-h,o-z)

!    real*8, dimension(5) :: bn
    real*8, pointer, dimension(:) :: x1, y1, y2, y3, y4, y5, y6, y7, y8, z1
    real*8, dimension(600,8) :: x,y,z
    real*8, pointer :: fn

    fn => fn2
    x1 => sa
    y1 => s
    y2 => r
    y3 => b
    y4 => p
    y5 => t
    y6 => e
    y7 => xc
    y8 => xhe
    z1 => sk

8   format(4i10/2f10.6/(6d13.7))
20  format(e15.8,f9.6,e15.8,f9.6,2e15.8)
21  format(2f9.3,f10.6)
153 format(e15.8,f9.6,1pe11.4,0pf10.6,f9.6,e14.7,2f8.6)
160 format(f6.3,f10.6,e14.7,e10.3,f5.2,2e9.2,i4,f9.6)
162 format(1h ,i2,f6.2,e13.6,6f8.3,7x,30hnn,sm,sg,p2,t2,ucent,rm,bl,tel)
163 format(1h ,'amx=',f6.3,'  amxo=',f6.3)
326 format(1h ,6h bnpl=f8.4,6h bnpa=f8.4,6h bnph=f8.4,6h bnrc=f8.4, &
         6h bnbr=f8.4) 
711 format(/7h   fcc=e12.4,6h,  s2=f8.3 ,7h, xdel=f12.6, &
         26h  max.c12 rate, x12-change)
712 format(7h    fn=e12.4,6h,  s2=f8.3 ,23h      max.neutrino rate) 
713 format(7h    bl=f12.4,6h, bcc=f12.4,6h,  bn=f12.4,7h,bgrav=f12.4, &
         8h  egrav=e14.5/)
724 format(2h  f7.3,e14.5,3e14.4,23h  m,r,v,a,shock at peak)
725 format(9h         e14.5,2e14.4,35h                r,v,a at last shell/)

    x(:,1) = x1
    y(:,1) = y1
    y(:,2) = y2
    y(:,3) = y3
    y(:,4) = y4
    y(:,5) = y5
    y(:,6) = y6
    y(:,7) = y7
    y(:,8) = y8
    z(:,1) = z1

! write header and model to tape9.

    ssg=ssg+10.**(gnew-7.499095)
    sn1=10.**sn1
    sc1=10.**sc1
    bl=dlog10(b(jb))
    blfin=bl
    rfin=rm
    bcc=dlog10(bcc)
    do i=1,5
       bn(i)=dlog10(bn(i))
    enddo
    bnt=dlog10(bnt)
    bax=dlog10(bax)

!*** grav. potential luminosity ***
    bgrav=dlog10(dabs(egrav-egrav1))-gx-33.591065 

!*** old expansion quantities. velocity and acceleration ***
    vexp=(10.**r(jb)-10.**rlast)*10.**(-gx)
    aexp=2.*(vexp-vexp1)/(10.**gx+10.**g1)
    spek=10.**speak
    rpek=10.**rpeak
    rlast=r(jb)
    rlas=10.**rlast
    if ( evoloutput ) then
       write(17, 711) fcc1,sc1,xdel 
       write(17, 712) fn1,sn1
       write(17, 326) (bn(i),i=1,5)
       write(17, 713) bl,bcc,bnt,bgrav,egrav 
       write(17, 724) spek,rpek,vpeak,apeak,shock
       write(17, 725) rlas,vexp,aexp
! dump out neutrino rates to separate file
       write(51, 333) (bn(i),i=1,5),bnt,bgrav,bl
    end if
333 format(8(1x,f8.4))


!*** dump out neutrino rates to separate file ***
!*** rename variables to those of previous model. ***
    egrav1=egrav
    vexp1=vexp
    vpeak1=vpeak
    rpeak1=rpeak
    g1=gx
    gx=g2
    g3=.5*max(bl,bcc)
    p2=p(1)
    t2=t(1)
    tel=bl/4.-rm/2.+9.18458 

    if ( ixswch .eq. 1 ) then
       pfreeze = phasec(t2) 
    elseif ( ixswch .eq. 2 ) then
       pfreeze = phaseo(t2) 
    elseif ( ixswch .eq. 3 ) then
       pfreeze = xc(1) * phasec(t2) + (1.-xc(1)) * phaseo(t2)
    elseif ( ixswch .eq. 4 ) then
       pfreeze = dmin1( phasec(t2), phaseo(t2) )
    else
       if (verbose) write(*,*) 'read2: ixswch not set. exiting.'
       stop
    endif
    xtest = y(1,4) - pfreeze
    amxc = -10.
    amxo = -10.
    if ( xtest .ge. 0. ) then
       call xbndry
    endif
    gold=gx

!*************************************************
!*** THIS IS WHERE THE OUTPUT MODEL IS WRITTEN ***
!*************************************************
    if (evoloutput) then
       if ((ixswch.ne.4).and.(tfit(ifit).lt.1.d-3)) then
!*** renew output model file for each write ***
          open(50,file=fname,status='unknown')
          write(50,161) modnr,ssg,p2,t2,ucent,rm,tel,bl,bnt,bax,10.**amxc
          smx=10.**(sm-33.298635)
       endif

161 format(i4,1p,e12.4,0p,f7.3,2f6.3,f7.3,f9.6,2f11.4,2f9.3)
171 format(i4,1p,e12.4,0p,f7.3,2f6.3,f7.3,f9.6,2f11.4,f9.3)

       if ( ixswch .eq. 4 ) then
          write(11,161) modnr,ssg,p2,t2,ucent,rm,tel,bl,bnt,10.**amxc,10.**amxo 
          write(9,161) modnr,ssg,p2,t2,ucent,rm,tel,bl,bnt,10.**amxc,10.**amxo 
       else
          write(11,171) modnr,ssg,p2,t2,ucent,rm,tel,bl,bnt,10.**amxc
          write(9,171) modnr,ssg,p2,t2,ucent,rm,tel,bl,bnt,10.**amxc
          write(20,171) modnr,ssg,p2,t2,ucent,rm,tel,bl,bnt,10.**amxc
       end if
       write(9,160) dg,gx,ssg,stpms,f,tmax1,tmax2,modnr,g3
261    format(f4.1,2i5,1p,3e9.2,0p,f8.0,f9.3,f8.3)
       write(9,20) speak,rpeak1,vpeak1,rlast,vexp1,egrav1
       write(9,21) bgrav,w,g1
       write(9,8) ks,ls,ms,kstart,rm,bm,u,v,ww
       write(17, 162)nu,sm,ssg,p2,t2,ucent,rm,bl,tel
       write(17, 163)10.**amxc,10.**amxo
    end if

876 format('mod',i4,'  age:',1pe12.2,'  te:',0pf8.0, &
         '  l:',f6.2,'  tc:',f6.2,'  written')
877 format('mod',i4,'  age:',1pe12.2,'  te:',0pf8.0, &
         '  l:',f6.2,'  tc:',f6.2,'  not written')
700 format(i10)

    if ( ip7 .gt. 0 ) then
       if(mod(m,md).eq.0) then
          if ( verbose ) write(*,876) modnr,ssg,10.**tel,bm,t2 
          if ( evoloutput ) then
             write(26,876) modnr,ssg,10.**tel,bm,t2 
             write(20,876) modnr,ssg,10.**tel,bm,t2 
             write(9,700) jb
             do j=1,jb
                write(9,153) y(j,1), y(j,2), y(j,3), y(j,4), y(j,5), &
                     y(j,6), y(j,7), 1.-y(j,7)-y(j,8)
             end do
          end if
       else
          if ( verbose ) write(*,877) modnr,ssg,10.**tel,bm,t2 
          if ( evoloutput ) then
             write(17,877) modnr,ssg,10.**tel,bm,t2 
             write(20,877) modnr,ssg,10.**tel,bm,t2   
          end if
       endif
       m=m-1
    endif
    
    nmod=nmod-1
!*** note: before 8 23 91 this was le 0. ***
    if ( nmod.le.0 )then
       nshint = jb
       call eprep
       ihotspot = 1
       return
    endif

! Probably not necessary, but just to be safe.
    x1 = x(:,1)
    y1 = y(:,1)
    y2 = y(:,2)
    y3 = y(:,3)
    y4 = y(:,4)
    y5 = y(:,5)
    y6 = y(:,6)
    y7 = y(:,7)
    y8 = y(:,8)
    z1 = z(:,1)
  
  end subroutine write3

!************************************************************************

  subroutine xbndry

! find the xtal boundary and print it's mass fraction

!Subroutines
    use phase_func

!Common blocks
    use acc
    use aoo
    use shells 
    use temp
    use contrl
    use xbnd
    use ixswchh
    use kill
    use flags, only: verbose, evoloutput

    implicit double precision(a-h,o-z)
    
    logical :: oncec, onceo
    real*8, dimension(600,8) :: y

! Fill y(:,:) with common block variables
    y(:,1) = s(:)
    y(:,2) = r(:)
    y(:,3) = b(:)
    y(:,4) = p(:)
    y(:,5) = t(:)
    y(:,6) = e(:)
    y(:,7) = xc(:)
    y(:,8) = xhe(:)

!   locate xtal boundary between two shells   ....

    jx = jb+2
    jxc = jx
    jxo = jx
    oncec = .true.
    onceo = .true.
    do i=1,jb
       temper  = y(i,5)
       xc2 = xc(i)
       xhe2 = xhe(i)
       xo2 = 1.0d0 - xc2 - xhe2
       if ( ixswch .eq. 1 ) then
          pfreeze = phasec(temper)
       elseif ( ixswch .eq. 2 ) then
          pfreeze = phaseo(temper)
       elseif ( ixswch .eq. 3 ) then
          pfrzc = phasec(temper)
          pfrzo = phaseo(temper)
          pfreeze = xc2 * pfrzc + xo2 * pfrzo
       elseif ( ixswch .eq. 4 ) then
          pfrzc = phasec(temper)
          pfrzo = phaseo(temper)
          pfreeze = max( pfrzc, pfrzo )
       else
          if (verbose) print *, 'ixswch not set. exiting.'
          if (evoloutput) write(20,*) 'ixswch not set. exiting.'
          stop
       endif
       press = y(i,4)
876    format('xbndry: i, xc, xo, t(i), p, px: ',i3,5f9.3)
       if ( ixswch .ne. 4 ) then
          if (press .le. pfreeze) then
             jx=i-1
             goto 3
          endif
       else
          if (press .le. pfrzc .and. oncec) then
             jxc = i-1
             oncec = .false.
          endif
          if (press .le. pfrzo  .and. onceo) then
             jxo = i-1
             onceo = .false.
          endif
          if (press .le. pfrzc .and. press .le. pfrzo ) goto 3
       endif
    enddo
      
3   continue
    if ( ixswch .ne. 4 ) then
       if ( evoloutput ) write(17, 91)pfreeze,jx,gold
       jxold=jx
       kx = jx + 1
       temp1 = y(jx,5)
       temp2 = y(kx,5)
       pres1 = y(jx,4)
       pres2 = y(kx,4)
       if ( ixswch .eq. 1 ) then
          px1 = phasec(temp1)
          px2 = phasec(temp2)
       elseif ( ixswch .eq. 2 ) then
          px1 = phaseo(temp1)
          px2 = phaseo(temp2)
       elseif ( ixswch .eq. 3 ) then
          pxc = phasec(temp1)
          pxo = phaseo(temp1)
          px1 = xc(jx) * pxc + (1.-xc(jx)) * pxo
          pxc = phasec(temp2)
          pxo = phaseo(temp2)
          px2 = xc(kx) * pxc + (1.-xc(kx)) * pxo
       endif
       sfrac = (px1 - pres1) / (pres2 - px2 + px1 - pres1 )
       am1 = 10.**s(jx)
       am2 = 10.**s(kx)
       amxc = dlog10( am1 + sfrac*(am2 - am1) ) 
       if ( evoloutput ) then
          write(17, 991)jx,pres1,pres2,px1,px2,amxc
          write(17,987) 10.**amxc, jx, 10.**y(jx,1), kx, 10.**y(kx,1)
       end if
       if (verbose) then
          write(*,987) 10.**amxc, jx, 10.**y(jx,1), kx, 10.**y(kx,1)
       end if
987    format(' amx: ',f7.3,2('  shell',i3,': ',f7.3))
991    format('xbndry: jx, pres1, pres2, px1, px2, amxc:',i3,5f7.2) 
    else
       if ( jxc .ne. jb + 2 .and. jxc .ne. 0 ) then
          jxcold = jxc
          pxc=phasec(y(jxc,5))
          kx=jxc+1
          dmdp=(y(kx,1)-y(jxc,1))/(y(kx,4)-y(jxc,4))
          amxc= y(jxc,1)+dmdp*(pxc-y(jxc,4))
          if ( evoloutput ) then
             write(17,988) 10.**amxc, jxc, 10.**y(jxc,1), kx,10.**y(kx,1)
          end if
          if ( verbose ) then
             write(*,988) 10.**amxc, jxc, 10.**y(jxc,1), kx, 10.**y(kx,1)
988          format(' amxc: ',f7.3,2('  shell',i3,': ',f7.3))
          end if
       endif
       if ( jxo .ne. jb + 2 .and. jxo .ne. 0 ) then
          jxoold = jxo
          pxo=phaseo(y(jxo,5))
          kx=jxo+1
          dmdp=(y(kx,1)-y(jxo,1))/(y(kx,4)-y(jxo,4))
          amxo= y(jxo,1)+dmdp*(pxo-y(jxo,4))
          if ( evoloutput ) then
             write(17,989) 10.**amxo, jxo, 10.**y(jxo,1), kx,10.**y(kx,1)
          end if
          if ( verbose ) then
             write(*,989) 10.**amxo, jxo, 10.**y(jxo,1), kx, 10.**y(kx,1)
          end if
       endif
    endif
989 format(' amxo: ',f7.3,2('  shell',i3,': ',f7.3))
      
!   find mx/mstar for this model              
      
91  format(1h ,'px,jx,g  ',e18.8,i5,e18.8)
    if ( jx .eq. jb-1 ) then
       ikill = ikill+1
       if ( ikill .gt. nxmods ) killit = .true.
    endif

    return

  end subroutine xbndry

!****************************************************************************
!  this routine is called once.  after reading in the desired
!  alph() constants (diffusion exponents), compute the
!  integrals 1 thru 4 to calibrate qmix1 and qmix2 and the upper
!  and lower boundaries of both transition zones.  integration
!  is by the trapazoidal rule.  the integration subroutines
!  are by press et al.
!       qmix1, qmix2: mass where concentrations of two species are equal
!       zone(1,1):    upper h transition zone boundary
!       zone(1,2):    lower h transition zone boundary
!       zone(2,1):    upper he transition zone boundary
!       zone(2,2):    lower he transition zone boundary
!  nominal values derived using 'trace' approximation and assuming
!  diffusive equilibrium:
!       alph(1): 5
!       alph(1): -5/4
!       alph(1): 2
!       alph(1): -2/3
!  MAW 2/14/88
!****************************************************************************
! As far as I can tell, this routine is no longer useful in current version
! of WDEC (where profiles are computed in a new way). Consider removing 
! altogether. For now, quantities never used are set to zero (and still written
! out) ABK 6/24/2016

  subroutine setprof
  
!Subroutines
    use eprep_subroutines, only: func1, func2, func3, func4
    use utils_subroutines, only: qsimp
    
!Common blocks
    use alpha
    use jee
    use comp
    use flags, only: verbose, evoloutput

    implicit double precision(a-h,o-z)

    zero = 0.
    cent = 100.
    one  = 1.
   
    if (verbose) write(*,*) 'alpha:', alph
    
!    call qsimp(func1,zero,   one, aint1)
!    call qsimp(func2, one,  cent, aint2)
!    call qsimp(func3,zero,   one, aint3)
!    call qsimp(func4, one,  cent, aint4)

    aint1 = 0.
    aint2 = 0.
    aint3 = 0.
    aint4 = 0.

10  format('Integral of Func',i1.1,' = ',f10.5)
    if (verbose) then
       write(*,*)
       write(*,10) 1,aint1
       write(*,10) 2,aint2
       write(*,10) 3,aint3
       write(*,10) 4,aint4
    end if
    if ( evoloutput ) then
       write(17,*)
       write(17,10) 1,aint1
       write(17,10) 2,aint2
       write(17,10) 3,aint3
       write(17,10) 4,aint4
    end if

!    qmix1 = amhyhe / ( aint1 + aint2 )
!    qmix2 = amheca / ( aint3 + aint4 )
    qmix1 = 0.
    qmix2 = 0.

    if ( evoloutput ) then
       write(17,30) 1,qmix1
       write(17,30) 2,qmix2
    end if
30  format('Qmix',i1.1,' = ',e12.5)

! Pushing the pure H boundary outward to accomodate a shallower He->H
! transition diffusion profile.

!    zone(1,1) = 0.001 * qmix1
!    zone(1,2) = 100.   * qmix1
!    zone(2,1) =   0.01 * qmix2
!    zone(2,2) = 400.   * qmix2

    zone(1,1) = 0.
    zone(1,2) = 0.
    zone(2,1) = 0.
    zone(2,2) = 0.

    if ( evoloutput ) then
       write(17,20) 1,1,zone(1,1)
       write(17,20) 1,2,zone(1,2)
       write(17,20) 2,1,zone(2,1)
       write(17,20) 2,2,zone(2,2)
    end if
20  format('Zone',i1.1,3x,i1.1,' = ',e12.5)

    return
  end subroutine setprof

!****************************************************************************

end module evol_subroutines
