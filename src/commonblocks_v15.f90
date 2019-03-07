module outfile

  character(15) :: fname
  
end module outfile

!****************************************************************************
! COMMON/cperiods/permin,permax,nper,lmax,calc_per,ncalc

module cperiods

  real*8 :: permin, permax
  real*8, dimension(2,100) :: calc_per
  integer :: nper, lmax
  integer, dimension(2) :: ncalc

end module cperiods

!****************************************************************************

module startmod
  
  real*8 :: Teff, M_h, M_env, M_he, SM_star, alph1, alph2, logg_init
  character(20) :: start_file

end module startmod

!****************************************************************************

module ax

  real*8 :: maxion

end module ax
  
!****************************************************************************

module flagtmp
  
  integer :: neuts

end module flagtmp

!****************************************************************************
! COMMON/crash/itcrashed,ittimedout,irestart,ihotspot,istpms,mode,
!              xheneg,xmaxtime,t0
module crash
  
  integer :: itcrashed,ittimedout,irestart,ihotspot,mode,xheneg,t0
!  integer :: istpms ! No longer used
  integer, parameter :: xmaxtime = 120 !seconds

end module crash

!****************************************************************************

module fracxtal

  real*8 :: fracm, fracm_new
  
end module fracxtal

!****************************************************************************
! common/contrl/ds,g,sm,wc,it,nite,ja,jb,j,k,l

module contrl
  
  real*8 :: ds, sm, wc
  real*8, target :: gx
  integer :: it, nite, ja, jb, j, k, l

end module contrl

!****************************************************************************

module writeitt

  logical :: writeit
  
end module writeitt

!****************************************************************************

module oldval
  
  real*8 :: bold, teold
  
end module oldval

!****************************************************************************

module ifaill

  integer :: ifail

end module ifaill

!****************************************************************************

module threshh
  
  real*8 :: thresh
  
end module threshh

!****************************************************************************

module deldgg

  real*8 :: deldg

end module deldgg

!****************************************************************************

module ixswchh

  integer :: ixswch

end module ixswchh

!****************************************************************************

module idiffuss

  integer :: idiffus

end module idiffuss

!****************************************************************************

module mixl

  real*8 :: aml
  integer :: mls

end module mixl

!****************************************************************************

module wprep

  integer :: iprep

end module wprep

!****************************************************************************
! common/dprep/aa(20,500), ecv(600), ext(600), exr(600)
! aa(1,:)  = radius
! aa(2,:)  = mr
! aa(3,:)  = lr
! aa(4,:)  = temperature
! aa(5,:)  = density
! aa(6,:)  = pressure
! aa(7,:)  = neutrino emission rate
! aa(8,:)  = cv
! aa(9,:)  = chr 
! aa(10,:) = cht 
! aa(11,:) = epst
! aa(12,:) = epsr
! aa(13,:) = kapr
! aa(14,:) = kapt
! aa(15,:) = del
! aa(16,:) = delad
! aa(17,:) = xhe
! aa(18,:) = kap (opacity)
! aa(19,:) = ledoux term
! aa(20,:) = oxygen abundance

module dprep

  real*8, dimension(20,650) :: aa
  real*8, dimension(650) :: ecv, ext, exr

end module dprep

!****************************************************************************

module kill

  logical :: killit
  integer :: nxmods, ikill
  
end module kill

!****************************************************************************
! common/tfit/tfit(100),told,ifit,itfit,nite1

module tfitt

  real*8, dimension(100) :: tfit
  real*8 :: told
  integer :: ifit, itfit, nite1

end module tfitt

!****************************************************************************
! common/temp/s1,r1,s2,r2,b2,p2,t2,ea2,xc2,xo2,fca2,f2,q2,w2,c
!t2 = logT
!xc2 = Carbon abundance
!xo2 = Oxygen abundance
!p2 = logP
module temp

  real*8 :: s1,r1,s2,r2,b2,p2,xc2,xo2,xhe2,xh2,fca2,f2,q2,w2,c
  real*8, target :: t2, ea2

end module temp

!****************************************************************************
! common/comptmp/amhyhe_tmp,amheca_tmp,XO_tmp,XOfm_tmp,sm_tmp,
!                hecdexp_tmp,w1_tmp,w3_tmp,h1frac_tmp,h2frac_tmp
!sm_tmp = Stellar mass in solar masses

module comptmp

  real*8 :: amhyhe_tmp,amheca_tmp,amheca2_tmp,XO_tmp,XOfm_tmp,sm_tmp
  real*8 :: hecdexp_tmp,w1_tmp,w3_tmp,h1frac_tmp,h2frac_tmp

end module comptmp

!****************************************************************************

module idone

  integer :: itconverged

end module idone

!****************************************************************************
! common/shells/ sa(600),ra(600),ba(600),pa(600),ta(600),
!                ea(600),xca(600),fca(600),s(600),r(600),b(600),
!                p(600),t(600),e(600),xc(600),sk(600),
!                rk(600),bk(600),pk(600),tk(600)
! s = log(Mr), r = log(r) , b = l/Lsun, p = log(P) , t =log(T) , e = entropy
! (sa, ra, ...) and (sk, rk, ...) copies of (s, r, ...)

module shells 

  integer :: nshint ! Number of shells in the core
  real*8, dimension(600) :: ba,pa,ta,ea,xca,xxhe,fca,bk,pk,tk
  real*8, target, dimension(600) :: s,r,b,p,t,e,xc,xhe,sa,ra,sk,rk

end module shells 

!****************************************************************************

module fin

  integer :: ifin

end module fin

!****************************************************************************
! Also includes
! common/comp2/xhebar,amheca2
! Comes from makedb. Used to make double decked helium layers

module comp 

  real*8 :: amhyhe, amheca, amheca2, xhebar
  real*8, dimension(4) :: x
  
end module comp

!****************************************************************************
! common/thermo/u2,up2,ut2,e2,ep2,et2,psi2,pg,o2,op2,ot2,fp,ft,
!               en2(5),fn,fcc,fa,ce,ci,cif,w,nu

module thermo 

  real*8 :: pg,o2,fp,ft,ce,ci,cif,w
  real*8, dimension(5) :: en2
  real*8, target :: u2  !d2 in istat1
  real*8, target :: up2,ut2,e2,ep2,et2,psi2,op2,ot2,fcc,fa,fn2
  integer :: nu

end module thermo 

!****************************************************************************
!common/cal/ut1,ot1,et1,op1,o1,ep1,qnp1,qnt1,p1,q1,f1,t1,w1,b1,e1,
!           qe1,ea1,qn1,up1,qe2,qn2,qnp2,qnt2

module cal

  real*8 :: ut1,ot1,et1,op1,o1,ep1,qnp1,qnt1,p1,q1,f1,t1,w1,b1,e1
  real*8 :: qe1,ea1,qn1,up1,qe2,qn2,qnp2,qnt2

end module cal

!****************************************************************************
! common/xx/sin,sout,smid,grid1,grid2,ucent,cste,kon,kom

module xxnew
  
  real*8 :: sin,sout,smid,grid1,grid2,ucent,cste
  integer :: kon,kom

end module xxnew

!****************************************************************************

module jee

  integer :: je
  
end module jee

!****************************************************************************
! common/corats/ams(10),delmass(10),corat(10),dcorat(10),
!               ncore, irdold, alph 

module corats

real*8, dimension(10) :: ams, delmass, corat, dcorat
integer :: ncore, irdold, alph

end module corats

!****************************************************************************

module bldxc

  real*8 :: bledc

end module bldxc

!****************************************************************************
! common/c/az(3),z(3),taux,ts,ps,sm,rs,err,te,gs,rho,
!          ol1x,rl1x,el1x,atgx,rtgx,tl,pl,xmass,acu,am,bk,an0,pi
! bk,an0,pi now in module misc_const
! rho = logrho
! sm = stellar mass
! tl = logT
! pl = logP
! smm = stellar mass in grams

module cc 

  real*8, dimension(3) :: az, z
  real*8 :: taux,ts,ps,rs,err,te,gs,rho 
  real*8 :: ol1x,rl1x,el1x,atgx,rtgx,tl,pl,acu,am
  real*8 :: xmass ! Mr in solar masses
  real*8, target :: smm  ! Stellar mass in solar masses
end module cc

!************************* Diffusion profile info **************************
! common/alpha/alph(4),qmix1,qmix2,zone(2,2)

module alpha

real*8, dimension(2,2) :: zone
real*8, dimension(4) :: alph
real*8 :: qmix1,qmix2

end module alpha

!****************************************************************************

module flags

  integer :: idif, alphachoice
  integer, parameter :: chemprofmode = 1 
!1 If using points to define a profile
!2 If reading in a ready made profile (not implemented in wdec yet)
  logical :: ios_new, kap_new, firstmod
  logical :: verbose
  logical :: evoloutput, makeplots, kernels, pulsoutput, pulsemodel
  logical, parameter :: newledoux = .TRUE.

end module flags

!****************************************************************************
! common/coeff/cfac1,cfac2,cfac3,qmx,qmy,qs,aa,bb

module coeff

  real*8 :: cfac1, cfac2, cfac3, qmx, qmy, qs, aa, bb
  
end module coeff

!****************************************************************************
! COMMON /hypg/ aa,bb,cc,z0,dz

module hypg
  
  complex*16 :: aa, bb, cc, z0, dz

end module hypg

!****************************************************************************
! COMMON /path/ kmax,kount,dxsav,xp,yp
! Also throwing in some integers, which define the size of the arrays

module path

  integer :: kmax, kount
  integer, parameter :: MAXSTP=10000,NMAX=50,KMAXX=200
  real*8, dimension(KMAXX) ::  xp
  real*8, dimension(NMAX,KMAXX) :: yp
  real*8 :: dxsav

end module path

!****************************************************************************
! common/coef/hrp(600),hrt(600),hra(600),hbp(600),hbt(600),hba(600),
!             hpp(600),hpt(600),hpa(600),htp(600),htt(600),hta(600)

module coef

  real*8, dimension(600) :: hrp,hrt,hra,hbp,hbt,hba,hpp,hpt,hpa,htp,htt,hta

end module coef

!****************************************************************************
! common/times/tmax1,tmax2,tmax3,time,time1,f,dg,dg1,g1,ip5,ip6,ip8

module times

  real*8 :: tmax1,tmax2,tmax3,time,time1,f,dg,dg1,g1
  integer :: ip5,ip6,ip8

end module times

!****************************************************************************

module vca

  integer :: modnr

end module vca

!****************************************************************************

module dvca

  real*8, target :: sg

end module dvca

!****************************************************************************
! common/surf/u(2,3),v(2,3),ww(2,3),rm,bm,is,ks,ls,ms,kstart,npass

module surf

  real*8, dimension(2,3) :: u,v, ww
  real*8, target :: rm,bm
  integer, target :: is,ks,ls,ms
  integer :: kstart,npass

end module surf

!****************************************************************************
! common/stuck/dgcavg,istuck

module stuck

  real*8 :: dgcavg
  integer :: istuck

end module stuck

!****************************************************************************

module gne

  real*8 :: gnew

end module gne

!****************************************************************************
! common/shells/r(800),m(800),p(800),t(800),el1(800),ol1(800),
!               rl1(800),atg(800),rtg(800),phs(800)

module shells2

  real*8, dimension(800) :: r,mm,p,t,el1,ol1,rl1,atg,rtg,phs

end module shells2

!*************************** EOS data ***************************************
! common/tenv/rhok(3,60,38),eta(3,60,38),pp(3,60,38),datg(3,60,38),
!             od(3,18,38),tk(3,60),uu(3,60,38),xtt(3,60,38),xrt(3,60,38)

module tenv

  real*8, dimension(3,60,38) :: rhok,eta,pp,datg,uu,xtt,xrt
  real*8, dimension(3,18,38) :: od
  real*8, dimension(3,60) :: tk

end module tenv

!****************************************************************************

module xcompp

  real*8, dimension(600,4) :: xcomp, xcomp_orig
! The following composition arrays are good to use in the pulsation code
  real*8, dimension(600) :: xoxyg,xheli,xcarb,xhydr
  real*8, dimension(:), allocatable :: ams_o, corat_o, ams_he, corat_he
! amr_hyhe is the input boundary between helium and hydrogen regions in
! fractional stellar mass (e.g. 0.99 if log(M_He) = -2) 
  real*8 :: amr_hyhe, ao, ahe, buffer_inner, buffer_outer
  integer :: ndimo, ndimhe

end module xcompp

!****************************************************************************

module slr

  real*8, dimension(11) :: sl

end module slr

!****************************************************************************
! common/d/ ml(3,60),nmax,kind,lmax,nsuff,nlim,nlim1,j,iter,lmam, 
!           nlum,ko,nhp
module d 

  integer, dimension(3,60) :: ml
  integer :: nmax,kind,lmax,nsuff,nlim,nlim1,iter,lmam, nlum,ko,nhp
  integer, target :: jj

end module d

!****************************************************************************
! common/grad/ drdtp(600),cp(600),ttg(600),drdtpx,cpx,pecx,ttgx,
!              deptx,dpdr,dtdr,dmdr,dr,xr,xt,cv,drdptx

module grad

  real*8, dimension(600) :: drdtp,cp,ttg
  real*8 :: drdtpx,cpx,pecx,ttgx,deptx,dpdr,dtdr,dmdr,dr,xr,xt,cv,drdptx

end module grad

!****************************************************************************

module terp

  real*8 :: stpms_orig,stpms,rat1

end module terp

!****************************************************************************

module wprop

  real*8, dimension(600) :: bn2,ac2

end module wprop

!****************************************************************************

module dfin

  real*8 :: rfin,blfin

end module dfin

!****************************************************************************

module bldx

  real*8 :: bled

end module bldx

!****************************************************************************

module kapder

  real*8, dimension(800) :: rkapr, rkapt

end module kapder

!****************************************************************************

module opacfin

  integer :: ifinal

end module opacfin

!****************************************************************************
! common/neut/en(5),fn,fa

module neut

  real*8, dimension(5) :: en
  real*8 :: fn, fax

end module neut

!****************************************************************************

module rhomooee

  real*8 :: rhomooe

end module rhomooee

!****************************************************************************
! common/phys/flam,gamma,degtemp,tfermi

module phys

  real*8 :: flam, gamma, degtemp, tfermi

end module phys

!****************************************************************************

module tablecc

  real*8, dimension(6,39,87) :: tablec
  
end module tablecc

!****************************************************************************

module tableoo
  
  real*8, dimension(6,36,87) :: tableo

end module tableoo

!****************************************************************************
! common/sizec/tminc,tmaxc,deltc,pminc,pmaxc,delpc,pmc(50),ipc(50),
!              iphasec(50),ioverc,itc

module sizec

  real*8 :: tminc,tmaxc,deltc,pminc,pmaxc,delpc
  real*8, dimension(50) :: pmc
  integer, dimension(50) :: ipc,iphasec
  integer :: ioverc,itc

end module sizec

! common/sizeo/tmino,tmaxo,delto,pmino,pmaxo,delpo,pmo(50),ipo(50),
!              iphaseo(50),iovero,ito

module sizeo
  
  real*8 :: tmino,tmaxo,delto,pmino,pmaxo,delpo
  real*8, dimension(50) :: pmo
  integer, dimension(50) :: ipo,iphaseo
  integer :: iovero,ito

end module sizeo

!****************************************************************************

module iiii

  integer :: iii
  
end module iiii

!****************************************************************************
! common/thermoc/dc,dpc,dtc,ec,epc,etc,psic

module thermoc

  real*8 :: dpc,dtc,ec,epc,etc,psic
  real*8, target :: dc

end module thermoc

!****************************************************************************
! common/thermoo/do,dpo,dto,eo,epo,eto,psio

module thermoo
  
  real*8 :: dpo,dto,eo,epo,eto,psio
  real*8, target :: do

end module thermoo

!****************************************************************************

module opcswch

  logical :: first

end module opcswch

!****************************************************************************
! common/names/t1,t6,t9,t50

module names

  character(30) :: t1,t6,t9,t50

end module names

!****************************************************************************
! common/opacs/cso(3,29,8),csr(29),cst(29)

module opacs

  real*8, dimension(3,29,8) :: cso
  real*8, dimension(29) :: csr, cst

end module opacs

!****************************************************************************
! common/dirac/ psi(141),f12(141),f32(141)

module dirac

  real*8, dimension(141) :: psi, f12, f32

end module dirac

!****************************************************************************
! melting coefficients for carbon.

module acc

  real*8, dimension(6) :: ac

  data ac/1.24605188e3,-9.32465995e2,2.776543312e2,-4.080201741e1, &
       2.979129728e0,-8.665569421e-2/

end module acc

!****************************************************************************
!  melting coefficients for oxygen from a fit to the output of
!  xtal for temps in the range 5.4 < log(t) < 8.0.

module aoo
  
  real*8, dimension(6) :: ao

  data ao/-2513.772484, 1823.347553, -528.588773, 76.834498, &
       -5.580515, 0.161850/

end module aoo

!****************************************************************************
! common/xbnd/amxc,amxo,gold

module xbnd

  real*8 :: amxc, amxo, gold

end module xbnd

!****************************************************************************
! New modules. Created to share variables among subroutines read2, write1, 
! write2, and write3 (formerly entry write1, entry write2 and entry write3
! in subroutine read2)

module rw

  real*8, dimension(600) :: cp, cv
  real*8, dimension(600,7) :: x,y

end module rw

module rw2
  
  real*8, dimension(5) :: bn
  real*8 :: g2, g3, bnt, rpeak, vpeak, bgrav, egrav, xh, yh, yyh
  real*8 :: fcc1, sc1, fn1, sn1, bcc, apeak, shock, bax

end module rw2

!****************************************************************************

module shells3

  real*8, dimension(8000) :: stash2

end module shells3

!****************************************************************************
! common /cdatl/ a(4), b(3), c(3), d(4), e(3), f(3)
! Eliminated. Replaced by cndlc_data in block_data.f90

!****************************************************************************
! common/cdats/ a(5), b(3), e(5), f(3), i(5), j(3), p(5), q(3),
!               alpha(4), beta(4), c, d, g, h, k, l, r, s
! Eliminated. Replaced by cndsc_data in block_data.f90

!****************************************************************************
! common /odatl/ a(4), b(3), c(3), d(4), e(3), f(3)
! Eliminated. Replaced by cndlo_data in block_data.f90

!****************************************************************************
! common/odats/ a(5), b(3), e(5), f(3), i(5), j(3), p(5), q(3),
!               alpha(4), beta(4), c, d, g, h, k, l, r, s
! Eliminated. Replaced by cndso_data in block_data.f90

!****************************************************************************
! common /hedatl/ a(4), b(3), c(3), d(4), e(3), f(3)
! Eliminated. Replaced by cndlhe_data in block_data.f90

!****************************************************************************
! common/hedats/ a(5), b(3), e(5), f(3), i(5), j(3), p(5), q(3),
!                alpha(4), beta(4), c, d, ggg, h, k, l, r, s
! Eliminated. Replaced by cndshe_data in block_data.f90

!****************************************************************************
! common /hydatl/ a(4), b(3), c(3), d(4), e(3), f(3)
! Eliminated. Replaced by cndlhy_data in block_data.f90

!****************************************************************************
! common/freq/partt(650),piece(650),acous(650),bvfreq(650),
!             bvfrq(650),tfreq(650)

module freq

  real*8, dimension(650) :: partt, piece, acous, bvfreq, bvfrq, tfreq

end module freq

!****************************************************************************

module element

  real*8, dimension(650) :: xhe

end module element

!****************************************************************************
! COMMON/modelp1/gnu0,per0,tdyn,mstar,model

module modelp1

  real*8 :: gnu0, per0, tdyn, mstar
  integer :: model

end module modelp1

!****************************************************************************
! COMMON/modelp2/age,llsun,rrsun,teff,np_tmp

module modelp2
  
  real*8 ::  age,llsun,rrsun,teff
  integer :: np_tmp

end module modelp2

!****************************************************************************
! COMMON/tape28a/brad,amass,nmod,nnsurf

module tape28a

  real*8 :: brad, amass
  integer :: nmod, nnsurf

end module tape28a

!****************************************************************************
! COMMON/tape28b/xi2,rn,ggrav,rrho,mr2

module tape28b

  real*8, dimension(650) :: xi2, rn, ggrav, rrho, mr2

end module tape28b

!****************************************************************************
! COMMON/tape29/y1,voga1,y3,u,y5

module tape29

  real*8, dimension(650) :: yy1, voga1, yy3, u, yy5

end module tape29

!****************************************************************************
! contains a copy of the quantities written out to evolved model (tape50)

module tape50

  real*8 :: tel,bl
  integer :: nshell, nl, nel

end module tape50

!****************************************************************************
! common g(650),x(650),rho(650),yliq(4,650),ray(4,650)
! no name common block in original code

module puls

  real*8, dimension(650) :: gg,x,rho
  real*8, dimension(4,650) :: yliq,rayy

  
end module puls

!****************************************************************************
! common/misc/l,lhat,lindex,nsurf

module misc

  real*8 :: l, lhat, lindex
  integer :: nsurf

end module misc

!****************************************************************************
! common/dmisc/period,grav,pi,pi4,p43,eps,verg,eig,eigt,y3i,y3t,amass
! Taking out pi. Leaving pi4 = 4pi and p43 = 4pi/3. pi4 and p43 calculated in 
! subroutine init.
! Any routine that required common block dmisc also requires 
! module misc_const (for pi).

module dmisc

  real*8 :: period,grav,pi4,p43,eps,verg,eig,eigt,y3i,y3t,amass

end module dmisc

!****************************************************************************

module ray
  
  integer :: iray

end module ray

!****************************************************************************

module rs

  real*8, dimension(650) :: r

end module rs

!****************************************************************************

module ekint
  
  real*8 :: ekin

end module ekint

!****************************************************************************

module setup

  integer :: isetup

end module setup

!****************************************************************************
! common/modect/y1m(1000),y2m(1000),nfine,nodes1,nodes2,modep

module modect

  real*8, dimension(1000) :: y1m, y2m
  integer :: nfine, nodes1, nodes2, modep

end module modect

!****************************************************************************
! common/eig1/rint(4,650),h(650),f(650),part(650)

module eig1

  real*8, dimension(4,650) :: rint
  real*8, dimension(650) :: h,f,part

end module eig1

!****************************************************************************

module ms
  
  real*8, dimension(650) :: mr

end module ms

!****************************************************************************
! common/rot1/rone(650),rpone(650),rptwo(650),t(650),clk,crone,crtwo

module rot1

  real*8, dimension(650) :: rone, rpone, rptwo, t
  real*8 :: clk,crone,crtwo

end module rot1

!****************************************************************************
! common/rot2/angfac(650),rpthr(650),rpfour(650),rpfive(650)

module rot2
  
  real*8, dimension(650) :: angfac, rpthr, rpfour, rpfive

end module rot2

!****************************************************************************
! common/perds/periods(2000),discr(2000),pguess(3,101),yguess(101)

module perds

  real*8, dimension(2000) :: periods, discr
  real*8, dimension(3,101) :: pguess
  real*8, dimension(101) :: yguess

end module perds

!****************************************************************************

module mfac

  integer :: m

end module mfac

!****************************************************************************
! common/yfacs/y1(650),y2(650),y3(650),y4(650),y5(650)

module yfacs

  real*8, dimension(650) :: y1,y2,y3,y4,y5

end module yfacs

!****************************************************************************

module splinq

  real*8, dimension(6) :: vq

end module splinq

!****************************************************************************

module savingl

  integer :: k, kq

end module savingl

!****************************************************************************

module read2input

  real*8 :: ssg, speak,rpeak1,vpeak1,rlast,vexp1,egrav1
  integer :: ip7,nmod,m,md,ip1,ip40

end module read2input

!****************************************************************************
! common/heat/cp(650),ga1(650),tthl(650),b4(650),c4(650)

module heat
  
  real*8, dimension(650) :: cp, ga1, tthl, b4, c4

end module heat

!****************************************************************************
! common/values/xx(21,650)
!!$       r(n)     = xx( 1,n) 
!!$       mr(n)    = xx( 2,n)
!!$       lr(n)    = xx( 3,n)
!!$       t(n)     = xx( 4,n)
!!$       rrho(n)   = xx( 5,n)
!!$       p(n)     = xx( 6,n)
!!$       eps(n)   = xx( 7,n)
!!$       kap(n)   = xx( 8,n)
!!$       cv(n)    = xx( 9,n)
!!$       chr(n)   = xx(10,n)
!!$       cht(n)   = xx(11,n)
!!$       epsr(n)  = xx(12,n)
!!$       epst(n)  = xx(13,n)
!!$       kapr(n)  = xx(14,n)
!!$       kapt(n)  = xx(15,n)
!!$       del(n)   = xx(16,n)
!!$       delad(n) = xx(17,n)
!!$       xhe(n)   = xx(18,n)
!!$       derro(n) = xx(19,n)
! xx(20,n) is the ledoux term
! xx(21,n) is the oxygen abundance

module values 

  real*8, dimension(21,650) :: xx
  real*8, dimension(650) :: bled_new
  
end module values

!****************************************************************************
! common/stuff/mr(650)

module stuff

  real*8, dimension(650) :: mr

end module stuff

!****************************************************************************
! common/lnrp/xi(650)

module lnrp

  real*8, dimension(650) :: xi

end module lnrp

!****************************************************************************
! common/bvorig/bvfreq3(650)

module bvorig

  real*8, dimension(650) :: bvfreq3

end module bvorig
!****************************************************************************
! common/res/res(650)
! added to store intermediate integral values in dosim (MHM)

module resvalue
  
  real*8, dimension(650) :: res
  
end module resvalue

!****************************************************************************
! common/trap/per0,perH,perHe,alpha2,ane,l
! Not creating a module for this one, it's only used by one subroutine now.
!****************************************************************************
! common/lager/trapez

module lager

  real*8 :: trapez

end module lager

!****************************************************************************
! common/delta/ beck,del

module delta

  real*8 :: beck, del

end module delta
!****************************************************************************
! common/gkern/gkernel(650),gkernel2(650),phi(650),
!              dphix(650),dxphi(650),partrev(650)

module gkern
  
  real*8, dimension(650) :: gkernel, gkernel2, phi, dphix, dxphi, partrev

end module gkern
!****************************************************************************
! common/if/ifirst

module if
  
  integer :: ifirst

end module if



!****************************************************************************

