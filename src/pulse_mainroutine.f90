  subroutine pulsate

!***********************************************************************
! This code is taken from PRPWDXDP.F the prep code for models produced *
! by the evolution code WDXDFIT.F last modified by Mike Montgomery. It *
! calls the subroutine pulse at the end which is taken from CJHANRO.F  *
! a pulsation code originally written by Carl Hansen.                  *
! Specialized for Metacomputer use by Travis Metcalfe, 1998/1999       *
!***********************************************************************
        
! Subroutines
    use pulse_subroutines
    use utils_subroutines
    use wd_eos_mod, only : wd_eos

! Common blocks
    use alpha, only : alph
    use bvorig
    use cc, only : te, smm, gs !temporary
    use comp, only : amhyhe, amheca
    use dprep, only : aa
    use element
    use flags
    use fracxtal
    use freq
    use heat
    use lnrp
    use mixl ! temporary
    use modelp1
    use misc_const
    use modelp2
    use outfile
    use read2input, only : ssg
    use rw2, only : bnt, bax
    use shells, only : nshint
    use startmod, only : logg_init
    use stuff
    use surf, only : rm
    use tape28a
    use tape28b
    use tape29
    use tape50
    use temp, only : p2, t2
    use terp
    use values
    use vca
    use xbnd, only : amxc 
    use xcompp
    use xxnew, only : ucent

    implicit double precision(a-h,o-z)
    
    real*8, dimension(650) :: r,lr,t,p, kap,cv,chr,cht,epsr,epst,kapt,del, &
         delad,eps,kapr,v,ra,sound,derro,chtor,ram,mu,alambda,akrdsig2, &
         dnum,dnom,dnom2,derdel,derdels,alfa,eta,b1,b2,b3,c1,w,wk,bvfreqmag, &
         cs,parr,fac1,fac2
    real*8, dimension(3) :: tab
    real*8, dimension(7) :: xavec
    real*8, dimension(11) :: res
    real*8 :: lnlsun,x,dens,y2,y4,glog
    real*8 :: m1, m2
    real*8 :: yp1,ypn,interp_value
    integer, dimension(2) :: iop
    integer, dimension(3) :: itab
    integer :: i,j,ndim,nsub,imin,imax
    integer, parameter :: nsm=5
    real*8, dimension(2*nsm) :: xsub, ysub, yder

!Open files

    call openfi
    
    model = modnr
    age = ssg
    pcen = p2
    tcen = t2
    ucen = ucent
    rstr = rm
    teff = tel
    llsun = bl
    lnlsun = bnt
    xtal = 10.**amxc
    
    rrsun=10.**rstr/6.96e+10
    llsun=10.**llsun
    teff=10.**teff
    
    np = nshell
    np_tmp = np

1000 format (i5)

! read in equilibrium quantities

    do i=1,20
       in=i
       if(i.eq.8) in=9
       if(i.eq.9) in=10
       if(i.eq.10) in=11
       if(i.eq.11) in=13
       if(i.eq.13) in=14
       if(i.eq.14) in=15
       if(i.eq.15) in=16
       if(i.eq.16) in=17
       if(i.eq.17) in=18
       if(i.eq.18) in=8
       if(i.eq.19) in=20
       if(i.eq.20) in=21
       xx(in,:) = aa(i,:)
    end do    

    if ( makeplots ) then
       mstar=xx(2,np)/amsun
       glog = 27386.8*mstar/rrsun
       glog = glog/rrsun
       glog = log10(glog)
    end if

! Take care of kink in quantities at the core-envelope
! boundary. The idea is to take a few shells before the boundary and a
! few shells after and right at the boundary (top shell of core), use
! an interpolated value of the quantity instead of the actual value.

    xx(2,np-1)=xx(2,np-2)
    xx(2,np)=xx(2,np-1)

    nsub = 2*nsm
    imin = nshint-nsm
    imax = nshint+1+nsm
    
    do i=1,nsm
       xsub(i) = real(nshint-1-nsm+i)
    end do
    do i=nsm+1,nsub
       xsub(i) = real(nshint-nsm+i) 
    end do

10 format (10F12.4)

    do j=1,21
       ysub(1:nsm) = xx(j,imin:nshint-1)
       ysub(nsm+1:nsub) = xx(j,nshint+1:imax)
       yp1 = (ysub(2)-ysub(1))/(xsub(2)-xsub(1))
       ypn = (ysub(nsub)-ysub(nsub-1))/(xsub(nsub)-xsub(nsub-1))
       call spline(xsub,ysub,nsub,yp1,ypn,yder)
       call splint(xsub,ysub,yder,nsub,real(nshint),interp_value)
       xx(j,nshint) = interp_value
    end do
       
1001  format(4e22.15)
    if ( makeplots ) then
       call dert(rrho,r,derro,np)
       write(18,1000) np 
       do i=1,19
          write (18,1001) (xx(i,1:np))
       end do
       write (18,1001) (xx(21,1:np))
    end if
      
! take 19 slices through 2d array

    do n=1,np

       r(n)     = xx( 1,n)
       mr(n)    = xx( 2,n)
       lr(n)    = xx( 3,n)
       t(n)     = xx( 4,n)
       rrho(n)   = xx( 5,n)
       p(n)     = xx( 6,n)
       eps(n)   = xx( 7,n)
       kap(n)   = xx( 8,n)
       cv(n)    = xx( 9,n)
       chr(n)   = xx(10,n)
       cht(n)   = xx(11,n)
       epsr(n)  = xx(12,n)
       epst(n)  = xx(13,n)
       kapr(n)  = xx(14,n)
       kapt(n)  = xx(15,n)
       del(n)   = xx(16,n)
       delad(n) = xx(17,n)
!       xhe(n)   = xx(18,n)
       xhe(n) = xheli(n)
       derro(n) = xx(19,n)
! xx(20,n) is the ledoux term
! xx(21,n) is the oxygen abundance
    end do

2000  format (51x,'model data for texas 60400 ml3 T=',f7.0,'K',/)
2001  format (64x,i3,'points')
2002  format(3x,'n',8x,'r',15x,'Mr',12x,'Lr',12x,'T',13x,'Rho',11x, &
           'P',13x,'Eps',11x,'Kap',11x,'Cv'/)
2003  format(i4,9(1pe14.6))
2004  format(3x,'n',7x,'Chr',11x,'Cht',11x,'Epsr',10x,'Epst',10x, &
           'Kapr',10x,'Kapt',10x,'Del',11x,'Delad',9x,'XHe'/)

    if ( makeplots ) then
       write(8,2000) teff
       write(8,2001) np 
       write(8,2002)
       do n=1,np
          write(8,2003) n,r(n),mr(n),lr(n),t(n),rrho(n),p(n),eps(n),kap(n),cv(n)
       end do
       write(8,2004)
       do n=1,np
          write(8,2003) n,chr(n),cht(n),epsr(n),epst(n),kapr(n),kapt(n), &
               del(n),delad(n),xhe(n)
       end do

! compute logarithmic derivatives
       kount = np/8
       do kj=1,kount
          krount = kj*8
          dnum(kj) = dlog10(delad(krount))
          dnom(kj) = dlog10(r(krount))
          dnom2(kj) = dlog10(t(krount))
       end do
       intt = 1
       iop(1) = 4
       iop(2) = 4
       call coeff1(kount,dnom,dnum,w,iop,intt,wk)
       itab(1) = 0
       itab(2) = 1
       itab(3) = 0
       do n=1,np-1
          ar1 = dlog10(r(n))
          call terp1(kount,dnom,dnum,w,ar1,intt,tab,itab)
          derdels(n) = tab(2)
       end do
! now interpolate the dln(delad)/dlnr points back on the model grid
       do ik=1,np-1
          diff = derdels(ik+1) - derdels(ik)
          derdel(ik) = derdels(ik) + diff/2.0
       end do
       call dert(lr,r,eta,np)

    end if

! start computing pulsation quantities
    pi4=4.d0*pi 
    np=np-1
    totr=r(np) 
    totm=mr(np)
    gorsr=sqrt(g*totm/totr**3)
    piece(2)=pi4*(r(2)**3*rrho(2) - r(1)**3*rrho(1))/3.
    call profgen(np)

    if (newledoux) then

! Beginning of "New Ledoux" calculation.
! Loop to calculate ledoux term using Montgomery's prescription,
! eqn. (7) of Paxton et al. (2013) 2013ApJS..208....4P (MHM June 2016)
   
       bled_new(:)=0.d0
       parr(:)=0.d0
       
       do n=2,np-1
! calculate pressure at nth shell -- should equal p(n)
          xavec(:)   = 0.d0
          xavec(1) = xhydr(n)
          xavec(2) = xheli(n)
          xavec(5) = xoxyg(n)
          xavec(3) = xcarb(n)
          call wd_eos(rrho(n),t(n),xavec,'d',res)
          pressn = res(1)
          chitn  = res(6)
          parr(n) = pressn

! calculate pressure at rrho(n), t(n), and xcomp(n+1,:)
       xavec(:)   = 0.d0
       xavec(1) = xhydr(n+1)
       xavec(2) = xheli(n+1)
       xavec(5) = xoxyg(n+1)
       xavec(3) = xcarb(n+1)
       call wd_eos(rrho(n),t(n),xavec,'d',res)
       pressp = res(1)

! calculate pressure at n+1st shell -- should equal p(n+1)
          xavec(:)   = 0.d0
          xavec(1) = xhydr(n+1)
          xavec(2) = xheli(n+1)
          xavec(5) = xoxyg(n+1)
          xavec(3) = xcarb(n+1)
          call wd_eos(rrho(n+1),t(n+1),xavec,'d',res)
          pressnp1 = res(1)
          chitp1   = res(6)
       
          oochitavg = 0.5*(1./chitn + 1./chitp1)

          bled_new(n) = -oochitavg*(dlog(pressp) - &
               dlog(pressn))/(dlog(pressnp1) - dlog(pressn))
       enddo
       bled_new(np)= bled_new(np-1)
       bled_new(1) = bled_new(2)

! do weighting so that the differences are centered on shell n
! Try calculations with and without this weighting to see if
! it does any good...or damage!
! fac1() is the weighting for bled(n+1/2) and fac2() is the weighting 
! for bled(n-1/2), if you know what I mean
    ! fac1(:) = (dlog(p(n) - dlog(p(n+1))/(p(n-1)-p(n+1))
       fac1(2:np-1) = (dlog(parr(2:np-1)) - dlog(parr(3:np)))/ &
            (dlog(parr(1:np-2))-dlog(parr(3:np)))
       fac1(np) = 0.
       fac2(2:np) = 1. - fac1(2:np)
       bled_new(2:np) = fac2(2:np)*bled_new(1:np-1) + fac1(2:np)*bled_new(2:np)
! Some NaNs close to surface
       bled_new(np-3:np) = 0.d0

! End of "New Ledoux" calculation
    end if

    do n=2,np
       xi(n)=dlog(r(n)/p(n))
       chtor(n)=cht(n)/chr(n)
       v(n)=g*rrho(n)*mr(n)/p(n)/r(n)
       ga1(n)=chr(n)+p(n)*cht(n)**2/rrho(n)/t(n)/cv(n)
       voga1(n)=v(n)/ga1(n)
       u(n)=pi4*rrho(n)*r(n)**3/mr(n)

!bv frequency calculated here
       if (newledoux) then
          ram(n)=v(n)*chtor(n)*(del(n)-delad(n)-bled_new(n))
       else
          ram(n)=v(n)*chtor(n)*(del(n)-delad(n)-xx(20,n))
       end if
     
       bvfrq(n)=-g*mr(n)*ram(n)/r(n)**3
!debugging
!       write(*,'(i4,10f10.5)') n,xhe(n),xcomp(np+2-n,1),xcomp(np+2-n,2),xcomp(np+2-n,3),xx(21,n),xx(20,n),bled_new(n)
!       write(*,'(i4,10f10.5)') n,xi(n),xcomp(np+2-n,1),xhe(n),xcomp(np+2-n,2),xx(21,n),xx(20,n),bled_new(n)

!       if (n .eq. 170) bvfrq(n) = bvfrq(n-1) !For testing purposes
!       if (n .gt. np-5) bvfrq(n) = -10.d0 !For testing purposes
! modified Ledoux case
       ra(n)=ram(n)
       bvfreq(n)=bvfrq(n)

       if ( makeplots .or. pulsoutput  ) then
! compute Schwarzschild version for plotting comparison purposes
          ram(n)=v(n)*chtor(n)*(del(n)-delad(n))
          bvfrq(n)=-g*mr(n)*ram(n)/r(n)**3
          bvfreq3(n)=bvfreq(n)
          llp1=2.
          acous(n)=(llp1*ga1(n)*p(n))/(r(n)*r(n)*rrho(n))
          soundd=p(n)*ga1(n)/rrho(n)
          sound(n)=dsqrt(soundd)
          cs(n)=sound(n)
          if (xhe(n).gt.xhe(n+1).and.xhe(n).gt.xhe(n-1)) ixheflag=1
          if (ixheflag.eq.0) then
             carb = 1.d0 - xhe(n) - xx(21,n)
             xh = 0.0
          else
             carb = 0.0
             xh = 1.d0 - xhe(n)
          endif
          aavg=16.*xx(21,n)+12.*carb+4.*xhe(n)+1.*xh
          zavg=8.*xx(21,n)+6.*carb+2.*xhe(n)+1.*xh
          ell=2.
          ellhat=ell*(ell+1.)
          anumdens=rrho(n)/(aavg*mp)
          asep=(4.*pi*anumdens/3.)**(-1./3.)
          mu(n)=0.37*anumdens*zavg**2*elec**2/asep
          alambda(n)=ga1(n)*p(n)-(2./3.)*mu(n)
          akrdsig2(n)=rrho(n)/mu(n)
          tfreq(n)=(ellhat-2.)*mu(n)/(r(n)**2*rrho(n))
       end if
       
    end do

    do n=2,np

       if ( makeplots ) then
          cp(n)=cv(n)*ga1(n)/chr(n)
          if(cp(n) .lt. 0.0)then
             write(*,*) 'Warning; Cp<0 ',n,cp(n),chr(n)
             cp(n) = -cp(n)
          endif         
          c4(n)=gorsr*cp(n)/a/c/kap(n)/t(n)**3
          alfa(n)=delad(n)/del(n) 
          bfac=pi4*r(n)**3*rrho(n)/lr(n)
          bfacep=bfac*eps(n)
          if (ra(n) .lt. 0.) then
             eta(n)=bfacep
          endif
          b1(n)=(delad(n)*(kapt(n)-4.)+kapr(n)/ga1(n)) *v(n)+(v(n)+ &
               derdel(n))*alfa(n)
          b2(n)=(epsr(n)*voga1(n)+epst(n)*v(n)*delad(n))*bfacep &
               -(eta(n)-bfacep)*voga1(n)
          b3(n)=(epst(n)-epsr(n)*chtor(n))*bfacep+(eta(n)-bfacep)*chtor(n)
          b4(n)=bfac*gorsr*t(n)*cp(n)
          adden1 = (r(n+1)**3 - r(n)**3)*(rrho(n+1)+rrho(n))/2.
          piece(n+1) = pi4*adden1/3.
       end if
       r(n)=r(n)/totr
       mr(n)=mr(n)/totm 
       if ( makeplots ) then
          c1(n)=r(n)**3/mr(n)
       end if
   
    end do

    if ( pulsoutput ) then
       call asymp(np,r,sound,totr,totm,akrdsig2,mu,mr,rrho)
    end if

!temporary
    do n=2,np
       if (bvfrq(n) .le. 0.0) then
          nconv = n
          exit
       end if
    end do
    
    if ( makeplots ) then
       
       if ( .not. pulsoutput ) then
          call asymp(np,r,sound,totr,totm,akrdsig2,mu,mr,rrho)
       end if
!Take logs of acoustic and BV frequencies
       do n=2,np
          if (abs(bvfreq(n)).ne.0.0) then
             bvfreqmag(n)=(bvfreq(n)/abs(bvfreq(n)))*dlog10(abs(bvfreq(n)))
          else
             bvfreqmag(n)=0.0
          endif
          if(acous(n).le.0.0)then 
             acous(n)=0.0 
          else
             acous(n)=dlog10(acous(n))
          endif
          if(tfreq(n).le.0.0)then 
             tfreq(n)=0.0 
          else
             tfreq(n)=dlog10(tfreq(n))
          endif
          if(bvfreq(n).le.0.0)then
             bvfreq(n)=-7.0
          else
             bvfreq(n)=dlog10(bvfreq(n))
          endif
          if(bvfrq(n).le.0.0)then
             bvfrq(n)=-7.0
          else
             bvfrq(n)=dlog10(bvfrq(n))
          endif
       enddo

! compute thermal timescale
       npm1=np-1
       tthl(np)=0.
       tthl(1)=0.
       tth=0.
       alstar=lr(np) 
       do i=2,npm1
          n=np-i+1
          dxi=xi(n+1)-xi(n)
          dtn=pi4*cv(n)*t(n)*rrho(n)*(r(n)*totr)**3
          dtd=alstar*(1.0+v(n))
          dtth=dtn/dtd
          tth=tth+(dtth*dxi)
          if(tth.le.0.0)then
             tth=-tth
          endif
          tthl(n)=dlog10(tth)
          if(tthl(n).lt.0)then
             tthl(n)=0.
          endif

       end do

332    format(F12.1,2F10.3)
333    format(F12.1,2F10.3,F6.2,2E24.16E2)
334    format(F8.1,2F8.3,F6.2,ES24.16E2,F24.16,3F24.16)
335    format(F8.1,2F8.3,F6.2,2F24.16)

!       write(8000,332) te, smm, dlog10(gs)
!       write(8000,333) te, smm, dlog10(gs), aml, 1.0d0 -  mr(nconv), r(nconv)
!       write(8000,334) te, smm, dlog10(gs), aml, 1.0d0 -  mr(nconv), r(nconv), &
!            dlog10(t(nconv)), dlog10(p(nconv)), tthl(nconv)
!       write(*,*) logg_init, dlog10(gs)
       write(8000,335) te,smm,dlog10(gs),aml,dlog10(t(nconv)),dlog10(p(nconv))
!       write(8000,333) te, dlog10(gs), tthl(nconv), aml, 1.0d0 -  mr(nconv), &
!            r(nconv)

! finished thermal timescale computations

       call plots(np,totm)
 
       write(8,3000)
       write(77,3010) 
3010   format('#  n     radius    -log(1-Mr/M*)        Mr         ln(r/p) &
                    U            V            Gamma1         A', &
            /'#  -  ------------  ------------  ------------  ------------ &
            ------------  ------------  ------------  ------------')
       do n=2,np
          write(8,3001) n,r(n),partt(n),xi(n),u(n),v(n),voga1(n),ra(n), &
               chtor(n),c1(n)
          write(77,3001) n,10**rstr*r(n),partt(n),mr(n),xi(n),u(n),v(n), &
               v(n)/voga1(n),ra(n)/(10**rstr*r(n))
       end do
       write(8,3002)
       do n=2,np
          write(8,3001) n,c4(n),b1(n),b2(n),b3(n),b4(n),alfa(n),tthl(n), &
               derdel(n),eta(n)
       end do
       
    end if

3000 format(60x,'computed quantities'/ 45x,'(defined in Saio and Cox, &
          Ap. J. 236, 558 (1980))'/ 53x,'center and surface points &
          deleted'/ 4x,'n',6x,'Fr',11x,'log q',8x,'ln(r/p)',7x,'U',13x, &
          'V',13x,'Voga1',9x,'Ra',12x,'chToR',9x,'C1'/)
3001 format(i4,9(1pe14.6))
3002 format(4x,'n',6x,'C4',12x,'B1',12x,'B2',12x,'B3',12x,'B4',12x, &
          'alfa',11x,'Tthl',9x,'derdelad',6x,'eta'/)
3005 format(f5.3,x,i4,x,e12.5,x,2(f14.10,x),f11.3,i4,x,f13.11)
3006 format(i4,x,8(e12.5,2x))
3007 format(2(1pe9.3,2x),0p,4(f7.4,x),x,1pe10.3)
9001 format(e22.15)
    
    if ( makeplots ) then
       do n=2, np
! write out log of rad. luminosity
          if(lr(n).gt.0.0)then
             radl = dlog10(lr(n))
          else
             radl = 0.0
          endif
          write(71,3006) n,partt(n),bvfreq(n),tthl(n),radl,bvfrq(n), &
               bvfreqmag(n),mr(n)/mr(np),xi(n)
       end do
       write(50,3005) mstar,model,age,llsun,rrsun,teff,np,glog
       write(50,3007) amhyhe,amheca,alph(1),alph(2),alph(3),alph(4),stpms
       npm1=np-1
       write (19,1000) npm1
       write (19,9001) totr
       write (19,9001) totm
       write (19,1001) tthl(2:np)
       write (19,1001) r(2:np)
       write (19,1001) xi(2:np)
       write (19,1001) u(2:np)
       write (19,1001) v(2:np)
       write (19,1001) voga1(2:np)
       write (19,1001) ra(2:np)
       write (19,1001) c1(2:np)
       write (19,1001) chtor(2:np)
       write (19,1001) eta(2:np)
       write (19,1001) c4(2:np)
       write (19,1001) b1(2:np)
       write (19,1001) b2(2:np)
       write (19,1001) b3(2:np)
       write (19,1001) b4(2:np)
       write (19,1001) del(2:np)
       write (19,1001) kapr(2:np)
       write (19,1001) kapt(2:np)
       write (19,1001) alfa(2:np)
       rewind 19
    end if

! below is everything needed for the pulse code

    nnsurf=np-1
    nmod=model
    amass=totm
    brad=dlog10(totr)

5007 format(f10.5,i5)
5008 format(e12.4)
5009 format(i4)
5010 format(i5)
5011 format(4e20.12)
5012 format(4e20.12)  
5013 format(5e20.12)
    if ( makeplots ) then
       write(28,5007) brad,nmod
       write(28,5008) amass
       write(28,5009) nnsurf
       write(29,5010) nnsurf
    end if
    do i=2,np
       x=xi(i)
       xi2(i-1) = xi(i)
       rn(i)=totr*r(i)
       mr(i)=mr(i)*totm
       mr2(i) = mr(i)
       ggrav(i)=g*mr(i)/rn(i)**2
       dens=rrho(i)
       yy1(i)=ggrav(i)/rn(i)
       y2=voga1(i)
       yy3(i)=-1.*ra(i)
       y4=u(i)
       yy5(i)=1./(1.+v(i))
       if ( makeplots ) then
          write(28,5013) x,rn(i),ggrav(i),dens,mr(i)
          write(29,5011) x,yy1(i),y2,yy3(i)
          write(29,5012) y4,yy5(i),mu(i),alambda(i)
!          if ( (yy3(i).gt.0) .and. (i.lt.(np-3)) ) then
          if (yy3(i) .gt. 0.) then
             sqrtn2=sqrt(yy3(i)*yy1(i))
          else
             sqrtn2=0.0
          endif
          gform=yy5(i)*sqrtn2
          pform=yy5(i)*rn(i)/cs(i)
          dNsum=dNsum+(sqrtn2/rn(i))*totr*(r(i+1)-r(i))
          if (partt(i-1).lt.18.) then
             write(72,5019) x,gform,partt(i),sqrtn2,cs(i),pform,dNsum
          endif
       end if
5019  format(7(e12.5,x))

    end do
    
    if ( pulsemodel ) then
       call pulse
    end if

    return

  end subroutine pulsate
