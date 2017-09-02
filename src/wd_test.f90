
program wd_test
  use wd_eos_mod, only : wd_eos, dp, log10_cr
  use kappa, only : opalz

  implicit none

  real(dp) :: rho, T, xh, xhe, xc, xo, P, rho_p, res3(11), xa(7), fkappa, opr, opt
  integer :: i
  character (len=1) :: ef1


  rho = 1.3519d2
!  T = 10**(8.15d0)
  T = 10**(6.15d0)
!  rho = 1.d-06
!  T = 10**(4.15d0)
  xh = 0.d0
  xhe = 0.2d0
  xc = 0.4d0
  xo = 0.4d0
  xa(:) = 0.0d0

  xa(1) = xh
  xa(2) = xhe
  xa(3) = xc
  xa(5) = xo

  ef1 = 'd'  ! 'd' for (Rho, T) EoS, 'p' for (P, T) EoS, 
  
  call wd_eos(rho,T,xa,'d',res3)
  call opalz(T, rho, xa, fkappa, opr, opt)
  write(*,'(12x,a,12x,a,10x,a,11x,a,12x,a,12x,a,10x,a,10x,a)') 'rho','T','Pgas','grad_ad','cp','kappa','dlnKdrho','dlnKdT'

  do i = 0,4

     T = 10**(4.d0 + real(i))
     rho = 10**(-8.d0 + 3.*real(i))
     call wd_eos(rho,T,xa,'d',res3)
     call opalz(T, rho, xa, fkappa, opr, opt)
     P = res3(1)
     write(*,'(a5,es16.9,es11.4,4es16.9,3es17.9)') '(d,T)',rho,T,P,res3(4),res3(2),fkappa,opr,opt
     
     call wd_eos(P,T,xa,'p',res3)
     rho_p = res3(1)
     call opalz(T, rho_p, xa, fkappa, opr, opt)
     write(*,'(a5,es16.9,es11.4,4es16.9,3es17.9)') '(P,T)',rho_p,T,P,res3(4),res3(2),fkappa,opr,opt

  enddo

end program wd_test

