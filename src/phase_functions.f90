module phase_func

contains

!****************************************************************************
  
  function phasec(t)

!  is it liquid metal, or is it xtal lattice?
!  ac(6) contains the melting coefficients for carbon.

    use acc

    implicit double precision(a-h,o-z)

    real*8 :: t
    
    phasec = ac(6)
    do i=1,5
       phasec = phasec * t + ac(6-i)
    enddo
    
    return
  end function phasec

!************************************************************************

  function phaseo(t)

!  is it liquid metal, or is it xtal lattice?
!  ao(6) contains the melting coefficients for oxygen.

    use aoo

    implicit double precision(a-h,o-z)
    
    real*8 :: t

    phaseo=ao(6)
    do i=1,5
       phaseo = phaseo * t + ao(6-i)
    enddo

    return

  end function phaseo

!**************************************************************************

end module phase_func
