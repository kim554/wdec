program rewritemods

implicit none

real*8, dimension(5) :: vector
integer :: i, nshells

open(unit=1,file="input")
open(unit=2,file="output")

do i=1,23
   read(1,*)
end do

read(1,*) nshells

10 format(0pE15.8,F9.6,ES11.4,F10.6,F9.6,0pE14.7,2F8.6)
20 format(E17.8,6F12.6)

do i=1,nshells
   read(1,10) vector
   write(2,20) vector, 1.0d0, 0.0d0
end do

end program rewritemods
