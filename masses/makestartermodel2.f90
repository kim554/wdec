program makestartermodel2

implicit none

real*8, allocatable, dimension(:,:) :: fid, quantities
real*8 :: logmr, b0, b, s
real*8, parameter :: en = 1.0d0, xc = 0.0d0
integer :: i, j, ndim, nshells, nsurf

10 format(0pE15.8,F9.6,ES11.4,F10.6,F9.6,0pE14.7,2F8.6)
20 format(E17.8,6F12.6)
30 format(I10)

open(unit=1,file="startermod.fid")
open(unit=2,file="startermod.raw")

do i=1,23
   read(1,*)
end do

read(1,*) ndim

call system('head -23 startermod.raw > startermod.new')
open(unit=3,file="startermod.new",position="append")
write(3,30) ndim

do i=1,23
   read(2,*)
end do

read(2,*) nshells

nsurf = ndim - nshells

allocate(fid(ndim,5),quantities(ndim,5))

! 1. logmr
! 2. logr
! 3. llsun
! 4. logp
! 5. logt

do i=1,ndim
   read(1,*) fid(i,:)
end do

do i=1,nshells
   read(2,*) quantities(i,:)
end do

quantities(nshells+1:nshells+nsurf,1) = fid(nshells+1:nshells+nsurf,1)

1 continue

do j=2,5

! "0" in the call below means that we are zero points away from nshells. 
! i.e. we use the points x(nshells), x(nshells-1) to get equation of line.
! b0 helps correct an offset between raw and fid model
! We handle luminosity differently, through straight linear interpolation
! based off constant slope from last two points of startmod.raw
   b0 = quantities(nshells,j) - fid(nshells,j)
   if (j .ne. 3) then
      call extrapolate(nshells,ndim,0,fid(:,1),fid(:,j),b0,b,s)
   else
      call extrapolate(nshells,ndim,0,quantities(:,1),quantities(:,j),0.d0,b,s)
   end if

   do i=1,nsurf
      quantities(nshells+i,j) = s * quantities(nshells+i,1) + b


      if (j .ne. 3) then
         call extrapolate(nshells,ndim,i,fid(:,1),fid(:,j),b0,b,s)
      else
         call extrapolate(nshells,ndim,0,quantities(:,1),quantities(:,j),0.d0,b,s)
      end if
   end do

end do
   
do i=1,ndim
   write(3,20) quantities(i,:), en, xc
end do
   
end program makestartermodel2

!*****************************************************************************
subroutine extrapolate(nshells,nmax,istart,x,y,b0,b,slope)

integer :: nshells, nmax, istart, i
real*8, dimension(nmax) :: x,y
real*8 :: slope, b, b0, dx, dy

! Figure out slope from last two points

i = nshells + istart

dy = y(i) - y(i-1)
dx = x(i) - x(i-1)
slope = dy/dx
b = y(i) - slope*x(i)

! Shift curve so it is correctly grafted at the end of what we have
b = b + b0

end subroutine extrapolate

!*****************************************************************************
subroutine extrapolate2(nshells,nmax,istart,x,y,b,slope0)

integer :: nshells, nmax, istart, i
real*8, dimension(nmax) :: x,y
real*8 :: slope0, slope1, slope, b, xx, yy, dx

! Take the last two slopes to extrapolate a new slope

i = nshells + istart

slope1 = (y(i-1) - y(i-2))/(x(i-1)-x(i-2))
slope0 = (y(i) - y(i-1))/(x(i)-x(i-1))

dx = x(i) - x(i-1)
slope = (slope0-slope1)/dx
b = slope0 - slope*x(i)

! Use that new slope to get an equation of the line that will help us place
! the next (r,rho) point
xx = x(i+1)
slope0 = slope*xx+b
b = y(i) - slope0*x(i)

end subroutine extrapolate2


