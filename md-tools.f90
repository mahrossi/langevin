! a couple of functions needed for dynamics, and shared variables
! obviously for any serious work better prng should be used

module md_tools
  use splines
  implicit none
  integer, parameter ::  ndim=1
  real*8 :: kt
  real*8 :: qmax
  real*8 :: wtw(1,1)
contains

real*8 function rang(idum)
  implicit none
  integer, intent(inout) :: idum
  integer iset
  real(8) gset,v1,v2,fac,rsq
  save iset,gset
  data iset/0/
  if (iset.eq.0) then
1    v1=2.d0*randu(idum)-1.d0
     v2=2.d0*randu(idum)-1.d0
     rsq=v1**2+v2**2  
     if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1 
     fac=sqrt(-2.d0*log(rsq)/rsq) 
     gset=v1 * fac
     rang=v2 * fac
     iset=1
  else 
     rang=gset
     iset=0 
  endif
  return
end function rang

 double precision function randu(idum)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
!
!     ------------------------------------------------------------------
!     Numerical Recipes ran2
!     ------------------------------------------------------------------
!
parameter (im1=2147483563,im2=2147483399,am=1.d0/im1,imm1=im1-1, &
    ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791, &
    ntab=32,ndiv=1+imm1/ntab,eps=1.d-15,rnmx=1.d0-eps)
    dimension iv(ntab)
    save iv,iy,idum2
    data iv/ntab*0/, iy/0/, idum2/123456789/

    if (idum .le. 0) then
       idum = max(-idum,1)
       idum2 = idum
       do j = ntab+8,1,-1
          k = idum/iq1
          idum = ia1*(idum-k*iq1)-k*ir1
          if (idum .lt. 0) idum = idum+im1
          if (j .le. ntab) iv(j) = idum
       enddo
       iy = iv(1)
    endif
    k = idum/iq1
    idum = ia1*(idum-k*iq1)-k*ir1
    if (idum .lt. 0) idum = idum+im1
    k = idum2/iq2
    idum2 = ia2*(idum2-k*iq2)-k*ir2
    if (idum2 .lt. 0) idum2 = idum2+im2
    j = 1+iy/ndiv
    iy = iv(j)-idum2
    iv(j) = idum
    if (iy .lt. 1) iy = iy+imm1
    randu = min(am*iy,rnmx)
    return
end function

!real*8 function ran2(idum)
!  implicit none
!  real *8 x
!  integer, intent(inout) :: idum
!  integer iseed,i
!  integer, allocatable :: seed(:)
!    
!  if(idum.le.0) then 
!     idum=-idum
!     call random_seed(size=iseed) 
!     allocate(seed(iseed))
!     do i=1,iseed  !ugly. once again, this is just a stub. you should use a GOOD prng!
!       seed(i)=idum+i
!     end do
!     call random_seed(put=seed)
!  endif
!  call random_number(harvest=x)
!  ran2=x
!  return
!end function ran2

! returns force in a 1d potential
subroutine force(q, potxdata, potydata, ndata, splinedy2, f)
  real*8, intent(inout) :: q
  real*8, intent(in) :: potxdata(ndata), potydata(ndata), splinedy2(ndata)
  integer, intent(in) :: ndata
  real*8, intent(out) :: f
  integer :: i
  real*8 :: dx, fp, fm

! dx=0.0001
!
!  q=q+dx
!  call splint(potxdata, potydata, splinedy2, ndata, q, fp)
!  q=q-2*dx
  call splintder(potxdata, potydata, splinedy2, ndata, q, f)
  f=-f
!  q=q+dx
!  f= -(fp-fm)/(dx+dx)
  if (q>qmax .or. q<-qmax) then
     f=f-2*1.098d-05*q
  endif

end subroutine

! Below force for a harmonic potential
!subroutine force(q, f)
!  real*8, intent(in) :: q(:)
!  real*8, intent(out) :: f(:)
!  integer i,j
!  f=0.d0
!  do i=1,ndim
!    do j=1,ndim
!      f(i)=f(i)-wtw(i,j)*q(j)
!    end do
!  end do
!end subroutine

subroutine pot(q, potxdata, potydata, ndata, splinedy2, v)
  real*8, intent(inout) :: q
  real*8, intent(in) :: potxdata(ndata), potydata(ndata), splinedy2(ndata)
  integer, intent(in) :: ndata
  real*8, intent(out) :: v
  integer :: i

  v=0.d0

  call splint(potxdata, potydata, splinedy2, ndata, q, v)
  if (q>qmax .or. q<-qmax) then
     v=v+0.0039+1.098d-05*q
  endif

end subroutine

! Below potential for a harmonic potential
!subroutine pot(q, v)
!  real*8, intent(in) :: q(:)
!  real*8, intent(out) :: v
!  integer i,j
!  v=0.d0
!  do i=1,ndim
!    do j=1,ndim
!      v=v+q(i)*wtw(i,j)*q(j)
!    end do
!  end do
!  
!  v=v*0.5
!end subroutine

subroutine kin(p, k, mass)
  real*8, intent(in) :: p, mass
  real*8, intent(out) :: k
  integer i
  k=0.d0

  k=p*p*0.5d0/mass
end subroutine

end module md_tools
