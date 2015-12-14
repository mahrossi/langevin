MODULE SPLINES
IMPLICIT NONE
CONTAINS

SUBROUTINE spline(x,y,n,y2) 
INTEGER n
REAL*8 x(n),y(n),y2(n) 
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi),
! with x1 < x2 < ... < xN, and given values yp1 and ypn for the first 
! derivative of the interpolating function at points 1 and n, 
! respectively, this routine returns an array y2(1:n) of length n 
! which contains the second derivatives of the interpolating function at 
! the tabulated points xi. We set yp1 and ypn to zero here in order
! to always have natural spline 

INTEGER i,k
REAL*8 p,qn,sig,un,u(n) 

y2(1)=0.
u(1)=0. 
qn=0.
un=0.

do i=2,n-1
   sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
   p=sig*y2(i-1)+2.
   y2(i)=(sig-1.)/p
   u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
enddo 


y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.) 

do k=n-1,1,-1
   y2(k)=y2(k)*y2(k+1)+u(k)
enddo

return 
END SUBROUTINE spline


SUBROUTINE splint(xa,ya,y2a,n,x,y) 
INTEGER n
REAL*8 x,y,xa(n),y2a(n),ya(n)
 ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate 
 ! a function (with the xai’s in order), and given the array y2a(1:n), 
 ! which is the output from spline above, and given a value of x, 
 ! this routine returns a cubic-spline interpolated value y.
INTEGER k,khi,klo 
REAL*8 a,b,h
klo=1
khi=n
! Here a better performance could be attained by storing the khi/klo in 
! the previous steps, if the x are in sequence
do while (khi-klo.gt.1)  
   k=(khi+klo)/2
   if(xa(k).gt.x)then 
      khi=k
   else 
      klo=k
   endif 
enddo

h=xa(khi)-xa(klo) 
!write(*,*) x, khi, klo, xa(khi), xa(klo), h 

a=(xa(khi)-x)/h  
b=(x-xa(klo))/h
y=a*ya(klo)+b*ya(khi)+ ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

return 

END SUBROUTINE splint

END MODULE SPLINES