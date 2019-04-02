subroutine thomas(n,a,b,c,r,x)
!-----------------------------
!Author: Vincent Emond
!-----------------------------
!Thomas Algorithm subroutine
!----------------------------
!input: n: problem size
!	a,b,c: arrays of legnth n containing coefficients of tridiagonal matrix A
!	r: array of length n containing right hand side
!output: x:solution of the linear system
!	l,d,u: coefficients for LU decomposition
!----------------------------------------------------------------------------------
implicit none
 integer, intent(in) :: n 
 real, dimension(n), intent(in)::a,b,c,r
 real dimension(n)::l,u,d
 real dimension(n), intent(out)::x

integer :: i

u(1)=c(1)
do i=2,n-1
 u(i)=c(i)
end do

c(0)=0
d(1)=a(1)-(b(1)*c(0))/d(0)
do i=2,n
d(i)=a(i)-(b(i)*c(i-1))/d(i-1)
enddo

l(2)=b(2)/a(1)-(b(1)*c(0))/d(0)
do i=3,n
l(i)=b(i)/a(i-1)-(b(i-1)*c(i-2))/d(i-2)
enddo

!Need forward and bacward sub in order to solve x
!Forward Lz=r
l(1)=0
z(1)=r(1)
do i=2,n
 z(i)=r(i)-z(i-1)*l(i)
enddo

!Backward Ux=z
x(n)=z(n)/d(n)
do i=n-1,1,-1
 x(i)=(z(i)-u(i)*x(i+1))/d(i)
enddo


endsubroutine
