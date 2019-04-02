program tridiag
!-------------------------------
!Author: Vincent Emond, 0746222
!-------------------------------
implicit none
 integer n
 parameter(n=9)
 real,dimension(n) :: a,b,c,r,x
 integer :: i

!setup the arrays a,b,c,r whhch correspond to the diagonal entries  which form matrix A, and the resultant vector r as given in the Assignment.
do i=1,n
a(i)=2
enddo

b(1)=0
do i=2,n
b(i)=-1
enddo

do i=1,n-1
c(i)=-1
enddo

r(1)=1.
do i=2,n
r(i)=0
enddo

!Calling the subroutine which solves the problem

call thomas(n,a,b,c,r,x)

do i=1,n
 write(*,*) i, x(i)
enddo

endprogram



subroutine thomas(n,a,b,c,r,x)
!-----------------------------
!Author: Vincent Emond
!-----------------------------
!Thomas Algorithm subroutine
!----------------------------
!input: n: problem size
!       a,b,c: arrays of legnth n containing coefficients of tridiagonal matrix A
!       r: array of length n containing right hand side
!output: x:solution of the linear system
!       l,d,u: coefficients for LU decomposition
!----------------------------------------------------------------------------------
implicit none
 integer,intent(in) :: n
 real,dimension(n),intent(in) :: a,b,c,r
 real,dimension(n) :: l,u,d,z
 real,dimension(n),intent(out) :: x

integer :: i
!convert the diagonal arrays of matrix A to corresponding coefficients l,u,d

do i=1,n-1
 u(i)=c(i)
enddo

d(1)=a(1)
do i=2,n
    d(i)=a(i)-(b(i)*c(i-1))/d(i-1)
enddo

l(1)=0
do i=2,n
l(i)=b(i)/d(i-1)
enddo

!Need forward and bacward sub in order to solve x
!Forward Lz=r

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

