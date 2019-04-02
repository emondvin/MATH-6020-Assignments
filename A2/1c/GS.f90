subroutine solveLSDIAG(n,ndiag,ioff,A,sol,rhs,nit,err1,err2,stopcrit)
!------------------------------------------------------------------------
! hje for M6020
!------------------------------------------------------------------------
! solves the linear system Ax=b in Diagonal format using the
! GAUSS-SEIDEL method
!------------------------------------------------------------
! input:  n       problem size
!         ndiag:  number of diagonals
!         ioff:   offsets (distance of sub diagonals to main diagonal)
!         A:      matrix values
!         sol:    initial guess for iteration (will be overwritten with result)
!         rhs:    the righ hand side of the linear system
!         nit:    maximum number of iterations to be carried out (will be overwritten)
!         err1:   tolerance for 1st stopping criterion (will be overwritten)
!         err2:   tolerance for 2nd stopping criterion (will be overwritten)
! output: sol:    solution of Ax=b
!         nit:    number of iterations taken
!         err1:   computed value for 1st stopping criterion
!         err2:   computed value for 2nd stopping criterion
!         stopcrit:  integer value that contains information which stopping criterions became active
!---------------------------------------------------------------------------------------------------
implicit none
 integer, intent(in)               :: n,ndiag
 real,dimension(n,ndiag),intent(in):: a
 integer, dimension(ndiag),intent(in)::ioff
 real,dimension(n),intent(in)      :: rhs
 real,dimension(n),intent(inout)   :: sol
 real,intent(inout)                :: err1,err2
 integer,intent(inout)             :: nit
 integer,intent(out)               :: stopcrit

 integer           :: maxit,j,i,id
 real              :: tol1,tol2,normx,normA,normb
 real,dimension(n) :: res,update,y,rdiag

open(10, file='GSr.txt')
open(11, file='GSx.txt')

 ! initialization: set tolerances, max number of iterations
 ! and rough estimates of matrix and rhs norms for stopping criteria
 maxit=nit; tol1=err1; tol2=err2
 normA=maxval(abs(a)); normb=maxval(abs(rhs))
 

 ! extract diagonal from matrix A
 ! and store its reciprocals in rdiag
 !-----------------------------------
 id=1
 do while (ioff(id)/=0)
  id=id+1
 enddo
 


 ! initialize iteration
 nit=0; stopcrit=0
 ! jacobi iteration loop
 do while(stopcrit==0)
   nit=nit+1
   
   ! the next four lines are the core of the Jacobi method
   y=sol
   call gsdiag(n,sol,rhs,A,ndiag,ioff,id)
   update=sol-y
   normx=sqrt(sum(y*y))/n
   ! compute certain norms for use in stopping criteria
   err1=1.
   err2=sqrt(sum(update*update))/n   

   ! test for convergence
   !---------------------
   if (nit>maxit) stopcrit=-1
   !if (err1<tol1*(norma*normx+normb)) stopcrit=stopcrit-10
   if (err2<tol2*normx) stopcrit=stopcrit-100

   ! uncomment the next line to monitor convergence progress
   write(*,'(I6,4(E14.7,X))') nit,err1,err2,maxval(sol),minval(sol)
   write(10,*) nit, err1
   write(11,*) nit, err2
 enddo
close(10)
close(11)
end subroutine


subroutine GSdiag (n,x,b,a,idiag,ioff,id) 
!----------------------------------------------------------------------------
! carries out the actual GAUSS-Seidel step
! see Barrett et al, Templates
!----------------------------------------------------------------------------
! input:   n:     problem size
!          idiag: number of diagonal
!          ioff:  offsets of diagonals
!          id:    index of main diagonal
!          b:     right hand side of system
!          A:     matrix values
!          x:     solution values of previous iteration (will be overwritten)
! output:  x:     new iterate
!---------------------------------------------------------------------------- 
implicit none
  integer, intent(in)::  n, idiag,id
  integer, intent(in),dimension(idiag) :: ioff
  real, dimension(n), intent(inout) :: x  
  real, dimension(n), intent(in) :: b  
  real, dimension(n,idiag), intent(in) :: a
  integer :: j, io, i1, i2, i,k
  
  do i=1,n
    x(i)=b(i)
    do j=1,idiag
       k=i+ioff(j)
       if (k>=1 .and. k<=n .and. j/=id) x(i)=x(i)-a(i,j)*x(k)
    enddo
    x(i)=x(i)/a(i,id)
  enddo

end subroutine

