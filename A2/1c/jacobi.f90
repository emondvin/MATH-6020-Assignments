subroutine solveLSDIAG(n,ndiag,ioff,A,sol,rhs,nit,err1,err2,stopcrit)
!--------------------------------------------------------------------------------------------------
! HJE FOR M6020
!--------------------------------------------------------------------------------------------------
! uses the Jacobi method to solve Ax=b in diagonal format
!--------------------------------------------------------------------------------------------------
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
! implements the method in the form
!
!     x_{k+1} =  x_k  + D^{-1} * (b -  A*x_k)
!
!  which can be broken down into
!           residual:= b-A*x_k
!           update  := D^{-1}*res
!           x_{k+1} := x_k + update
!--------------------------------------------------------------------------------------------------
implicit none
 integer, intent(in)               :: n,ndiag
 real,dimension(n,ndiag),intent(in):: a
 integer, dimension(ndiag),intent(in)::ioff
 real,dimension(n),intent(in)      :: rhs
 real,dimension(n),intent(inout)   :: sol
 real,intent(inout)                :: err1,err2
 integer,intent(inout)             :: nit
 integer,intent(out)               :: stopcrit

 integer           :: maxit,j,i
 real              :: tol1,tol2,normx,normA,normb
 real,dimension(n) :: res,update,y,rdiag

open(10, file='Jr.txt')
open(11, file='Jx.txt')

 ! initialization: set tolerances, max number of iterations
 ! and rough estimates of matrix and rhs norms for stopping criteria
 maxit=nit; tol1=err1; tol2=err2
 normA=maxval(abs(a)); normb=maxval(abs(rhs))
 

 ! extract diagonal from matrix A
 ! and store its reciprocals in rdiag
 !-----------------------------------
 i=1
 do while (ioff(i)/=0)
  i=i+1
 enddo
 rdiag=1./a(:,i)


 ! initialize iteration
 nit=0; stopcrit=0
 ! jacobi iteration loop
 do while(stopcrit==0)
   nit=nit+1
   
   ! the next four lines are the core of the Jacobi method
   call amuxd(n,sol,y,A,ndiag,ioff)
   res=(rhs-y)
   update=rdiag*res
   sol=sol+update
   
   ! compute certain norms for use in stopping criteria
   normx=sqrt(sum(sol*sol))/n
   err1=sqrt(sum(res*res))/n
   err2=sqrt(sum(update*update))/n   

   ! test for convergence
   !---------------------
   if (nit>maxit) stopcrit=-1
   if (err1<tol1*(norma*normx+normb)) stopcrit=stopcrit-10
   if (err2<tol2*normx) stopcrit=stopcrit-100

   ! uncomment the next line to monitor convergence progress
   write(*,'(I6,4(E14.7,X))') nit,err1,err2,maxval(sol),minval(sol)
   write(10,*) nit,err1
   write(11,*) nit,err2
 enddo
end subroutine




