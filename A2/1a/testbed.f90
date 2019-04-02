program test_linear_solvers_in_diagonal_format
!----------------------------------------------------------------------------
! hje, for M6020
!----------------------------------------------------------------------------
! a test bed for solvers of linear systems Ax=b 
! where A is stored in sparse diagonal format
!----------------------------------------------------------------------------
! in the main program
! (1) a test problem is set up:  subroutine genDIAG
! (2) problem is solved: subroutine solveLSDIAG
! (3) report on results is given
!----------------------------------------------------------------------------
! the linear solver is supplied externally, regardless of the method used
! but they all should use the same interface to exchange data with the
! calling program; this is explicitly defined in an interface
!----------------------------------------------------------------------------
! parameters
!----------------------------------------------------------------------------
implicit none

INTERFACE
subroutine solveLSDIAG(n,ndiag,ioff,A,sol,rhs,nit,err1,err2,stopcrit)
!---------------------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------------------
implicit none
 integer, intent(in)               :: n,ndiag
 real,dimension(n,ndiag),intent(in):: a
 integer, dimension(ndiag),intent(in)::ioff
 real,dimension(n),intent(in)      :: rhs
 real,dimension(n),intent(inout)   :: sol
 real,intent(inout)                :: err1,err2
 integer,intent(inout)             :: nit
 integer,intent(out)               :: stopcrit
end subroutine
END INTERFACE
!----------------------------------------------------
! this is where the variable declaration of he actual
! program begins
!----------------------------------------------------
 integer :: n,ndiag
 parameter (n=1000000, ndiag=9)
 real,dimension(n,ndiag) ::A 
 integer,dimension(ndiag) :: ioff
 real,dimension(n) :: rhs,sol,res,y,x,defect
 integer :: nit,stopcrit,i
 real :: err1,err2

 integer :: tstart,tfinish,clock_rate
 real :: elapsed_time


 ! initialise: set tolerances, max no iterations
 ! and initial guess
 !----------------------------------------------
 err1=1e-12;  err2=1e-12;  nit=n*100;  x=0.


 !(1) setup a test case: prescribe solution sol
 ! and generate a matrix A and rhs, such that A*sol=rhs
 !-----------------------------------------------------
 call genDIAG(n,ndiag,ioff,A,rhs,sol)


 !(2) call linear solver
 !-----------------------------------------------------------
 call system_clock(tstart)
 call solveLSDIAG(n,ndiag,ioff,A,x,rhs,nit,err1,err2,stopcrit)
 call system_clock(tfinish, clock_rate)
 elapsed_time = float(tfinish-tstart) / float(clock_rate)


 !(3) report results: 
 !-------------------------------
 call amuxd (n,x,y,a,ndiag,ioff)
 res=y-rhs    ! compute residual
 defect=sol-x ! compute defect
 write(*,'(A30,I6,X,E14.7,X,E14.7,XI6)') 'done: it/err1/err2/stopcrit',nit,err1,err2,stopcrit
 write(*,'(A30,2(E14.7,X))')'max defect, max residual: ',maxval(abs(defect)),maxval(abs(res))
 write(*,'(A30,E14.7)') 'CPU time in seconds: ',elapsed_time
end program




subroutine genDIAG(n,ndiag,ioff,A,rhs,x)
!---------------------------------------
! author: hje for M6020
!---------------------------------------
! generates a matrix A in diagonal format
! and determines the rhs to go with known
! solution x
!----------------------------------------
! input
!   n:     problem size
!   ndiag: number of off-diagonals
! output
!   ioff:  diagonal offsets
!   A:     matrix values as n x ndiag array
!   rhs:   rhs for matrix A and solution x
!   x:     solution for Ax=b
!-------------------------------------------
implicit none
 integer, intent(in) :: n,ndiag
 integer,intent(out),dimension(ndiag) :: ioff
 real, dimension(n),intent(out) :: rhs
 real, dimension(n,ndiag),intent(out) :: a
 real,dimension(n),intent(out) :: x

 integer :: i,j,sw,nd
 real :: rn


sw =1 ! sets type of matrix to be created:
!---------------------------------------------------
! (0) random, diagonally dominant
! (1) symmetric; in each off-diagonal, all entries are the same
!---------------------------------------------------

write(*,'(A23,I9,A10,I3,A17,I3)') 'genDIAG creates matrix with ',n,' rows and ',ndiag,' diagonals, type:',sw


! some preparatory work that is the same for all times of matrices
!----------------------------------------------------------------

! set the location of off-diagonals
!-----------------------------------
! make sure that
! (+) all offsets are <n
! (+) no diagonal index is used more than once
! (+) order offsets from smallest to largest (not required but recommended)
! (+) the main diagonal is in the matrix (not required but recommended)
!--------------------------------------------------------------------------- 
i=floor(sqrt(n*1.))
ioff=(/ -i,-100,-40,-1,0,1,40,100,i /)


! automatically detect position of main diagonal
!-----------------------------------------------
do i=1,ndiag
  if (ioff(i)==0) nd=i
enddo

! some tests
!--------------------------------------------------------
if (nd>n) then
  write(*,'(A39)') 'genDIAG: main diagonal not in matrix'
else
  write(*,'(A28,I9)') 'genDIAG has main diagonal at',nd
endif


if (maxval(abs(ioff))>n) then
  write(*,*) 'offdiagonal outside matrix'
  stop
endif 


select case(sw)
case(0)
  write(*,'(A30,I9,A10,I3,A10)') 'genDIAG creates random matrix with ',n,' rows and ',ndiag,' diagonals'

  ! fill matrix with random numbers between -1 and 1 
  ! note: also elements of the array that are not part of
  ! the matrix are filled here, for simplicity;
  ! these are not accessed by the mat-vec product,
  ! so we should be good
  do i=1,n
     do j=1,ndiag
        call random_number(rn)
        a(i,j)=rn
     enddo
  enddo

  ! this next line assures diagonal dominance
  a(:,nd)=a(:,nd)+ndiag*1.

case(1)
  write(*,'(A30,I9,A10,I3,A10)') 'genDIAG creates fixed matrix with ',n,' rows and ',ndiag,' diagonals'

  a(:,1) =-1.;  a(:,9) =a(:,1)+0.01
  a(:,2) =-2.;  a(:,8) =a(:,2)
  a(:,3) =-3.;  a(:,7) =a(:,3)
  a(:,4) =-4.;  a(:,6) =a(:,4)
  a(:,5) =25
  
end select



  ! define solution vector x
  x=1.

  ! compute rhs as rhs=A*x
  write(*,'(A31)') 'genDIAG creates right hand side'
  call amuxd(n,x,rhs,A,ndiag,ioff)
  write(*,'(A15)') 'genDIAG is done'


end subroutine



subroutine amuxd (n,x,y,diag,idiag,ioff) 
!-----------------------------------------------------------------------
!        A times a vector in Diagonal storage format (DIA) 
!        f90/f95 version of the sparskit f77 subroutine
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector when the original matrix is stored 
! in the diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! diag   = real array containing the diagonals stored of A.
! idiag  = number of diagonals in matrix.
! diag   = real array of size (ndiag x idiag) containing the diagonals
! ioff   = integer array of length idiag, containing the offsets of the
!   	   diagonals of the matrix:
!          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=A*x
!
!-----------------------------------------------------------------------
implicit none
  integer, intent(in)::  n, idiag
  integer, intent(in),dimension(idiag) :: ioff
  real, dimension(n), intent(in) :: x    
  real, dimension(n,idiag), intent(in) :: diag
  real, dimension(n), intent(out) :: y
  integer :: j, io, i1, i2, i       

      y=0.

      do j=1, idiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)
         do i=i1,i2
           y(i) = y(i)+diag(i,j)*x(i+io)
         enddo
      enddo  
end subroutine amuxd

