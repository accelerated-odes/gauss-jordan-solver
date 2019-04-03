! Copyright (c) 2016, Donald E. Willcox
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! * Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
! * Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
! * Neither the name of gauss-jordan-solver nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


subroutine get_nonzero_dp(x)
  implicit none
  
  double precision :: x
  x = 0.0d0
  do while ( x == 0.0d0)
     call RANDOM_NUMBER(x)
  end do
end subroutine get_nonzero_dp

subroutine random_solve(which_iters, which_solver, N, niter, Asparse, aresid)
  use gauss_jordan_module, only : gauss_jordan_solve
  implicit none

  integer, intent(in) :: N, which_solver, which_iters, niter
  integer, dimension(N,N), intent(in) :: Asparse
  double precision, dimension(N), intent(out) :: aresid
  double precision, dimension(N,N) :: A, Acopy
  double precision, dimension(N) :: xsolv, xtrue, b
  integer, dimension(N) :: IPIV
  integer :: i, j, INFO

  do i = 1, N
     call get_nonzero_dp(xtrue(i))
     do j = 1, N
        if ( Asparse(i, j) == 1 ) then
           call get_nonzero_dp(A(i, j))
        else
           A(i, j) = 0.0d0
        end if
     end do
  end do
  b = MATMUL(A, xtrue)
  if ( which_solver == 1 ) then
     !     write(*,*) 'Using gauss-jordan solver'
     if ( which_iters == 2 ) then
        do i = 1, niter
           Acopy = A
           call gauss_jordan_solve(Acopy, xsolv, b)
        end do
     elseif ( which_iters == 1 ) then
        call gauss_jordan_solve(A, xsolv, b)
     end if
  elseif ( which_solver == 2) then
     !     write(*,*) 'Using LAPACK DGESV'
     if ( which_iters == 2 ) then
        do i = 1, niter
           xsolv = b
           Acopy = A
           call DGESV(N, 1, Acopy, N, IPIV, xsolv, N, INFO)
        end do
     elseif ( which_iters == 1 ) then
        xsolv = b
        call DGESV(N, 1, A, N, IPIV, xsolv, N, INFO)
     end if
     if ( INFO /= 0 ) then
        write(*,*) 'LAPACK could not find a solution!'
        stop
     end if
  else
     write(*,*) 'NO SOLVER SPECIFIED. ENTER 1 or 2.'
     stop
  end if
  aresid = abs(xtrue - xsolv)
end subroutine random_solve

program test_gj_solve
  ! Takes 4 command-line arguments:
  ! 1: Name of sparsity structure file to read
  ! 2: Number of random matrix equations Ax=b to solve
  ! 3: Which solver to use (1 = gauss jordan, 2 = LAPACK DGESV)
  ! 4: Matrix generation (1 = resample matrix every iteration,
  !                       2 = randomly sample matrix only once)
  
  implicit none
  
  character(LEN=50) :: sparsity_file
  character(LEN=50) :: str_num_iters, str_solver_choice, str_iters_choice
  integer :: N, unit, niter, which_solver, which_iters
  integer, dimension(:,:), allocatable :: Asparsity
  double precision, dimension(:,:), allocatable :: aresid
  double precision, dimension(:), allocatable :: arave, arstd
  double precision :: scratch
  integer :: i, j
  logical :: invalid_arguments

  ! Initialize the random number generator
  call init_random_seed()

  invalid_arguments = .false.
  
  ! Read the first command-line argument as the sparsity file name
  call GET_COMMAND_ARGUMENT(1, sparsity_file)

  if (len(trim(sparsity_file)) .eq. 0) then
     invalid_arguments = .true.
  end if

  ! Read the second command-line argument as the number of random solutions to run
  call GET_COMMAND_ARGUMENT(2, str_num_iters)
  read(str_num_iters, '(I50)') niter
  if ( niter <= 0 ) then
     invalid_arguments = .true.
  end if
  ! Read the third command-line argument as the selection for which solver to use
  ! if which_solver = 1, then use the nonzero gauss-jordan solver
  ! if which_solver = 2, then use LAPACK DGESV
  call GET_COMMAND_ARGUMENT(3, str_solver_choice)
  read(str_solver_choice, '(I50)') which_solver
  if ( .not. ( which_solver == 1 .or. which_solver == 2 )) then
     invalid_arguments = .true.
  end if
  if ( which_solver == 1 ) then
     write(*,*) 'Using Nonzero Gauss-Jordan Solver'
  elseif ( which_solver == 2 ) then
     write(*,*) 'Using LAPACK DGESV Solver'
  end if
  ! Read the fourth command-line argument as the selection for which randomization to use
  ! if which_iters = 1, then solve different randomized matrix systems
  ! if which_iters = 2, then solve the same randomized matrix system over and over
  call GET_COMMAND_ARGUMENT(4, str_iters_choice)
  read(str_iters_choice, '(I50)') which_iters
  if ( .not. ( which_iters == 1 .or. which_iters == 2 )) then
     invalid_arguments = .true.
  end if

  if (invalid_arguments) then
     print *, "test_gj_solve takes 4 command-line arguments:"
     print *, "  1: Name of sparsity structure file to read"
     print *, "  2: Number of iterations of matrix equations Ax=b to solve"
     print *, "  3: Which solver to use (1 = gauss jordan, 2 = LAPACK DGESV)"
     print *, "  4: Matrix generation (1 = resample random matrix every iteration,"
     print *, "                        2 = randomly sample matrix only once)"
     stop
  end if

  ! Read the sparsity file and set up sparsity matrix
  open(newunit=unit, file=sparsity_file)
  ! Get size of matrix
  read(unit,*) N
  allocate( Asparsity(N, N) )
  ! Read sparsity matrix in row-major order
  read(unit,*) Asparsity
  ! Transpose matrix into column-major format for Fortran
  Asparsity = TRANSPOSE(Asparsity)
  close(unit=unit)
  
  ! Test solver with a random matrix equation using this sparsity
  allocate( arave(N) )
  allocate( arstd(N) )
  if ( which_iters == 1 ) then
     allocate( aresid(N, niter) )
     do i = 1, niter
        call random_solve(which_iters, which_solver, N, niter, Asparsity, aresid(:, i))
     end do
     ! Get statistics for residuals
     do i = 1, N
        arave(i) = SUM(aresid(i, :))/dble(niter)
        scratch = 0.0d0
        do j = 1, niter
           scratch = scratch + ((aresid(i, j) - arave(i))**2)/dble(niter)
        end do
        arstd(i) = sqrt(scratch)
     end do
     ! Print statistics for absolute values of residuals
     write(*,*) niter, ' Random Systems Solved'
     write(*,*) 'Residual Average +/- StdDev for solutions of the system Ax=b.'
     do i = 1, N
        write(*,*) 'x(', i, '): ', arave(i), ' +/- ', arstd(i)
     end do
  elseif ( which_iters == 2 ) then
     allocate( aresid(N, 1) )
     call random_solve(which_iters, which_solver, N, niter, Asparsity, aresid(:, 1))
     ! Print residuals
     write(*,*) niter, ' Solution Iterations of a Random System'
     write(*,*) 'Residuals for solution of the system Ax=b.'
     do i = 1, N
        write(*,*) 'x(', i, '): ', aresid(i, 1)
     end do
  end if

  ! Cleanup memory
  deallocate( aresid )
  deallocate( arave )
  deallocate( arstd )
  deallocate( Asparsity )
end program test_gj_solve
