subroutine get_nonzero_dp(x)
  implicit none
  
  double precision :: x
  x = 0.0d0
  do while ( x == 0.0d0)
     call RANDOM_NUMBER(x)
  end do
end subroutine get_nonzero_dp

subroutine random_solve(N, Asparse, aresid)
  use gauss_jordan, only : gauss_jordan_solve
  implicit none

  integer, intent(in) :: N
  integer, dimension(N,N), intent(in) :: Asparse
  double precision, dimension(N), intent(out) :: aresid
  double precision, dimension(N,N) :: A
  double precision, dimension(N) :: xsolv, xtrue, b
  integer :: i, j

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
  call gauss_jordan_solve(A, xsolv, b)
  aresid = xtrue - xsolv
end subroutine random_solve

program test_gj_solve
  ! Takes 2 command-line arguments:
  ! 1: Name of sparsity structure file to read
  ! 2: Number of random matrix equations Ax=b to solve
  
  implicit none
  
  character(LEN=50) :: sparsity_file
  character(LEN=50) :: str_num_iters
  integer :: N, unit, niter
  integer, dimension(:,:), allocatable :: Asparsity
  double precision, dimension(:,:), allocatable :: aresid
  double precision, dimension(:), allocatable :: arave, arstd
  double precision :: scratch
  integer :: i, j

  ! Initialize the random number generator
  call init_random_seed()
  
  ! Read the first command-line argument as the sparsity file name
  call GET_COMMAND_ARGUMENT(1, sparsity_file)
  ! Read the second command-line argument as the number of random solutions to run
  call GET_COMMAND_ARGUMENT(2, str_num_iters)
  read(str_num_iters, '(I50)') niter
  
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
  allocate( aresid(N, niter) )
  allocate( arave(N) )
  allocate( arstd(N) )
  do i = 1, niter
     call random_solve(N, Asparsity, aresid(:, i))
  end do
  
  ! Get statistics for residuals
  do i = 1, N
     arave(i) = SUM(aresid(i, :))/max(1, niter)
     scratch = 0.0d0
     do j = 1, niter
        scratch = scratch + ((aresid(i, j) - arave(i))**2)/max(1, niter)
     end do
     arstd(i) = sqrt(scratch)
  end do
  
  ! Print statistics for absolute values of residuals
  write(*,*) 'Residual Average +/- StdDev for solutions of the system Ax=b.'
  do i = 1, N
     write(*,*) 'x(', i, '): ', arave(i), ' +/- ', arstd(i)
  end do

  ! Cleanup memory
  deallocate( aresid )
  deallocate( arave )
  deallocate( arstd )
  deallocate( Asparsity )
end program test_gj_solve
