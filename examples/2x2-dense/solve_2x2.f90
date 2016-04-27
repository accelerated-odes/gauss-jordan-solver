module gauss_jordan
  implicit none

contains

  subroutine gauss_jordan_solve(A, x, b)
    double precision, dimension(2,2), intent(in) :: A
    double precision, dimension(2), intent(out) :: x
    double precision, dimension(2), intent(in) :: b

    x(1) = -A(1,2)*(b(2) - A(2,1)*b(1)/A(1,1))/(A(1,1)*(A(2,2) - A(1,2)*A(2,1)/ &
      A(1,1))) + b(1)/A(1,1)
    x(2) = (b(2) - A(2,1)*b(1)/A(1,1))/(A(2,2) - A(1,2)*A(2,1)/A(1,1))
  end subroutine gauss_jordan_solve
end module gauss_jordan
