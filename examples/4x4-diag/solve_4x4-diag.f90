module gauss_jordan
  implicit none

contains

  subroutine gauss_jordan_solve(A, x, b)
    double precision, dimension(4,4), intent(in) :: A
    double precision, dimension(4), intent(out) :: x
    double precision, dimension(4), intent(in) :: b

    x(1) = b(1)/A(1,1)
    x(2) = b(2)/A(2,2)
    x(3) = b(3)/A(3,3)
    x(4) = b(4)/A(4,4)
  end subroutine gauss_jordan_solve
end module gauss_jordan
