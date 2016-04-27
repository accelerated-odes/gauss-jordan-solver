module gauss_jordan
  implicit none

contains

  subroutine gauss_jordan_solve(A, x, b)
    double precision, dimension(4,4), intent(in) :: A
    double precision, dimension(4), intent(out) :: x
    double precision, dimension(4), intent(in) :: b

    x(1) = A(1,3)*A(2,4)*(b(4) - A(4,2)*b(3)/A(3,2) + A(1,3)*A(4,1)*b(2)/(A(1,1)* &
      A(2,3)) - A(4,1)*b(1)/A(1,1))/(A(1,1)*A(2,3)*(-A(3,4)*A(4,2)/ &
      A(3,2) + A(1,3)*A(2,4)*A(4,1)/(A(1,1)*A(2,3)))) - A(1,3)*b(2)/( &
      A(1,1)*A(2,3)) + b(1)/A(1,1)
    x(2) = -A(3,4)*(b(4) - A(4,2)*b(3)/A(3,2) + A(1,3)*A(4,1)*b(2)/(A(1,1)*A(2,3)) &
      - A(4,1)*b(1)/A(1,1))/(A(3,2)*(-A(3,4)*A(4,2)/A(3,2) + A(1,3)* &
      A(2,4)*A(4,1)/(A(1,1)*A(2,3)))) + b(3)/A(3,2)
    x(3) = -A(2,4)*(b(4) - A(4,2)*b(3)/A(3,2) + A(1,3)*A(4,1)*b(2)/(A(1,1)*A(2,3)) &
      - A(4,1)*b(1)/A(1,1))/(A(2,3)*(-A(3,4)*A(4,2)/A(3,2) + A(1,3)* &
      A(2,4)*A(4,1)/(A(1,1)*A(2,3)))) + b(2)/A(2,3)
    x(4) = (b(4) - A(4,2)*b(3)/A(3,2) + A(1,3)*A(4,1)*b(2)/(A(1,1)*A(2,3)) - A(4,1) &
      *b(1)/A(1,1))/(-A(3,4)*A(4,2)/A(3,2) + A(1,3)*A(2,4)*A(4,1)/( &
      A(1,1)*A(2,3)))
  end subroutine gauss_jordan_solve
end module gauss_jordan
