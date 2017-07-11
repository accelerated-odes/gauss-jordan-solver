module gauss_jordan_module
  implicit none

contains

  subroutine gauss_jordan_solve(A, x, b)
    double precision, dimension(10,10), intent(in) :: A
    double precision, dimension(10), intent(out) :: x
    double precision, dimension(10), intent(in) :: b
    double precision :: scratch_0
    double precision :: scratch_1
    double precision :: scratch_2
    double precision :: scratch_3
    double precision :: scratch_4
    double precision :: scratch_5
    double precision :: scratch_6
    double precision :: scratch_7
    double precision :: scratch_8
    double precision :: scratch_9
    double precision :: scratch_10
    double precision :: scratch_11
    double precision :: scratch_12
    double precision :: scratch_13
    double precision :: scratch_14
    double precision :: scratch_15
    double precision :: scratch_16
    double precision :: scratch_17
    double precision :: scratch_18
    double precision :: scratch_19
    double precision :: scratch_20
    double precision :: scratch_21
    double precision :: scratch_22
    double precision :: scratch_23
    double precision :: scratch_24
    double precision :: scratch_25
    double precision :: scratch_26
    double precision :: scratch_27
    double precision :: scratch_28
    double precision :: scratch_29
    double precision :: scratch_30

    scratch_0 = 1.0/A(1,1)
    scratch_1 = b(1)*scratch_0
    scratch_2 = A(1,9)*scratch_0
    scratch_3 = 1.0/A(3,3)
    scratch_4 = A(3,4)*scratch_3
    scratch_5 = 1.0/(-A(4,3)*scratch_4 + A(4,4))
    scratch_6 = A(3,9)*scratch_3
    scratch_7 = scratch_5*(-A(4,3)*scratch_6 + A(4,9))
    scratch_8 = A(1,4)*scratch_0
    scratch_9 = 1.0/A(2,2)
    scratch_10 = scratch_9*(-A(2,1)*scratch_2 + A(2,9))
    scratch_11 = 1.0/A(5,5)
    scratch_12 = -A(5,3)*scratch_4 + A(5,4)
    scratch_13 = scratch_11*(-A(5,3)*scratch_6 + A(5,9) - scratch_12*scratch_7)
    scratch_14 = 1.0/A(6,6)
    scratch_15 = scratch_14*(-A(6,4)*scratch_7 + A(6,9))
    scratch_16 = 1.0/A(7,7)
    scratch_17 = scratch_16*(-A(7,4)*scratch_7 + A(7,9))
    scratch_18 = 1.0/A(8,8)
    scratch_19 = scratch_18*(-A(8,4)*scratch_7 + A(8,9))
    scratch_20 = scratch_9*(-A(2,1)*scratch_8 + A(2,4))
    scratch_21 = -A(9,1)*scratch_8 - A(9,2)*scratch_20 - A(9,3)*scratch_4 + A(9,4)
    scratch_22 = scratch_9*(-A(2,1)*scratch_1 + b(2))
    scratch_23 = b(3)*scratch_3
    scratch_24 = scratch_5*(-A(4,3)*scratch_23 + b(4))
    scratch_25 = scratch_11*(-A(5,3)*scratch_23 + b(5) - scratch_12*scratch_24)
    scratch_26 = scratch_14*(-A(6,4)*scratch_24 + b(6))
    scratch_27 = scratch_16*(-A(7,4)*scratch_24 + b(7))
    scratch_28 = scratch_18*(-A(8,4)*scratch_24 + b(8))
    scratch_29 = (-A(9,1)*scratch_1 - A(9,2)*scratch_22 - A(9,3)*scratch_23 - A(9,5)* &
      scratch_25 - A(9,6)*scratch_26 - A(9,7)*scratch_27 - A(9,8)* &
      scratch_28 + b(9) - scratch_21*scratch_24)/(-A(9,1)*scratch_2 - &
      A(9,2)*scratch_10 - A(9,3)*scratch_6 - A(9,5)*scratch_13 - A(9,6) &
      *scratch_15 - A(9,7)*scratch_17 - A(9,8)*scratch_19 + A(9,9) - &
      scratch_21*scratch_7)
    scratch_30 = -A(10,1)*scratch_8 - A(10,2)*scratch_20 - A(10,3)*scratch_4 + A(10,4)

    x(1) = scratch_1 - scratch_24*scratch_8 - scratch_29*(scratch_2 - scratch_7* &
      scratch_8)
    x(2) = -scratch_20*scratch_24 + scratch_22 - scratch_29*(scratch_10 - &
      scratch_20*scratch_7)
    x(3) = scratch_23 - scratch_24*scratch_4 - scratch_29*(-scratch_4*scratch_7 + &
      scratch_6)
    x(4) = scratch_24 - scratch_29*scratch_7
    x(5) = -scratch_13*scratch_29 + scratch_25
    x(6) = -scratch_15*scratch_29 + scratch_26
    x(7) = -scratch_17*scratch_29 + scratch_27
    x(8) = -scratch_19*scratch_29 + scratch_28
    x(9) = scratch_29
    x(10) = (-A(10,1)*scratch_1 - A(10,2)*scratch_22 - A(10,3)*scratch_23 - A(10,5)* &
      scratch_25 - A(10,6)*scratch_26 - A(10,7)*scratch_27 - A(10,8)* &
      scratch_28 + b(10) - scratch_24*scratch_30 - scratch_29*(-A(10,1) &
      *scratch_2 - A(10,2)*scratch_10 - A(10,3)*scratch_6 - A(10,5)* &
      scratch_13 - A(10,6)*scratch_15 - A(10,7)*scratch_17 - A(10,8)* &
      scratch_19 + A(10,9) - scratch_30*scratch_7))/A(10,10)
  end subroutine gauss_jordan_solve
end module gauss_jordan_module
