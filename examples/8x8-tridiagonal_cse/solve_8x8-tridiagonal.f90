module gauss_jordan
  implicit none

contains

  subroutine gauss_jordan_solve(A, x, b)
    double precision, dimension(8,8), intent(in) :: A
    double precision, dimension(8), intent(out) :: x
    double precision, dimension(8), intent(in) :: b
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
    double precision :: scratch_31
    double precision :: scratch_32
    double precision :: scratch_33
    double precision :: scratch_34
    double precision :: scratch_35
    double precision :: scratch_36
    double precision :: scratch_37
    double precision :: scratch_38
    double precision :: scratch_39
    double precision :: scratch_40
    double precision :: scratch_41
    double precision :: scratch_42

    scratch_0 = 1.0/A(1,1)
    scratch_1 = b(1)*scratch_0
    scratch_2 = A(1,2)*scratch_0
    scratch_3 = 1.0/(-A(2,1)*scratch_2 + A(2,2))
    scratch_4 = scratch_3*(-A(2,1)*scratch_1 + b(2))
    scratch_5 = A(2,3)*scratch_3
    scratch_6 = 1.0/(-A(3,2)*scratch_5 + A(3,3))
    scratch_7 = scratch_6*(-A(3,2)*scratch_4 + b(3))
    scratch_8 = scratch_5*scratch_7
    scratch_9 = A(3,4)*scratch_6
    scratch_10 = 1.0/(-A(4,3)*scratch_9 + A(4,4))
    scratch_11 = scratch_10*(-A(4,3)*scratch_7 + b(4))
    scratch_12 = scratch_11*scratch_9
    scratch_13 = scratch_12*scratch_5
    scratch_14 = A(4,5)*scratch_10
    scratch_15 = 1.0/(-A(5,4)*scratch_14 + A(5,5))
    scratch_16 = scratch_15*(-A(5,4)*scratch_11 + b(5))
    scratch_17 = scratch_14*scratch_16
    scratch_18 = scratch_17*scratch_9
    scratch_19 = scratch_18*scratch_5
    scratch_20 = A(5,6)*scratch_15
    scratch_21 = 1.0/(-A(6,5)*scratch_20 + A(6,6))
    scratch_22 = scratch_21*(-A(6,5)*scratch_16 + b(6))
    scratch_23 = scratch_20*scratch_22
    scratch_24 = scratch_14*scratch_23
    scratch_25 = scratch_24*scratch_9
    scratch_26 = scratch_25*scratch_5
    scratch_27 = A(6,7)*scratch_21
    scratch_28 = 1.0/(-A(7,6)*scratch_27 + A(7,7))
    scratch_29 = scratch_28*(-A(7,6)*scratch_22 + b(7))
    scratch_30 = scratch_27*scratch_29
    scratch_31 = scratch_20*scratch_30
    scratch_32 = scratch_14*scratch_31
    scratch_33 = scratch_32*scratch_9
    scratch_34 = scratch_33*scratch_5
    scratch_35 = A(7,8)*scratch_28
    scratch_36 = (-A(8,7)*scratch_29 + b(8))/(-A(8,7)*scratch_35 + A(8,8))
    scratch_37 = scratch_35*scratch_36
    scratch_38 = scratch_27*scratch_37
    scratch_39 = scratch_20*scratch_38
    scratch_40 = scratch_14*scratch_39
    scratch_41 = scratch_40*scratch_9
    scratch_42 = scratch_41*scratch_5

    x(1) = scratch_1 - scratch_13*scratch_2 + scratch_19*scratch_2 - scratch_2* &
      scratch_26 + scratch_2*scratch_34 - scratch_2*scratch_4 - &
      scratch_2*scratch_42 + scratch_2*scratch_8
    x(2) = scratch_13 - scratch_19 + scratch_26 - scratch_34 + scratch_4 + &
      scratch_42 - scratch_8
    x(3) = -scratch_12 + scratch_18 - scratch_25 + scratch_33 - scratch_41 + &
      scratch_7
    x(4) = scratch_11 - scratch_17 + scratch_24 - scratch_32 + scratch_40
    x(5) = scratch_16 - scratch_23 + scratch_31 - scratch_39
    x(6) = scratch_22 - scratch_30 + scratch_38
    x(7) = scratch_29 - scratch_37
    x(8) = scratch_36
  end subroutine gauss_jordan_solve
end module gauss_jordan
