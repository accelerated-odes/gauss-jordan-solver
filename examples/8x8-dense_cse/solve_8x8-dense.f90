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
    double precision :: scratch_43
    double precision :: scratch_44
    double precision :: scratch_45
    double precision :: scratch_46
    double precision :: scratch_47
    double precision :: scratch_48
    double precision :: scratch_49
    double precision :: scratch_50
    double precision :: scratch_51
    double precision :: scratch_52
    double precision :: scratch_53
    double precision :: scratch_54
    double precision :: scratch_55
    double precision :: scratch_56
    double precision :: scratch_57
    double precision :: scratch_58
    double precision :: scratch_59
    double precision :: scratch_60
    double precision :: scratch_61
    double precision :: scratch_62
    double precision :: scratch_63
    double precision :: scratch_64
    double precision :: scratch_65
    double precision :: scratch_66
    double precision :: scratch_67
    double precision :: scratch_68
    double precision :: scratch_69
    double precision :: scratch_70
    double precision :: scratch_71
    double precision :: scratch_72
    double precision :: scratch_73
    double precision :: scratch_74
    double precision :: scratch_75
    double precision :: scratch_76
    double precision :: scratch_77
    double precision :: scratch_78

    scratch_0 = 1.0/A(1,1)
    scratch_1 = b(1)*scratch_0
    scratch_2 = A(1,3)*scratch_0
    scratch_3 = A(1,2)*scratch_0
    scratch_4 = 1.0/(-A(2,1)*scratch_3 + A(2,2))
    scratch_5 = scratch_4*(-A(2,1)*scratch_2 + A(2,3))
    scratch_6 = scratch_2 - scratch_3*scratch_5
    scratch_7 = -A(3,1)*scratch_3 + A(3,2)
    scratch_8 = 1.0/(-A(3,1)*scratch_2 + A(3,3) - scratch_5*scratch_7)
    scratch_9 = scratch_4*(-A(2,1)*scratch_1 + b(2))
    scratch_10 = scratch_8*(-A(3,1)*scratch_1 + b(3) - scratch_7*scratch_9)
    scratch_11 = A(1,4)*scratch_0
    scratch_12 = scratch_4*(-A(2,1)*scratch_11 + A(2,4))
    scratch_13 = scratch_8*(-A(3,1)*scratch_11 + A(3,4) - scratch_12*scratch_7)
    scratch_14 = scratch_11 - scratch_12*scratch_3 - scratch_13*scratch_6
    scratch_15 = -A(4,1)*scratch_3 + A(4,2)
    scratch_16 = -A(4,1)*scratch_2 + A(4,3) - scratch_15*scratch_5
    scratch_17 = 1.0/(-A(4,1)*scratch_11 + A(4,4) - scratch_12*scratch_15 - scratch_13* &
      scratch_16)
    scratch_18 = scratch_17*(-A(4,1)*scratch_1 + b(4) - scratch_10*scratch_16 - &
      scratch_15*scratch_9)
    scratch_19 = A(1,5)*scratch_0
    scratch_20 = scratch_4*(-A(2,1)*scratch_19 + A(2,5))
    scratch_21 = scratch_8*(-A(3,1)*scratch_19 + A(3,5) - scratch_20*scratch_7)
    scratch_22 = scratch_17*(-A(4,1)*scratch_19 + A(4,5) - scratch_15*scratch_20 - &
      scratch_16*scratch_21)
    scratch_23 = -scratch_14*scratch_22 + scratch_19 - scratch_20*scratch_3 - scratch_21* &
      scratch_6
    scratch_24 = -A(5,1)*scratch_3 + A(5,2)
    scratch_25 = -A(5,1)*scratch_2 + A(5,3) - scratch_24*scratch_5
    scratch_26 = -A(5,1)*scratch_11 + A(5,4) - scratch_12*scratch_24 - scratch_13* &
      scratch_25
    scratch_27 = 1.0/(-A(5,1)*scratch_19 + A(5,5) - scratch_20*scratch_24 - scratch_21* &
      scratch_25 - scratch_22*scratch_26)
    scratch_28 = scratch_27*(-A(5,1)*scratch_1 + b(5) - scratch_10*scratch_25 - &
      scratch_18*scratch_26 - scratch_24*scratch_9)
    scratch_29 = A(1,6)*scratch_0
    scratch_30 = scratch_4*(-A(2,1)*scratch_29 + A(2,6))
    scratch_31 = scratch_8*(-A(3,1)*scratch_29 + A(3,6) - scratch_30*scratch_7)
    scratch_32 = scratch_17*(-A(4,1)*scratch_29 + A(4,6) - scratch_15*scratch_30 - &
      scratch_16*scratch_31)
    scratch_33 = scratch_27*(-A(5,1)*scratch_29 + A(5,6) - scratch_24*scratch_30 - &
      scratch_25*scratch_31 - scratch_26*scratch_32)
    scratch_34 = -scratch_14*scratch_32 - scratch_23*scratch_33 + scratch_29 - scratch_3* &
      scratch_30 - scratch_31*scratch_6
    scratch_35 = -A(6,1)*scratch_3 + A(6,2)
    scratch_36 = -A(6,1)*scratch_2 + A(6,3) - scratch_35*scratch_5
    scratch_37 = -A(6,1)*scratch_11 + A(6,4) - scratch_12*scratch_35 - scratch_13* &
      scratch_36
    scratch_38 = -A(6,1)*scratch_19 + A(6,5) - scratch_20*scratch_35 - scratch_21* &
      scratch_36 - scratch_22*scratch_37
    scratch_39 = 1.0/(-A(6,1)*scratch_29 + A(6,6) - scratch_30*scratch_35 - scratch_31* &
      scratch_36 - scratch_32*scratch_37 - scratch_33*scratch_38)
    scratch_40 = scratch_39*(-A(6,1)*scratch_1 + b(6) - scratch_10*scratch_36 - &
      scratch_18*scratch_37 - scratch_28*scratch_38 - scratch_35* &
      scratch_9)
    scratch_41 = A(1,7)*scratch_0
    scratch_42 = scratch_4*(-A(2,1)*scratch_41 + A(2,7))
    scratch_43 = scratch_8*(-A(3,1)*scratch_41 + A(3,7) - scratch_42*scratch_7)
    scratch_44 = scratch_17*(-A(4,1)*scratch_41 + A(4,7) - scratch_15*scratch_42 - &
      scratch_16*scratch_43)
    scratch_45 = scratch_27*(-A(5,1)*scratch_41 + A(5,7) - scratch_24*scratch_42 - &
      scratch_25*scratch_43 - scratch_26*scratch_44)
    scratch_46 = scratch_39*(-A(6,1)*scratch_41 + A(6,7) - scratch_35*scratch_42 - &
      scratch_36*scratch_43 - scratch_37*scratch_44 - scratch_38* &
      scratch_45)
    scratch_47 = -scratch_14*scratch_44 - scratch_23*scratch_45 - scratch_3*scratch_42 - &
      scratch_34*scratch_46 + scratch_41 - scratch_43*scratch_6
    scratch_48 = -A(7,1)*scratch_3 + A(7,2)
    scratch_49 = -A(7,1)*scratch_2 + A(7,3) - scratch_48*scratch_5
    scratch_50 = -A(7,1)*scratch_11 + A(7,4) - scratch_12*scratch_48 - scratch_13* &
      scratch_49
    scratch_51 = -A(7,1)*scratch_19 + A(7,5) - scratch_20*scratch_48 - scratch_21* &
      scratch_49 - scratch_22*scratch_50
    scratch_52 = -A(7,1)*scratch_29 + A(7,6) - scratch_30*scratch_48 - scratch_31* &
      scratch_49 - scratch_32*scratch_50 - scratch_33*scratch_51
    scratch_53 = 1.0/(-A(7,1)*scratch_41 + A(7,7) - scratch_42*scratch_48 - scratch_43* &
      scratch_49 - scratch_44*scratch_50 - scratch_45*scratch_51 - &
      scratch_46*scratch_52)
    scratch_54 = scratch_53*(-A(7,1)*scratch_1 + b(7) - scratch_10*scratch_49 - &
      scratch_18*scratch_50 - scratch_28*scratch_51 - scratch_40* &
      scratch_52 - scratch_48*scratch_9)
    scratch_55 = A(1,8)*scratch_0
    scratch_56 = scratch_4*(-A(2,1)*scratch_55 + A(2,8))
    scratch_57 = scratch_8*(-A(3,1)*scratch_55 + A(3,8) - scratch_56*scratch_7)
    scratch_58 = scratch_17*(-A(4,1)*scratch_55 + A(4,8) - scratch_15*scratch_56 - &
      scratch_16*scratch_57)
    scratch_59 = scratch_27*(-A(5,1)*scratch_55 + A(5,8) - scratch_24*scratch_56 - &
      scratch_25*scratch_57 - scratch_26*scratch_58)
    scratch_60 = scratch_39*(-A(6,1)*scratch_55 + A(6,8) - scratch_35*scratch_56 - &
      scratch_36*scratch_57 - scratch_37*scratch_58 - scratch_38* &
      scratch_59)
    scratch_61 = scratch_53*(-A(7,1)*scratch_55 + A(7,8) - scratch_48*scratch_56 - &
      scratch_49*scratch_57 - scratch_50*scratch_58 - scratch_51* &
      scratch_59 - scratch_52*scratch_60)
    scratch_62 = -A(8,1)*scratch_3 + A(8,2)
    scratch_63 = -A(8,1)*scratch_2 + A(8,3) - scratch_5*scratch_62
    scratch_64 = -A(8,1)*scratch_11 + A(8,4) - scratch_12*scratch_62 - scratch_13* &
      scratch_63
    scratch_65 = -A(8,1)*scratch_19 + A(8,5) - scratch_20*scratch_62 - scratch_21* &
      scratch_63 - scratch_22*scratch_64
    scratch_66 = -A(8,1)*scratch_29 + A(8,6) - scratch_30*scratch_62 - scratch_31* &
      scratch_63 - scratch_32*scratch_64 - scratch_33*scratch_65
    scratch_67 = -A(8,1)*scratch_41 + A(8,7) - scratch_42*scratch_62 - scratch_43* &
      scratch_63 - scratch_44*scratch_64 - scratch_45*scratch_65 - &
      scratch_46*scratch_66
    scratch_68 = (-A(8,1)*scratch_1 + b(8) - scratch_10*scratch_63 - scratch_18* &
      scratch_64 - scratch_28*scratch_65 - scratch_40*scratch_66 - &
      scratch_54*scratch_67 - scratch_62*scratch_9)/(-A(8,1)*scratch_55 &
      + A(8,8) - scratch_56*scratch_62 - scratch_57*scratch_63 - &
      scratch_58*scratch_64 - scratch_59*scratch_65 - scratch_60* &
      scratch_66 - scratch_61*scratch_67)
    scratch_69 = scratch_12 - scratch_13*scratch_5
    scratch_70 = scratch_20 - scratch_21*scratch_5 - scratch_22*scratch_69
    scratch_71 = scratch_30 - scratch_31*scratch_5 - scratch_32*scratch_69 - scratch_33* &
      scratch_70
    scratch_72 = scratch_42 - scratch_43*scratch_5 - scratch_44*scratch_69 - scratch_45* &
      scratch_70 - scratch_46*scratch_71
    scratch_73 = -scratch_13*scratch_22 + scratch_21
    scratch_74 = -scratch_13*scratch_32 + scratch_31 - scratch_33*scratch_73
    scratch_75 = -scratch_13*scratch_44 + scratch_43 - scratch_45*scratch_73 - scratch_46 &
      *scratch_74
    scratch_76 = -scratch_22*scratch_33 + scratch_32
    scratch_77 = -scratch_22*scratch_45 + scratch_44 - scratch_46*scratch_76
    scratch_78 = -scratch_33*scratch_46 + scratch_45

    x(1) = scratch_1 - scratch_10*scratch_6 - scratch_14*scratch_18 - scratch_23* &
      scratch_28 - scratch_3*scratch_9 - scratch_34*scratch_40 - &
      scratch_47*scratch_54 - scratch_68*(-scratch_14*scratch_58 - &
      scratch_23*scratch_59 - scratch_3*scratch_56 - scratch_34* &
      scratch_60 - scratch_47*scratch_61 + scratch_55 - scratch_57* &
      scratch_6)
    x(2) = -scratch_10*scratch_5 - scratch_18*scratch_69 - scratch_28*scratch_70 - &
      scratch_40*scratch_71 - scratch_54*scratch_72 - scratch_68*( &
      -scratch_5*scratch_57 + scratch_56 - scratch_58*scratch_69 - &
      scratch_59*scratch_70 - scratch_60*scratch_71 - scratch_61* &
      scratch_72) + scratch_9
    x(3) = scratch_10 - scratch_13*scratch_18 - scratch_28*scratch_73 - scratch_40* &
      scratch_74 - scratch_54*scratch_75 - scratch_68*(-scratch_13* &
      scratch_58 + scratch_57 - scratch_59*scratch_73 - scratch_60* &
      scratch_74 - scratch_61*scratch_75)
    x(4) = scratch_18 - scratch_22*scratch_28 - scratch_40*scratch_76 - scratch_54* &
      scratch_77 - scratch_68*(-scratch_22*scratch_59 + scratch_58 - &
      scratch_60*scratch_76 - scratch_61*scratch_77)
    x(5) = scratch_28 - scratch_33*scratch_40 - scratch_54*scratch_78 - scratch_68* &
      (-scratch_33*scratch_60 + scratch_59 - scratch_61*scratch_78)
    x(6) = scratch_40 - scratch_46*scratch_54 - scratch_68*(-scratch_46*scratch_61 &
      + scratch_60)
    x(7) = scratch_54 - scratch_61*scratch_68
    x(8) = scratch_68
  end subroutine gauss_jordan_solve
end module gauss_jordan
