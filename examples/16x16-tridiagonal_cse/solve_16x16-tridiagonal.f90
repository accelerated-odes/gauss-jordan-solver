module gauss_jordan
  implicit none

contains

  subroutine gauss_jordan_solve(A, x, b)
    double precision, dimension(16,16), intent(in) :: A
    double precision, dimension(16), intent(out) :: x
    double precision, dimension(16), intent(in) :: b
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
    double precision :: scratch_79
    double precision :: scratch_80
    double precision :: scratch_81
    double precision :: scratch_82
    double precision :: scratch_83
    double precision :: scratch_84
    double precision :: scratch_85
    double precision :: scratch_86
    double precision :: scratch_87
    double precision :: scratch_88
    double precision :: scratch_89
    double precision :: scratch_90
    double precision :: scratch_91
    double precision :: scratch_92
    double precision :: scratch_93
    double precision :: scratch_94
    double precision :: scratch_95
    double precision :: scratch_96
    double precision :: scratch_97
    double precision :: scratch_98
    double precision :: scratch_99
    double precision :: scratch_100
    double precision :: scratch_101
    double precision :: scratch_102
    double precision :: scratch_103
    double precision :: scratch_104
    double precision :: scratch_105
    double precision :: scratch_106
    double precision :: scratch_107
    double precision :: scratch_108
    double precision :: scratch_109
    double precision :: scratch_110
    double precision :: scratch_111
    double precision :: scratch_112
    double precision :: scratch_113
    double precision :: scratch_114
    double precision :: scratch_115
    double precision :: scratch_116
    double precision :: scratch_117
    double precision :: scratch_118
    double precision :: scratch_119
    double precision :: scratch_120
    double precision :: scratch_121
    double precision :: scratch_122
    double precision :: scratch_123
    double precision :: scratch_124
    double precision :: scratch_125
    double precision :: scratch_126
    double precision :: scratch_127
    double precision :: scratch_128
    double precision :: scratch_129
    double precision :: scratch_130
    double precision :: scratch_131
    double precision :: scratch_132
    double precision :: scratch_133
    double precision :: scratch_134
    double precision :: scratch_135
    double precision :: scratch_136
    double precision :: scratch_137
    double precision :: scratch_138
    double precision :: scratch_139
    double precision :: scratch_140
    double precision :: scratch_141
    double precision :: scratch_142
    double precision :: scratch_143
    double precision :: scratch_144
    double precision :: scratch_145
    double precision :: scratch_146
    double precision :: scratch_147
    double precision :: scratch_148
    double precision :: scratch_149
    double precision :: scratch_150

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
    scratch_36 = 1.0/(-A(8,7)*scratch_35 + A(8,8))
    scratch_37 = scratch_36*(-A(8,7)*scratch_29 + b(8))
    scratch_38 = scratch_35*scratch_37
    scratch_39 = scratch_27*scratch_38
    scratch_40 = scratch_20*scratch_39
    scratch_41 = scratch_14*scratch_40
    scratch_42 = scratch_41*scratch_9
    scratch_43 = scratch_42*scratch_5
    scratch_44 = A(8,9)*scratch_36
    scratch_45 = 1.0/(-A(9,8)*scratch_44 + A(9,9))
    scratch_46 = scratch_45*(-A(9,8)*scratch_37 + b(9))
    scratch_47 = scratch_44*scratch_46
    scratch_48 = scratch_35*scratch_47
    scratch_49 = scratch_27*scratch_48
    scratch_50 = scratch_20*scratch_49
    scratch_51 = scratch_14*scratch_50
    scratch_52 = scratch_51*scratch_9
    scratch_53 = scratch_5*scratch_52
    scratch_54 = A(9,10)*scratch_45
    scratch_55 = 1.0/(A(10,10) - A(10,9)*scratch_54)
    scratch_56 = scratch_55*(-A(10,9)*scratch_46 + b(10))
    scratch_57 = scratch_54*scratch_56
    scratch_58 = scratch_44*scratch_57
    scratch_59 = scratch_35*scratch_58
    scratch_60 = scratch_27*scratch_59
    scratch_61 = scratch_20*scratch_60
    scratch_62 = scratch_14*scratch_61
    scratch_63 = scratch_62*scratch_9
    scratch_64 = scratch_5*scratch_63
    scratch_65 = A(10,11)*scratch_55
    scratch_66 = 1.0/(-A(11,10)*scratch_65 + A(11,11))
    scratch_67 = scratch_66*(-A(11,10)*scratch_56 + b(11))
    scratch_68 = scratch_65*scratch_67
    scratch_69 = scratch_54*scratch_68
    scratch_70 = scratch_44*scratch_69
    scratch_71 = scratch_35*scratch_70
    scratch_72 = scratch_27*scratch_71
    scratch_73 = scratch_20*scratch_72
    scratch_74 = scratch_14*scratch_73
    scratch_75 = scratch_74*scratch_9
    scratch_76 = scratch_5*scratch_75
    scratch_77 = A(11,12)*scratch_66
    scratch_78 = 1.0/(-A(12,11)*scratch_77 + A(12,12))
    scratch_79 = scratch_78*(-A(12,11)*scratch_67 + b(12))
    scratch_80 = scratch_77*scratch_79
    scratch_81 = scratch_65*scratch_80
    scratch_82 = scratch_54*scratch_81
    scratch_83 = scratch_44*scratch_82
    scratch_84 = scratch_35*scratch_83
    scratch_85 = scratch_27*scratch_84
    scratch_86 = scratch_20*scratch_85
    scratch_87 = scratch_14*scratch_86
    scratch_88 = scratch_87*scratch_9
    scratch_89 = scratch_5*scratch_88
    scratch_90 = A(12,13)*scratch_78
    scratch_91 = 1.0/(-A(13,12)*scratch_90 + A(13,13))
    scratch_92 = scratch_91*(-A(13,12)*scratch_79 + b(13))
    scratch_93 = scratch_90*scratch_92
    scratch_94 = scratch_77*scratch_93
    scratch_95 = scratch_65*scratch_94
    scratch_96 = scratch_54*scratch_95
    scratch_97 = scratch_44*scratch_96
    scratch_98 = scratch_35*scratch_97
    scratch_99 = scratch_27*scratch_98
    scratch_100 = scratch_20*scratch_99
    scratch_101 = scratch_100*scratch_14
    scratch_102 = scratch_101*scratch_9
    scratch_103 = scratch_102*scratch_5
    scratch_104 = A(13,14)*scratch_91
    scratch_105 = 1.0/(-A(14,13)*scratch_104 + A(14,14))
    scratch_106 = scratch_105*(-A(14,13)*scratch_92 + b(14))
    scratch_107 = scratch_104*scratch_106
    scratch_108 = scratch_107*scratch_90
    scratch_109 = scratch_108*scratch_77
    scratch_110 = scratch_109*scratch_65
    scratch_111 = scratch_110*scratch_54
    scratch_112 = scratch_111*scratch_44
    scratch_113 = scratch_112*scratch_35
    scratch_114 = scratch_113*scratch_27
    scratch_115 = scratch_114*scratch_20
    scratch_116 = scratch_115*scratch_14
    scratch_117 = scratch_116*scratch_9
    scratch_118 = scratch_117*scratch_5
    scratch_119 = A(14,15)*scratch_105
    scratch_120 = 1.0/(-A(15,14)*scratch_119 + A(15,15))
    scratch_121 = scratch_120*(-A(15,14)*scratch_106 + b(15))
    scratch_122 = scratch_119*scratch_121
    scratch_123 = scratch_104*scratch_122
    scratch_124 = scratch_123*scratch_90
    scratch_125 = scratch_124*scratch_77
    scratch_126 = scratch_125*scratch_65
    scratch_127 = scratch_126*scratch_54
    scratch_128 = scratch_127*scratch_44
    scratch_129 = scratch_128*scratch_35
    scratch_130 = scratch_129*scratch_27
    scratch_131 = scratch_130*scratch_20
    scratch_132 = scratch_131*scratch_14
    scratch_133 = scratch_132*scratch_9
    scratch_134 = scratch_133*scratch_5
    scratch_135 = A(15,16)*scratch_120
    scratch_136 = (-A(16,15)*scratch_121 + b(16))/(-A(16,15)*scratch_135 + A(16,16))
    scratch_137 = scratch_135*scratch_136
    scratch_138 = scratch_119*scratch_137
    scratch_139 = scratch_104*scratch_138
    scratch_140 = scratch_139*scratch_90
    scratch_141 = scratch_140*scratch_77
    scratch_142 = scratch_141*scratch_65
    scratch_143 = scratch_142*scratch_54
    scratch_144 = scratch_143*scratch_44
    scratch_145 = scratch_144*scratch_35
    scratch_146 = scratch_145*scratch_27
    scratch_147 = scratch_146*scratch_20
    scratch_148 = scratch_14*scratch_147
    scratch_149 = scratch_148*scratch_9
    scratch_150 = scratch_149*scratch_5

    x(1) = scratch_1 + scratch_103*scratch_2 - scratch_118*scratch_2 - scratch_13* &
      scratch_2 + scratch_134*scratch_2 - scratch_150*scratch_2 + &
      scratch_19*scratch_2 - scratch_2*scratch_26 + scratch_2* &
      scratch_34 - scratch_2*scratch_4 - scratch_2*scratch_43 + &
      scratch_2*scratch_53 - scratch_2*scratch_64 + scratch_2* &
      scratch_76 + scratch_2*scratch_8 - scratch_2*scratch_89
    x(2) = -scratch_103 + scratch_118 + scratch_13 - scratch_134 + scratch_150 - &
      scratch_19 + scratch_26 - scratch_34 + scratch_4 + scratch_43 - &
      scratch_53 + scratch_64 - scratch_76 - scratch_8 + scratch_89
    x(3) = scratch_102 - scratch_117 - scratch_12 + scratch_133 - scratch_149 + &
      scratch_18 - scratch_25 + scratch_33 - scratch_42 + scratch_52 - &
      scratch_63 + scratch_7 + scratch_75 - scratch_88
    x(4) = -scratch_101 + scratch_11 + scratch_116 - scratch_132 + scratch_148 - &
      scratch_17 + scratch_24 - scratch_32 + scratch_41 - scratch_51 + &
      scratch_62 - scratch_74 + scratch_87
    x(5) = scratch_100 - scratch_115 + scratch_131 - scratch_147 + scratch_16 - &
      scratch_23 + scratch_31 - scratch_40 + scratch_50 - scratch_61 + &
      scratch_73 - scratch_86
    x(6) = scratch_114 - scratch_130 + scratch_146 + scratch_22 - scratch_30 + &
      scratch_39 - scratch_49 + scratch_60 - scratch_72 + scratch_85 - &
      scratch_99
    x(7) = -scratch_113 + scratch_129 - scratch_145 + scratch_29 - scratch_38 + &
      scratch_48 - scratch_59 + scratch_71 - scratch_84 + scratch_98
    x(8) = scratch_112 - scratch_128 + scratch_144 + scratch_37 - scratch_47 + &
      scratch_58 - scratch_70 + scratch_83 - scratch_97
    x(9) = -scratch_111 + scratch_127 - scratch_143 + scratch_46 - scratch_57 + &
      scratch_69 - scratch_82 + scratch_96
    x(10) = scratch_110 - scratch_126 + scratch_142 + scratch_56 - scratch_68 + &
      scratch_81 - scratch_95
    x(11) = -scratch_109 + scratch_125 - scratch_141 + scratch_67 - scratch_80 + &
      scratch_94
    x(12) = scratch_108 - scratch_124 + scratch_140 + scratch_79 - scratch_93
    x(13) = -scratch_107 + scratch_123 - scratch_139 + scratch_92
    x(14) = scratch_106 - scratch_122 + scratch_138
    x(15) = scratch_121 - scratch_137
    x(16) = scratch_136
  end subroutine gauss_jordan_solve
end module gauss_jordan
