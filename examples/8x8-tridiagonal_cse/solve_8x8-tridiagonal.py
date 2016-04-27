import numpy as np

def gauss_jordan_solve(A, b):
    x = np.empty(8, dtype=np.float64)

    scratch_0 = 1/A[0][0]
    scratch_1 = b[0]*scratch_0
    scratch_2 = A[0][1]*scratch_0
    scratch_3 = 1/(-A[1][0]*scratch_2 + A[1][1])
    scratch_4 = scratch_3*(-A[1][0]*scratch_1 + b[1])
    scratch_5 = A[1][2]*scratch_3
    scratch_6 = 1/(-A[2][1]*scratch_5 + A[2][2])
    scratch_7 = scratch_6*(-A[2][1]*scratch_4 + b[2])
    scratch_8 = scratch_5*scratch_7
    scratch_9 = A[2][3]*scratch_6
    scratch_10 = 1/(-A[3][2]*scratch_9 + A[3][3])
    scratch_11 = scratch_10*(-A[3][2]*scratch_7 + b[3])
    scratch_12 = scratch_11*scratch_9
    scratch_13 = scratch_12*scratch_5
    scratch_14 = A[3][4]*scratch_10
    scratch_15 = 1/(-A[4][3]*scratch_14 + A[4][4])
    scratch_16 = scratch_15*(-A[4][3]*scratch_11 + b[4])
    scratch_17 = scratch_14*scratch_16
    scratch_18 = scratch_17*scratch_9
    scratch_19 = scratch_18*scratch_5
    scratch_20 = A[4][5]*scratch_15
    scratch_21 = 1/(-A[5][4]*scratch_20 + A[5][5])
    scratch_22 = scratch_21*(-A[5][4]*scratch_16 + b[5])
    scratch_23 = scratch_20*scratch_22
    scratch_24 = scratch_14*scratch_23
    scratch_25 = scratch_24*scratch_9
    scratch_26 = scratch_25*scratch_5
    scratch_27 = A[5][6]*scratch_21
    scratch_28 = 1/(-A[6][5]*scratch_27 + A[6][6])
    scratch_29 = scratch_28*(-A[6][5]*scratch_22 + b[6])
    scratch_30 = scratch_27*scratch_29
    scratch_31 = scratch_20*scratch_30
    scratch_32 = scratch_14*scratch_31
    scratch_33 = scratch_32*scratch_9
    scratch_34 = scratch_33*scratch_5
    scratch_35 = A[6][7]*scratch_28
    scratch_36 = (-A[7][6]*scratch_29 + b[7])/(-A[7][6]*scratch_35 + A[7][7])
    scratch_37 = scratch_35*scratch_36
    scratch_38 = scratch_27*scratch_37
    scratch_39 = scratch_20*scratch_38
    scratch_40 = scratch_14*scratch_39
    scratch_41 = scratch_40*scratch_9
    scratch_42 = scratch_41*scratch_5

    x[0] = scratch_1 - scratch_13*scratch_2 + scratch_19*scratch_2 - scratch_2*scratch_26 + scratch_2*scratch_34 - scratch_2*scratch_4 - scratch_2*scratch_42 + scratch_2*scratch_8
    x[1] = scratch_13 - scratch_19 + scratch_26 - scratch_34 + scratch_4 + scratch_42 - scratch_8
    x[2] = -scratch_12 + scratch_18 - scratch_25 + scratch_33 - scratch_41 + scratch_7
    x[3] = scratch_11 - scratch_17 + scratch_24 - scratch_32 + scratch_40
    x[4] = scratch_16 - scratch_23 + scratch_31 - scratch_39
    x[5] = scratch_22 - scratch_30 + scratch_38
    x[6] = scratch_29 - scratch_37
    x[7] = scratch_36

    return x
