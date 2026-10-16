prefix = r"/Users/shawn/Desktop/Ground_represent/"

NameList = [
    r"A/20241128_3_7keV_H_A_21_600_150_100_TH_0_BT_P5_Z_P55_Y_N3",
    #r"A/20241112_3keV_H_A_4_400_150_100",
    r"A/20241204_10keV_O_A_22_600_150_100_TH_0_BT_P5_Z_P62_Y_0",
    ]

InputList = [prefix + name for name in NameList]
Parse_data = 1
Parse_parameter = 1
PartialOutput = 1
L0 = 0
L1 = 0
#triggered one
L1A = 1
L1B = 0
L1C = 0
L1D = 0
L2 = 1
L2A = 1
L2B = 1


#ns limit need divided by 2
A_X_LEFT_NS  = -50.44742
A_X_RIGHT_NS = 54.91726

A_Y_LEFT_NS  = -34.71570
A_Y_RIGHT_NS = 23.59747
 
B_X_LEFT_NS  = -51.18360
B_X_RIGHT_NS = 53.40479

#SN_ID 1 orbit 2 ground
SN_ID = 2
if SN_ID == 1:
    #B_Y_LEFT_NS  = -16.68598
    #B_Y_RIGHT_NS = 41.62719
    B_Y_LEFT_NS  = -34.71570
    B_Y_RIGHT_NS = 23.59747

if SN_ID == 2:
    B_Y_LEFT_NS  = -34.71570
    B_Y_RIGHT_NS = 23.59747


