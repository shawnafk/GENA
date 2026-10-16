prefix = r"/Users/shawn/Desktop/0224/"

NameList = [
    r"ENA_202502240913_202503011545"    
]

InputList = [prefix + name for name in NameList]
Parse_data = 1
Parse_parameter = 1
PartialOutput = 1
L0 = 0
L1 = 0
#triggered one
L1A = 0
L1B = 0
L1C = 0
L1D = 1
L2 = 0
L2A = 0
L2B = 0


#ns limit need divided by 2
A_X_LEFT_NS  = -50.44742
A_X_RIGHT_NS = 54.91726

A_Y_LEFT_NS  = -34.71570
A_Y_RIGHT_NS = 23.59747
 
B_X_LEFT_NS  = -51.18360
B_X_RIGHT_NS = 53.40479

#SN_ID 1 orbit 2 ground
SN_ID = 1
if SN_ID == 1:
    #B_Y_LEFT_NS  = -16.68598
    #B_Y_RIGHT_NS = 41.62719
    B_Y_LEFT_NS  = -34.71570
    B_Y_RIGHT_NS = 23.59747

if SN_ID == 2:
    B_Y_LEFT_NS  = -34.71570
    B_Y_RIGHT_NS = 23.59747


