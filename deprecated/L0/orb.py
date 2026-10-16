#prefix = r"/Volumes/ASM2464_ZJS/GENA-RAW/"
#prefix = r"/Users/shawn/Desktop/deepsky/"
#prefix = r"/mnt/e/GENA-RAW/"
prefix = r"/home/zhengjs/"

NameList = [
    #r"CE7ZJ_1_202407-202410.raw",
    #r"CE7ZJ_2_202410-202411_01.raw",
    #r"CE7ZJ_2_202410-202411_02.raw",
    #r"CE7ZJ_2_202410-202411_03.raw",
    #r"CE7ZJ_3_202412-.raw",
    #r"CE7ZJ_4_20250221T120324-20250221T124727.raw",
    #r"CE7ZJ_5_20250223T115438-20250223T121250.raw",
    #r"CE7ZJ_6_ENA_202502240913_202503011545.raw",
    r"CE7ZJ_8_20250318T005939-20250318T041103.raw",
    #r"CE7ZJ-ZP-PLD-CMENAVCID-LVDS-20250402T104139-20250402T165225.raw",
    #r"CE7ZJ-ZP-PLD-CMENAVCID-LVDS-20250402T170723-20250403T004915.raw",
    #r"CE7ZJ-ZP-PLD-CMENAVCID-LVDS-20250427T184534-20250427T191206.raw"
    ]

InputList = [prefix + name for name in NameList]
Parse_data = 1
Parse_parameter = 0
PartialOutput = 1
L0 = 0
L1 = 0
#triggered one
L1A = 0
L1B = 0
L1C = 0
L1D = 0
L2 = 1
L2A = 1
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


